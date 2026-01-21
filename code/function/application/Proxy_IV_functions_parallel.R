suppressPackageStartupMessages({
  library(dplyr); library(data.table); library(Matrix); library(glmnet); library(tidyr)
  library(qs); library(sandwich); library(lmtest); library(splines)
  library(onlineFDR); library(MASS); library(igraph); library(ggplot2); library(future)
})

# ---------------- helpers ----------------
calc_intervention_score <- function(nodes, exclude, desc_sets) {
  if (is.null(nodes) || length(nodes) == 0) return(list(nu=0,mu=0))
  gset <- names(desc_sets); which_gene <- which(gset %in% nodes)
  scores <- lapply(desc_sets[which_gene], function(d) {
    if (length(d) == 0) return(list(nu=0,mu=0))
    list(nu=sum(!names(d) %in% exclude), mu=sum(d[!names(d) %in% exclude]))
  })
  names(scores) <- nodes; scores
}

p_adjust_generic <- function(p, method = "BH") {
  if (!length(p)) return(p)
  base_methods <- c("BH","BY","bonferroni","holm","hochberg","hommel","fdr")
  if (method %in% base_methods) return(p.adjust(p, method = method))
  if (method == "qvalue" && requireNamespace("qvalue", quietly = TRUE))
    return(qvalue::qvalue(p)$qvalues)
  p.adjust(p, method = "BH")
}

safe_glm_nb <- function(formula, data, off_vec, x=TRUE, y=TRUE, maxit1=50, maxit2=100) {
  off_name <- ".__off__"; while (off_name %in% names(data)) off_name <- paste0(off_name, "_")
  data[[off_name]] <- as.numeric(off_vec)
  form2 <- stats::update(formula, paste0("~ . + offset(", off_name, ")"))
  fit <- tryCatch(MASS::glm.nb(form2, data=data, control=glm.control(maxit=maxit1), x=x, y=y), error=function(e) NULL)
  bad <- is.null(fit) || !is.finite(as.numeric(fit$theta)) || as.numeric(fit$theta) <= 1e-6
  if (bad) fit <- tryCatch(MASS::glm.nb(form2, data=data, control=glm.control(maxit=maxit2), x=x, y=y), error=function(e) NULL)
  bad <- is.null(fit) || !is.finite(as.numeric(fit$theta)) || as.numeric(fit$theta) <= 1e-6
  if (bad) { fit <- glm(form2, family=poisson(), data=data, x=x, y=y); attr(fit,"poisson_fallback") <- TRUE; attr(fit,"theta") <- Inf }
  fit
}

.clamp <- function(x, lo=1e-12, hi=1e12) pmax(lo, pmin(hi, x))
.inv   <- function(A) tryCatch(solve(A), error=function(e) MASS::ginv(A))

# ---------------- Proxy IV (Stage-1 uses mclapply + shared X for Poisson) ----------------
proxy_IV <- function(target_node, ancestors, data, D, depth,
                     ridge = 1e-8, family = c("poisson","nb"),
                     cc_df = NULL, cc_cols = NULL,
                     cores = 1) {
  
  family <- match.arg(family)
  tgt <- if (is.numeric(target_node)) colnames(data)[target_node] else target_node
  anc <- ancestors; if (is.numeric(anc)) anc <- colnames(data)[anc]; anc <- as.character(anc)
  stopifnot(length(tgt) == 1L)
  n <- nrow(data)
  if (length(anc) == 0L) {
    return(list(coefs_all=setNames(numeric(0),character(0)),
                se_all=setNames(numeric(0),character(0)),
                p_all=setNames(numeric(0),character(0)),
                theta=setNames(numeric(0),character(0))))
  }
  
  yj      <- as.numeric(data[, tgt])
  off_log <- log(pmax(as.numeric(depth), 1e-12))
  Dj      <- as.numeric(D[, tgt, drop=FALSE])
  
  use_cc <- !is.null(cc_df) && !is.null(cc_cols) && length(cc_cols) > 0
  Z <- if (use_cc) {
    mf_cc <- model.frame(~ . - 1, data = cc_df[, cc_cols, drop = FALSE], na.action = na.pass)
    model.matrix(~ . - 1, data = mf_cc)
  } else {
    NULL
  }
  
  
  ## ---------- Stage-1 design (SHARED, single copy in memory) ----------
  
  X_shared    <- cbind(as.matrix(D[, anc, drop=FALSE]), Z)
  colnames(X_shared) <- c(anc, if (use_cc) colnames(Z) else NULL)
  X1_shared   <- cbind(`(Intercept)` = 1, X_shared)
  
  ## ---------- Stage 1 ----------
  if (family == "poisson") {
    stg1 <- parallel::mclapply(seq_along(anc), function(idx) {
      k  <- anc[idx]
      yk <- as.numeric(data[, k])
      
  
      fit1 <- glm.fit(x = X_shared, y = yk, family = poisson(),
                     offset = off_log, intercept = TRUE)
      
      mu1 <- pmax(1e-12, pmin(1e12, fit1$fitted.values))
      r1  <- yk - mu1
      
      A1  <- crossprod(X1_shared, X1_shared * mu1) + diag(ridge, ncol(X1_shared))
      B1  <- crossprod(X1_shared, X1_shared * (r1^2))
      iA1 <- tryCatch(solve(A1), error = function(e) MASS::ginv(A1))
      Vpi <- iA1 %*% B1 %*% iA1
      
      list(k = k,
           S_hat_k = as.numeric(fit1$linear.predictors) - off_log,
           Vpi = Vpi)   # <-- no X1 per-k returned
    }, mc.cores = max(1, min(cores, length(anc))), mc.preschedule = FALSE)
    
  } else {
    # NB path unchanged (sequential); keeps memory predictable
    stg1 <- lapply(anc, function(k) {
      df1 <- if (use_cc) data.frame(y = as.numeric(data[, k]),
                                    as.data.frame(D[, anc, drop=FALSE]),
                                    as.data.frame(Z)) else
                                      data.frame(y = as.numeric(data[, k]),
                                                 as.data.frame(D[, anc, drop=FALSE]))
      rhs1 <- setdiff(colnames(df1), "y")
      fml1 <- as.formula(paste0("y ~ ", paste(rhs1, collapse = " + ")))
      fit1 <- safe_glm_nb(fml1, data = df1, off_vec = off_log, x = TRUE, y = TRUE)
      
      eta_k <- as.numeric(stats::predict(fit1, type="link"))
      S_hat_k <- eta_k - off_log
      
      X1  <- fit1$x
      mu1 <- pmax(1e-12, pmin(1e12, as.numeric(stats::fitted(fit1))))
      if (isTRUE(attr(fit1, "poisson_fallback"))) {
        r1  <- fit1$y - mu1
        A1  <- crossprod(X1, X1 * mu1) + diag(ridge, ncol(X1))
        B1  <- crossprod(X1, X1 * (r1^2))
        iA1 <- tryCatch(solve(A1), error = function(e) MASS::ginv(A1))
        Vpi <- iA1 %*% B1 %*% iA1
      } else {
        th1 <- as.numeric(fit1$theta); ws1 <- 1/(1 + mu1/th1); wI1 <- mu1*ws1
        U1  <- X1 * as.numeric((fit1$y - mu1)*ws1)
        A1  <- crossprod(X1, X1 * wI1) + diag(ridge, ncol(X1))
        B1  <- crossprod(U1); iA1 <- tryCatch(solve(A1), error=function(e) MASS::ginv(A1))
        Vpi <- iA1 %*% B1 %*% iA1
      }
      list(k = k, S_hat_k = S_hat_k, Vpi = Vpi, X1_nb = X1) # NB keeps per-k X1
    })
  }
  
  ## Collect Stage-1
  S_hat_list <- setNames(vector("list", length(anc)), anc)
  Vpi_list   <- setNames(vector("list", length(anc)), anc)
  for (res in stg1) {
    S_hat_list[[res$k]] <- res$S_hat_k
    Vpi_list[[res$k]]   <- res$Vpi
  }
  S_hat <- as.data.frame(S_hat_list, check.names = FALSE)
  
  ## ---------- Stage 2 (unchanged interface) ----------
  if (use_cc) {
    df2  <- cbind.data.frame(y = yj, D_j = Dj, S_hat, Z)
    rhs2 <- c("D_j", colnames(S_hat), colnames(Z))
  } else {
    df2  <- cbind.data.frame(y = yj, D_j = Dj, S_hat)
    rhs2 <- c("D_j", colnames(S_hat))
  }
  fml2 <- as.formula(paste0("y ~ ", paste(rhs2, collapse = " + ")))
  
  if (family == "poisson") {
    fit2 <- glm(fml2, family = poisson(), offset = off_log, data = df2, x = TRUE, y = TRUE)

    beta_hat <- stats::coef(fit2)
    X2  <- fit2$x
    mu2 <- pmax(1e-12, pmin(1e12, as.numeric(stats::fitted(fit2))))
    r2  <- fit2$y - mu2
    A2  <- crossprod(X2, X2 * mu2) + diag(ridge, ncol(X2))
    B2  <- crossprod(X2, X2 * (r2^2))
    
    ## Use SHARED X1 for all k (no per-k copies)
    C_blocks <- lapply(anc, function(k) {
      w <- mu2 * as.numeric(beta_hat[k])
      - crossprod(X2, X1_shared * w)
    })
    
  } else {
    fit2 <- safe_glm_nb(fml2, data = df2, off_vec = off_log, x = TRUE, y = TRUE)
    beta_hat <- stats::coef(fit2)
    X2  <- fit2$x
    mu2 <- pmax(1e-12, pmin(1e12, as.numeric(stats::fitted(fit2))))
    if (isTRUE(attr(fit2,"poisson_fallback"))) {
      r2 <- fit2$y - mu2
      A2 <- crossprod(X2, X2 * mu2) + diag(ridge, ncol(X2))
      B2 <- crossprod(X2, X2 * (r2^2))
      # Prefer shared design if Poisson-fallback makes them compatible
      if (!is.null(attr(stg1[[1]], "X1_nb"))) {
        C_blocks <- lapply(anc, function(k) {
          w <- mu2 * as.numeric(beta_hat[k])
          - crossprod(X2, stg1[[1]]$X1_nb * w)  # best-effort; NB path can keep per-k X1 if needed
        })
      } else {
        C_blocks <- lapply(anc, function(k) {
          w <- mu2 * as.numeric(beta_hat[k])
          - crossprod(X2, X1_shared * w)
        })
      }
    } else {
      th2 <- as.numeric(fit2$theta); w_info <- (mu2 / (1 + mu2/th2))
      A2  <- crossprod(X2, X2 * (mu2 * (1/(1 + mu2/th2)))) + diag(ridge, ncol(X2))
      U2  <- X2 * as.numeric((fit2$y - mu2) * (1/(1 + mu2/th2)))
      B2  <- crossprod(U2)
      C_blocks <- lapply(anc, function(k) {
        - crossprod(X2, X1_shared * (w_info * as.numeric(beta_hat[k])))
      })
    }
  }
  
  Cmat <- if (length(C_blocks) == 1) C_blocks[[1]] else do.call(cbind, C_blocks)
  Vpi_block <- if (length(Vpi_list) == 1) Vpi_list[[1]] else Matrix::bdiag(Vpi_list)
  iA2 <- tryCatch(solve(A2), error=function(e) MASS::ginv(A2))
  V_MT <- iA2 %*% (B2 + as.matrix(Cmat %*% Vpi_block %*% t(Cmat))) %*% iA2
  se   <- sqrt(pmax(diag(V_MT), 0) + 1e-16)
  
  beta_names <- names(beta_hat)
  anc_idx <- match(anc, beta_names)
  coefs <- beta_hat[anc_idx]; ses <- se[anc_idx]
  zval <- coefs / ses; pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  
  list(coefs_all=setNames(as.numeric(coefs), anc),
       se_all=setNames(as.numeric(ses),   anc),
       p_all=setNames(as.numeric(pval),   anc),
       theta=beta_hat)
}

# ---------------- descendant screening + closure ----------------
build_descendants_cov_adjust <- function(
    X_raw, D, control_ind, alpha_desc = 0.05, depth = NULL, nodes,
    padj_method = "BH", cc_cols = NULL, cc_df = NULL,
    ncores = 1, engine = c("sequential","multisession","multicore"),
    blas_threads = 1
) {
  stopifnot(!is.null(depth), length(depth) == nrow(X_raw))
  engine <- match.arg(engine)
  p <- length(nodes)
  Pmat <- matrix(NA_real_, p, p, dimnames = list(nodes, nodes))
  
  ## --- clamp BLAS/OMP threads on workers to avoid nested parallelism ---
  if (blas_threads == 1) {
    Sys.setenv(
      OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
      BLIS_NUM_THREADS = "1", GOTO_NUM_THREADS = "1"
    )
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1)
      RhpcBLASctl::omp_set_num_threads(1)
    }
    if (requireNamespace("data.table", quietly = TRUE)) {
      data.table::setDTthreads(1)
    }
  }
  
  ## --- choose plan temporarily ---
  if (requireNamespace("future", quietly = TRUE) &&
      requireNamespace("future.apply", quietly = TRUE)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    if (ncores > 1) {
      if (engine == "multicore") future::plan(future::multicore, workers = ncores)
      else if (engine == "multisession") future::plan(future::multisession, workers = ncores)
      else future::plan(future::sequential)
    } else {
      future::plan(future::sequential)
    }
    
    rows <- future.apply::future_lapply(seq_len(p), function(i) {
      Pi <- rep(NA_real_, p)
      perturb_ind <- D[, i] == 1
      for (k in seq_len(p)) if (k != i) {
        x1 <- X_raw[perturb_ind, k]
        x0 <- X_raw[control_ind, k]
        y <- c(x1, x0)
        d <- c(D[perturb_ind, i], D[control_ind, i])
        depth_i <- c(depth[perturb_ind], depth[control_ind])
        
        has_cc <- !is.null(cc_df) && !is.null(cc_cols) && length(cc_cols) > 0
        if (has_cc) {
          Xc <- rbind(cc_df[perturb_ind, cc_cols, drop = FALSE],
                      cc_df[control_ind, cc_cols, drop = FALSE])
          df <- data.frame(y = y, d = d, Xc, check.names = TRUE)
          fml <- y ~ .
        } else {
          df <- data.frame(y = y, d = d)
          fml <- y ~ d
        }
        
        pval <- NA_real_
        fit <- try(glm(fml,
                       offset = log(pmax(depth_i, 1e-12)),
                       family = poisson(), data = df), silent = TRUE)
        if (!inherits(fit, "try-error")) {
          V <- try(sandwich::vcovHC(fit, type = "HC0"), silent = TRUE)
          if (!inherits(V, "try-error")) {
            ct <- try(lmtest::coeftest(fit, vcov = V), silent = TRUE)
            if (!inherits(ct, "try-error") && "d" %in% rownames(ct)) {
              pval <- suppressWarnings(ct["d", "Pr(>|z|)"])
            }
          }
        }
        Pi[k] <- pval
      }
      Pi
    })
    Pmat <- do.call(rbind, rows)
    rownames(Pmat) <- nodes
    colnames(Pmat) <- nodes
    
  } else {
    ## Fallback: sequential loop if future not available
    for (i in seq_len(p)) {
      perturb_ind <- D[, i] == 1
      for (k in seq_len(p)) if (k != i) {
        x1 <- X_raw[perturb_ind, k]
        x0 <- X_raw[control_ind, k]
        y <- c(x1, x0)
        d <- c(D[perturb_ind, i], D[control_ind, i])
        depth_i <- c(depth[perturb_ind], depth[control_ind])
        
        has_cc <- !is.null(cc_df) && !is.null(cc_cols) && length(cc_cols) > 0
        if (has_cc) {
          Xc <- rbind(cc_df[perturb_ind, cc_cols, drop = FALSE],
                      cc_df[control_ind, cc_cols, drop = FALSE])
          df <- data.frame(y = y, d = d, Xc, check.names = TRUE)
          fml <- y ~ .
        } else {
          df <- data.frame(y = y, d = d)
          fml <- y ~ d
        }
        
        pval <- NA_real_
        fit <- try(glm(fml,
                       offset = log(pmax(depth_i, 1e-12)),
                       family = poisson(), data = df), silent = TRUE)
        if (!inherits(fit, "try-error")) {
          V <- try(sandwich::vcovHC(fit, type = "HC0"), silent = TRUE)
          if (!inherits(V, "try-error")) {
            ct <- try(lmtest::coeftest(fit, vcov = V), silent = TRUE)
            if (!inherits(ct, "try-error") && "d" %in% rownames(ct)) {
              pval <- suppressWarnings(ct["d", "Pr(>|z|)"])
            }
          }
        }
        Pmat[i, k] <- pval
      }
    }
  }
  
  ## p-adjust
  pv_idx <- which(is.finite(Pmat) & !is.na(Pmat))
  pvec <- Pmat[pv_idx]
  P_adj <- Pmat
  if (length(pvec)) {
    P_adj[pv_idx] <- p_adjust_generic(pvec, method = padj_method)
  }
  
  ## descendant sets + scores
  descendant_set <- setNames(vector("list", p), nodes)
  descendant_set_score <- setNames(vector("list", p), nodes)
  for (i in seq_len(p)) {
    j_all <- setdiff(seq_len(p), i)
    valid_adj <- j_all[is.finite(P_adj[i, j_all]) & !is.na(P_adj[i, j_all])]
    keep <- valid_adj[P_adj[i, valid_adj] < alpha_desc]
    descendant_set[[nodes[i]]] <- if (length(keep)) colnames(P_adj)[keep] else character(0)
    srow <- as.numeric(Pmat[i, j_all]); srow[!is.finite(srow) | is.na(srow)] <- 1
    names(srow) <- colnames(P_adj)[j_all]
    descendant_set_score[[nodes[i]]] <- -log10(srow + 1e-16)
  }
  
  get_full_descendants <- function(j_name, adj_list) {
    seen <- j_name; frontier <- setdiff(unique(adj_list[[j_name]]), seen); full <- character(0)
    while (length(frontier) > 0) {
      full <- union(full, frontier); seen <- union(seen, frontier)
      next_frontier <- unique(unlist(adj_list[intersect(frontier, names(adj_list))], use.names = FALSE))
      frontier <- setdiff(next_frontier, seen)
    }
    full
  }
  
  descendant_set_closure <- setNames(vector("list", length(descendant_set)), names(descendant_set))
  for (j in names(descendant_set)) descendant_set_closure[[j]] <- get_full_descendants(j, descendant_set)
  
  descendant_set_score_closure <- mapply(function(s, d) s[d], descendant_set_score,
                                         descendant_set_closure, SIMPLIFY = FALSE)
  
  list(
    P_raw = Pmat, P_adj = P_adj,
    descendant_set = descendant_set,
    descendant_set_closure = descendant_set_closure,
    descendant_set_score_closure = descendant_set_score_closure,
    padj_method = padj_method
  )
}


# ---------------- main DAG search (per-call cores) ----------------
dag_search_perturbseq <- function(
    desc_sets, data_raw, D, depth,
    family = c("poisson","nb"),
    cc_df = NULL, cc_cols = NULL,
    alpha_parent = 0.05,
    padj_method_parent = 'onlineBH',
    cores = NULL,blas_threads=1  # used by proxy_IV Stage-1 Poisson via mclapply
) {
  
  
  if (blas_threads == 1) {
    Sys.setenv(
      OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
      BLIS_NUM_THREADS = "1", GOTO_NUM_THREADS = "1"
    )
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1)
      RhpcBLASctl::omp_set_num_threads(1)
    }
    if (requireNamespace("data.table", quietly = TRUE)) {
      data.table::setDTthreads(1)
    }
  }
  
  
  
  
  
  family <- match.arg(family)
  p <- ncol(data_raw); gene_set <- colnames(data_raw)
  
  score <- calc_intervention_score(gene_set, exclude=NULL, desc_sets=desc_sets)
  nu <- sapply(score, function(sc) sc$nu); mu <- sapply(score, function(sc) sc$mu)
  pi_hat <- rep(NA_character_, p)
  if (sum(nu == min(nu)) == 1) pi_hat[p] <- names(which.min(nu)) else {
    ties <- names(nu)[nu == min(nu)]; pi_hat[p] <- names(which.min(mu[ties]))
  }
  
  E_hat <- vector("list", p); E_all <- vector("list", p)
  cum_dat <- NULL; bb <- 0
  for (j in seq(p - 1, 1)) {
    cat(paste0('Node ',j,' done','\n'))
    bb <- bb + 1
    already <- pi_hat[(j + 1):p]; target <- pi_hat[j + 1]
    anc_all <- gene_set[sapply(desc_sets, function(s) target %in% names(s))]
    anc <- setdiff(anc_all, already)
    
    if (length(anc) == 0L) {
      empty_df <- data.frame(from=character(0), to=character(0), coef=numeric(0), p_raw=numeric(0))
      E_hat[[j + 1]] <- empty_df; E_all[[j + 1]] <- empty_df
    } else {
      iv_res <- proxy_IV(target_node=target, ancestors=anc, data=data_raw, D=D, depth=depth,
                         family=family, cc_df=cc_df, cc_cols=cc_cols, cores=cores)
      if (padj_method_parent == "onlineBH") {
        pvec_all <- iv_res$p_all[anc]; coef_all <- iv_res$coefs_all[anc]
        E_all[[j + 1]] <- data.frame(from=anc, to=rep(target, length(anc)),
                                     coef=as.numeric(coef_all), p_raw=as.numeric(pvec_all),
                                     row.names=NULL, check.names=FALSE)
        batch <- data.frame(id=anc, batch=bb, pval=pvec_all)
        cum_dat <- rbind(cum_dat, batch)
        online_res <- onlineFDR::BatchBH(cum_dat, alpha = alpha_parent)
        parents <- (online_res %>% dplyr::filter(batch==bb & R==1))$id
        if (length(parents) > 0L) {
          E_hat[[j + 1]] <- data.frame(from=parents, to=rep(target, length(parents)),
                                       coef=as.numeric(coef_all[parents]), p_raw=as.numeric(pvec_all[parents]),
                                       row.names=NULL, check.names=FALSE)
        } else {
          E_hat[[j + 1]] <- data.frame(from=character(0), to=character(0), coef=numeric(0), p_raw=numeric(0))
        }
      } else {
        pvec_all <- iv_res$p_all[anc]; coef_all <- iv_res$coefs_all[anc]
        padj_all <- p_adjust_generic(pvec_all, method = padj_method_parent)
        E_all[[j + 1]] <- data.frame(from=anc, to=rep(target, length(anc)),
                                     coef=as.numeric(coef_all), p_raw=as.numeric(pvec_all),
                                     p_adj=as.numeric(padj_all), row.names=NULL, check.names=FALSE)
        pa_mask <- is.finite(padj_all) & !is.na(padj_all) & (padj_all < alpha_parent)
        parents <- names(padj_all)[pa_mask]
        if (length(parents) > 0L) {
          E_hat[[j + 1]] <- data.frame(from=parents, to=rep(target, length(parents)),
                                       coef=as.numeric(coef_all[parents]), p_raw=as.numeric(pvec_all[parents]),
                                       p_adj=as.numeric(padj_all[parents]), row.names=NULL, check.names=FALSE)
        } else {
          E_hat[[j + 1]] <- data.frame(from=character(0), to=character(0),
                                       coef=numeric(0), p_raw=numeric(0), p_adj=numeric(0))
        }
      }
    }
    
    remaining_nodes <- setdiff(gene_set, already)
    score_minus <- calc_intervention_score(remaining_nodes, exclude=already, desc_sets=desc_sets)
    nu_minus <- sapply(score_minus, function(sc) sc$nu); mu_minus <- sapply(score_minus, function(sc) sc$mu)
    if (sum(nu_minus == min(nu_minus)) == 1) pi_hat[j] <- names(which.min(nu_minus)) else {
      ties <- names(nu_minus)[nu_minus == min(nu_minus)]
      pi_hat[j] <- names(which.min(mu_minus[ties]))
    }
    
    if (((p - j) %% 3) == 0) gc()
  }
  
  list(causal_order=pi_hat, edge_structure=E_hat, edge_structure_all=E_all)
}

# ---------------- summarizing ----------------
edges_from_dagsearch <- function(dag_obj) {
  out <- lapply(dag_obj$edge_structure, function(df) if (is.null(df) || !nrow(df)) NULL else df[complete.cases(df), , drop=FALSE])
  out <- do.call(rbind, out)
  if (is.null(out) || !nrow(out)) data.frame(from=character(0), to=character(0), coef=numeric(0)) else unique(out)
}

compute_ancestors_map <- function(nodes, edges_df) {
  parents_map <- setNames(lapply(nodes, function(nm) character(0)), nodes)
  if (nrow(edges_df)) { split_by_to <- split(edges_df$from, edges_df$to); for (to in names(split_by_to)) parents_map[[to]] <- unique(split_by_to[[to]]) }
  ancestors <- setNames(lapply(nodes, function(nm) character(0)), nodes)
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    for (v in nodes) {
      cur <- ancestors[[v]]
      new <- unique(c(parents_map[[v]], unlist(ancestors[parents_map[[v]]], use.names=FALSE)))
      new <- setdiff(new, v)
      if (length(setdiff(new, cur)) > 0) { ancestors[[v]] <- unique(c(cur, new)); changed <- TRUE }
    }
  }
  list(parents_map=parents_map, ancestors_map=ancestors)
}
