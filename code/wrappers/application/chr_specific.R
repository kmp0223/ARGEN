DIR='/Users/kpark243/Library/CloudStorage/Box-Box/Penn/Research/ARGEN/'
ncores=5
alpha_desc=0.05
alpha_IV = 0.05
chr=23




source(paste0(DIR,'/code/function/application/Proxy_IV_functions_parallel.R'))
chrlist=paste0('chr',c(1:22,"X"))
options(future.globals.maxSize = 10 * 1024^3) 

  cat(paste0('Running codes for ',chrlist[chr],'\n'))
  set.seed(1)
  
  # -------------------------------
  # Paths / output
  # -------------------------------
  paste0(DIR,'data/',chrlist[chr])
  DAG_DIR <- paste0(DIR,'results/application/',chrlist[chr])
  if (!dir.exists(DAG_DIR)) dir.create(DAG_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # -------------------------------
  # Load data
  # -------------------------------
  X_raw  <- qs::qread(paste0(DIR,'data/',chrlist[chr],"/X_raw_",chrlist[chr],'.qs'))
  meta_data <- qs::qread(paste0(DIR,'data/',chrlist[chr],"/meta_data_",chrlist[chr],'.qs'))
  

  perturbation_unique <- unique(meta_data$gene_id)
  meta_data$gene_id <- factor(meta_data$gene_id,
                                            levels = unique(meta_data$gene_id))
  meta_data$gem_group<-factor(meta_data$gem_group)
  
  D <- model.matrix(~ gene_id - 1, data = meta_data)
  colnames(D) <- sub("^gene_id", "", colnames(D))
  D <- D[, -c(length(perturbation_unique))]
  control_ind <- which(rowSums(D) == 0)
  
  
  C <- meta_data$guide_UMI_sum / mean(meta_data$guide_UMI_sum)
  
  
  stopifnot(nrow(X_raw) == nrow(D))
  nodes <- colnames(X_raw)
  p <- length(nodes)
  
  cat(paste0('Estimating Ancestors and Descendants...','\n'))
  
  ds_out <- build_descendants_cov_adjust(
    X_raw = X_raw,
    D = D, 
    control_ind = control_ind,
    alpha_desc = alpha_desc, nodes = nodes,
    padj_method = "BH",
    cc_cols=c('mitopercent', "gem_group"),
    cc_df=meta_data %>% as.data.frame,depth=C,
    ncores = ncores, engine = "multisession", blas_threads = 1
  )
  
  setwd(DAG_DIR)
  qs::qsave(ds_out,'ds_out_covadjust.qs')

  
  cat(paste0('Estimating Parents...','\n'))
  dag_obj <- dag_search_perturbseq(
    desc_sets = ds_out$descendant_set_score_closure,
    data_raw  = X_raw ,
    D = D ,depth =C,family = 'poisson',cc_df=data.frame(meta_data),
    cc_cols = c('mitopercent', "gem_group"),alpha_parent=alpha_IV,
    padj_method_parent='onlineBH',
    cores = ncores, blas_threads = 1
  )
  
  edges_k <- edges_from_dagsearch(dag_obj)
  maps_k  <- compute_ancestors_map(nodes, edges_k)
  dag_artifact <- list(
    dag_object    = dag_obj,
    dag_edges     = edges_k,
    parents_map   = maps_k$parents_map,
    ancestors_map = maps_k$ancestors_map
  )
  
  
  saveRDS(dag_artifact, file = file.path(DAG_DIR,  'dag_0.05BH_0.05onlineBH.rds'))
  cat(paste0(chrlist[chr],' complete!','\n'))



