# ARGEN : Causal Network Recovery in Perturb-seq Experiments Using Proxy and Instrumental Variables
![ARGEN diagram](/figures/main_summary.png)

## ARGEN Usage
### Preparation
First of all, before cloning the ARGEN github package, go to the right directory that you would like to implement ARGEN. In the cmd terminal, do
```
cd {ARGEN directory}
```
then go on to the next step.
#### 1. Repository clone
For cloning the github repository, again on the cmd terminal, run the linux code 
```
git clone git@github.com:kmp0223/ARGEN.git
```
For R, we need the requirements as below : 
-   R: [R installation](https://www.r-project.org)  (>=4.5.3)

#### 2. Install and load required R packages
In R, run the following commands to install the packages required for running **ARGEN** (skip any packages that are already installed):
```r
# ============================================================
# 1. Bootstrap installers
# ============================================================
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages('dplyr')

# ============================================================
# 2. Install everything at current/native versions. This also
#    pulls in AnnotationDbi/RSQLite/DBI/etc. as dependencies.
# ============================================================
pacman::p_load(
  dplyr, data.table, Matrix, glmnet, tidyr,
  sandwich, lmtest, MASS, igraph, ggplot2, future,
  future.apply, purrr, rlang, scales, visNetwork,
  htmlwidgets, tibble, ggraph, tidygraph,
  patchwork, here, RColorBrewer, Gviz, GenomicRanges, TxDb.Hsapiens.UCSC.hg38.knownGene,
  clusterProfiler, org.Hs.eg.db, enrichplot, rtracklayer, onlineFDR, ggpubr, caret,
  hexbin, ggtext, grid, stringr, RhpcBLASctl, Rcpp
)

```
Details about implementing codes and reproducing Figures in the manuscript can be found in the [Tutorial](https://github.com/kmp0223/ARGEN/tree/main/tutorials) of this github. There are nine tutorials provided, each reproducing a specific set of figures. Several notebooks read cached results written by earlier ones, so it's best to run them in order, Tutorial 1 through Tutorial 9:
- Tutorial 1: [Simulation study](https://github.com/kmp0223/ARGEN/tree/main/tutorials/Tutorial1-Fig_2.ipynb) — comparing ARGEN, Naive, and INSPRE (Brown et al., 2025, Nature Communications) on a synthetic DAG (Fig. 2)
- Tutorial 2: [Single-chromosome application](https://github.com/kmp0223/ARGEN/tree/main/tutorials/Tutorial2-Fig_3a_b_c.ipynb) — the MCM3 sub-DAG and intervention effects (Fig. 3a–c)
- Tutorial 3: [Across-chromosome application](https://github.com/kmp0223/ARGEN/tree/main/tutorials/Tutorial3-Fig_3d_e_f_g.ipynb) — ARGEN and INSPRE (Brown et al., 2025, Nature Communications) across all chromosomes (Fig. 3d–g)
- Tutorial 4: [Confounding diagnostics](https://github.com/kmp0223/ARGEN/tree/main/tutorials/Tutorial4-Fig_4a_b.ipynb) — checking the perturbation indicator against GEM group / mitochondrial percentage (Fig. 4a–b)
- Tutorial 5: [Robustness to GEM group](https://github.com/kmp0223/ARGEN/tree/main/tutorials/Tutorial5-Fig_4c.ipynb) — re-fitting ARGEN with GEM group as a covariate (Fig. 4c)
- Tutorial 6: [Robustness to dropping genes](https://github.com/kmp0223/ARGEN/tree/main/tutorials/Tutorial6-Fig_4d_e_f.ipynb) — re-fitting ARGEN after dropping the confounding gene MED12 (Fig. 4d–f)
- Tutorial 7: [3D genomics-based exploration](https://github.com/kmp0223/ARGEN/tree/main/tutorials/Tutorial7-Fig_5a_b.ipynb) — edge strength vs. genomic distance and A/B compartment (Fig. 5a–b)
- Tutorial 8: [External validation](https://github.com/kmp0223/ARGEN/tree/main/tutorials/Tutorial8-Fig_5c_d.ipynb) — validating ARGEN edges against TAD co-membership and ChIP-seq overlap (Fig. 5c–d)
- Tutorial 9: [Inter-chromosomal network](https://github.com/kmp0223/ARGEN/tree/main/tutorials/Tutorial9-Fig_6_7.ipynb) — the pooled, genome-wide ARGEN fit, module/GO characterization, and hub visualization (Fig. 6–7)
  
The processed data from Replogle et al. (2022, Cell), along with auxiliary data, are available via Zenodo at: [Zenodo link](https://doi.org/10.5281/zenodo.22052266). Once you download `data.zip` from Zenodo, put the extracted directory right under your `{ARGEN directory}`. The raw data is available at [the link provided by Replogle et al.](https://gwps.wi.mit.edu/).

## References
Replogle, J. M., Saunders, R. A., Pogson, A. N., Hussmann, J. A., Lenail, A., Guna, A., Mascibroda, L., Wagner, E. J., Adelman, K., Lithwick-Yanai, G., Iremadze, N., Oberstrass, F., Lipson, D., Bonnar, J. L., Jost, M., Norman, T. M., & Weissman, J. S. (2022). *Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq*. **Cell, 185**(14), 2559–2575.e28. https://doi.org/10.1016/j.cell.2022.05.013

Brown, B. C., Tokolyi, A., Morris, J. A., Lappalainen, T., & Knowles, D. A. (2025). *Large-scale causal discovery using interventional data sheds light on gene network structure in K562 cells*. **Nature Communications, 16**, 9628. https://doi.org/10.1038/s41467-025-64353-7
