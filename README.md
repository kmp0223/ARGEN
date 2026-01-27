# ARGEN : Causal Network Recovery in Perturb-seq Experiments Using Proxy and Instrumental Variables
![ARGEN diagram](/figures/intro.png)

## ARGEN Usage

### Preparation


First of all, before cloning the ARGEN github package, go to the right directory that you would like to implement ARGEN. In the cmd terminal, do

```
cd {ARGEN directory}
```

then go on to the next step. {ARGEN directory} could be for example /kwangmoon/ARGEN.


#### 1. Repository clone

For cloning the github repository, again on the cmd terminal, run the linux code 

```
git clone git@github.com:kmp0223/ARGEN.git
```

For R, we need the requirements as below : 


-   R: [R installation](https://www.r-project.org)  (>=4.2.2)

#### 2. Install and load required R packages

In R, run the following commands to install the packages required for running **ARGEN** (skip any packages that are already installed):

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages(c(
  "dplyr", "data.table", "Matrix", "glmnet", "tidyr",
  "qs", "sandwich", "lmtest", "MASS",
  "igraph", "ggplot2", "future",
  "purrr", "rlang", "scales", "visNetwork",
  "ggpubr", "tibble", "ggraph", "patchwork"
))

BiocManager::install(c("GenomicFeatures", "onlineFDR"))


```

Details about implementing codes and reproducing Figures in the manuscript can be found in the [Tutorial](https://github.com/kmp0223/ARGEN/tree/main/code/tutorials) of this github. There are mainly three tutorials provided: 


- Tutorial 1: [K562 cell line application](https://github.com/kmp0223/ARGEN/tree/main/code/tutorials/Tutorial.ipynb)
- Tutorial 2: [K562 cell line application across chromosomes](https://github.com/kmp0223/ARGEN/tree/main/code/tutorials/Tutorial2.ipynb)
- Tutorial 3: [Simulation](https://github.com/kmp0223/ARGEN/tree/main/code/tutorials/Tutorial3.ipynb)

  
The processed data from Replogle et al. (2022, Cell), along with auxiliary data, are available via Google Drive at: [Google Drive link](https://drive.google.com/drive/u/0/folders/1luFVh6laubmLKERc0wWDYh3k8_gCCVv0). The raw data is available at [the link provided by Replogle et al.](https://gwps.wi.mit.edu/).
## Reference

Replogle, J. M., Saunders, R. A., Pogson, A. N., Hussmann, J. A., Lenail, A., Guna, A., Mascibroda, L., Wagner, E. J., Adelman, K., Lithwick-Yanai, G., Iremadze, N., Oberstrass, F., Lipson, D., Bonnar, J. L., Jost, M., Norman, T. M., & Weissman, J. S. (2022). *Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq*. **Cell, 185**(14), 2559–2575.e28. https://doi.org/10.1016/j.cell.2022.05.013

