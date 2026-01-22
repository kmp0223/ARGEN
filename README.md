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

#### 2. Install/load required R packages

In R, run those codes that download the required packages for running ARGEN.

```
install.packages('MASS')
```

Details about implementing codes and reproducing Figures in the manuscript can be found in the [Tutorial](https://github.com/kmp0223/ARGEN/tree/main/code/tutorials) of this github. There are mainly two tutorials provided: 

https://github.com/kmp0223/ARGEN/tree/main/code/tutorials/Tutorial.ipynb
- Tutorial 1: [K562 cell line application](https://github.com/kmp0223/ARGEN/tree/main/code/tutorials/Tutorial.ipynb)
- Tutorial 2: [Simulation](https://github.com/kmp0223/ARGEN/tree/main/code/tutorials/Tutorial.ipynb)
  
The required datasets can be found in the [Google Drive link](link).

