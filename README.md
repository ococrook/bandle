[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)

# The `bandle` package


The [**bandle**](https://ococrook.github.io/bandle/) (Bayesian Analysis of Differential Localisation Experiments) package is R/Bioconductor package for analysis differential localisation experiments, include storage, computation, statistics and visulisations. 


# Installation requirements

Users will require a working version of R, currently at least version >4. It is recommended to use [RStudio](https://www.rstudio.com). The package can then be installed using the Bioconductor or `devtools` package. The package should take a few minutes to install on a regular desktop or laptop. The package will need to be loaded using `library(bandle)`

To install the stable Bioconductor release (recommended):

```{r,}
## unless BiocManager is already installed
install.packages("BiocManager")
## then
BiocManager::install("bandle")
```


Unstable/development package version install:

```{r,}
devtools::install_github("ococrook/bandle")
```

We do not advise you install the unstable version unless you know what you are doing as not all pre-release features may be tested or documented.


# Basic ideas and concepts

The [`bandle` package](https://bioconductor.org/packages/bandle) implements the [BANDLE method](https://www.biorxiv.org/content/10.1101/2021.01.04.425239v3) for the analysis of comparative/dynamic mass spectrometry based proteomics experiments. Data from form such experiments most commonly yield a matrix of measurements where we have proteins/peptides/peptide spectrum matches (PSMs) along the rows, and samples/fractions along the columns. In comparative/dynamic experiments where we expect re-localisation upon some stimulus to sub-cellular environment, the data analysis is more challenging The BANDLE method takes two (replicated) datasets as input and uses these data to compute the probability that a protein differentially localises upon cellular perturbation, as well quantifying the uncertainty in these estimates.

To use bandle the data must be stored as a *list* of *MSnSet* instances, as implemented in the Bioconductor [`MSnbase`](https://bioconductor.org/packages/release/bioc/html/MSnbase.html) package. Please see the relevant vignettes in MSnbase for constructing these data containers.

# Vignettes

There are currently two vignettes that accompany the package. The first vignette [v01-getting-started](https://ococrook.github.io/bandle/articles/v01-getting-started.html) provides an introduction to the `bandle` package and follows a short theortical example of how to perform differential localisation analysis of quantitative proteomics data using the BANDLE model (Crook et al. 2022 doi: https://doi.org/10.1101/2021.01.04.425239). Explanation and general recommendations of the input parameters are provided in this vignette.

The second vignette [v02-workflow](https://ococrook.github.io/bandle/articles/v02-workflow.html) is a more comprehensive workflow which follows a real-life use case applying the BANDLE methods and workflow to the analysis of data generated from spatial proteomics of a human THP-1 monocyte system.

These vignettes can be found through the Bioconductor landing page for [`bandle`](https://bioconductor.org/packages/bandle), here in the repo and also the *Articles* tab of the accompanying [web page](https://ococrook.github.io/bandle).

Notes on run time: A small dataset can take around an hour to run; for large dataset we recommend a a compute server. The longest the analysis has taken has been a couple of hours on a single compute node. The demo take a few minutes to run.

# Documentation

Documentation to run the main functions can be found in the vignette or by typing `?bandle` in the console after loading the package. 

We recommend reading our other workflow manuscripts:

Basic processing and machine learning:

https://f1000research.com/articles/5-2926

Bayesian analysis:

https://f1000research.com/articles/8-446

The BANDLE manusript is currently on biorxiv:

https://www.biorxiv.org/content/10.1101/2021.01.04.425239v3


For manuscripts that apply bandle, see:

https://www.biorxiv.org/content/10.1101/2022.01.24.477541v1
https://linkinghub.elsevier.com/retrieve/pii/S1535-9476(22)00002-0

# Contribution

Contributions are welcome, please open issue so we can discuss any contribution in advance.

# Feature requests

This package is actively being developed and maintained, please open Github issue if you would like to request or discuss a particular feature.
