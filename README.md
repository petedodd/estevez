# estevez



## License

This work is licensed under a Creative Commons Attribution 4.0 International License

http://creativecommons.org/licenses/by/4.0/

![http://creativecommons.org/licenses/by/4.0/](https://i.creativecommons.org/l/by/4.0/88x31.png)


## Download

The repository is relatively *large* due to some of the output data.
It might therefore be preferable to download using

```
git clone --depth 1 git@github.com:petedodd/estevez.git
```

which will only keep the most recent file versions, and should allow everything to run. Alternatively, GitHub allows download as a zip file, but this will be substantially larger since it includes previous versions of output data.


## Directory structure

```
.
├── data
├── indata
├── plots
├── R
└── scripts
```

- *scripts* : shell scripts (bash) for running analyses on PCs/HPCs
- *R* : analysis R code
- *indata* : analysis input data required, including fit data
- *data* :  data generated during the analyses
- *plots* : graphical output from analyses



### Analysis order ###

1. IRR analysis
    - scripts/IRRanalyses.sh runs -> R/IRRanalysis.R
    - depends {R/AIMdynInc.R, R/getCountryParmsInc.R}
    - creates data/irrAOdata.Rdata
2. Inference
    - scripts/jobarray.sh runs -> R/inference.R
    - uses data/irrAOdata.Rdata from above
    - depends {R/modelprep.R, R/plotters.R, R/getCountryAimParmsInc.R, R/AIMdynInc.R, R/dataLL.R}
3. Computation of results
    - scripts/MRallPar.sh runs -> R/multirunner.R
    - uses data/sampe*.Rdata from above
    - depends {R/modelprep.R, R/plotters.R, R/getCountryAimParmsInc.R, R/AIMdynInc.R, R/dataLL.R}
4. Analysis of results
    - R/resultplots.R
    - uses {data/MR*.Rdata, data/OD*.Rdata} from above
    - depends {R/modelprep.R, R/plotters.R, R/dataplots.R, R/dataLL.R}
    - scripts/corners.sh -> R/cornerplots.R
    - uses data/sampe*.Rdata from above


### Required packages ###

The following R packages are required:

here, data.table, ggplot2, lhs, glue, ggthemes, scales, ggpubr, ggrepel, grid, patchwork, odin, GGally, metafor, reticulate

which may themselves have dependencies. In addition, the python package zeus is required.

