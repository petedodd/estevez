# estevez



## License

This work is licensed under a Creative Commons Attribution 4.0 International License

http://creativecommons.org/licenses/by/4.0/

![http://creativecommons.org/licenses/by/4.0/](https://i.creativecommons.org/l/by/4.0/88x31.png)


## Download

The repository is relatively large due to some of the output data.
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
    - uses data/sampe*R.data from above
    - depends {TODO}
4. Analysis of results
    - R/resultplots.R
    - uses {data/MR*.Rdata, data/OD*.Rdata} from above
    - depends {TODO}

TODO list of libraries needed

TODO

- check whont
- ARIS & ARIA missing from MR output
- HIV effects check
- TODO make sure any missing directories are created by scripts
- continue plotting/result work
