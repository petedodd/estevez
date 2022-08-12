# estevez



## License

This work is licensed under a Creative Commons Attribution 4.0 International License

http://creativecommons.org/licenses/by/4.0/

![http://creativecommons.org/licenses/by/4.0/](https://i.creativecommons.org/l/by/4.0/88x31.png)

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
    - uses data/irrAOdata.Rdata from above
    - depends {R/modelprep.R 
      {R/plotters.R,R/getCountryAimParmsInc.R,R/AIMdynInc.R},
      R/dataLL.R}
3. Computation of results


TODO

- check whont
- SWZ bug see out tail
- see multirun end bug 
