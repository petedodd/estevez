#!/bin/bash
# for multicore machine NOTE edit to batch by < number of available cores
R --slave --vanilla --args < ../R/multirunnerInc.R ETH  & R --slave --vanilla --args < ../R/multirunnerInc.R KEN  & R --slave --vanilla --args < ../R/multirunnerInc.R LSO  & R --slave --vanilla --args < ../R/multirunnerInc.R MOZ  & R --slave --vanilla --args < ../R/multirunnerInc.R MWI  & R --slave --vanilla --args < ../R/multirunnerInc.R NGA & R --slave --vanilla --args < ../R/multirunnerInc.R SWZ  & R --slave --vanilla --args < ../R/multirunnerInc.R TZA  & R --slave --vanilla --args < ../R/multirunnerInc.R UGA  & R --slave --vanilla --args < ../R/multirunnerInc.R ZAF  & R --slave --vanilla --args < ../R/multirunnerInc.R ZMB  & R --slave --vanilla --args < ../R/multirunnerInc.R ZWE
