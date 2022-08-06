#!/bin/bash
# for multicore machine
R --slave --vanilla --args < MRmakerInc.R ETH  & R --slave --vanilla --args < MRmakerInc.R KEN  & R --slave --vanilla --args < MRmakerInc.R LSO  & R --slave --vanilla --args < MRmakerInc.R MOZ  & R --slave --vanilla --args < MRmakerInc.R MWI  & R --slave --vanilla --args < MRmakerInc.R NGA

R --slave --vanilla --args < MRmakerInc.R SWZ  & R --slave --vanilla --args < MRmakerInc.R TZA  & R --slave --vanilla --args < MRmakerInc.R UGA  & R --slave --vanilla --args < MRmakerInc.R ZAF  & R --slave --vanilla --args < MRmakerInc.R ZMB  & R --slave --vanilla --args < MRmakerInc.R ZWE
