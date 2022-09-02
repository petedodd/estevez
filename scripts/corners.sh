#!/bin/bash
# NOTE these plots take a surprisingly long time
# for multicore machine NOTE edit to batch by < number of available cores
R --slave --vanilla --args < ../R/cornerplots.R ETH  & R --slave --vanilla --args < ../R/cornerplots.R KEN  & R --slave --vanilla --args < ../R/cornerplots.R LSO  & R --slave --vanilla --args < ../R/cornerplots.R MOZ  & R --slave --vanilla --args < ../R/cornerplots.R MWI  & R --slave --vanilla --args < ../R/cornerplots.R NGA & R --slave --vanilla --args < ../R/cornerplots.R SWZ  & R --slave --vanilla --args < ../R/cornerplots.R TZA  & R --slave --vanilla --args < ../R/cornerplots.R UGA  & R --slave --vanilla --args < ../R/cornerplots.R ZAF  & R --slave --vanilla --args < ../R/cornerplots.R ZMB  & R --slave --vanilla --args < ../R/cornerplots.R ZWE

cd ../plots/corner/
for p in *.pdf
do
   pdftoppm "$p" "PNG/$p" -png
done
