#!/bin/bash

make
./cholera > output
Rscript plot.R
