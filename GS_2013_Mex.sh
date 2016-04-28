#!/bin/sh

R --vanilla --slave --args pheno_mex2013_ver0.3_g 1 < GS_2013_Mex.R &
R --vanilla --slave --args pheno_mex2013_ver0.3_g 2 < GS_2013_Mex.R &
R --vanilla --slave --args pheno_mex2013_ver0.3_g 3 < GS_2013_Mex.R &
R --vanilla --slave --args pheno_mex2013_ver0.3_g 4 < GS_2013_Mex.R &
R --vanilla --slave --args pheno_mex2013_ver0.3_g 5 < GS_2013_Mex.R &
