#!/bin/sh

R --vanilla --slave --args pheno_mex2013_ver0.3_g lodging < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex2013_ver0.3_g culm.num < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex2013_ver0.3_g panicle.length < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex2013_ver0.3_g plant.height < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex2013_ver0.3_g culm.length < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex2013_ver0.3_g leaf.culm.weight < heritability_sorghum.R &
wait

R --vanilla --slave --args pheno_mex2013_ver0.3_g juicy < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex2013_ver0.3_g brix < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex2013_ver0.3_g log.leaf.culm.weight < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2014_inbred_ABEF lodging < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2014_inbred_ABEF culmnum < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2014_inbred_ABEF juice < heritability_sorghum.R &
wait

R --vanilla --slave --args pheno_mex_2014_inbred_ABEF brix < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2014_inbred_ABEF weight < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2014_inbred_ABEF panicle.length < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2014_inbred_ABEF plant.height < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2014_inbred_ABEF culm.length < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2014_inbred_ABEF culm.diameter.mean < heritability_sorghum.R &
wait

R --vanilla --slave --args pheno_mex_2014_inbred_ABEF log.weight < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2015_A-K culmnum < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2015_A-K insects < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2015_A-K chemical < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2015_A-K plant.height < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2015_A-K panicle.length < heritability_sorghum.R &
wait

R --vanilla --slave --args pheno_mex_2015_A-K culm.length < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2015_A-K weight < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2015_A-K juice < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2015_A-K brix < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2015_A-K culm.diameter.mean < heritability_sorghum.R &
R --vanilla --slave --args pheno_mex_2015_A-K log.weight < heritability_sorghum.R &
