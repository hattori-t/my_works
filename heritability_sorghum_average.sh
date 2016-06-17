#!/bin/sh

R --vanilla --slave --args Mexico2013_average juice < heritability_sorghum.R &
R --vanilla --slave --args Mexico2013_average brix < heritability_sorghum.R &
R --vanilla --slave --args Mexico2013_average leaf.culm.weight < heritability_sorghum.R &
R --vanilla --slave --args Mexico2013_average log.leaf.culm.weight < heritability_sorghum.R &
R --vanilla --slave --args Mexico2013_average plant.height < heritability_sorghum.R &
R --vanilla --slave --args Mexico2013_average panicle.length < heritability_sorghum.R &
R --vanilla --slave --args Mexico2013_average culm.length < heritability_sorghum.R &
R --vanilla --slave --args Mexico2013_average culm.number < heritability_sorghum.R &
R --vanilla --slave --args Mexico2013_average lodging < heritability_sorghum.R &
wait

R --vanilla --slave --args Mexico2014_average juice < heritability_sorghum.R &
R --vanilla --slave --args Mexico2014_average brix < heritability_sorghum.R &
R --vanilla --slave --args Mexico2014_average total.weight < heritability_sorghum.R &
R --vanilla --slave --args Mexico2014_average log.total.weight < heritability_sorghum.R &
R --vanilla --slave --args Mexico2014_average plant.height < heritability_sorghum.R &
R --vanilla --slave --args Mexico2014_average panicle.length < heritability_sorghum.R &
wait
R --vanilla --slave --args Mexico2014_average culm.length < heritability_sorghum.R &
R --vanilla --slave --args Mexico2014_average culm.number < heritability_sorghum.R &
R --vanilla --slave --args Mexico2014_average culm.diameter.1 < heritability_sorghum.R &
R --vanilla --slave --args Mexico2014_average culm.diameter.2 < heritability_sorghum.R &
R --vanilla --slave --args Mexico2014_average culm.diameter.mean < heritability_sorghum.R &
R --vanilla --slave --args Mexico2014_average culm.area < heritability_sorghum.R &
R --vanilla --slave --args Mexico2014_average culm.volume < heritability_sorghum.R &
R --vanilla --slave --args Mexico2014_average lodging < heritability_sorghum.R &
wait

R --vanilla --slave --args Mexico2015_average juice < heritability_sorghum.R &
R --vanilla --slave --args Mexico2015_average brix < heritability_sorghum.R &
R --vanilla --slave --args Mexico2015_average total.weight < heritability_sorghum.R &
R --vanilla --slave --args Mexico2015_average log.total.weight < heritability_sorghum.R &
R --vanilla --slave --args Mexico2015_average plant.height < heritability_sorghum.R &
R --vanilla --slave --args Mexico2015_average panicle.length < heritability_sorghum.R &
R --vanilla --slave --args Mexico2015_average culm.length < heritability_sorghum.R &
wait
R --vanilla --slave --args Mexico2015_average culm.number < heritability_sorghum.R &
R --vanilla --slave --args Mexico2015_average culm.diameter.1 < heritability_sorghum.R &
R --vanilla --slave --args Mexico2015_average culm.diameter.2 < heritability_sorghum.R &
R --vanilla --slave --args Mexico2015_average culm.diameter.mean < heritability_sorghum.R &
R --vanilla --slave --args Mexico2015_average culm.area < heritability_sorghum.R &
R --vanilla --slave --args Mexico2015_average culm.volume < heritability_sorghum.R &
R --vanilla --slave --args Mexico2015_average insects < heritability_sorghum.R &
R --vanilla --slave --args Mexico2015_average chemical < heritability_sorghum.R &
