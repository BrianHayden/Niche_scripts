# Niche_scripts
R scripts analysing and assessing geographic variation in dietary niche width of niches.

The purpose of this repository is to enable the reproduction of the figures in the [/figs](./figs) directory using the scripts in [/R](./R) and the data in [/data](./data). 

Roughly, the recipe for reproducing the figures is:

 1. get species interaction data from [GloBI](http://globalbioticinteractions.org) using the [rglobi](https://github.com/ropensci/rglobi) package. (see [R/GeoDietNiche.R](./R/GeoDietNiche.R))
 1. get species occurrence data from [GBIF](http://gbif.org) using the [rgbif](https://github.com/ropensci/rgbif) package. (see [R/GeoDietNiche.R](./R/GeoDietNiche.R))
 1. integrate fresh water and marine biodiversity and environmental parameters (see [R/Freshwater.R[('./R/Freshwater') and [R/Marine.R](./R/Marine.R)
 1. calculate fresh water and marine dietary niche models (see [R/Freshwater.R](./R/Freshwater.R) and [R/Marine.R](./R/Marine.R))
 1. overlay results on a world map and plot figures (see [R/GeoDietNiche.R](.R/GeoDietNiche.R))

Please open an [issue](https://github.com/BrianHayden/issues/new) if need help reproducing the result.
