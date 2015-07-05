# Niche_scripts
R scripts analysing and assessing geographic variation in dietary niche width of niches.

The purpose of this repository is to enable the reproduction of the figures in the ```/figs``` directory using the scripts in ```/R``` and the data in ```/data```. 

Please open an [issue](https://github.com/BrianHayden/issues/new) if need help reproducing the result.

Roughly, the recipe for reproducing the figure are:

1. get species interaction data from [GloBI](http://globalbioticinteractions.org) using the [rglobi](https://github.com/ropensci/rglobi) package. (see R/CompileGeoDiet.R)
2. get species occurrence data from [GBIF](http://gbif.org) using the [rgbif](https://github.com/ropensci/rgbif) package. (see [R/CompileGeoDiet.R](./blob/master/R/CompileGeoDiet.R))
3. integrate fresh water and marine biodiversity and environmental parameters (see R/Freshwater.R and R/Marine.R)
4. calculate fresh water and marine dietary niche models (see R/Freshwater.R and R/Marine.R)
5. overlay results on a world map and plot figures (see R/CompileGeoDiet.R)


