# an example of overlaying a lat lng raster onto a map of the world.

# rasterize locations, retrieve stats for bounding boxes, plot graphs
if (!'package:rgdal' %in% search()) {
  # at time of writing, rgdal no binary was available for mavericks / R v 3.1.2 so it has to be compiled from source
  # to install rgdal, you also have to install "gdal" (e.g. using http://brew.sh/ and [brew install gdal])
  install.packages('rgdal', repos='http://cran.us.r-project.org', type='source')
	install.packages(c('mapproj','maps','classInt', 'RColorBrewer','rworldmap'), repos='http://cran.us.r-project.org')
}

library('mapproj')
library('maps')
library('rgdal')
library('classInt')
library('RColorBrewer')
library('rworldmap')

plotMap <- function(dataFilename, divmap) {
  cat(paste(dataFilename, 'loading...'))
  diet.matrix <- read.csv(dataFilename)
  cat(paste(dataFilename, 'loaded.'))
  divmap <- mergeNicheWidthDiversityMap(diet.matrix, divmap)
  plotNicheMap(diet.matrix, divmap)
}
# transform to Mollweide equal area world projection
equiAreaProj4 <- function() {
  CRS("+proj=moll +lon_0=0")
  squareProj4()
}

squareProj4 <- function() {
  CRS("+proj=longlat +datum=WGS84 +lon_0=0")
}

plotWorldTisseuil <- function() {
  map(interior=FALSE, boundary=FALSE, add=FALSE, fill=TRUE, border='grey', col='grey')
}

plotWorldTittensor <- function() {
  map(interior=FALSE, boundary=FALSE, add=FALSE, fill=TRUE, border='white', col='white')
}

# calculated diet matrix (including niche.width) and diversity map (e.g. Tittensor)
plotNicheMap <- function(divmap, param.name = 'niche.width', main.title = 'niche width', add = TRUE) {
  vals <- divmap[[param.name]]

  summary(vals)
    ( cuts <- classIntervals(vals, style="fixed", dataPrecision=2,
                            fixedBreaks=c(summary(na.omit(vals)))) )
  plotcolors <- rev(brewer.pal(5, "Spectral"))
  divmap.colors <- findColours(cuts, plotcolors)
  
  #divmap.mollweide <- spTransform(divmap, equiAreaProj4())
  plot(divmap, col=divmap.colors, border=NA, add=add)
  title(main=main.title)
  legend("bottomleft", horiz=FALSE, legend=names(attr(divmap.colors, "table")), fill=attr(divmap.colors, "palette"), bty="n")
}

mergeNicheWidthDiversityMapMarine <- function(diet.matrix, divmap, aggfun = mean) {
  # merge niche and div map after aggregating by gridcode using aggregator function aggfun
  diet.matrix.meanPerGridCode <- aggregate(niche.width ~ GRIDCODE, data = diet.matrix, aggfun)
  divmap <- merge(divmap, diet.matrix.meanPerGridCode, by='GRIDCODE')
}

mergeNicheWidthDiversityMapFresh <- function(diet.matrix, divmap, aggfun = mean) {
  # merge niche and div map after aggregating by gridcode using aggregator function aggfun
  divmap <- merge(divmap, diet.matrix.meanPerGridCode, by='ID_JF')
}

# example using a diet matrix data file
samplePlot <- function(divmap = ReadBiodiversityShapes()$diversity.map) { 
  plotMap('data/hayden2014_Actinopterygii_GBIF_tittensor2010_2014-11-05-153153.csv', divmap)
}

setup_device <- function(filename, format = 'png') {
  if ('pdf' == format) {
    pdf(file=paste(filename, '.pdf', sep=''))  
  } else {
    png(file=paste(filename, '.png', sep=''), width=640, height=320)  
  }
}

graph_formats <- c('png', 'pdf')

plotFreshwaterMaps <- function(biodiversityShapes = ReadBiodiversityShapes(GetTisseuilShapes()), format='pdf') {
  fresh.basins <- biodiversityShapes
  fresh.niche.data <- read.csv('scripts/fresh.data.csv')
  fresh.trophic.data <- read.csv('scripts/fresh.TL.data.csv')
  fresh.basins <- merge(fresh.basins$diversity.map, fresh.niche.data, by='ID_JF')
  fresh.basins <- merge(fresh.basins, fresh.trophic.data, by='ID_JF')

  # print plots
  lapply(graph_formats, function(format) {
    setup_device(filename = 'scripts/fresh.niche.median', format=format)  
    map(interior=FALSE, boundary=FALSE, add=FALSE, fill=TRUE, border='grey', col='grey')
    plotNicheMap(fresh.basins, param.name='niche.med', main.title='Median DNW')
    dev.off()
  })
  
  lapply(graph_formats, function(format) {
    setup_device(filename='scripts/fresh.niche.mean', format=format)  
    map(interior=FALSE, boundary=FALSE, add=FALSE, fill=TRUE, border='grey', col='grey')
    plotNicheMap(fresh.basins, param.name='niche.mean', main.title='Mean DNW')
    dev.off()
  })

  lapply(graph_formats, function(format) {
    setup_device(filename='scripts/fresh.trophic.mean', format=format)  
    map(interior=FALSE, boundary=FALSE, add=FALSE, fill=TRUE, border='grey', col='grey')
    plotNicheMap(fresh.basins, param.name='TL.mean', main.title='Mean TL')
    dev.off()
  })
}

plotMarineMaps <- function(biodiversityShapes = ReadBiodiversityShapes()) {
  marine.squares <- biodiversityShapes
  marine.niche.data <- read.csv('scripts/marine.sumry.csv')
  marine.trophic.data <- read.csv('scripts/marine.TL.sumry.csv')
  marine.squares <- merge(marine.squares$diversity.map, marine.niche.data, by='GRIDCODE')
  marine.squares <- merge(marine.squares, marine.trophic.data, by='GRIDCODE')

  # print plots

  lapply(graph_formats, function(format) {
    setup_device(file='scripts/marine.niche.median', format= format)  
    plotNicheMap(marine.squares, param.name='niche.med', main.title='Median DNW', add = FALSE)
    map(interior=FALSE, boundary=FALSE, add=TRUE, fill=TRUE, border='gray', col='gray')
    dev.off()
  });
  
  lapply(graph_formats, function(format) {
    setup_device(filename='scripts/marine.niche.mean', format = format)  
    plotNicheMap(marine.squares, param.name='niche.mean', main.title='Mean DNW', add = FALSE)
    map(interior=FALSE, boundary=FALSE, add=TRUE, fill=TRUE, border='gray', col='gray')
    dev.off()
  })

  lapply(graph_formats, function(format) {
    setup_device(filename='scripts/marine.trophic.mean', format = format)  
    plotNicheMap(marine.squares, param.name='TL.mean', main.title='Mean TL', add = FALSE)
    map(interior=FALSE, boundary=FALSE, add=TRUE, fill=TRUE, border='gray', col='gray')
    dev.off()
  })
}



