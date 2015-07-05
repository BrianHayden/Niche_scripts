# Compiles geospatial diet matrices for species and associates environmental (e.g. sea surface temperature) and ecological (e.g. biodiversity). 

# Uses GloBI (globalbioticinteractions.org) and GBIF (gbif.org) for species interactions and occurrence data. 

# To compile a diet for ray-finned fishes (Actinopterygii)
#   CompileGeoDiet()
# 
# For any other taxon (e.g. Elasmobranchii)
#   CompileGeoDiet('Elasmobranchii')

# Please report issues running this script by opening a new issues at https://github.com/BrianHayden/Niche_scripts

install.packages(c('rglobi','rgbif', 'rfishbase', 'maptools','data.table', 'RCurl'), repos='http://cran.us.r-project.org')


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

# Generate Diet Matrix using rglobi (see https://github.com/ropensci/rglobi)
#
# To install: 
# (a) install devtools using install.packages('devtools') 
# (b) install rglobi using devtools::install_github('rglobi', 'ropensci')
# (c) load the contents of this gist into R using devtools::source_url('https://gist.githubusercontent.com/jhpoelen/81b4eced3963633b96cb/raw/dietMatrix.R')
#
# Examples
#
#  Create a diet matrix for humans (Homo sapiens) and old world rat (Rattus rattus) using prey categories Aves and Mammalia:
#  CreateDietMatrix(c('Homo sapiens', 'Rattus rattus'), c('Aves', 'Mammalia'))
#
#  Example 1. Create a diet matrix of first 25 predator species taxa using Haydens Diet Category list (see defaults in GetPreyCategories()):
#
#  diet.matrix <- CreateDietMatrix(GetPredatorTaxa(limit=25))
#
#
#  Example 2. Same as example 1. except that port [1234] is used to query. This can be used to run queries again something
#    other (e.g. Dark GloBI) than the default, open, GloBI instance. Note that port [1234] will have to be replaced with a 
#    valid port value in order for it to work.
# 
#  opts <- list(port=1234)
#  predator.taxa <- GetPredatorTaxa(limit=25, opts=opts)
#  diet.matrix <- CreateDietMatrix(predator.taxa, opts = opts)
#
# Result is a data frame with predator taxon name and food categories as columns and predator taxon occurrences as rows.
#

predator.match.clause <- "preyTaxon<-[:CLASSIFIED_AS]-prey<-[:ATE|PREYS_UPON]-predator-[:CLASSIFIED_AS]->predatorTaxon"

ReportProgress <- function() {
  message('.', appendLF=FALSE)
}

# Limited to Sharks (Elasmobranchii) for now.
GetPredatorTaxa <- function(rank = "Species", predator.taxa = list('Elasmobranchii'), skip = 0, limit = 25, opts = list(port = 7474)) {
  luceneQuery <- paste('path:', predator.taxa, ' ', sep='', collapse='')
  cypher <- paste("START predatorTaxon = node:taxonPaths('", luceneQuery , "') MATCH ", predator.match.clause, " WHERE has(predatorTaxon.rank) AND predatorTaxon.rank = '", rank, "' RETURN distinct(predatorTaxon.name) as `predator.taxon.name`", sep="")
  if (!is.null(limit) && !is.null(skip)) {
    cypher <- paste(cypher, 'skip', skip, 'limit', limit)
  }
  rglobi::query(cypher, opts = opts)$predator.taxon.name
}

EscapeRegex <- function(string) {
  gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", string)
}

# OtherCat
#
# Creates a regex expression that helps to identify taxa with a parent that is not part of the provided child taxa list
#
# Example:
#   >expr <- OtherCat('Mollusca', c('Bivalvia','Gastropoda'))
#   >grep(expr, 'Mollusca | somethingelsa', perl=TRUE) == 1 
#   TRUE
#   >grep(expr, 'Mollusca | Bivalvia', perl=TRUE) == 1
#   FALSE
#
OtherCat <- function(parent.taxon.name, excluded.child.taxon.names) {
  excl.terms <- paste('(', excluded.child.taxon.names, ')', sep='', collapse='|')
  paste('(',parent.taxon.name,') \\| (?!',excl.terms,')', sep='')
}

# Creates individual categories to cover all child taxa in the form (child1, child2, ..., other)
#
#
CreateCats <- function(parent.taxon.name, child.taxon.names) {
  c(child.taxon.names, OtherCat(parent.taxon.name, child.taxon.names))
}

GetPreyCategories <- function() {
  c(
     CreateCats('Cnidaria',c('Anthozoa', 'Medusozoa')),
	 CreateCats('Annelida',c('Polychaeta', 'Clitellata', 'Haplodrili')),
	 CreateCats('Mollusca',c('Polyplacophora', 'Bivalvia', 'Gastropoda', 'Cephalopoda')),
	 CreateCats('Arthropoda',c('Ostracoda', 'Maxiliopoda', 'Malacostraca', 'Crustacea', 'Insecta', 'Brachiopoda')),
	 CreateCats('Echinodermata', c('Ophiuroidea', 'Echinoidea', 'Holothuroidea')),
	 CreateCats('Chordata', c('Actinopterygii','Chondrichthyes', 'Amphibia','Reptilia', 'Aves','Mammalia')),
	 OtherCat('Animalia', c('Cnidaria','Annelida','Mollusca','Arthropoda','Echinodermata','Chordata')))
}

# Retrieves diet items of given predator and classifies them by matching the prey categories against 
# both taxon hierarchy of prey and the name that was originally used to describe the prey.
UniquePreyTaxaOfPredator <- function(predator.taxon.name, prey.categories, opts = list(port = 7474)) {
  cypher <- paste("START predatorTaxon = node:taxons(name='", predator.taxon.name, "') MATCH ", predator.match.clause, ", prey-[:ORIGINALLY_DESCRIBED_AS]->preyTaxonOrig WHERE has(preyTaxon.path) RETURN distinct(preyTaxonOrig.name) as `prey.taxon.name.orig`, preyTaxon.path as `prey.taxon.path`", sep="")
  result <- rglobi::query(cypher, opts = opts)
  ReportProgress()
  all.taxa.paths <- Reduce(function(accum, path) paste(accum, path), paste('{',result$prey.taxon.path,' | }', sep=''))
  has.prey.category <- lapply(prey.categories, function(prey.category) {
    match <- grep(prey.category, all.taxa.paths, perl=TRUE)
    ifelse(length(match) > 0, 1, 0)
  })
  df <- data.frame(t(c(predator.taxon.name, has.prey.category)))
  names(df) <- c('predator.taxon.name', prey.categories)
  df
}

CreatePredatorPreyCategories <- function(predator.taxon.name, prey.categories = GetPreyCategories(), opts = list()) {
  UniquePreyTaxaOfPredator(predator.taxon.name, prey.categories, opts)
}

CreateDietMatrix <- function(predator.taxon.names = NULL, prey.categories = GetPreyCategories(), opts = list()) {
  message(paste('building diet matrix for [', length(predator.taxon.names), '] predators', sep=""), appendLF=FALSE) 
  if (is.null(predator.taxon.names)) {
    predator.taxon.names <- GetPredatorTaxa(opts = opts)
  }
  diet.matrix <- Reduce(function(accum, predator.taxon.name) rbind(accum, CreatePredatorPreyCategories(predator.taxon.name, opts = opts)), predator.taxon.names, init=data.frame())
  diet.matrix$niche.width <- apply(diet.matrix[2:length(diet.matrix)], 1, function(row) sum(unlist(row)))
  # normalize using number of prey categories, assuming number of prey categories > 0
  diet.matrix$niche.width <- diet.matrix$niche.width / length(prey.categories)
  message(' done.')
  diet.matrix
}
# GetBiodiversityAtPoint
#
# Retrieve biodiversity information for latitude longitude. 
#
#
# Example:
#  # get a point in SpatialPoints class
#` points <- CreateSpatialPoint(37.7584754, -122.235328)
#  shapes <- ReadBiodiversityShapes()
#  GetBiodiversityAtPoints(points, shapes)
#
# Output:
#  GRIDCODE X_COORD  Y_COORD  AllTaxa   AllNorm   CoastNorm  OceanNorm
#1 98       -120     40.03754 404.984   0.3326864 0.1891753  0.500116


library(maptools)
library(rgbif)
library(data.table)
require(RCurl)

DownloadFile <- function(url, file) {
  f = CFILE(file, mode='wb')
  a = curlPerform(url = url, writedata = f@ref, noprogress=FALSE)
  close(f)
  a 
}

GetTittensorProjection <- function() {
  CRS("+proj=longlat +datum=WGS84")
}

GetTisseuilProjection <- function() {
  CRS("+proj=longlat +datum=WGS84 +no_defs")
}

# see Tittensor et al. 2010 in nature 
# http://dx.doi.org/10.1038/nature09329
GetTittensorShapes <- function() {
  list(url='http://www.mathstat.dal.ca/~derekt/publications/global_patterns_and_predictors_of_marine_biodiversity_across_taxa.zip', file = 'Global_patterns_predictors_marine_biodiversity_across_taxa.shp', diversity.metric = 'AllNorm', name='tittensor2010', proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs")) 
}

# see 2013 Tisseuil et al. paper in Journal of Animal Ecology http://dx.doi.org/10.1111/1365-2656.12018
# Only Fish shapefiles stored in S3 repository for convenience and ability to reproduce results.
GetTisseuilShapes <- function() {
  shapefile.name <- '/commondata/data0/bv_clemANR.shp'
  shapefile.url <- 'https://globi.s3.amazonaws.com/datasets/org/eol/globi/geo/tisseuil2013/0.1/tisseuil2013-0.1-Fish_TotNativSpec.zip'
  list(url=shapefile.url, file=shapefile.name, diversity.metric='fish_endem', name='tisseuil2013', proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs")) 
}

# Retrieve and read biodiversity maps, uses tittensor by default
ReadBiodiversityShapes <- function(shapes = GetTittensorShapes()) {
  td <- tempdir()
  tf <- tempfile(tmpdir=td, fileext=".zip")
  DownloadFile(shapes$url, tf)
  unzip(tf, exdir=td)
  shapefile.tmppath <- file.path(td, shapes$file)
  div.map <- readShapeSpatial(shapefile.tmppath, shapes$proj4string)
  unlink(td)
  shapes$diversity.map <- div.map
  shapes
}

CreateSpatialPoint <- function(lat, lng) {
  SpatialPoints(data.frame(c(lng), c(lat)), proj4string = GetProjection())
}

# Retrieve biodiversity information at specified points
GetBiodiversityAtPoints <- function(points, shapes = ReadBiodiversityShapes()) {
   points %over% shapes$diversity.map
}

# retrieve occurrence count (see https://github.com/ropensci/rgbif/issues/94)
GetOccurrenceCount <- function(taxonKey) {
  resp <- httr::GET(paste('http://api.gbif.org/v1/occurrence/count?taxonKey=', taxonKey, '&isGeoreferenced=true', sep=''))
  as.numeric(httr::content(resp, as='text', encoding='UTF-8'))
}

SelectUniqueRegions <- function(occurrence.data, div.map) {
  zones <- NA
  if ('decimalLongitude' %in% names(occurrence.data)) {
	  points <- SpatialPoints(cbind(occurrence.data$decimalLongitude, occurrence.data$decimalLatitude), proj4string = GetProjection())
	  if (length(points) > 0) {
      zones <- unique(points %over% div.map)
		  names(zones) <- names(div.map)
	  }
  }
	zones
}

# Retrieve known locations for species with name
GetLocationsForSpeciesGBIF <- function(scientificName, opts = list(occurrence.limit = 20)) {
  key <- name_backbone(scientificName, rank='SPECIES', strict=TRUE)
  if (key$matchType == 'EXACT' && ('speciesKey' %in% names(key))) {
    result <- occ_search(taxonKey = key$speciesKey, hasCoordinate=TRUE, spatialIssues=FALSE, limit = opts$occurrence.limit)
    result$data
  } else {
    NA
  }
}

GetLocationsForSpeciesGLOBI <- function(scientificName, opts = list(occurrence.limit = 20, port = '7474')) {
  rglobi::query(paste("START taxon = node:taxons(name='", scientificName, "') MATCH loc<-[:COLLECTED_AT]-specimen-[:CLASSIFIED_AS]->taxon RETURN loc.latitude as `decimalLatitude`, loc.longitude as `decimalLongitude` LIMIT ", opts$occurrence.limit, sep=''), opts) 
}

GetLocationsForSpecies <- function(scientificName, div.map, opts = list(occurrence.source = 'GBIF', occurrence.limit = 20, port = '7474')) {
  factory <- list('GBIF' = GetLocationsForSpeciesGBIF, 'GLOBI' = GetLocationsForSpeciesGLOBI)
  zones <- NA
  source <- opts$occurrence.source
  if (source %in% names(factory)) {
    locations <- factory[[source]](scientificName, opts)
    zones <- SelectUniqueRegions(locations, div.map)
  } 
  zones
}

# Appends Tittensor square info for each predator item in the diet matrix. 
# 
# Uses occurrence lat/lng pairs from specified source. Supported Sources are 'GBIF' and 'GLOBI'. By default GBIF is selected.
#
# Example:
#	DiversityVsNicheWidth(diet.matrix, opts = list(occurrence.source = 'GLOBI', occurrence.limit = 1024, port = 7474))
#
# The example selects up to 1024 occurrence locations from GloBI on port 7474. Also, because the diversity map is not provided, a map is downloaded. 
#

DiversityVsNicheWidth <- function(diet.matrix, shapes = ReadBiodiversityShapes(), opts = list(occurrence.source = 'GBIF', occurrence.limit = 20, port = '7474')) {
  message('aggregating diet.matrix with provided diversity map', appendLF=FALSE)
  
  rows <- apply(diet.matrix, MARGIN=1, function(row) {
    ReportProgress()
	  zones <- GetLocationsForSpecies(row$predator.taxon.name, shapes$diversity.map, opts)
	  if (length(zones) > 1) {
	    row.new <- cbind(row, zones)
	    names(row.new) <- c(names(row), names(zones))
	    row.new
	  } else {
      NULL
	  }
  })
  message(' done.')
  # exclude rows with NAs
  na.omit(rbindlist(rows))
}

ExportToFile <- function(diet.matrix, file = paste('diet-matrix.csv', format(Sys.time(), '%Y-%m-%d-%H%M%S'), sep='-')) {
  out <- data.frame(lapply(diet.matrix, function(x) factor(unlist(x))))
  write.csv(out, file = file)
  message(paste('wrote file [', file, ']', sep=''))
}

CreateAndPlotNicheWidthModel <- function(diet.matrix, title='Elasmobranchii', diversity.metric = 'AllNorm') {
  niche.diversity.model <- lm(diet.matrix$niche.width~diet.matrix[[diversity.metric]])
  
  # plot data
  plot(diet.matrix[[diversity.metric]], diet.matrix$niche.width, main = title)

  # plot model
  abline(niche.diversity.model, col='red')

  # show confidence intervals of model
  niche.diversity.model.pred <- predict(niche.diversity.model, interval='confidence')
  lines(diet.matrix[complete.cases(diet.matrix)][[diversity.metric]], niche.diversity.model.pred[,2], lty=2, col='red')
  lines(diet.matrix[complete.cases(diet.matrix)][[diversity.metric]], niche.diversity.model.pred[,3], lty=2, col='red')

  niche.diversity.model
}


# Recipe to calculate diet niche and associated diversity for species.
#
#
CalculateDietDiversity <- function(predator.taxa, shapes, occurrence.source, opts, name.filter = function (x) {x}) {
  predator.taxon.names <- name.filter(GetPredatorTaxa(predator.taxa=predator.taxa, limit=opts$predator.species.limit, opts = opts))
  diet.matrix <- CreateDietMatrix(predator.taxon.names = predator.taxon.names, opts = opts)
  diet.matrix <- DiversityVsNicheWidth(diet.matrix, shapes = shapes, opts = c(opts, list(occurrence.source = occurrence.source, occurrence.limit = opts$occurrence.limit)))
  predator.title <-  paste(predator.taxa, collapse='_')
  diet.filename <- paste('hayden2014', predator.title, occurrence.source, shapes$name, format(Sys.time(), '%Y-%m-%d-%H%M%S'), sep='_')
  ExportToFile(diet.matrix, file = paste(diet.filename, '.csv', sep=''))
  list(diet.matrix=diet.matrix, occurrence.source=occurrence.source, opts=opts, diversity.shapes=shapes) 
}

library(rfishbase)
data(fishbase)

SpeciesNamesForHabitat <- function(habitat) {
  species_ids <- which_fish(habitat, "habitat", fish.data)
  fish_names(fish.data[species_ids])
}

IncludeMarineButNotFreshwaterSpecies <- function(speciesNames) {
  speciesNames[speciesNames %in% SpeciesNamesForHabitat('marine') 
    & !speciesNames %in% SpeciesNamesForHabitat('freshwater')]
}

IncludeFreshwaterButNotMarineSpecies <- function(speciesNames) {
  speciesNames[!speciesNames %in% SpeciesNamesForHabitat('marine') 
    & speciesNames %in% SpeciesNamesForHabitat('freshwater')]
}


# Configuration options (change for access to "Dark" GloBI)
opts.dark <- list(port=80, predator.species.limit = 128, occurrence.limit = 32)

CompileGeoDiet <- function (consumers = list('Actinopterygii'), opts = opts.dark) {
  results <- list()
  tittensor.shapes <- ReadBiodiversityShapes(shapes = GetTittensorShapes())
  results$tittensor <- CalculateDietDiversity(consumers, tittensor.shapes, 'GBIF', opts = opts, name.filter = IncludeMarineButNotFreshwaterSpecies)
  
  tisseuil.shapes <- ReadBiodiversityShapes(shapes = GetTisseuilShapes()) 
  results$tisseuil <- CalculateDietDiversity(consumers, tisseuil.shapes, 'GBIF', opts = opts, name.filter=IncludeFreshwaterButNotMarineSpecies)
  results
}
# Plot Global Dietary Niche Maps by overlaying a lat lng raster onto a map of the world.


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

GetData <- function(datafile) {
  paste('data/', datafile, sep='')
}

GetFigure <- function(datafile) {
  paste('figs/', datafile, sep='')
}

PlotFreshwaterMaps <- function(biodiversityShapes = ReadBiodiversityShapes(GetTisseuilShapes()), format='pdf') {
  fresh.basins <- biodiversityShapes
  fresh.niche.data <- read.csv(GetData('fresh.data.csv'))
  fresh.trophic.data <- read.csv(GetData('fresh.TL.data.csv'))
  fresh.basins <- merge(fresh.basins$diversity.map, fresh.niche.data, by='ID_JF')
  fresh.basins <- merge(fresh.basins, fresh.trophic.data, by='ID_JF')

  # print plots
  lapply(graph_formats, function(format) {
    setup_device(filename = GetFigure('fresh.niche.median'), format=format)  
    map(interior=FALSE, boundary=FALSE, add=FALSE, fill=TRUE, border='grey', col='grey')
    plotNicheMap(fresh.basins, param.name='niche.med', main.title='Median DNW')
    dev.off()
  })
  
  lapply(graph_formats, function(format) {
    setup_device(filename= GetFigure('fresh.niche.mean'), format=format)  
    map(interior=FALSE, boundary=FALSE, add=FALSE, fill=TRUE, border='grey', col='grey')
    plotNicheMap(fresh.basins, param.name='niche.mean', main.title='Mean DNW')
    dev.off()
  })

  lapply(graph_formats, function(format) {
    setup_device(filename= GetFigure('fresh.trophic.mean'), format=format)  
    map(interior=FALSE, boundary=FALSE, add=FALSE, fill=TRUE, border='grey', col='grey')
    plotNicheMap(fresh.basins, param.name='TL.mean', main.title='Mean TL')
    dev.off()
  })
}

PlotMarineMaps <- function(biodiversityShapes = ReadBiodiversityShapes()) {
  marine.squares <- biodiversityShapes
  marine.niche.data <- read.csv(GetData('marine.sumry.csv'))
  marine.trophic.data <- read.csv(GetData('marine.TL.sumry.csv'))
  marine.squares <- merge(marine.squares$diversity.map, marine.niche.data, by='GRIDCODE')
  marine.squares <- merge(marine.squares, marine.trophic.data, by='GRIDCODE')

  # print plots

  lapply(graph_formats, function(format) {
    setup_device(file=GetFigure('marine.niche.median'), format= format)  
    plotNicheMap(marine.squares, param.name='niche.med', main.title='Median DNW', add = FALSE)
    map(interior=FALSE, boundary=FALSE, add=TRUE, fill=TRUE, border='gray', col='gray')
    dev.off()
  });
  
  lapply(graph_formats, function(format) {
    setup_device(filename=GetFigure('marine.niche.mean'), format = format)  
    plotNicheMap(marine.squares, param.name='niche.mean', main.title='Mean DNW', add = FALSE)
    map(interior=FALSE, boundary=FALSE, add=TRUE, fill=TRUE, border='gray', col='gray')
    dev.off()
  })

  lapply(graph_formats, function(format) {
    setup_device(filename=GetFigure('marine.trophic.mean'), format = format)  
    plotNicheMap(marine.squares, param.name='TL.mean', main.title='Mean TL', add = FALSE)
    map(interior=FALSE, boundary=FALSE, add=TRUE, fill=TRUE, border='gray', col='gray')
    dev.off()
  })
}

message('PlotMarineMaps() : plot marine diet niche maps');
message('PlotFreshwaterMaps() : plot fresh water diet niche maps')
message('CompileGeoDiet(): compile geo diet for fishes (Actinopterygii)')
