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

# Retrieve first known location for species with name
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
