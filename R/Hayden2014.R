# Recipe to calculate diet niche and associated diversity for species.
# 
# Steps:
# 1. install dependencies using install.packages(c('devtools','maptools','rgbif','data.table','RCurl'))
# 2. install rglobi using devtools::install_github('ropensci/rglobi')
# 3. install rfishbase using devtools::install_github('ropensci/rfishbase')
# 4. devtools::source_url('https://gist.githubusercontent.com/jhpoelen/81b4eced3963633b96cb/raw/dietMatrix.R')
# 5. devtools::source_url('https://gist.githubusercontent.com/jhpoelen/81b4eced3963633b96cb/raw/GetBiodiversityAtPoints.R')
# 6. devtools::source_url('https://gist.githubusercontent.com/jhpoelen/81b4eced3963633b96cb/raw/Hayden2014.R')
#
#
# Configuration options
opts.dark <- list(port=7574, predator.species.limit = 128, occurrence.limit = 32)

CalculateDietDiversity <- function(predator.taxa, shapes, occurrence.source, opts, name.filter = function (x) {x}) {
  predator.taxon.names <- name.filter(GetPredatorTaxa(predator.taxa=predator.taxa, limit=opts$predator.species.limit, opts = opts))
  diet.matrix <- CreateDietMatrix(predator.taxon.names = predator.taxon.names, opts = opts)
  diet.matrix <- DiversityVsNicheWidth(diet.matrix, shapes = shapes, opts = c(opts, list(occurrence.source = occurrence.source, occurrence.limit = opts$occurrence.limit)))
  predator.title <-  paste(predator.taxa, collapse='_')
  diet.filename <- paste('hayden2014', predator.title, occurrence.source, shapes$name, format(Sys.time(), '%Y-%m-%d-%H%M%S'), sep='_')
  ExportToFile(diet.matrix, file = paste(diet.filename, '.csv', sep=''))
  plot.title <- paste(predator.title, '\n<=', opts$occurrence.limit, ' lat/lng per species\n', occurrence.source, '  occurrences;', shapes$name, ' diversity map', sep='')
  
  figure.filename <- paste(diet.filename, '.png', sep='')
  png(file=figure.filename)
  niche.model <- CreateAndPlotNicheWidthModel(diet.matrix[diet.matrix[[shapes$diversity.metric]]>0], title=plot.title, diversity.metric=shapes$diversity.metric)
  dev.off()
  message(paste('wrote file [', figure.filename, ']', sep='')) 
  list(diet.matrix=diet.matrix, niche.model=niche.model, occurrence.source=occurrence.source, opts=opts, diversity.shapes=shapes) 
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

RunCalculations <- function (consumers = list('Actinopterygii')) {
  results <- list()
  tittensor.shapes <- ReadBiodiversityShapes(shapes = GetTittensorShapes())
  results$tittensor <- CalculateDietDiversity(consumers, tittensor.shapes, 'GBIF', opts = opts.dark, name.filter = IncludeMarineButNotFreshwaterSpecies)
  
  tisseuil.shapes <- ReadBiodiversityShapes(shapes = GetTisseuilShapes()) 
  results$tisseuil <- CalculateDietDiversity(consumers, tisseuil.shapes, 'GBIF', opts = opts.dark, name.filter=IncludeFreshwaterButNotMarineSpecies)
  results
}
