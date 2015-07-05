# Generate Diet Matrices using https://github.com/ropensci/rglobi 
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
