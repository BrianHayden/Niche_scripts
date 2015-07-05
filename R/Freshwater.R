setwd("~/Dropbox/Projects/Niche/Freshwater/Analysis")

library(foreign)
library(rgdal)
library(maptools)
library(plyr)
library(plotrix)
library(fields)
library(mgcv)
library(car)
library('raster')
library('maps')

#Get basin centroids and data
basindata<-read.dbf("Basin.dbf")
SPRICH<-read.csv("FISH.SPRICH.csv")
basins<-merge(basindata, SPRICH, by="BASIN")


##Test DNW
# Read diet matrix
fresh <- read.csv("~/Dropbox/Projects/Niche/Freshwater/Analysis/hayden2014_Actinopterygii_GBIF_tisseuil2013_2015-02-22-002548.csv", header=T)
fresh <-subset(fresh, niche.width > 0)

#Remove basins with <30 species
fresh.working<-ddply(fresh, c("ID_JF"), summarise,
                   N=length(niche.width))
fresh.working<-subset(fresh.working, N > 20)                   
fresh.working<-merge(fresh.working, fresh, by = "ID_JF")
IDs<-read.csv("~/Dropbox/Projects/Niche/Freshwater/Analysis/ID.csv")
fresh.working<-merge(IDs, fresh.working, by = "ID_JF")
                   
#Resample to 1000 datapoints per basin
fresh.sp <- split(fresh.working, fresh.working$ID)
fresh.samples <- lapply(fresh.sp, function(x) x[sample(1:nrow(x), 1000, TRUE),])
fresh.new <- do.call(rbind, fresh.samples)

fresh.sumry<-ddply(fresh.new, c("ID_JF"), summarise,
                   N=length(niche.width),
                   niche.mean=mean(niche.width),
                   niche.med=median(niche.width),
                   niche.sd=sd(niche.width),
                   bio=mean(fish_endem))

fresh.species<-ddply(fresh.working, c("predator.taxon.name"), summarise,
                     N=length(niche.width))

#read and add envirmontal data
enviros<-read.csv("~/Dropbox/Projects/Niche/Freshwater/Analysis/environmentals.csv")
fresh.data<-merge(fresh.sumry, enviros, by = "ID_JF")
fresh.data<-fresh.data[-c(20,26),]
fresh.data$Alt<-as.numeric(fresh.data$Alt)
fresh.data$Area<-as.numeric(fresh.data$Area)

pairs(~ niche.mean + bio + long + lat +Temp + Alt + Area, data=fresh.data, panel = panel.smooth)
pairs(~ niche.med + bio + long + lat +Temp + Alt + Area, data=fresh.data, panel = panel.smooth)
pairs(~ niche.sd + bio + long + lat +Temp + Alt + Area, data=fresh.data, panel = panel.smooth)

## GAM to match marine analysis
# Interaction term to account for spatial auriocorrelation
#Mean niche width
gamF1<-gam(niche.mean ~ s(bio, k=5) + s(lat, k=5) + s(Temp, k=5) + s(Alt, k=5) + s(Area, k=5) + s(long, lat, k = 30), data=fresh.data)
summary(gamF1)
plot(gamF1, pages=1)
gam.check(gamF1)
AIC(gamF1)

# Remove non significant terms, starting with Area 
gamF1.1<-gam(niche.mean ~ s(bio, k=5) + s(lat, k=5) + s(Temp, k=5) + s(Alt, k=5) + s(long, lat, k = 30), data=fresh.data)
summary(gamF1.1)
plot(gamF1.1, pages=1)
gam.check(gamF1.1)
anova(gamF1, gamF1.1, test="F")

#Remove temp
gamF1.2<-gam(niche.mean ~ s(bio, k=5) + s(lat, k=5) + s(Alt, k=5) + s(long, lat, k = 30), data=fresh.data)
summary(gamF1.2)
plot(gamF1.2, pages=1)
gam.check(gamF1.2)
anova(gamF1.1, gamF1.2, test="F")

#Remove lat
gamF1.3<-gam(niche.mean ~ s(bio, k=5) + s(Alt, k=5) + s(long, lat, k = 30), data=fresh.data)
summary(gamF1.3)
plot(gamF1.3, pages=1)
gam.check(gamF1.3)
anova(gamF1.1, gamF1.2, test="F")

#############      #Median Niche width
gamF2<-gam(niche.med ~ s(bio, k=5) + s(lat, k=5) + s(Temp, k=10) + s(Alt, k=5) + s(Area, k=5) + s(long, lat, k = 30), data=fresh.data)
summary(gamF2)
plot(gamF2, pages=1)
gam.check(gamF2)

#Remove non-significant terms, starting with Area
gamF2.1<-gam(niche.med ~ s(bio, k=5) + s(lat, k=5) + s(Temp, k=5) + s(Alt, k=10) + s(long, lat, k = 30), data=fresh.data)
summary(gamF2.1)
plot(gamF2.1, pages=1)
gam.check(gamF2.1)
anova(gamF2, gamF2.1, test="F")

#Remove lat
gamF2.2<-gam(niche.med ~ s(bio, k=5) + s(Temp, k=5) + s(Alt, k=10) + s(long, lat, k = 30), data=fresh.data)
summary(gamF2.2)
plot(gamF2.2, pages=1)
gam.check(gamF2.2)
anova(gamF2.2, gamF2.1, test="F")

#Remove temp
gamF2.3<-gam(niche.med ~ s(bio, k=5) + s(Alt, k=10) + s(long, lat, k = 30), data=fresh.data)
summary(gamF2.3)
plot(gamF2.3, pages=1)
gam.check(gamF2.3)
anova(gamF2.2, gamF2.3, test="F")

#########   SD of mean niche width    ##########################

gamF3<-gam(niche.sd ~ s(bio, k=3) + s(lat, k=3) + s(Temp, k=3) + s(Alt, k=3) + s(Area, k=3) + s(long, lat, k = 30), data=fresh.data)
summary(gamF3)
plot(gamF3, pages=1)
gam.check(gamF3)

#Remove non-significant terms, starting with Area
gamF3.1<-gam(niche.sd ~ s(bio, k=3) + s(lat, k=3) + s(Temp, k=3) + s(Alt, k=3) + s(long, lat, k = 30), data=fresh.data)
summary(gamF3.1)
plot(gamF3.1, pages=1)
gam.check(gamF3.1)
anova(gamF3, gamF3.1, test="F")

#Remove temp
gamF3.2<-gam(niche.sd ~ s(bio, k=3) + s(lat, k=3) + s(Alt, k=3) + s(long, lat, k = 30), data=fresh.data)
summary(gamF3.2)
plot(gamF3.2, pages=1)
gam.check(gamF3.2)
anova(gamF3.2, gamF3.1, test="F")

#Remove lat
gamF3.3<-gam(niche.sd ~ s(bio, k=3) + s(Alt, k=3) + s(long, lat, k = 30), data=fresh.data)
summary(gamF3.3)
plot(gamF3.3, pages=1)
gam.check(gamF3.3)
anova(gamF3.2, gamF3.3, test="F")

#Remove lat
gamF3.4<-gam(niche.sd ~ s(Alt, k=3) + s(long, lat, k = 30), data=fresh.data)
summary(gamF3.4)
plot(gamF3.4, pages=1)
gam.check(gamF3.4)
anova(gamF3.4, gamF3.3, test="F")

#############         TEST  TROPHIC LEVEL       ###################

## Read freshwater diet matrix
fresh <- read.csv("~/Dropbox/Projects/Niche/Freshwater/Analysis/hayden2014_Actinopterygii_GBIF_tisseuil2013_2015-02-22-002548.csv", header=T)
fresh <-subset(fresh, niche.width > 0)

# Read trophic levels
TL <- read.csv("~/Dropbox/Projects/Niche/Trophic.csv", header=T)
fresh.TL <-merge(fresh, TL, by = "predator.taxon.name")
fresh.TL <- subset(fresh.TL, foodTroph > 1)

#Remove basins with <30 species
fresh.TL.working<-ddply(fresh.TL, c("ID_JF"), summarise,
                     N=length(niche.width),
                     TL=mean(foodTroph))
fresh.TL.working<-subset(fresh.TL.working, N > 20)                   
fresh.TL.working<-merge(fresh.TL.working, fresh.TL, by = "ID_JF")
IDs<-read.csv("ID.csv")
fresh.TL.working<-merge(IDs, fresh.TL.working, by = "ID_JF")

#Resample to 1000 datapoints per basin
fresh.TL.sp <- split(fresh.TL.working, fresh.TL.working$ID)
fresh.TL.samples <- lapply(fresh.TL.sp, function(x) x[sample(1:nrow(x), 1000, TRUE),])
fresh.TL.new <- do.call(rbind, fresh.TL.samples)

fresh.TL.sumry<-ddply(fresh.TL.new, c("ID_JF"), summarise,
                   N=length(niche.width),
                   TL.mean=mean(foodTroph),
                   TL.med=median(foodTroph),
                   bio=mean(fish_endem))

## Read environmental characteristic of basins
enviros<-read.csv("environmentals.csv")
fresh.TL.data<-merge(fresh.TL.sumry, enviros, by = "ID_JF")
fresh.TL.data<-fresh.TL.data[-c(22,26),]
fresh.TL.data$Area<-as.numeric(fresh.TL.data$Area)
fresh.TL.data$Alt<-as.numeric(fresh.TL.data$Alt)
write.csv(fresh.TL.data, "fresh.TL.data.csv")

#Pairs plots to eyeball data
pairs(~ TL.mean + bio + long + lat + Temp + Area + Alt, data=fresh.TL.data, panel = panel.smooth)
pairs(~ TL.med + bio + long + lat + Temp + Area + Alt, data=fresh.TL.data, panel = panel.smooth)

#GAM of mean trophic level
gamTLF1<-gam(TL.mean ~ s(bio, k=5) + s(lat, k=5) + s(Temp, k=5) + s(Alt, k=5) + s(Area, k=10) + s(long, lat, k = 30), data=fresh.TL.data)
summary(gamTLF1)
plot(gamTLF1, pages=1)
gam.check(gamTLF1)

#Drop Area
gamTLF1.1<-gam(TL.mean ~ s(bio, k=5) + s(lat, k=5) + s(Temp, k=5) + s(Alt, k=5) + s(long, lat, k = 30), data=fresh.TL.data)
summary(gamTLF1.1)
plot(gamTLF1.1, pages=1)
gam.check(gamTLF1.1)
anova(gamTLF1, gamTLF1.1, test="F")

#Drop Lat
gamTLF1.2<-gam(TL.mean ~ s(bio, k=5) + s(Temp, k=5) + s(Alt, k=5) + s(long, lat, k = 30), data=fresh.TL.data)
summary(gamTLF1.2)
plot(gamTLF1.2, pages=1)
gam.check(gamTLF1.2)
anova(gamTLF1.1, gamTLF1.2, test="F")

#Drop Temp
gamTLF1.3<-gam(TL.mean ~ s(bio, k=5) + s(Alt, k=5) + s(long, lat, k = 30), data=fresh.TL.data)
summary(gamTLF1.3)
plot(gamTLF1.3, pages=1)
gam.check(gamTLF1.3)
anova(gamTLF1.3, gamTLF1.2, test="F")


#########################################################################################################
#GAM of median trophic level
gamTLF2<-gam(TL.med ~ s(bio, k=5) + s(lat, k=5) + s(Temp, k=5) + s(Alt, k=5) + s(Area, k=10) + s(long, lat, k = 30), data=fresh.TL.data)
summary(gamTLF2)
plot(gamTLF2, pages=1)
gam.check(gamTLF2)

#Drop Area
gamTLF2.1<-gam(TL.med ~ s(bio, k=5) + s(lat, k=5) + s(Temp, k=5) + s(Alt, k=5) + s(long, lat, k = 30), data=fresh.TL.data)
summary(gamTLF2.1)
plot(gamTLF2.1, pages=1)
gam.check(gamTLF2.1)
anova(gamTLF2, gamTLF2.1, test="F")

#Drop Lat
gamTLF2.2<-gam(TL.med ~ s(bio, k=5) + s(Temp, k=5) + s(Alt, k=5) + s(long, lat, k = 30), data=fresh.TL.data)
summary(gamTLF2.2)
plot(gamTLF2.2, pages=1)
gam.check(gamTLF2.2)
anova(gamTLF2, gamTLF2.2, test="F")


# Rough working plots
plot(fresh.data$bio, fresh.data$niche.mean, col=fresh.data$Brealm)
plot(fresh.sumry$INSIDE_X, fresh.sumry$INSIDE_Y)

plot(fresh.sumry$bio, fresh.sumry$niche.mean)
plot(fresh.data$Lonmed, fresh.data$Latmed)
plot(basindata$INSIDE_X, basindata$INSIDE_Y)


## Data plots

dev.off()
dev.new()
#Plot bio - DNW by region

par(mfrow=c(2,2))
plot(NULL, NULL, ylim=c(0.1, 0.3), xlim=c(0, 1), main ="a",
     ylab = "Mean niche width", xlab = "Biodiversity (endemic fish richness")
points(fresh.plots$Palearctic$bio, fresh.plots$Palearctic$niche.mean, pch=20, cex=1.5, col="deepskyblue4")
points(fresh.plots$Nearctic$bio, fresh.plots$Nearctic$niche.mean, pch=20, cex=1.5, col="deepskyblue2")
points(fresh.plots$Afrotropical$bio, fresh.plots$Afrotropical$niche.mean, cex=1.5, pch=20, col="mediumaquamarine")
points(fresh.plots$Neotropical$bio, fresh.plots$Neotropical$niche.mean, cex=1.5, pch=20, col="lightgoldenrod2")
points(fresh.plots$Australian$bio, fresh.plots$Australian$niche.mean, cex=1.5, pch=20, col="salmon3")
points(fresh.plots$Oriental$bio, fresh.plots$Oriental$niche.mean, pch=20, cex=1.5, col="red")
legend("topright", pch=20, cex = 0.5,
       col=c("deepskyblue4", "deepskyblue2", "mediumaquamarine", "lightgoldenrod2", "salmon3", "red"),
       legend=c("Palearctic", "Nearctic", "Afrotropical", "Neotropical", "Australian", "Oriental"))


### DNW plots
plot.gam(gamF1.3, select=1, ylim=c(-0.075, 0.075), seWithMean=T,scheme =1, all.terms=T, 
         residuals=T, cex=0.4, pch=20,
         rug=T, col="lightgrey", lwd=1,shade.col="lightblue",
         main="b",
         ylab= "Mean niche width (GAM)",
         xlab= "Biodiversity (endemic fish richness)")

plot.gam(gamF3.3, select=1, ylim=c(-0.03, 0.03), seWithMean=T,scheme =1, all.terms=T, 
         residuals=T, cex=0.4, pch=20,
         rug=T, col="lightgrey", lwd=1,shade.col="lightblue",
         main="c",
         ylab= "Standard deviation of mean niche width (GAM)",
         xlab= "Biodiversity (endemic fish richness)")


plot.gam(gamTLF1.3, select=1, ylim=c(-0.1, 0.2), seWithMean=T,scheme =1, all.terms=T, 
         residuals=T, cex=0.4, pch=20,
         rug=T, col="lightgrey", lwd=1,shade.col="lightblue",
         main="d",
         ylab= "Mean trophic level (GAM)",
         xlab= "Biodiversity (endemic fish richness)")



#### TL Plots
plot.gam(gamTLF1.3, select=1, ylim=c(-0.1, 0.2), seWithMean=T,scheme =1, all.terms=T, 
         residuals=T, cex=0.4, pch=20,
         rug=T, col="lightgrey", lwd=1,shade.col="lightblue",
         main="c",
         ylab= "Mean trophic level (GAM)",
         xlab= "Biodiversity (endemic fish richness)")

boxplot(niche.mean ~ lat, data=fresh.data,
        main="b",
        ylab="Mean DNW (0 - 1)",
        xlab="Latitude (oN)",
        outline=F,
        col="lightblue")

bplot.xy(fresh.data$niche.mean, fresh.data$lat, N = 16)



#Plot global heat map
if (!'package:raster' %in% search()) {
  install.packages(c('raster','maps'), repos='http://cran.us.r-project.org')
  library('raster')
  library('maps')
}

##Plot DNW map
plotMap <- function(dataFilename) {
  # create raster with cells 5' x 5'
  r <- raster()
  nrow(r) <- 36
  ncol(r) <- 72
  
  cat(paste(dataFilename, 'loading...'))
  niche <- read.csv(dataFilename)
  cat(paste(dataFilename, 'loaded.'))
  
  # extract lat lng with their values (niche width)
  latLng <- cbind(niche$long, niche$lat)
  vals <- niche$niche.mean
  rasterWithData <- rasterize(latLng, r, vals, fun=median)
  mapFilename = 'DNWmap_fresh_>20sp_mean.png'
  png(file=mapFilename, width=640, height=320)
  world <- map(interior=FALSE)
  plot(rasterWithData, main='Median DNW')
  # add world contours for reference
  lines(world, col=rgb(0,0,1,1.0))
  dev.off()
  cat(paste('wrote distribution graph to', mapFilename, '\n'))
} 
plotMap("fresh.data.csv")

##Plot TL map
plotMap <- function(dataFilename) {
  # create raster with cells 5' x 5'
  r <- raster()
  nrow(r) <- 36
  ncol(r) <- 72
  
  cat(paste(dataFilename, 'loading...'))
  niche <- read.csv(dataFilename)
  cat(paste(dataFilename, 'loaded.'))
  
  # extract lat lng with their values (niche width)
  latLng <- cbind(niche$long, niche$lat)
  vals <- niche$TL.med
  rasterWithData <- rasterize(latLng, r, vals, fun=mean)
  mapFilename = 'TLmap_fresh_>20sp_median.png'
  png(file=mapFilename, width=640, height=320)
  world <- map(interior=FALSE)
  plot(rasterWithData, main='Median TL')
  # add world contours for reference
  lines(world, col=rgb(0,0,1,1.0))
  dev.off()
  cat(paste('wrote distribution graph to', mapFilename, '\n'))
} 
plotMap("fresh.TL.data.csv")
