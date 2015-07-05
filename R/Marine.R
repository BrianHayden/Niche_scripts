##Load relevant packages
library(plyr)
library(plotrix)
library(fields)
library(mgcv)
library(car)

##Load diet matrix
setwd("~/Dropbox/Projects/Niche")
marine<-read.delim("hayden2014_Actinopterygii_GBIF_tittensor2010_4500sp.100k.csv", header=T, sep=",")

##Load summary data including SST, Productivity, Dist to Coast (created seperately)
marine.data<-read.delim("marine.data.csv", sep=",", header=T)

##Load trophic levels
trophic<-read.delim("trophic.csv", header=T, sep=",")


##Merge daatsets
marine<-merge(marine, trophic, by="predator.taxon.name")
marine<-merge(marine, marine.data, by="GRIDCODE")


marine<-subset(marine, select = c(predator.taxon.name, GRIDCODE, 
                                   niche.width, lat, lat.t, long, SST, Prod, 
                                   Dist_Coast, AllNorm, foodTroph))
marine<-arrange(marine, GRIDCODE)

marine<-subset(marine, AllNorm > 0)
marine<-subset(marine, niche.width > 0)



#Check sample size
marine.sum<-ddply(marine, c("GRIDCODE"), summarise,
                       N=length(niche.width))
marine.sum<-arrange(marine.sum, N)

##Remove gridsqaures with <30 species (done in excel.....)
write.csv(marine, file="marinexx.csv")
marine.new<-read.csv("~/Dropbox/Projects/Niche/marinexy.csv", sep=",", header = T)

##Remove 0'S
marine.new<-subset(marine.new, AllNorm > 0)
marine.new<-subset(marine.new, niche.width > 0)

##Remove distance outlier gridsqaure 502
marine.new<-subset(marine.new, Dist_Coast < 4000)

##Remove 3 productivity outliers
marine.new<-subset(marine.new, Prod < 1)

#Resample to 1000 rows for each gridsquare
sp <- split(marine.new, marine.new$GRIDCODE)
samples <- lapply(sp, function(x) x[sample(1:nrow(x), 1000, TRUE),])
marine.new <- do.call(rbind, samples)
write.csv(marine.sumry, file="marine.sumry.csv")


##Summarise by gridsquare

marine.sumry<-ddply(marine.new, c("GRIDCODE"), summarise,
                       lat=min(lat),
                      lat.t=min(lat.t),
                       long=min(long),
                       N=length(niche.width),
                       niche.mean=mean(niche.width),
                       niche.med=median(niche.width),
                       niche.sd=sd(niche.width),
                       niche.min=min(niche.width),
                       niche.max=max(niche.width),
                       niche.range=niche.max - niche.min,
                       SST=min(SST),
                       Prod=min(Prod),
                       Dist=min(Dist_Coast),
                       bio=mean(AllNorm),
                       TL.mean=mean(foodTroph),
                       TL.med=median(foodTroph))


##Visualise data with pairs plots
plot(marine.new$AllNorm, marine.new$niche.width)
pairs(~ niche.mean + bio + Dist+ lat + SST + Prod, data = marine.sumry, panel = panel.smooth)
pairs(~ niche.mean + bio + Dist+ lat + SST + Prod, data = marine.sumry, 
      panel = function(x,y){points(x,y);lines(lowess(x,y))},
      )


pairs(~ niche.mean + bio + Dist + lat + SST + Prod, data = marine.sumry, panel = panel.smooth)
pairs(~ niche.sd + bio + Dist + N + lat + SST + Prod + TL.mean, data = marine.sumry, panel = panel.smooth)

pairs(~ TL.mean + bio.mean + Dist_Coast + N + lat + SST + Prod, data = marine.data.TL, panel = panel.smooth)
pairs(~ TL.med + bio.mean + Dist_Coast + N + lat + SST + Prod, data = marine.data.TL, panel = panel.smooth)

## Model niche and TL independently

## Non linear relationsips (despite many trials with transformations etc) so uisng a GAM rather than GLM
#Niche width
  
gam1<-gam(niche.mean ~ s(bio, k=5) + s(lat, k=5) + s(Dist, k=5) + s(Prod, k=5) + s(SST, k=8) + s(lat,long, k=70), data=marine.sumry)
summary(gam1)
plot(gam1)
gam.check(gam1)

# Prod non significant so try removing it from the model
gam1.1<-gam(niche.mean ~ s(bio, k=5) + s(lat, k=5) + s(Dist, k=5) + s(SST, k=8) + s(lat,long, k=70), data=marine.sumry)
summary(gam1.1)
plot(gam1.1)
gam.check(gam1.1)
# looks good but are the models significantly different
anova(gam1, gam1.1, test="F")  #No, therefore drop productivity

#Dist_Coast non-significant sop try removing it also
gam1.2<-gam(niche.mean ~ s(bio, k=8) + s(lat, k=5) + s(SST, k=8) + s(lat,long, k=70), data=marine.sumry)
summary(gam1.2)
plot(gam1.2, pages=1)
gam.check(gam1.2)
# looks good but are the models significantly different
anova(gam1.2, gam1.1, test="F")  #No, therefore drop Dist_Coast
## gam1.2 is the simplest model explaining variation in mean niche wdith.



## Next we look at variation in median ncihe width - k term for SST is increased slightly due to higher edf
gam2<-gam(niche.med ~ s(bio, k=5) + s(lat, k=5) + s(Dist, k=5) + s(Prod, k=5) + s(SST, k=8), data=marine.sumry)
summary(gam2)
plot(gam2)
gam.check(gam2)

#Dist_Coast non-significant sop try removing it
gam2.1<-gam(niche.med ~ s(bio, k=5) + s(lat, k=5) + s(Prod, k=5) + s(SST, k=8) + s(lat,long, k=70), data=marine.sumry)
summary(gam2.1)
plot(gam2.1, pages=1)
gam.check(gam2.1)
# looks good but are the models significantly different
anova(gam2, gam2.1, test="F")  #No, therefore remove Dist_Coast from the main model

#Prod non-significant so try removing it
gam2.2<-gam(niche.med ~ s(bio, k=5) + s(lat, k=5) + s(SST, k=8) + s(lat,long, k=70), data=marine.sumry)
summary(gam2.2)
plot(gam2.2, pages=1)
gam.check(gam2.2)
# looks good but are the models significantly different
anova(gam2.2, gam2.1, test="F")  #No, therefore remove Dist_Coast from the main model
#gam2.2 simplest model explaining variation in median niche width


## Next up, varition instandard deviation of niche width
gam3<-gam(niche.sd ~ s(bio, k=8) + s(lat, k=5) + s(Dist, k=5) + s(Prod, k=5) + s(SST, k=8) + s(lat,long, k=70), data=marine.sumry)
summary(gam3)
gam.check(gam3)
#bio, lat and dist are non significant so trying removing these terms

#Remove lat
gam3.1<-gam(niche.sd ~ s(bio, k=8) + s(Dist, k=5) + s(Prod, k=5) + s(SST, k=8) + s(lat,long, k=70), data=marine.sumry)
summary(gam3.1)
gam.check(gam3.1)
anova(gam3, gam3.1, test="F")
# No effect removing latitude, so drop it from the model

#Remove Dist
gam3.2<-gam(niche.sd ~ s(bio, k=8) + s(Prod, k=5) + s(SST, k=8) + s(lat,long, k=70), data=marine.sumry)
summary(gam3.2)
gam.check(gam3.2)
anova(gam3.1, gam3.2, test="F")
# No effect removing Distance, so drop it from the model

#Remove bio
gam3.3<-gam(niche.sd ~ s(Prod, k=5) + s(SST, k=8) + s(lat,long, k=70), data=marine.sumry)
summary(gam3.3)
gam.check(gam3.3)
anova(gam3.2, gam3.3, test="F")
plot(gam3.3)
# No effect removing bio, so drop it from the model

#Remove Prod
gam3.4<-gam(niche.sd ~  s(SST, k=8) + s(lat,long, k=70), data=marine.sumry)
summary(gam3.4)
gam.check(gam3.4)
anova(gam3.3, gam3.4, test="F")
plot(gam3.4)
# No effect removing bio, so drop it from the model
## gam3.4 is the simplest model explaining variation in standard deviation in niche width

### Split oceanic and coastal regions



### Trophic level

marine.new.tl<-subset(marine.new, foodTroph > 1)
write.csv(marine.new.tl, file="marine.new.TL.csv")

marine.tl.sumry<-ddply(marine.new.tl, c("GRIDCODE"), summarise,
                    lat=min(lat),
                    long=min(long),
                    SST=min(SST),
                    Prod=min(Prod),
                    Dist=min(Dist_Coast),
                    bio=mean(AllNorm),
                    TL.mean=mean(foodTroph),
                    TL.med=median(foodTroph))

gamTL1<-gam(TL.mean ~ s(bio, k=5) + s(lat, k=5) + s(Dist, k=5) + s(Prod, k=5) + s(SST, k=5) + s(lat,long, k=60), data=marine.tl.sumry)
summary(gamTL1)
plot(gamTL1, pages=1)
gam.check(gamTL1)

#Drop non significant terms

# SST
gamTL1.1<-gam(TL.mean ~ s(bio, k=5) + s(Dist, k=5) + s(Prod, k=5) + s(lat, k=5) + s(lat,long, k=60), data=marine.tl.sumry)
summary(gamTL1.1)
plot(gamTL1.1, pages=1)
gam.check(gamTL1.1)
anova(gamTL1, gamTL1.1, test="F")

# lat
gamTL1.2<-gam(TL.mean ~ s(bio, k=5) + s(Dist, k=5) + s(Prod, k=5) + s(lat,long, k=60), data=marine.tl.sumry)
summary(gamTL1.2)
plot(gamTL1.2, pages=1)
gam.check(gamTL1.2)
anova(gamTL1.2, gamTL1.1, test="F")






##Median Niche width
gamTL2<-gam(TL.med ~ s(bio, k=5) + s(lat, k=5) + s(Dist, k=5) + s(Prod, k=5) + s(SST, k=5) + s(lat,long, k=60), data=marine.tl.sumry)
summary(gamTL2)
gam.check(gamTL2)
plot(gamTL2, pages=1)
#Lat:Long interaction is the only significant term so try dropping the others

#Drop Bio
gamTL2.1<-gam(TL.med ~ s(lat, k=5) + s(Dist, k=5) + s(Prod, k=5) + s(SST, k=5) + s(lat,long, k=60), data=marine.tl.sumry)
summary(gamTL2.1)
gam.check(gamTL2.1)
anova(gamTL2, gamTL2.1, test="F")

#Drop Dist
gamTL2.2<-gam(TL.med ~ s(lat, k=5) + s(Prod, k=5) + s(SST, k=5) + s(lat,long, k=60), data=marine.tl.sumry)
summary(gamTL2.2)
gam.check(gamTL2.2)
anova(gamTL2.2, gamTL2.1, test="F")

#Drop Prod
gamTL2.3<-gam(TL.med ~ s(lat, k=5) + s(SST, k=5) + s(lat,long, k=60), data=marine.tl.sumry)
summary(gamTL2.3)
gam.check(gamTL2.3)
anova(gamTL2.2, gamTL2.3, test="F")

#Drop SST
gamTL2.4<-gam(TL.med ~  s(lat, k=5) + s(lat,long, k=60), data=marine.tl.sumry)
summary(gamTL2.4)
gam.check(gamTL2.4)
anova(gamTL2.4, gamTL2.3, test="F")

#Drop Lat
gamTL2.5<-gam(TL.med ~ s(lat,long, k=60), data=marine.tl.sumry)
summary(gamTL2.5)
gam.check(gamTL2.5)
anova(gamTL2.4, gamTL2.5, test="F")








#Plots

par(mfrow=c(3,2))


plot(marine.sumry$bio, marine.sumry$niche.mean,
     pch=16, col="lightblue",
     ylim=c(0.1, 0.35),
     xlab="Biodiversity (0 - 1)",
     main="a",
     ylab="Mean DNW (0 - 1)")

boxplot(niche.mean ~ lat, data=marine.sumry,
        main="b",
        ylab="Mean DNW (0 - 1)",
        xlab="Latitude (oN)",
        outline=F,
        col="lightblue")

plot.gam(gam1.2, select=1, ylim=c(-0.075, 0.075), seWithMean=T,scheme =1, all.terms=T, 
         residuals=T, cex=0.4, pch=20,
         rug=T, col="lightgrey", lwd=1, shade.col="lightblue",
         main="c",
         ylab= "Mean DNW (GAM)",
         xlab= "Biodiversity (0 - 1)")

plot.gam(gam1.2, select=3, ylim=c(-0.075, 0.125), seWithMean=T,scheme =1, all.terms=T, 
         residuals=T, cex=0.4, pch=20,
         rug=T, col="lightgrey", lwd=1,shade.col="lightblue",
         main="d",
         ylab= "Mean DNW (GAM)",
         xlab= "SST (oC)")

plot.gam(gam3.4, select=1, ylim=c(-0.075, 0.125), seWithMean=T,scheme =1, all.terms=T, 
         residuals=T, cex=0.4, pch=20,
         rug=T, col="lightgrey", lwd=1, shade.col="lightblue",
         main="e",
         ylab= "SD of mean DNW (GAM)",
         xlab= "SST (oC)")

plot.gam(gamTL1.2, select=2, ylim=c(-0.1, 0.15), seWithMean=T,scheme =1, all.terms=T, 
         residuals=T, cex=0.4, pch=20,
         rug=T, col="lightgrey", lwd=1, shade.col="lightblue",
         main="f",
         ylab= "Mean TL (GAM)",
         xlab= "Distance to shore (km)")





plot.gam(gam2.2, select=2, ylim=c(-0.075, 0.125), seWithMean=T,scheme =1, all.terms=T, 
         residuals=T, cex=0.4, pch=20,
         rug=T, col="lightgrey", lwd=1, shade.col="lightblue",
         main="b",
         ylab= "Median niche width (GAM)",
         xlab= "latitude")

plot.gam(gam2.2, select=3, ylim=c(-0.075, 0.125), seWithMean=T,scheme =1, all.terms=T, 
         residuals=T, cex=0.4, pch=20,
         rug=T, col="lightgrey", lwd=1,shade.col="lightblue",
         main="c",
         ylab= "Mean niche width (GAM)",
         xlab= "SST")

boxplot(niche.mean ~ lat, data=marine.sumry,
        main="b",
        ylab="Mean niche width (0 - 1)",
        xlab="Latitude (oN) of ceneterpoint of gridsquare",
        outline=F,
        col="lightblue")

plot(marine.sumry$bio, marine.sumry$niche.mean,
     pch=16, col="lightblue",
     ylim=c(0.05, 0.35),
     xlab="Biodiversity",
     ylab="Mean niche width")
lm.bio.niche<-lm(niche.mean ~ bio, data=marine.sumry)
abline(lm.bio.niche, col="red")
lmbio.niche.pred <- predict(lm.bio.niche, interval='confidence')
lines(marine.sumry$bio, lmbio.niche.pred[,2], lty=2, col='red')
lines(marine.sumry$bio, lmbio.niche.pred[,3], lty=2, col='red')


scatter.smooth(marine.sumry$bio, marine.sumry$niche.mean, family="gaussian",
               span=0.75, pch=16, col="lightblue",
               lpars=c(lwd=2, col="red"),
               xlab="Biodiversity",
               ylab="Mean niche width")

boxplot(niche.sd ~ lat, data=marine.sumry,
        main="b",
        ylab="Mean niche width (0 - 1)",
        xlab="Latitude (oN) of ceneterpoint of gridsquare",
        outline=F,
        col="lightblue")




Lat.SST<-vis.gam(gam1, view=c("lat", "SST"),
                 plot.type="persp", color='heat',
                 theta = 205, phi = 30,
                 main="c",
                 xlab="Latitude",
                 ylab="SST",
                 zlab="Median niche width")

