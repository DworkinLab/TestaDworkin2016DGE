#setwd("~/Documents/Documents/Lab/SSD/Side Projects/Dworkin-Gibson Reanalysis")

# Note from Ian - This overwrites your previous use of the setwd() in the other file
#setwd("~/Dropbox/NickTesta (1)/SShD_9LM/") # terrible name for folder (1)?? 


# Note from Ian - but you call geomorph from the main script as well. What are you trying to do? It only needs to be called the once.
library(geomorph)

#################################################
#################################################
##       read in and define data               ##
##                                             ##   

#wings.TPS<-readland.tps("~/Drive/Documents/Lab/SSD/Side Projects/Dworkin-Gibson Reanalysis/Dworkin_Gibson_2006_LM.TPS")


metadata<-read.csv("data/metadata.csv")
wings1<-read.csv("data/Dworkin_Gibson_2006_PC.csv")
wings1<-subset(wings1, LINE!="10911" & LINE!="11843" & LINE!="14953" & LINE!="15028" & LINE!="10478" & LINE!="10514" & LINE!="12862" & LINE!="14986")
wings2<-wings1[,-c(1:5, 24:25)]
wings<-arrayspecs(wings2, 9, 2) #turns data into 3D array

##run Generalized Procrustes Analysis/plot out raw data
gpa.wings<-gpagen(wings, ShowPlot=F)

##turns data back to 2D
wings2d<-two.d.array(gpa.wings$coords)

#defining variables
sex<-wings1[,5]
geno<-wings1[,4]
csize<-wings1[,24]
backgr<-wings1[,2]
line<-wings1[,3]
gps<-as.factor(paste(sex, backgr, geno, line))

###figure out which landmarks are which
wing.link<-matrix(c(6,5,5,4,5,3,3,2,4,2,2,1,5,7,3,8,2,9),ncol=2,byrow=T)

master<-data.frame(csize, sex, backgr, line, geno)
A<-gpa.wings$coords
Y<-two.d.array(A)
#Y<-gpa.both$coords[][][1:3705]
#shape.model.1 <- lm( prcomp(Y)$x[,1:14] ~ csize * sex * geno)

master <- data.frame( master, prcomp( Y )$x[,1:14] )

######################
#
# Defining variables that I will need later
#
######################

#loop variables
line.names<-unique(wings1$LINE)
genotype.names<- c("w","p")
background.names<- c("Ore", "Sam")
sex.names<-c("M","F")
background<-c("Ore", "Sam")
genotype<-c("w", "p")

#data variables

wings1$Csize<- wings1$Centroid.Size #this just makes centroid size easier to type/reference
wingsm <- aggregate( wings1[,6:24], c( wings1["LINE"], wings1["SEX"], wings1["BACKGROUND"], wings1["GENOTYPE"] ), mean )


wings<- wings1 #commented out for allometry removed analysis... uncomment when finished.
CSresid <- lm( as.matrix( wings1[,6:23] ) ~ wings1$Centroid.Size )$resid
#wings[,6:23] <- wings[,27:44] #use this line to replace normal wing values with allometry removed values. Comment this out and rerun when finished.
wingsN <- aggregate( wings1[,1], c( wings1["LINE"], wings1["SEX"], wings1["BACKGROUND"], wings1["GENOTYPE"] ), length )  #this yields sample sizes for each treatment

PCmultiplied <- master[,6:19] * 10000
colnames(PCmultiplied) <- sub("PC", "PCm", colnames(PCmultiplied))
wingMaster <- data.frame(master, PCmultiplied)
wings$GENOTYPE <- relevel(wings$GENOTYPE, "w")

#stuff and junk for correlating gene names to lines for accurate labeling 
junk<-vector()
stuff<-metadata$gene[1]
stuff[1]<-metadata$gene[1]
gene.reference<-unique(wings1$LINE)
for (i in 1:42)
{
  for (j in 1:50)
  {
    if (gene.reference[i] == metadata$stock[j])
    {
      stuff[i]<-metadata$gene[j]
      junk[i]<- metadata$stock[j]
    }
  }
}
stuff
junk



