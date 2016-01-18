library(geomorph)
library(car)
library(extrafont)
loadfonts()
#font_import()

#for ian
setwd("../") # 

#for nick
setwd("~/Dropbox/NickTesta (1)/SShD_9LM/")

#pulls in required scripts first, then runs dataread script (which reads in data)
source("scripts/WRP_FUNCTIONS.R")
source("scripts/NDT_functions.R")
source("scripts/dataread.R")

#to restore any alterations to par()
#par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)


#####################
#
# Analysis of data
#
#####################


##overall tests

#wings <- cbind(wings, CSresid)   #use this modification to account for allometry
full_model <- data.frame()
for (l in line.names)
{  
  temp <- subset(wings, LINE==l)
  shape <- two.d.array(gpagen(arrayspecs(temp[,6:23], 9,2), ShowPlot=F)$coords) #for regular, non-correcting-for-allometry model
  #shape <- two.d.array(gpagen(arrayspecs(temp[,27:44], 9,2), ShowPlot=F)$coords) #use when correcting for allometry
  sex <- temp[,5]
  size <- temp[,24]
  genotype <- temp[,4]
  genotype <- relevel(genotype, "w")
  background <- temp[,2]
  peas <- t(procD.lm(shape ~ size*genotype*sex*background, RRPP=TRUE)[1:15,7])
  temp.line <- cbind(l, peas)
  full_model <- rbind(full_model, temp.line)
}
colnames(full_model) <- c("line", "size", "genotype", "sex", "background", "size_genotype", "size_sex", "genotype_sex", "size_background", "genotype_background", "sex_background", "size_genotype_sex", "size_genotype_background", "size_sex_background", "genotype_sex_background", "size_genotype_sex_background")
for (i in 2:length(full_model))
{
  full_model[,i]<-as.numeric(as.character(full_model[,i]))
}

full_model_size <- data.frame()
for (l in line.names)
{  
  temp <- subset(wings, LINE==l)
  sex <- temp[,5]
  size <- temp[,24]
  genotype <- temp[,4]
  genotype <- relevel(genotype, "w")
  background <- temp[,2]
  peas <- t(summary(lm(size ~ genotype*sex*background))$coefficients[2:8,4])
  temp.line <- cbind(l, peas)
  full_model_size <- rbind(full_model_size, temp.line)
}
colnames(full_model_size) <- c("line", "genotype", "sex", "background", "genotype_sex", "genotype_background", "sex_background", "genotype_sex_background")
for (i in 2:length(full_model_size))
{
  full_model_size[,i]<-as.numeric(as.character(full_model_size[,i]))
}

#overall test for effects on control population... just to see
wt_model <- data.frame()
  temp <- subset(wings, GENOTYPE=='w')
  shape <- two.d.array(gpagen(arrayspecs(temp[,6:23], 9,2), ShowPlot=F)$coords)
  shape <- as.matrix(temp[,6:23])
  sex <- temp[,5]
  size <- temp[,24]
  background <- temp[,2]
  line <- temp[,3]
  peas <- t(procD.lm(shape ~ line*sex*background, RRPP=TRUE)[1:15,7])
  temp.line <- cbind(l, peas)
  full_model <- rbind(full_model, temp.line)
colnames(full_model) <- c("line", "size", "genotype", "sex", "background", "size_genotype", "size_sex", "genotype_sex", "size_background", "genotype_background", "sex_background", "size_genotype_sex", "size_genotype_background", "size_sex_background", "genotype_sex_background", "size_genotype_sex_background")
for (i in 2:length(full_model))
{
  full_model[,i]<-as.numeric(as.character(full_model[,i]))
}










#####################
#
# Figure 1a
#  Natural variation in backgrounds 
#
#####################

###########
#
# this block calculates values for figures 1 & 2
#
###########

shape.sig.values <-  data.frame()
for (b in background.names)
{
  for (l in line.names)
  {
    temp.wings <- subset(wings, BACKGROUND==b & LINE==l)
    shape <- two.d.array(gpagen(arrayspecs(temp.wings[,6:23], 9,2), ShowPlot=F)$coords)
    sex <- temp.wings[,5]
    size <- temp.wings[,24]
    genotype <- temp.wings[,4]
    genotype <- relevel(genotype, "w")
    background <- temp.wings[,2]
    
    sex.genotype.effect <- advanced.procD.lm(shape ~ sex + genotype, ~ sex*genotype, groups= ~genotype + sex)$anova.table[2,6]
    
    data.shape.line <- cbind(b, l, 
                             round(sex.genotype.effect,3))
    shape.sig.values <- rbind(shape.sig.values, data.shape.line)
  }
}

SSD.values <- data.frame()
for (b in background.names)
{
  for (l in line.names)
  {
    temp.wings <- subset(wings, BACKGROUND==b & LINE==l)
    
    SSD.wt.boot <- na.omit(replicate(1000, SSD_wt_func()))
    SSD.mutant.boot <- na.omit(replicate(1000, SSD_mutant_func()))
    
    dat_temp <- temp.wings
    size_lm <- lm(Centroid.Size ~ GENOTYPE*SEX, data=dat_temp)
    SSD.wt <- as.numeric((coef(size_lm)[1]/(coef(size_lm)[1] + coef(size_lm)[3])) -1)
    SSD.mutant <- as.numeric(((coef(size_lm)[1] + coef(size_lm)[2])/(sum(coef(size_lm)))) -1 )

    data.size.line <- cbind(b, l, 
                            round(SSD.wt, 3), 
                            round(as.numeric(quantile(SSD.wt.boot, probs=0.025)), 3), 
                            round(as.numeric(quantile(SSD.wt.boot, probs=0.975)), 3),  
                            round(SSD.mutant, 3), 
                            round(as.numeric(quantile(SSD.mutant.boot, probs=0.025)), 3), 
                            round(as.numeric(quantile(SSD.mutant.boot, probs=0.975)), 3) )
    SSD.values <- rbind(SSD.values, data.size.line)
  }
}
#fix column names
colnames(shape.sig.values) <- c("background", "line", "interaction")
colnames(SSD.values) <- c("background", "line", "wt", "wt.low", "wt.high", "mutant", "mutant.low", "mutant.high")

#coerce data into numeric format (from factor... why does it always do this!?)
for (i in 3:length(SSD.values))
{
  SSD.values[,i]<-as.numeric(as.character(SSD.values[,i]))
}
for (i in 3:length(shape.sig.values))
{
  shape.sig.values[,i] <- as.numeric(as.character(shape.sig.values[,i]))
}

#calculate pvalues for SSD
size.sig.values <- data.frame()
for (b in background.names)
{
  for (l in line.names)
  {
    temp.wings <- subset(wings, BACKGROUND==b & LINE==l)
    size.sig <- summary(lm(Csize ~ GENOTYPE*SEX, data=temp.wings))$coefficients[4,4]
    temp.line <- cbind(b, l, round(size.sig, 3))
    size.sig.values <- rbind( size.sig.values, temp.line)
  }
}
#fix colnames
colnames(size.sig.values) <- c("background", "line", "interaction")
#coerce data into proper format
for (i in 3:length(size.sig.values))
{
  size.sig.values[,i]<-as.numeric(as.character(size.sig.values[,i]))
}

#calculate procrustes distance between each treatment group
SShD.values <- data.frame()
for (b in background.names)
{
  for (l in line.names)
  {
    temp.master <- subset(master, backgr==b & line==l)
    
    #get SShD values
    SShD.wt.boot <- replicate(1000, SShD_wt_func())
    SShD.mutant.boot <- na.omit(replicate(1000, SShD_mutant_func()))
    
    dat_temp <- temp.master
    shape_PCs <- as.matrix(dat_temp[,6:19])
    PC_lm <- lm(shape_PCs ~ geno*sex, data=dat_temp)
    wt_f_shape <- coef(PC_lm)[1,]
    wt_m_shape <- coef(PC_lm)[1,] + coef(PC_lm)[3,]
    wt_SShD <- wt_f_shape - wt_m_shape
    SShD.wt <- PD(wt_SShD) # This corresponds to the distance you see in the geomorph distance for wild type SShD.
    
    p_f_shape <- coef(PC_lm)[1,] + coef(PC_lm)[2,]
    p_m_shape <- coef(PC_lm)[1,] + coef(PC_lm)[2,] + coef(PC_lm)[3,] + coef(PC_lm)[4,]
    mt_SShD <- p_f_shape - p_m_shape
    SShD.mutant <- PD(mt_SShD)
    
    data.line <- cbind(b, l, 
                       round(SShD.wt, 4), 
                       round(as.numeric(quantile(SShD.wt.boot, probs=0.025)), 4), 
                       round(as.numeric(quantile(SShD.wt.boot, probs=0.975)), 4),  
                       round(SShD.mutant, 4), 
                       round(as.numeric(quantile(SShD.mutant.boot, probs=0.025)), 4), 
                       round(as.numeric(quantile(SShD.mutant.boot, probs=0.975)), 4) )
    SShD.values <- rbind(SShD.values, data.line)
  }
}
colnames(SShD.values) <- c("background", "line", "mutant", "mutant.low", "mutant.high", "wt", "wt.low", "wt.high" )
for (i in 3:length(SShD.values))
{
  SShD.values[,i]<-as.numeric(as.character(SShD.values[,i]))
}

##########
#
# stuff just for figure 1
#
###########
wingsm <- aggregate( wings[,6:24], c( wings["LINE"], wings["SEX"], wings["BACKGROUND"], wings["GENOTYPE"] ), mean )
wingsm.Sam<-subset(wingsm, BACKGROUND=="Sam" & GENOTYPE=="w")
wingsm.Ore<-subset(wingsm, BACKGROUND=="Ore" & GENOTYPE=="w")
###

ssd.Sam <- rep( NA, nlevels( factor(wingsm$LINE) ) )
sshd.Sam <- rep( NA, nlevels( factor(wingsm$LINE) ) )
for (i in 1:nlevels( factor(wingsm.Sam$LINE) ) )
{
  l <- unique(SShD.values$line)[i]
  data1 <- wingsm.Sam[ wingsm.Sam$LINE == levels(factor(wingsm.Sam$LINE))[i], ]
  sam_lm <- lm(Centroid.Size ~ SEX, data=data1)
  ssd.Sam[i] <- subset(SSD.values, background=="Sam" & line==l)$wt
  sshd.Sam[i] <- subset(SShD.values, background=="Sam" & line==l)$wt
}
sexdim.Sam <- data.frame( levels(factor(wingsm$LINE)), ssd.Sam, sshd.Sam )
## 
ssd.Ore <- rep( NA, nlevels( factor(wingsm$LINE) ) )
sshd.Ore <- rep( NA, nlevels( factor(wingsm$LINE) ) )
for (i in 1:nlevels( factor(wingsm.Ore$LINE) ) )
{
  l <- unique(SShD.values$line)[i]
  data2 <- wingsm.Ore[ wingsm.Ore$LINE == levels(factor(wingsm.Ore$LINE))[i], ]
  ore_lm <- lm(Centroid.Size ~ SEX, data=data2)
  ssd.Ore[i] <- subset(SSD.values, background=="Ore" & line==l)$wt
  sshd.Ore[i] <- subset(SShD.values, background=="Ore" & line==l)$wt  
}
sexdim.Ore <- data.frame( levels(factor(wingsm$LINE)), ssd.Ore, sshd.Ore)

##
sshd.Sam<-na.omit(sshd.Sam)
ssd.Sam<- na.omit(ssd.Sam )
sshd.Ore<-na.omit(sshd.Ore)
ssd.Ore<- na.omit(ssd.Ore )

####
x<-as.matrix(cbind(sexdim.Sam$ssd, sexdim.Sam$sshd))
x2<-na.omit(as.matrix(cbind(sexdim.Ore$ssd, sexdim.Ore$sshd)))
lhist<-20
num.dnorm<-100
## histogram (for barplot-ting the density)
xhist  <- hist(x[,1],  plot=FALSE, breaks=seq(from= min(x[,1]), to=max(x2[,1]),  length.out=lhist))
xhist2 <- hist(x2[,1], plot=FALSE, breaks=seq(from= min(x[,1]), to=max(x2[,1]), length.out=lhist))  

yhist  <- hist(x[,2],  plot=FALSE, breaks=seq(from=min(x[,2]), to=max(x2[,2]), length.out=lhist)) # note: this uses probability=TRUE
yhist2 <- hist(x2[,2], plot=FALSE, breaks=seq(from=min(x[,2]), to=max(x2[,2]), length.out=lhist))
## determine the plot range and all the stuff needed for the barplots and lines
xx <-  seq(min(x[,1]), max(x2[,1]),  length.out=num.dnorm) 
xx2 <- seq(min(x[,1]), max(x2[,1]),  length.out=num.dnorm) 
xy <-  dnorm(xx,  mean=mean(x[,1]),  sd=sd(x[,1])) # density points
xy2 <- dnorm(xx2, mean=mean(x2[,1]), sd=sd(x2[,1])) # density points

yx <-  seq(min(x[,2]), max(x2[,2]), length.out=num.dnorm)
yx2 <- seq(min(x[,2]), max(x2[,2]), length.out=num.dnorm)  

yy <-  dnorm(yx, mean=mean(x[,2]), sd=sd(x[,2]))
yy2 <- dnorm(yx2, mean=mean(x2[,2]), sd=sd(x2[,2]))
## barplot and line for x (top)

#############
#
# plot figure 1B
#
############

pdf("Fig.1a Natural Variation.pdf", family="Helvetica", width=11, height=9.5)

#this section adjusts layout to accommodate multiple plots
layMat <- matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
layout(layMat, widths=c(5/7, 2/7), heights=c(2/7, 5/7))
ospc <- 3 # outer space
pext <- 5 # par extension down and to the left
bspc <- 0 # space between scatter plot and bar plots
par. <- par(mar=c(pext, pext, bspc, bspc),
            oma=rep(ospc, 4)) # plot parameters

#initial plot of variation. x[,2] is sam shape, x[,1] is sam size
plot(x[,2] ~ x[,1], xlim=c(min(x[,1]), max(x2[,1])), ylim=c(min(x[,2]), max(x2[,2])), xlab="SSD", ylab="SShD", pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.5)#, ...)  #cex.lab= for text size

#x2[,2] is ore shape, x2[,1] is ore size
points(x2[,2] ~ x2[,1], col="black", pch=21, bg="gray", cex=1.5)
legend("topleft", c("Sam", "Ore"), pt.bg=c("black", "gray"), pch=c(21,21), bty='n', cex=1.5)

me.sam.shape <- qt(0.95, length(x[,2])-1) * (sd(x[,2])/sqrt(length(x[,2])))
me.sam.size <-  qt(0.95, length(x[,1])-1) * (sd(x[,1])/sqrt(length(x[,1])))
lines(c(mean(x[,1]), mean(x[,1])), c(mean(x[,2])-me.sam.shape, mean(x[,2])+me.sam.shape), lwd=2)
lines(c(mean(x[,1])-me.sam.size, mean(x[,1])+me.sam.size), c(mean(x[,2]), mean(x[,2])), lwd=2)
points(mean(x[,2]) ~ mean(x[,1]), pch=22, cex=2, col="black", bg="white")

me.ore.shape <- qt(0.95, length(x2[,2])-1) * (sd(x2[,2])/sqrt(length(x2[,2])))
me.ore.size <-  qt(0.95, length(x2[,1])-1) * (sd(x2[,1])/sqrt(length(x2[,1])))
lines(c(mean(x2[,1]), mean(x2[,1])), c(mean(x2[,2])-me.ore.shape, mean(x2[,2])+me.ore.shape), col="gray", lwd=2)
lines( c(mean(x2[,1])-me.ore.size, mean(x2[,1])+me.ore.size)  , c(mean(x2[,2]), mean(x2[,2])), col="gray", lwd=2)
points(mean(x2[,2]) ~ mean(x2[,1]), pch=22, col="gray", cex=2, bg="white")

#mtext("Natural variation in control populations", side=3, line=12, adj=-0.25, cex=1.5)

par(mar=c(0, pext, 0, 0))
barplot(xhist$density,  axes=FALSE, ylim=c(0, max(xhist2$density, xy2)), col="black", space=0) # barplot
barplot(xhist2$density, axes=FALSE, ylim=c(0, max(xhist2$density, xy2)), space=0, col="lightgray", add=TRUE)
lines(seq(from=0, to=lhist-1, length.out=num.dnorm), xy, col="black", lwd=2) # line
lines(seq(from=0, to=lhist-1, length.out=num.dnorm), xy2, col="darkgray", lwd=2) # line

## barplot and line for y (right)
par(mar=c(pext, 0, 0, 0))
barplot(yhist$density,  axes=FALSE, col="black",     xlim=c(0, max(yhist$density, yy)), 
        space=0, horiz=TRUE) # barplot
barplot(yhist2$density, axes=FALSE, col="lightgray", xlim=c(0, max(yhist$density, yy2)), 
        space=0, horiz=TRUE, add=TRUE) # barplot
lines(yy, seq(from=0, to=lhist-1, length.out=num.dnorm), col="black", lwd=2) # line
lines(yy2, seq(from=0, to=lhist-1, length.out=num.dnorm), col="darkgray", lwd=2) # line


dev.off()
####

par(mfrow=c(1,1))


#####################
#
# Figure 1b
#
####################
###show vector of change on actual wing landmarks
w.M.ore.omb<-gpagen(arrayspecs(subset(wings1, BACKGROUND=="Ore" & GENOTYPE=="w" & SEX=="M" & LINE=="3045")[,6:23], 9, 2), ShowPlot=F)
w.F.ore.omb<-gpagen(arrayspecs(subset(wings1, BACKGROUND=="Ore" & GENOTYPE=="w" & SEX=="F" & LINE=="3045")[,6:23], 9, 2), ShowPlot=F)
w.M.sam.omb<-gpagen(arrayspecs(subset(wings1, BACKGROUND=="Sam" & GENOTYPE=="w" & SEX=="M" & LINE=="3045")[,6:23], 9,2), ShowPlot=F)
w.F.sam.omb<-gpagen(arrayspecs(subset(wings1, BACKGROUND=="Sam" & GENOTYPE=="w" & SEX=="F" & LINE=="3045")[,6:23], 9,2), ShowPlot=F)

par(mfrow=c(1,2))
plotRefToTarget(w.F.sam.omb$coords[,,1], w.M.sam.omb$coords[,,1], method="vector", mag=5, links=wing.link)
mtext("Sam dimorphism")
plotRefToTarget(w.F.ore.omb$coords[,,1], w.F.ore.omb$coords[,,1], method="vector", mag=5, links=wing.link)
mtext("Ore dimorphism")


#####################
#
# Figure 2
#  Mutant lines SSD vs SShD
#
#####################


###########
#
# figure 2 plots
#
##########


pdf("Fig.2 SSD vs SShD.pdf", family="Helvetica", width=11, height=8.5)

layMat <- matrix(c(1,2,3,4), ncol=2, byrow=TRUE)
layout(layMat, widths=c(1/2, 1/2))
ospc <- 3 # outer space
pext <- 1#4 # par extension down and to the left
bspc <- 4 # space between scatter plot and bar plots
par. <- par(mar=c(pext, pext, bspc, bspc),
            oma=rep(ospc, 4)) # plot parameters

for (p in 1:2)
{
for (b2 in background.names)
{
	SSD <- subset(SSD.values, background==b2)
	SShD <- subset(SShD.values, background==b2)
	SShD.pvals <- subset(shape.sig.values, background==b2)
	SSD.pvals <- subset(size.sig.values, background==b2)
	
  if (p ==1)
  {
    
  
  if (b2 == "Ore")
  {
    plot(0,0, col="white", 
             xlim=c(min(SSD.values$mutant, SSD.values$wt), 
                    max(SSD.values$mutant, SSD.values$wt) ), 
             ylim=c(min(SShD.values$mutant, SShD.values$wt),
                    max(SShD.values$mutant, SShD.values$wt)+0.001),  
             xlab="SSD", cex.lab=1.5, cex.axis=1.5,
             ylab="SShD" )
    mtext(paste("Ore"), 3, cex=1.5)
    mtext("SSD", 1, line=2.3, cex=1.33)
    mtext("SShD", 2, line=2.3, cex=1.33)
    legend("topleft", c("wt", "mutant", "SSD", "SShD", "SSD + SShD"), pch=c(21,21,21,21,21), col=c("black", "black", "red", "blue", "darkviolet"), pt.bg=c("black", "gray", "red", "blue", "darkviolet"), bty='n')
    lines(c(0.075, 0.075), c(0.012, 0.025), lwd=1)
    lines(c(0.2, 0.2), c(0.012, 0.025), lwd=1)
    lines(c(0.075, 0.2), c(0.012, 0.012), lwd=1)
    lines(c(0.075, 0.2), c(0.025, 0.025), lwd=1)
    par(xpd=NA)
    lines(c(0.075, 0.0475), c(0.025, -0.0058), lty=2)
    lines(c(0.2, 0.2225), c(0.025, -0.0058), lty=2)
    lines(c(0.075, 0.0675), c(0.012, -0.0058), lty=2)
    lines(c(0.2, 0.205), c(0.012, -0.0058), lty=2)
    par(xpd=F)
  }
  else
  {
    plot(0,0, col="white", 
         xlim=c(min(SSD.values$mutant, SSD.values$wt), 
                max(SSD.values$mutant, SSD.values$wt) ), 
         ylim=c(min(SShD.values$mutant, SShD.values$wt),
                max(SShD.values$mutant, SShD.values$wt)+0.001),  
         xlab="SSD", cex.lab=1.5, cex.axis=1.5,
         ylab="SShD",
         yaxt='n')
    mtext(paste("Sam"), 3, cex=1.5) 
    
    #mtext(paste("SSD vs SShD in:"), 3, cex=1, line=2, adj=-0.4, font=2)  
    lines(c(0.075, 0.075), c(0.012, 0.025), lwd=1)
    lines(c(0.2, 0.2), c(0.012, 0.025), lwd=1)
    lines(c(0.075, 0.2), c(0.012, 0.012), lwd=1)
    lines(c(0.075, 0.2), c(0.025, 0.025), lwd=1)
    par(xpd=NA)
    mtext("SSD", 1, line=2.3, cex=1.33)
    lines(c(0.075, 0.0475), c(0.025, -0.0058), lty=2)
    lines(c(0.2, 0.2225), c(0.025, -0.0058), lty=2)
    lines(c(0.075, 0.0675), c(0.012, -0.0058), lty=2)
    lines(c(0.2, 0.205), c(0.012, -0.0058), lty=2)
    par(xpd=F)
  }
  }
  else
  {
    if (b2 == "Ore")
    {
      plot(0,0, col="white", 
           xlim=c(0.075, 
                  0.2), 
           ylim=c(0.012,
                  0.025),  
           xlab="SSD", 
           ylab="SShD" )
      #mtext(paste("Ore"), 3, cex=1)
      #legend()
    }
    else
    {
      plot(0,0, col="white", 
           xlim=c(0.075, 
                  0.2), 
           ylim=c(0.012,
                  0.025),  
           xlab="SSD", 
           ylab="SShD",
           yaxt='n')
    }
  }
	
	for (i in 1:length(SSD[,1]))
	{
		if (SShD.pvals$interaction[i] > 0.05 & SSD.pvals$interaction[i] > 0.05) #neither
		{
			lines( c(SSD$mutant[i], SSD$wt[i]), c(SShD$mutant[i], SShD$wt[i]), col="gray")
			points( SShD$wt[i] 	   ~ SSD$wt[i], col="gray", pch=19)
			points( SShD$mutant[i] ~ SSD$mutant[i], pch=21, col="gray")
		}
	}	
	for (i in 1:length(SSD[,1]))
	{		
		if (SSD.pvals$interaction[i] < 0.05 & SShD.pvals$interaction[i] > 0.05) #size
		{
			lines( c(SSD$wt.low[i], SSD$wt.high[i]), 
				   c(SShD$wt[i],    SShD$wt[i]),
				   col="pink")
			lines( c(SSD$wt[i],       SSD$wt[i]), 
				   c(SShD$wt.high[i], SShD$wt.low[i]),
				   col="pink")
			lines( c(SSD$mutant.low[i], SSD$mutant.high[i]), 
				   c(SShD$mutant[i],    SShD$mutant[i]),
				   col="pink")
			lines( c(SSD$mutant[i],       SSD$mutant[i]), 
				   c(SShD$mutant.high[i], SShD$mutant.low[i]),
				   col="pink")
		}
		
		if (SShD.pvals$interaction[i] < 0.05 & SSD.pvals$interaction[i] > 0.05) #shape
		{
			lines( c(SSD$wt.low[i], SSD$wt.high[i]), 
				   c(SShD$wt[i],    SShD$wt[i]),
				   col="lightblue")
			lines( c(SSD$wt[i],       SSD$wt[i]), 
				   c(SShD$wt.high[i], SShD$wt.low[i]),
				   col="lightblue")
			lines( c(SSD$mutant.low[i], SSD$mutant.high[i]), 
				   c(SShD$mutant[i],    SShD$mutant[i]),
				   col="lightblue")
			lines( c(SSD$mutant[i],       SSD$mutant[i]), 
				   c(SShD$mutant.high[i], SShD$mutant.low[i]),
				   col="lightblue")
		}
		
		if (SSD.pvals$interaction[i] < 0.05 & SShD.pvals$interaction[i] < 0.05) #both
		{				
			lines( c(SSD$wt.low[i], SSD$wt.high[i]), 
				   c(SShD$wt[i],    SShD$wt[i]),
				   col="plum2")
			lines( c(SSD$wt[i],       SSD$wt[i]), 
				   c(SShD$wt.high[i], SShD$wt.low[i]),
				   col="plum2")
			lines( c(SSD$mutant.low[i], SSD$mutant.high[i]), 
				   c(SShD$mutant[i],    SShD$mutant[i]),
				   col="plum2")
			lines( c(SSD$mutant[i],       SSD$mutant[i]), 
				   c(SShD$mutant.high[i], SShD$mutant.low[i]),
				   col="plum2")
			lines( c(SSD$mutant[i], SSD$wt[i]), 
			       c(SShD$mutant[i], SShD$wt[i]), 
			       col="darkviolet")
		}
	}
	for (i in 1:length(SSD[,1]))
	{	
		if (SSD.pvals$interaction[i] < 0.05 & SShD.pvals$interaction[i] > 0.05) #size
		{
			lines( c(SSD$mutant[i], SSD$wt[i]), c(SShD$mutant[i], SShD$wt[i]), col="red", lwd=2)
			points( SShD$wt[i] 	   ~ SSD$wt[i], col="red", pch=19)
			points( SShD$mutant[i] ~ SSD$mutant[i], pch=21, col="red", bg="pink")
			name <- stuff[i]
      text(   SSD$mutant[i], SShD$mutant[i], name, font=3, cex=0.75, pos=3, col="red")
		}
		
		if (SShD.pvals$interaction[i] < 0.05 & SSD.pvals$interaction[i] > 0.05) #shape
		{
			lines( c(SSD$mutant[i], SSD$wt[i]), c(SShD$mutant[i], SShD$wt[i]), col="blue", lwd=2)
			points( SShD$wt[i] 	   ~ SSD$wt[i], col="blue", pch=19)
			points( SShD$mutant[i] ~ SSD$mutant[i], pch=21, col="blue", bg="lightblue")
      name <- stuff[i]
			text(   SSD$mutant[i], SShD$mutant[i], name, font=3, cex=0.75, pos=3, col="blue")
		}
		
		if (SSD.pvals$interaction[i] < 0.05 & SShD.pvals$interaction[i] < 0.05) #both
		{				
			lines( c(SSD$mutant[i], SSD$wt[i]), 
			       c(SShD$mutant[i], SShD$wt[i]), 
			       col="darkviolet", lwd=2)
			points( SShD$wt[i] 	   ~ SSD$wt[i], col="darkviolet", pch=19)
			points( SShD$mutant[i] ~ SSD$mutant[i], pch=21, col="darkviolet", bg="plum2")
      name <- stuff[i]
			text(   SSD$mutant[i], SShD$mutant[i], name, font=3, cex=0.75, pos=3, col="darkviolet")
		}
	}
}
}
dev.off()


#####################
#
# Figure 3
#  vector correlations
#
#####################

sex<-wings[,5]
size<- wings[,24]
genotype <-wings[,4]
genotype <- relevel(genotype, "w")
background <- wings[,2]
line <- wings$LINE
shape <- two.d.array(gpagen(arrayspecs(wings[,6:23], 9,2), ShowPlot=F)$coords)
shape_PCs <- prcomp(shape, center=TRUE, scale.=FALSE, retx=TRUE)$x[,1:(ncol(shape)-4)]
shape_PCs <- wingMaster[,6:19]
data_together <- data.frame(shape_PCs, sex, genotype, size, background)

Shape.Vec.wt <- data.frame()

dat_together <- subset(data_together, genotype=="w")
SShD.Vec <- na.omit(t(replicate(1000, boot_SShD_VC_wt_test())))
Vec.cor <- mean(SShD.Vec[,1])
Vec.angle <- mean(SShD.Vec[,2])
Vec.cor.low <- quantile(SShD.Vec[,1], probs=0.025)
Vec.cor.high <- quantile(SShD.Vec[,1], probs=0.975)
Vec.angle.low <- quantile(SShD.Vec[,2], probs=0.025)  
Vec.angle.high <- quantile(SShD.Vec[,2], probs=0.975)
temp.line <- cbind(b, l, Vec.cor, Vec.angle, Vec.cor.low, Vec.cor.high, Vec.angle.low, Vec.angle.high)

Shape.Vec.wt <- rbind(Shape.Vec.wt, temp.line)
for (i in 3:length(Shape.Vec))
{
  Shape.Vec[,i] <- as.numeric(as.character(Shape.Vec[,i]))
}

sex<-wings[,5]
size<- wings[,24]
genotype <-wings[,4]
genotype <- relevel(genotype, "w")
background <- wings[,2]
line <- wings$LINE
shape <- two.d.array(gpagen(arrayspecs(wings[,6:23], 9,2), ShowPlot=F)$coords)
shape_PCs <- prcomp(shape, center=TRUE, scale.=FALSE, retx=TRUE)$x[,1:(ncol(shape)-4)]
shape_PCs <- wingMaster[,6:19]
data_together <- data.frame(shape_PCs, sex, genotype, size, background)

Shape.Vec <- data.frame()
for (b in background.names)
{
  for (l in line.names)
  {
    dat_together <- subset(data_together, background==b & line==l)
    SShD.Vec <- na.omit(t(replicate(1000, boot_SShD_VC_test())))
    Vec.cor <- mean(SShD.Vec[,1])
    Vec.angle <- mean(SShD.Vec[,2])
    Vec.cor.low <- quantile(SShD.Vec[,1], probs=0.025)
    Vec.cor.high <- quantile(SShD.Vec[,1], probs=0.975)
    Vec.angle.low <- quantile(SShD.Vec[,2], probs=0.025)
    Vec.angle.high <- quantile(SShD.Vec[,2], probs=0.975)
    temp.line <- cbind(b, l, Vec.cor, Vec.angle, Vec.cor.low, Vec.cor.high, Vec.angle.low, Vec.angle.high)
    Shape.Vec <- rbind(Shape.Vec, temp.line)
  }
}
for (i in 3:length(Shape.Vec))
{
  Shape.Vec[,i] <- as.numeric(as.character(Shape.Vec[,i]))
}


######
#
# plot figure 3
# 
######


par(mfrow=c(1,1))
pdf("Fig.3 Vector Correlations.pdf", family="Helvetica", width=14, height=7)

plot(0,0, col="white", xlim=c(1, 42), ylim=c(0,1), xlab="", ylab="Vector Correlation", xaxt='n', cex.axis=1.5, cex.lab=1.5)
axis(1, at=c(1:42), labels=stuff, las=2, cex.axis=1.2, font=3)

for (i in 1:42)
{
  j<- i+42
  lines(c((i+0.33),(i+0.33)), c(Shape.Vec$Vec.cor.low[j], Shape.Vec$Vec.cor.high[j]), lwd=1.2)
  points( (i+0.33), Shape.Vec$Vec.cor[j], pch=19, cex=1.2)
  lines(c(i,i), c(Shape.Vec$Vec.cor.low[i], Shape.Vec$Vec.cor.high[i]), lwd=1.2, col="gray")
  points(i, Shape.Vec$Vec.cor[i], pch=21, bg="gray", cex=1.2)
  
}

legend("bottomleft", c("Sam", "Ore"), pch=c(21,21), col=c("black", "black"), pt.bg=c("black", "gray"), bty='n', cex=1.5)

lines(c(0,50), c(0.825,0.825), col="white")

dev.off()





#####################
#
# Figure 4
#  allometry 
#
#####################

Allometry_pvals <- data.frame()
for (b in background.names)
{
  for (l in line.names)
  {
    temp_dat <- subset(wingMaster, backgr==b & line==l)
    shape <- temp_dat[,6:19] 
    size <- log(temp_dat[,1])
    sex <- temp_dat[,2]
    genotype <- temp_dat[,5]
    Allometry_analysis <- advanced.procD.lm(shape 
    ~ genotype + sex + size + genotype:size + size:sex + genotype:sex + genotype:sex:size, 
    ~ genotype + sex + size + genotype:size + size:sex + genotype:sex, 
      groups= ~sex + genotype, 
      slope = ~size)

    all_pval <- Allometry_analysis$anova.table[2,6]
    mutant_pval <- Allometry_analysis$Prob.Slopes.dist[3,1]
    wt_pval <- Allometry_analysis$Prob.Slopes.dist[4,2]
    m_pval <- Allometry_analysis$Prob.Slopes.dist[4,3]
    f_pval <- Allometry_analysis$Prob.Slopes.dist[2,1]
    temp.line <- cbind(b, l, round(all_pval, 3), round(mutant_pval, 3), round(wt_pval, 3), round(m_pval, 3), round(f_pval, 3))
    Allometry_pvals <- rbind(Allometry_pvals, temp.line)
    
  }
}
colnames(Allometry_pvals) <- c("background", "line", "allP", "mutant", "wt", "male", "female")
for (i in 3:length(Allometry_pvals))
{
  Allometry_pvals[,i]<-as.numeric(as.character(Allometry_pvals[,i]))
}



background.allometry <- data.frame()
g <- "w"
for (b in background.names)
{
  for (s in sex.names)
  {
      temp.master <- subset(wingMaster, backgr==b & geno==g & sex==s)
      slopes <- replicate(1000, allometry_func())
     
      dat_temp <- temp.master
      shape <- as.matrix(dat_temp[,6:19])
      size <- log(dat_temp[,1])   #csize
      sex <- dat_temp[,2]
      genotype <- dat_temp[,5]
      
      allometry.model <- lm(shape ~ size)
      Beta <- allometry.model$coefficients[2,] 
      allometry.shape <- ShapeScore(t(Beta), shape)
      shape <- allometry.shape
      allometry <- lm(shape ~ size)
      slope <- coef(allometry)[[2]]
      
      #slope <- mean(slopes)
      slope.low <- quantile(slopes, probs=0.025)
      slope.high <- quantile(slopes, probs=0.975)
      temp.line <- cbind(b, s, g, slope, slope.low, slope.high)
      background.allometry <- rbind(background.allometry, temp.line)
    
  }
}
#fixes errors of R making numbers into factors
for (i in 3:length(background.allometry))
{
  background.allometry[,i]<-as.numeric(as.character(background.allometry[,i]))
}
colnames(background.allometry) <- c("background", "sex", "slope", "slope.low", "slope.high")


Allometry <-data.frame()
for (i in 1:length(line.names)) #length of line.names
{
  for (s in sex.names) 
  {
    for (b in background.names) 
    {
      for (g in genotype.names)
      {    
        l<-line.names[i]
        
        temp.master <- subset(wingMaster, backgr==b & line==l & geno==g & sex==s)
        
        dat_temp <- temp.master
        shape <- as.matrix(dat_temp[,6:19])
        size <- log(dat_temp[,1])   #csize
        sex <- dat_temp[,2]
        genotype <- dat_temp[,5]

        allometry.model <- lm(shape ~ size)
        Beta <- allometry.model$coefficients[2,] 
        allometry.shape <- ShapeScore(t(Beta), shape)
        shape <- allometry.shape
        allometry <- lm(shape ~ size)
        slope <- coef(allometry)[[2]]
    
        slopes <- replicate(1000, allometry_func())
        slope.low <- quantile(slopes, probs=0.025)
        slope.high <- quantile(slopes, probs=0.975)
        temp.line <- cbind(i, l, s, b, g, slope, slope.low, slope.high)
        Allometry <- rbind(Allometry, temp.line)
      }
    }
  }
}
#fixes errors of R making numbers into factors
for (i in 6:length(Allometry))
{
  Allometry[,i]<-as.numeric(as.character(Allometry[,i]))
}


Allometry_size_log <- Allometry
Allometry_size_log_pvals <- Allometry_pvals




#######
#
# plot figure 4
#
#######

par(xpd=FALSE)

pdf("Fig.4 Allometry.pdf", family="Helvetica", width=14, height=7)

#this section adjusts layout to accommodate multiple plots
layMat <- matrix(c(1,2), ncol=1, byrow=TRUE)
layout(layMat, heights=c(1/2, 1/2))
ospc <- 4 # outer space
pext <- 1#4 # par extension down and to the left
bspc <- 0 # space between scatter plot and bar plots
par. <- par(mar=c(pext, pext, bspc, bspc),
            oma=rep(ospc, 4)) # plot parameters

k<-0
for (b2 in background.names)
{
  j<-1
  for (l2 in line.names) 
  {
    for (g2 in genotype.names) 
    {
      for (s2 in sex.names) 
      {    
        #make temporary subset based on position in loop
        wings.temp <- subset(Allometry, l==l2 & b==b2 & s==s2 & g==g2)
        wings.temp[,6:8] <- wings.temp[,6:8]
        pval.temp  <- subset(Allometry_pvals, line==l2 & background==b2)
        #on very first itteration, make the plot
        if (j==1)
        {
          k<-k+1
          #par(mgp=c(2,0.5,0))
          #par(mar=c(6.5,4,3,2))
          plot(0, 0, col="white", 
               xlim=c(1,42), ylim=c(0.05,0.3), 
               xlab="", ylab="", xaxt='n', cex.axis=1.25)
          mtext(b2, side=4, cex=1.5, line=0.5)
        }
        
        
        if (s2=="F" & g2=="p")
        {
          if (pval.temp$allP < 0.05)
          {
            lines(c(j,j), c(wings.temp[,7], wings.temp[,8]), col="red", lty=1)
            points(wings.temp$slope ~ j, col="red", pch=21, bg="pink")
          }
          else
          {
            lines(c(j,j), c(wings.temp[,7], wings.temp[,8]), col="pink", lty=1)
            points(wings.temp$slope ~ j, col="pink", pch=21, bg="white")
          }
          j<-j+0.4
        }
        
        if (s2=="F" & g2=="w")
        {
          if (pval.temp$allP < 0.05)
          {
            lines(c(j,j), c(wings.temp[,7], wings.temp[,8]), col="red", lty=1)
            points(wings.temp$slope ~ j, col="red", pch=19)
          }
          else
          {          
            lines(c(j,j), c(wings.temp[,7], wings.temp[,8]), col="pink", lty=1)
            points(wings.temp$slope ~ j, col="pink", pch=19)
          }          
          j<-j+0.2
        }
        
        if (s2=="M" & g2=="p")
        {
          if (pval.temp$allP < 0.05)
          {
            lines(c(j,j), c(wings.temp[,7], wings.temp[,8]), col="blue", lty=1)
            points(wings.temp$slope ~ j, col="blue", pch=21, bg="lightblue")
          }
          else
          { 
            lines(c(j,j), c(wings.temp[,7], wings.temp[,8]), col="lightblue", lty=1)
            points(wings.temp$slope ~ j, col="lightblue", pch=21, bg="white")
          }
          
          j<-j+0.2
        }
        
        if (s2=="M" & g2=="w")
        {
          if (pval.temp$allP < 0.05)
          {
            lines(c(j,j), c(wings.temp[,7], wings.temp[,8]), col="blue", lty=1)
            points(wings.temp$slope ~ j, col="blue", pch=19)
          }
          else
          {
            lines(c(j,j), c(wings.temp[,7], wings.temp[,8]), col="lightblue", lty=1)
            points(wings.temp$slope ~ j, col="lightblue", pch=19)
          }
          j<-j+0.2
        }
        
      }
    }
  }
  par(xpd=TRUE)
  if (k ==1)
  {
#    mtext("Allometry by sex and genotype", cex=1)
    legend("topleft", c("wt Male", "wt Female", "mutant Male","mutant Female"), pch=c(21, 21, 21, 21), col=c("blue","red", "blue","red"), pt.bg=c("blue", "red", "lightblue", "pink"), cex=1, bty='n')
  }
  else 
  {
#    mtext("Gene of interest", side=1, line=5, cex=0.75)
    axis(1, at=c(1:42), labels=stuff[1:42], las=2, cex.axis=1.2, font=3)
    mtext("Slope of Allometry", side=2, line=2.5, adj=0.6, cex=1.5, font=3)
  }
  par(xpd=FALSE)
}
dev.off()




#####################################

