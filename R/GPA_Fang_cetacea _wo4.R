setwd("P:/My Documents/MANU/Hubei/Hupehsuchus_Filter/R")

library(geomorph)
library(abind)
library(StereoMorph)
library(Morpho)
library(calibrate)
source(".\\Utility.R")
#library(concaveman)

###
### Reading landmarks
###

chars <- paste0("F-",seq(1,9))
DAT <- read.csv("Fang_EA_2023_Hupehsuchus_Filter.csv")
#CET <- DAT[DAT$Clade1 == "Cetacea" | DAT$Clade1 == "Hupehsuchia",]
CET <- DAT[DAT$Clade1 == "Cetacea",]
 CET <- CET[CET$Genus != "Physeter",]
 CET <- CET[CET$Genus != "Kogia",]
 CET <- CET[CET$Species != "eueu",]
 rownames(CET) <- fnames <- paste0(CET$Genus,"_",CET$Species)
nlevels(as.factor(CET$Family))
nlevels(as.factor(CET$Genus))
nrow(CET)
cet <- CET[,sapply(CET, FUN = is.numeric)]

LMcoords <- arrayspecs(cet[,],9,2)
dimnames(LMcoords) <- list(chars,c("X","Y"),fnames)

Cb_9St <- read.csv("Cb_9St.csv",row.names=1)
LMcoords <- abind(LMcoords,Cb_9St[,c(3,1)],along=3)
Cb_9Bo <- read.csv("Cb_9Bo.csv",row.names=1)
LMcoords <- abind(LMcoords,Cb_9Bo[,c(3,1)],along=3)
Hn_9St <- read.csv("Hn_9St.csv",row.names=1)
LMcoords <- abind(LMcoords,Hn_9St[,c(3,1)],along=3)
Hn_9Bo <- read.csv("Hn_9Bo.csv",row.names=1)
LMcoords <- abind(LMcoords,Hn_9Bo[,c(3,1)],along=3)
dimnames(LMcoords) <- list(chars,c("X","Y"),c(fnames,"Cb_9St","Cb_9Bo","Hn_9St","Hn_9Bo"))

fam <- c(CET$Family,"Cb_9St","Cb_9Bo","Hn_9St","Hn_9Bo")
cl2 <- c(CET$Clade2,"Ichthyosauromorpha","Ichthyosauromorpha","Ichthyosauromorpha","Ichthyosauromorpha")

###
### GPA + PCA
###

##
## 2D
##

LMcoords2D <- -LMcoords
# LMcoords2D[,1,] <- -LMcoords2D[,1,]
# LMcoords2D[,2,] <- -LMcoords2D[,2,]
gpa2 <- gpagen(LMcoords2D)
shape.data <- gpa2$coords
pdf(h=4.5, w=4, file="AllLM_Fang_Cetacea_wo4.pdf")
  plotAllSpecimens(shape.data)
  for(i in 1:9) add.ellipse.one(shape.data[i,1,],shape.data[i,2,],series=i)
dev.off()

plotLD(shape.data,file="LM_Fang_Cetacea_wo4_asis.pdf")

meanLM <- mshape(shape.data)
pca.res <- gm.prcomp(shape.data)
summary(pca.res)
pca.coord <- pca.res$x
pdf(h=9,w=8,file="PCA_Fang_Cetacea_wo4.pdf")
layout(matrix(seq(1,4),2,2))
  chplot(-pca.coord[,1],pca.coord[,2],cl2,plotlegend=F,xlab="-PC 1",ylab="PC2")
  chplot(-pca.coord[,1],pca.coord[,2],cl2,plotlegend=T,xlab="-PC 1",ylab="PC2")
  chplot(-pca.coord[,1],pca.coord[,2],fam,plotlegend=F,xlab="-PC 1",ylab="PC2")
  chplot(-pca.coord[,1],pca.coord[,2],fam,plotlegend=T,xlab="-PC 1",ylab="PC2")
dev.off()
#
# Background tangential plot
# Code quoated and slightly modified from:
# https://aaronolsen.github.io/tutorials/morphometrics/backtransform.html
#

polyseq <- chars

lm_array <- LMcoords2D
# Get generalized Procrustes coordinates
gpa_array <- gpagen(lm_array)$coords

# Convert array to matrix for PCA
gpa_mat <- t(apply(gpa_array, 3, function(y) matrix(t(y), 1)))

# Perform non-phylogenetic PCA
resEig <- eigen(cov(gpa_mat))

# Get PC scores
scores <- gpa_mat %*% resEig$vectors

# Get percent variance explained along each axis
per_var <- (resEig$values / sum(resEig$values))*100

# Define function to draw shape
plot_beak_lateral <- function(xy, coor, size=1, col='black'){
  
  # If 3D, rotate points about x-axis using 3D rotation matrix
  if(ncol(coor) == 3){
    coor <- coor %*% matrix(c(1,0,0, 0,cos(-pi/2),sin(-pi/2), 
                              0,-sin(-pi/2),cos(-pi/2)), nrow=3, ncol=3)
  }
  
  # Get just x,y coordinates (orthographic projection into xy-plane)
  coor <- coor
  
  # Get plot aspect ratio
  w <- par('pin')[1]/diff(par('usr')[1:2])
  h <- par('pin')[2]/diff(par('usr')[3:4])
  asp <- w/h
  
  # Correct for plot aspect ratio not necessarily being 1:1
  coor[, 1] <- coor[, 1] * (1/asp)
  
  # Scale points and place back in position
  coor <- coor*size
  
  # Center about zero based on range of coordinates
  coor <- coor - matrix(colMeans(apply(coor, 2, range)), 
                        nrow=nrow(coor), ncol=ncol(coor), byrow=TRUE)
  
  # Move shape to PC score
  coor <- coor + matrix(xy, nrow(coor), ncol(coor), byrow=TRUE)
  
  # Set order in which to draw points to create polygon
  # polygon_order <- c(which(grepl('upper_bill_culmen', rownames(coor))),
  #                    rev(which(grepl('upper_bill_tomium_L', rownames(coor)))))
  
  polygon_order <- polyseq
  
  # Create filled polygon
  polygon((coor[polygon_order, ]), col=col, border=col)
}

# Set PCs to plot
pcs <- 1:2

# Open PDF graphics device
pdf('Bt_Fang_Cetacea_wo4.pdf', width=9, height=6.5)

# Create plot box with axes and axis labels
plot(scores[, pcs], type='n', main='Backtransform morphospace',
     xlab=paste0('PC', pcs[1], ' (', round(per_var[pcs[1]]), '%)'),
     ylab=paste0('PC', pcs[2], ' (', round(per_var[pcs[2]]), '%)'))

# Plot backtransform shapes
btShapes(scores=scores, vectors=resEig$vectors, fcn=plot_beak_lateral, 
         pcs=pcs, n=c(5,7), m=dim(lm_array)[2], row.names=dimnames(lm_array)[[1]], 
         pc.margin=c(0.06,0.05), size=0.1, col=gray(0.7))

# Plot points for each species
points(scores[, pcs])

# Add text labels
text(scores[, pcs], labels=substr(rownames(scores), 0, 10), cex=0.8, 
     pos=1, offset=0.3)

# Close the PDF graphics device
dev.off()

#
# End quotation
#

