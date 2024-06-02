setwd("P:/My Documents/Github/hupehsuchus/R")


library(geomorph)
library(abind)
library(StereoMorph)
library(Morpho)
library(calibrate)
source(".\\Utility.R")

###
### Reading landmarks
###

fs <- list.files("../3D_Cetacea_Principal/@o_fcsv",pattern=".fcsv")
chars <- paste0("F-",seq(1,9))

LMcoords <- NULL
for(i in 1:length(fs)){
  tmpdat <- read.csv(paste("../3D_Cetacea_Principal/@o_fcsv/",fs[i],sep=""),skip=3,header=F)
  rownames(tmpdat) <- tmpdat$V12
  tmpdat2 <- tmpdat[chars,2:4]
  LMcoords <- abind(LMcoords,tmpdat2,along=3)
}
dimnames(LMcoords) <- list(chars,c("X","Y","Z"),substr(fs,3,nchar(fs)-5))

## Adding Ichthyosauromorpha
Cb_9St <- read.csv("Cb_9St.csv",row.names=1)
LMcoords <- abind(LMcoords,Cb_9St,along=3)
Cb_9Bo <- read.csv("Cb_9Bo.csv",row.names=1)
LMcoords <- abind(LMcoords,Cb_9Bo,along=3)
Hn_9St <- read.csv("Hn_9St.csv",row.names=1)
LMcoords <- abind(LMcoords,Hn_9St,along=3)
Hn_9Bo <- read.csv("Hn_9Bo.csv",row.names=1)
LMcoords <- abind(LMcoords,Hn_9Bo,along=3)
dimnames(LMcoords) <- list(chars,c("X","Y","Z"),c(substr(fs,3,nchar(fs)-5),"Cb_9St","Cb_9Bo","Hn_9St","Hn_9Bo"))

#LMcoords <- LMcoords[,,-10]

###
### Reading taxonomy
###

taxa <- read.csv("Cetacea.csv")
taxa <- taxa[,-ncol(taxa)]
rownames(taxa) <- paste0(taxa$Genus,"_",taxa$Species)
tnew <- taxa[dimnames(LMcoords)[[3]],]

## Adding Ichthyosauromorpha
tnew$Clade1[66:69] <- c("Ichthyosauromorpha","Ichthyosauromorpha","Ichthyosauromorpha","Ichthyosauromorpha")
tnew$Clade2[66:69] <- c("Ichthyosauromorpha","Ichthyosauromorpha","Ichthyosauromorpha","Ichthyosauromorpha")
tnew$Family[66:69] <- c("Chaohusauridae","Chaohusauridae","Hupehsuchidae","Hupehsuchidae")
tnew$Genus[66:69] <- c("Chaohusaurus","Chaohusaurus","Hupehsuchus","Hupehsuchus")
tnew$Species[66:69] <- c("St","Bo","St","Bo")
cbind(dimnames(LMcoords)[[3]],rownames(tnew),sort(rownames(taxa)))
dimnames(LMcoords) <- list(chars,c("X","Y","Z"),paste0(substr(tnew$Genus,1,5),"_",substr(tnew$Species,1,5)))


###
### GPA + PCA
###

##
## 3D
##

gpa1 <- gpagen(LMcoords)
shape.data <- gpa1$coords
plotAllSpecimens(shape.data)

meanLM <- mshape(shape.data)
pca.res <- gm.prcomp(shape.data)
summary(pca.res)
pca.coord <- pca.res$x
pdf(h=9,w=8,file="PCA_9_3D.pdf")
layout(matrix(seq(1,4),2,2))
  chplot(pca.coord[,1],pca.coord[,2],tnew$Clade2,plotlegend=F,xlab="PC 1",ylab="PC2")
  chplot(pca.coord[,1],pca.coord[,2],tnew$Clade2,plotlegend=T,xlab="PC 1",ylab="PC2")
  chplot(pca.coord[,1],pca.coord[,2],tnew$Family,plotlegend=F,xlab="PC 1",ylab="PC2")
  chplot(pca.coord[,1],pca.coord[,2],tnew$Family,plotlegend=T,xlab="PC 1",ylab="PC2")
dev.off()

x11(h=8.5,w=11)
plotRefToTarget(meanLM,pca.res$shapes$shapes.comp1$min)
pdf("Cetacea_9_PCA_TPS_3D.pdf",h=8.5,w=11)
  layout(matrix(seq(1,4),2,2))
  plotRefToTarget(meanLM,pca.res$shapes$shapes.comp1$min)
  plotRefToTarget(meanLM,pca.res$shapes$shapes.comp1$max)
  plotRefToTarget(meanLM,pca.res$shapes$shapes.comp1$min)
  plotRefToTarget(meanLM,pca.res$shapes$shapes.comp2$max)
dev.off()


### Saving deformed meshes
# MeshIn <- readOBJ("../3D_Cetacea_Principal/Orcinus_orca_AM_M_22839.obj")
# plot3d(MeshIn,type="shade",col="lightgreen",asp=F)
# MeshOutMean_9 <- warpRefMesh(MeshIn,LMcoords[,,"Orcinus_orca"],meanLM)
# mesh2obj(MeshOutMean_9,"Orcinus_orca_meanLM_9.obj")
# 
# MeshOutPC1min_9 <- plotRefToTarget(meanLM,pca.res$shapes$shapes.comp1$min,method="surface",mesh=MeshOut_9,mesh.coord=meanLM)
# mesh2obj(MeshOutPC1min_9,filename="Orcinus_orca_C1min_9.obj")
# MeshOutPC1max_9 <- plotRefToTarget(meanLM,pca.res$shapes$shapes.comp1$max,method="surface",mesh=MeshOut_9,mesh.coord=meanLM)
# mesh2obj(MeshOutPC1max_9,filename="Orcinus_orca_C1max_9.obj")
# MeshOutPC2min_9 <- plotRefToTarget(meanLM,pca.res$shapes$shapes.comp2$min,method="surface",mesh=MeshOut_9,mesh.coord=meanLM)
# mesh2obj(MeshOutPC2min_9,filename="Orcinus_orca_C2min_9.obj")
# MeshOutPC2max_9 <- plotRefToTarget(meanLM,pca.res$shapes$shapes.comp2$max,method="surface",mesh=MeshOut_9,mesh.coord=meanLM)
# mesh2obj(MeshOutPC2max_9,filename="Orcinus_orca_C2max_9.obj")


#
# Background tangential plot
# Code quoated and slightly modified from:
# https://aaronolsen.github.io/tutorials/morphometrics/backtransform.html
#

polyseq <- chars

lm_array <- LMcoords
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
  coor <- coor[, c(1,3)]
  
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
pdf('Bt_9_3D.pdf', width=9, height=6.5)

# Create plot box with axes and axis labels
plot(scores[, pcs], type='n', main='Backtransform morphospace',
     xlab=paste0('PC', pcs[1], ' (', round(per_var[pcs[1]]), '%)'),
     ylab=paste0('PC', pcs[2], ' (', round(per_var[pcs[2]]), '%)'))

# Plot backtransform shapes
btShapes(scores=scores, vectors=resEig$vectors, fcn=plot_beak_lateral, 
         pcs=pcs, n=c(5,7), m=dim(lm_array)[2], row.names=dimnames(lm_array)[[1]], 
         pc.margin=c(0.06,0.05), size=0.11, col=gray(0.7))

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


##
## 2D
##

LMcoords2D <- LMcoords[,c(3,1),]
gpa2 <- gpagen(LMcoords2D)
shape.data <- gpa2$coords
pdf(h=4.5, w=4, file="AllLM_9.pdf")
  plotAllSpecimens(shape.data)
  for(i in 1:9) add.ellipse.one(shape.data[i,1,],shape.data[i,2,],series=i)
dev.off()

plotLD(shape.data,file="LM_9.pdf")


meanLM <- mshape(shape.data)
pca.res <- gm.prcomp(shape.data)
summary(pca.res)
pca.coord <- pca.res$x
pdf(h=9,w=8,file="PCA_9_2D.pdf")
layout(matrix(seq(1,4),2,2))
  chplot(pca.coord[,1],pca.coord[,2],tnew$Clade2,plotlegend=F,xlab="PC 1",ylab="PC2")
  chplot(pca.coord[,1],pca.coord[,2],tnew$Clade2,plotlegend=T,xlab="PC 1",ylab="PC2")
  chplot(pca.coord[,1],pca.coord[,2],tnew$Family,plotlegend=F,xlab="PC 1",ylab="PC2")
  chplot(pca.coord[,1],pca.coord[,2],tnew$Family,plotlegend=T,xlab="PC 1",ylab="PC2")
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
pdf('Bt_9_2D.pdf', width=9, height=6.5)

# Create plot box with axes and axis labels
plot(scores[, pcs], type='n', main='Backtransform morphospace',
     xlab=paste0('PC', pcs[1], ' (', round(per_var[pcs[1]]), '%)'),
     ylab=paste0('PC', pcs[2], ' (', round(per_var[pcs[2]]), '%)'))

# Plot backtransform shapes
btShapes(scores=scores, vectors=resEig$vectors, fcn=plot_beak_lateral, 
         pcs=pcs, n=c(5,7), m=dim(lm_array)[2], row.names=dimnames(lm_array)[[1]], 
         pc.margin=c(0.06,0.05), size=0.08, col=gray(0.7))

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
