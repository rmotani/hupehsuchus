setwd("P:/My Documents/MANU/Hubei/Hupehsuchus_Filter/R")
library(geomorph)
library(abind)
source(".\\Utility.R")

#DAT <- read.csv("Fang_EA_2023_Hupehsuchus_Filter.csv")
DAT <- read.csv("Fang_EA_2023.csv")
rownames(DAT) <- fnames <- paste0(DAT$Genus,"_",DAT$Species)
dat <- DAT[,sapply(DAT,FUN=is.numeric)]
LMcoords <- arrayspecs(dat,9,2)
chars <- paste0("F-",seq(1,9))
dimnames(LMcoords) <- list(chars,c("X","Y"),fnames)
plotLD(LMcoords,file="LM_Fang__asis.pdf")



###
### Species list from Fang et al. (2023)
###

CET <- DAT[DAT$Clade1 == "Cetacea",]
#CET <- CET[CET$Genus != "Physeter",]
#CET <- CET[CET$Genus != "Kogia",]
rownames(CET) <- fnames <- paste0(CET$Genus,"_",CET$Species)
nlevels(as.factor(CET$Family))
nlevels(as.factor(CET$Genus))
nrow(CET)



###
### Reading taxonomy for 3D data
###

taxa <- read.csv("Cetacea.csv")
taxa <- taxa[,-ncol(taxa)]
rownames(taxa) <- tnames <- paste0(taxa$Genus,"_",taxa$Species)
shared <- intersect(tnames,fnames)
nlevels(as.factor(taxa$Family))
nlevels(as.factor(taxa$Genus))
nrow(taxa)
setdiff(taxa$Genus,CET$Genus)
setdiff(CET$Genus,taxa$Genus)

###
### Reading 3D landmarks --9
###

fs.9 <- list.files("../3D_Cetacea_Principal/@o_fcsv",pattern=".fcsv")
#chars.9 <-seq(1,9)
chars.9 <- paste0("F-",seq(1,9))

LMcoords.9 <- NULL
for(i in 1:length(fs.9)){
  tmpdat <- read.csv(paste("../3D_Cetacea_Principal/@o_fcsv/",fs.9[i],sep=""),skip=3,header=F)
  rownames(tmpdat) <- tmpdat$V12
  tmpdat2 <- tmpdat[chars.9,2:4]
  LMcoords.9 <- abind(LMcoords.9,tmpdat2,along=3)
}
dimnames(LMcoords.9) <- list(chars.9,c("X","Y","Z"),substr(fs.9,3,nchar(fs.9)-5))
LMcoords.9 <- LMcoords.9[,,shared]
LMcoords.92D <- LMcoords.9[,c(3,1),]
gpa2.9 <- gpagen(LMcoords.92D)
shape.data.9 <- gpa2.9$coords


###
### Reading 3D landmarks --15
###
fs.15 <- list.files("../3D_Cetacea_Principal/@fcsv",pattern=".fcsv")
#chars.15 <-c("1","2a","2b","3a","3b","3c","4","5c","5b","5a","6a","6b","7","8","9")
chars.15 <- paste0("F-",seq(1,15))

LMcoords.15 <- NULL
for(i in 1:length(fs.15)){
  tmpdat <- read.csv(paste("../3D_Cetacea_Principal/@fcsv/",fs.15[i],sep=""),skip=3,header=F)
  rownames(tmpdat) <- tmpdat$V12
  tmpdat2 <- tmpdat[chars.15,2:4]
  LMcoords.15 <- abind(LMcoords.15,tmpdat2,along=3)
}
dimnames(LMcoords.15) <- list(chars.15,c("X","Y","Z"),substr(fs.15,1,nchar(fs.15)-5))
LMcoords.15 <- LMcoords.15[,,shared]
LMcoords.152D <- LMcoords.15[,c(3,1),]
gpa2.15 <- gpagen(LMcoords.152D)
shape.data.15 <- gpa2.15$coords

###
### Landmarks from Fang et al. (2023)
###
cet <- CET[,sapply(CET, FUN = is.numeric)]
FangLM <- arrayspecs(cet,9,2)
dimnames(FangLM) <- list(chars,c("X","Y"),fnames)
FangLM <- FangLM[,,shared]

###
### Plotting side by side
###

pdf(h=7,w=8,file="CompareLM_Shared.pdf")
layout(t(matrix(seq(1,6),3,2)))
ntax <- dim(shape.data.9)[3]
  if(is.na(ntax)) ntax <-1
  nchar <- dim(shape.data.9)[1]
  naxes <- dim(shape.data.9)[2]
  chars4print <- seq(1,9)
  chars4print.15 <-c("1","2a","2b","3a","3b","3c","4","5c","5b","5a","6a","6b","7","8","9")
  if(naxes > 2) shape.data.9 <- shape.data.9[,1:2,]
  psequence <- c(seq(1,nchar-1),nchar,nchar-1,1)
  psequence.15 <- c(seq(1,14),15,14,1)
  
  for(i in 1: ntax){
    shape.tmp.15 <- shape.data.15[,,i]
      shape.tmp.15 <- shape.tmp.15[,c(2,1)]
      shape.tmp.15[,1] <- -shape.tmp.15[,1] 
      plot(shape.tmp.15,type="n",asp=1,main=dimnames(shape.data.15)[[3]][i])
      polygon(shape.tmp.15[psequence.15,],border="black",col="gray")
      points(shape.tmp.15,type="p",pch=19,cex=2)
      textxy(shape.tmp.15[,1],shape.tmp.15[,2],chars4print.15,cex=2)
    shape.tmp.9 <- shape.data.9[,,i]
      shape.tmp.9 <- shape.tmp.9[,c(2,1)]
      shape.tmp.9[,1] <- -shape.tmp.9[,1] 
      plot(shape.tmp.9,type="n",asp=1,main=dimnames(shape.data.9)[[3]][i])
      polygon(shape.tmp.9[psequence,],border="black",col="gray")
      points(shape.tmp.9,type="p",pch=19,cex=2)
      textxy(shape.tmp.9[,1],shape.tmp.9[,2],chars4print,cex=2)
    Fang.tmp <- FangLM[,,i]
      plot(Fang.tmp,type="n",asp=1,main=dimnames(shape.data.9)[[3]][i])
      polygon(Fang.tmp[psequence,],border="black",col="gray")
      points(Fang.tmp,type="p",pch=19,cex=2)
      textxy(Fang.tmp[,1],Fang.tmp[,2],chars4print,cex=2)
  }
dev.off()
