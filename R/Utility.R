#setwd("P:/My Documents/Github/hupehsuchus/R")

library(calibrate)
library(scales)

pointcolor <- rep(c("royalblue4","salmon4","seagreen4","gold4","plum4","steelblue4","turquoise4","mistyrose4","lightgoldenrod4","palevioletred4"),3)
fillcolor <- scales::alpha(pointcolor,0.3) 
bgcolor <- "white"
pointsymbol <- rep(c(15,17,19,4,3,8,0,1,2,5,6,7),3)
#pointsymbol <-rep(c(21:25,3,4,8),5)


# 2D plot of landmarks in 3D array from geomorph

plotLD <- function(shape,file="landmark.pdf",w=7,h=11,nw=2,nh=3){
  ntax <- dim(shape)[3]
  if(is.na(ntax)) ntax <-1
  nchar <- dim(shape)[1]
  naxes <- dim(shape)[2]
  chars <- dimnames(shape)[[1]]
  if(naxes > 2) shape <- shape[,1:2,]
  psequence <- c(seq(1,nchar-1),nchar,nchar-1,1)
  pdf(w=w,h=h,file=file)
  layout(matrix(seq(1,nw*nh),nh,nw))
  for(i in 1: ntax){
    shape.tmp <- shape[,,i]
    plot(shape.tmp,type="n",asp=1,main=dimnames(shape)[[3]][i])
    polygon(shape.tmp[psequence,],border="black",col="gray")
    points(shape.tmp,type="p",pch=19,cex=2)
    textxy(shape.tmp[,1],shape.tmp[,2],chars,cex=2)
  }
  dev.off()
}


### Add a group of points from X and Y data
# add.point.one <- function(x,y,series=2,colorfill=TRUE,fcolor=fillcolor,
#                           pcolor=pointcolor,psymbol=pointsymbol,cex=1)
# {
#   points(x,y,pch=psymbol[series],col=pcolor[series],cex=cex,bg=pcolor[series])
# }

### Add a confidence ellipse to X and Y data
add.ellipse.one <- function(x,y,series=2,pval=0.95,num=30,colorfill=TRUE,lwd=2,
                          pcolor=pointcolor,fcolor=fillcolor)
{
  acc <- num
  trans <- 1-pval
  vx <- var(x,na.rm=T)
  vy <- var(y,na.rm=T)
  vxy <- var(x,y,na.rm=T)
  lambda <- eigen(var(na.omit(cbind(x,y))))$values
  a <- sqrt(vxy^2/((lambda[2]-vx)^2+vxy^2))
  b <- (lambda[2]-vx)*a/vxy
  theta <- atan(a/b)
  k <- sqrt(-2*log(trans))
  l1 <- sqrt(lambda[1])*k
  l2 <- sqrt(lambda[2])*k
  pvec <- 0:num
  x2right <- sin((pi*pvec)/(num*2))*l1
  x2 <- c(-rev(x2right),x2right )
  tmp <- 1-x2^2/l1^2
  y2 <- l2*sqrt(ifelse(tmp < 0,0,tmp))
  x2 <- c(x2,rev(x2))
  y2 <- c(y2,-rev(y2))
  s0 <- sin(theta)
  c0 <- cos(theta)
  xx <- c0*x2+s0*y2+mean(x,na.rm=T)
  yy <- -s0*x2+c0*y2+mean(y,na.rm=T)
  if(colorfill) {polygon(xx,yy,border=pcolor[series],col=fillcolor[series],lwd=lwd)}else{
    polygon(xx,yy,border=pcolor[series],lwd=lwd)}
  epp <- cbind(xx,yy)
  return(epp)
}

### Add a convex hull to X and Y data
add.ch.one <- function(x,y,series=2,xlim=range(x),ylim=range(y),main="",
                       xlab="x",ylab="y",colorfill=TRUE,fcolor=fillcolor,pcolor=pointcolor,
                       lwdv=lwd,trans=0.5){
  ch.list <- chull(x,y)
  ch.points <- cbind(x[ch.list],y[ch.list])
  if(colorfill) {polygon(ch.points,border=pcolor[series],col=scales::alpha(fcolor[series],trans),lwd=lwdv)
  }else{polygon(ch.points,border=pcolor[series],lwd=lwdv)}
  invisible(ch.points)
}

### Add a bag hull to X and Y data
# add.bag.one <- function(x,y,series=2,colorfill=TRUE,fcolor=fillcolor,
#                         pcolor=pointcolor,lwd=lwd,trans=0.5,bagfactor=3)
# {
#   bag.one <- compute.bagplot(x,y,bagfactor)
#   print("1")
#   bag.points <- bag.one$hull.loop
#   print(bag.points)
#   if(colorfill) {polygon(bag.points,border=pcolor[series],col=scales::alpha(fcolor[series],trans),lwd=lwd)
#   }else{polygon(bag.points,border=pcolor[series],lwd=lwd)}
#   print("3")
#   invisible(bag.points)
# }

### Biplot with convex hulls
chplot <- function(x,y,grp=rep("g1",length(x)),lwd=2,lwd.axis=1.5,
                   cex=1,xlab="X",ylab="Y",main="",plotlegend=T,
                   plothull=T,plotellipse=F,legend.pos="topleft"){
  grp <- as.factor(grp)
  gnames <- levels(grp)
  gnums <- table(grp)
  plot(x,y,type="n",main=main,xlab=xlab,ylab=ylab)
  for(i in 1:nlevels(grp)){
    X <- x[grp==gnames[i]]
    Y <- y[grp==gnames[i]]
    if(gnums[i]>2) if(plotellipse) add.ellipse.one(X,Y,series=i,lwd=lwd)
    if(gnums[i]>1) if(plothull) add.ch.one(X,Y,series=i,lwd=lwd)
  }
  for(i in 1:nlevels(grp)){
    X <- x[grp==gnames[i]]
    Y <- y[grp==gnames[i]]
    points(X,Y,pch=pointsymbol[i],col=pointcolor[i])
  }
  if(plotlegend) legend(legend.pos,gnames,pch=pointsymbol,col=pointcolor)
}

# chplot(pca.coord[,1],pca.coord[,2],tnew$Family)
# chplot(pca.coord[,1],pca.coord[,2],tnew$Family,plothull=F,plotellipse=T)

# add.cipi.one <- function(x,y,series=2,pval=0.95,num=30,sfig=4,colorfill=TRUE,
#                        lwd=2,orgn=FALSE,pcolor=pointcolor,fcolor=fillcolor,olscipi=c("CIPI","CI","PI"),
#                        X=x,trans=0.5)
# {
#   require(stats)
#   olscipi=match.arg(olscipi)
#   new <- data.frame(x = seq(min(c(x,X),na.rm=T),max(c(x,X),na.rm=T),by=dist(range(c(x,X),na.rm=T))/num))
#   if(orgn){
#     pred.w.plim <- predict(lm(y ~ x-1),new,interval="prediction",level=pval)[,-1]
#     pred.w.clim <- predict(lm(y ~ x-1),new,interval="confidence",level=pval)[,-1]
#   }else{
#     pred.w.plim <- predict(lm(y ~ x),new,interval="prediction",level=pval)[,-1]
#     pred.w.clim <- predict(lm(y ~ x),new,interval="confidence",level=pval)[,-1]
#   }
#   cipipnts <- cbind(new$x,pred.w.clim,pred.w.plim)
#   pipnts <- cbind(c(cipipnts[,1],rev(cipipnts[,1])),c(cipipnts[,4],rev(cipipnts[,5])))
#   cipnts <- cbind(c(cipipnts[,1],rev(cipipnts[,1])),c(cipipnts[,2],rev(cipipnts[,3])))
#   if(olscipi=="CIPI"){
#     if(colorfill) {polygon(pipnts[,1],pipnts[,2],col=scales::alpha(fcolor[series],trans),border=NA)}
#     matplot(new$x,cbind(pred.w.clim,pred.w.plim[,]),lty=1,lwd=c(lwd/2,lwd/2,lwd/2,lwd/2),type="l",ylab="predicted y",add=TRUE,col=pcolor[series])
#   }else if(olscipi=="CI"){
#     if(colorfill) {polygon(cipnts[,1],cipnts[,2],col=scales::alpha(fcolor[series],trans),border=NA)}
#     matplot(new$x,pred.w.clim,lty=1,lwd=c(lwd/2,lwd/2),type="l",ylab="predicted y",add=TRUE,col=pcolor[series])
#   }else{
#     if(colorfill) {polygon(pipnts[,1],pipnts[,2],col=scales::alpha(fcolor[series],trans),border=NA)}
#     matplot(new$x,pred.w.plim,lty=1,lwd=c(lwd/2,lwd/2),type="l",ylab="predicted y",add=TRUE,col=pcolor[series])
#   }
#   return(cipipnts)
# }
