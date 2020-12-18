options(stringsAsFactors=FALSE)

## moving average
ma <- function(x,n=10,circular=FALSE)
    stats::filter(x, rep(1/n, n), sides = 2, circular = circular)
## diameter -> volume
d2v <- function(x) (x/2)^3*pi *4/3


expid <- "20201127_MD_plasmidserie_CoilHack"
PATH <- "/mnt/synmibi/Studierende/DATA/CASY"
if ( Sys.info()["nodename"]=="intron" )
    PATH <- "/data/synmibi/CASY"
in.path <- file.path(PATH,expid)
out.path <- file.path("~/work/CoilHack/experiments", expid)


## experiment parameters
expid <- "+rha" #"-rha" #  
dil <- 3000

## analysis&plot parameters
normalize <- TRUE
max.cnt <- 15.83e6 # max cell count for color scheme breaks
## size filter
min.counts <- 1.5
max.counts <- 5
min.norm <- min.counts
max.norm <- max.counts

## TODO
## *calculate median and peak sizes/volumes
## *calculate total cell volume


for ( expid in c("+rha","-rha") ) {
sizes <- counts <- matrix(NA, nrow=1024, ncol=9)
for ( i in 0:8 ) {
    file.name <- paste0("20201127_plasmid_",i,expid,".TXT")
    if ( i==0 ) file.name <- sub(expid,"",file.name,fixed=TRUE)
    data <- read.delim(file.path(in.path, file.name),header=FALSE)
    
    from <- which(data[,1]=="Size Channel")+1
    to <- which(data[,1]=="Counts Repeat 1")-1
    dat <- data[from:to,]

    counts[,i+1] <- as.numeric(dat[,2])*dil
    sizes[,i+1] <- as.numeric(gsub(" ","",dat[,1]))
    #plot(dat, type="l")
}

## check that all sizes are the same
unique(apply(sizes,1,function(x) length(unique(x))))
size <- sizes[,1]

## size filters
size.counts <- size>=min.counts & size<=max.counts
size.norm   <- size>=min.norm & size<=max.norm

## add half of difference to duplicate x (size) values
dups <- which(duplicated(size))
size[dups] <- size[dups] + median(diff(size))/2

## total cell number in size range W/O debris
idx <- size.counts

## total cell number, cells/mL
total <- apply(counts[idx,],2,sum)
## total cell volume, uL/mL - TODO: numbers seem too high?
volume <- apply(counts[idx,],2, function(x) sum(x*d2v(size[idx]))/1e9)

## plot to size range W/O debris
idx <- size.norm

## colors and breaks
cols <- viridis::viridis(100)
#cols <- grey.colors(100,start=1,end=0)
brks <- seq(0,max.cnt,length.out=101) #1

## normalize counts
counts.nrm <- counts[idx,] # apply(counts[idx,],2,function(x) x/max(x))
if ( normalize ) {
    counts.nrm <-apply(counts[idx,],2,function(x) x/max(x))
    brks <- brks/max(brks)
}

## moving average of counts
counts.nrm <- apply(counts.nrm,2,ma)


png(file.path(out.path,paste0("20201127_casy_",expid,".png")),
    width=400, height=200)
par(mai=c(.5,1,.1,.5),mgp=c(1.3,.4,0),tcl=-.25,xaxs="i",yaxs="i")
image(y=d2v(size[idx]),x=0:8,z=t(counts.nrm), col=cols,breaks=brks,
      ylab=expression("cell volume, "*fL), xlab="sample",ylim=c(0,30))
axis(1, at=0:8)
axis(2)
box()
par(new=TRUE)
plot(0:8, total/1e8, type="p",col=2, axes=FALSE,xlab=NA,ylab=NA,ylim=c(0,11),
     xlim=par("usr")[1:2],pch=19)
lines(0:8, total/1e8,col=2)
axis(4,col=2,col.axis=2)
mtext("1e8 cells/mL",4, par("mgp")[1],col=2)
par(new=TRUE)
plot(0:8, volume, type="p",col="white",
     axes=FALSE,xlab=NA,ylab=NA,ylim=c(0,6.5),
     xlim=par("usr")[1:2],pch=3,lwd=2)
lines(0:8, volume,col="white")
axis(2,col="black",col.axis="black",line=par("mgp")[1]*2)
mtext(expression("total cell volume, "*mu*L/mL),2, 3*par("mgp")[1],col="black")
legend("top", expid, text.col="white", text.font=2,cex=1.2,bty="n",x.intersp=0)
legend("topleft",c("cell count","total cell volume"),pch=c(19,3),
       col=c("red","white"),bty="n",text.col="white",pt.lwd=c(1,2))
dev.off()

sample.cols <- rev(viridis::viridis(ncol(counts)))

png(file.path(out.path,paste0("20201127_casy_raw_",expid,".png")),
    width=400, height=200)
par(mai=c(.5,.5,.1,.1),mgp=c(1.3,.4,0),tcl=-.25,xaxs="i",yaxs="i")
matplot(d2v(size), counts,type="l",lty=1,xlim=c(0,30),col=sample.cols,xlab=expression("cell volume, "*fL),ylim=c(0,1.75e7))
legend("topright", paste(0:8), title="Samples:",col=sample.cols, lty=1,y.intersp=.7)
legend("top", expid, text.col="black", text.font=2, cex=1.2,bty="n",x.intersp = 0)
dev.off()
}
