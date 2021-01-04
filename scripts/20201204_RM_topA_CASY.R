options(stringsAsFactors=FALSE)

## moving average
ma <- function(x,n=10,circular=FALSE)
    stats::filter(x, rep(1/n, n), sides = 2, circular = circular)
## diameter -> volume
d2v <- function(x) (x/2)^3*pi *4/3


expid <- "20201204_RM_topA"
PATH <- "/mnt/synmibi/Studierende/DATA/CASY/"
if ( Sys.info()["nodename"]=="intron" )
    PATH <- "/data/synmibi/CASY"

in.path <- file.path(PATH,expid)
out.path <- file.path("~/work/CoilHack/experiments/reactor/pcc6803/",
                      expid,"analysis")

cat(paste("PARSING CASY DATA:", date(), "\n"))

## SAMPLE NOTES AND SETTINGS

## SAMPLE_27: original sample clogged CASY, cells were washed,
## OD decreased from 0,64 to 0,60 -> use this OD change as a factor
## to correct cell counts?
## 

## experiment parameters
dil <- 4000
undil <- 1000 # "UNDILUTED"

## analysis&plot parameters
normalize <- TRUE
max.cnt <- 10e6 # max cell count for color scheme breaks
## size filter
min.counts <- 1.5
max.counts <- 5
min.norm <- min.counts
max.norm <- max.counts

## TODO
## *calculate median and peak sizes/volumes
## *calculate total cell volume

files <- list.files(path=in.path, pattern="measurement_.*\\.TXT$")

sizes <- counts <- matrix(NA, nrow=1024, ncol=length(files))
sampleIDs <- rep(NA, length(files))
for ( i in seq_along(files) ) {
    file.name <- files[i]
    data <- try(read.delim(file.path(in.path, file.name),header=FALSE))

    if ( inherits(data, 'try-error') ) next

    ## comment contains sample ID
    comment <- toupper(data[which(data[,1]=="Comment 1"),2])

    cat(paste("parsing", comment, "\n"))
    sampleIDs[i] <- comment

    DIL <- dil
    if ( length(grep("UNDILUTED",sampleIDs[i]))>0 )
        DIL <- undil
    
    ## number of cycles: total cell counts need to be divided by this
    ## cell count correction factor
    ## TODO: urgently check how Cycles, Sample Volume and Volume Correction
    ##       are interpreted.
    cycles <- as.numeric(data[which(data[,1]=="Cycles"),2])
    volume <- as.numeric(data[which(data[,1]=="Sample Volume (\xb5l)"),2])
    volcor <- as.numeric(data[which(data[,1]=="Volume Correction"),2])

    corr <- DIL /(cycles*volume*volcor/1000)

    from <- which(data[,1]=="Size Channel")+1
    to <- which(data[,1]=="Counts Repeat 1")-1
    dat <- data[from:to,]

    counts[,i] <- as.numeric(dat[,2]) * corr
    sizes[,i] <- as.numeric(gsub(" ","",dat[,1]))
    #plot(dat, type="l")
}

## filter data
filter <- grep("SAMPLE",sampleIDs)
#filter <- filter[grep("_[A-Z]", sampleIDs[filter], invert=TRUE)]
counts <- counts[,filter]
sizes <- sizes[,filter]
sampleIDs <- sampleIDs[filter]

sampleLabels <- sub("_[A-Z].*","",sub(" .*","",sub("SAMPLE_","",sampleIDs)))

## check that all sizes are the same
if ( unique(apply(sizes,1,function(x) length(unique(x))))!=1 )
    stop("different size vectors!")
size <- sizes[,1]

## filters
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
    counts.all <-apply(counts,2,function(x) x/max(x))
}

## moving average of counts
counts.nrm <- apply(counts.nrm,2,ma)
counts.all <- apply(counts.all,2,ma)


png(file.path(out.path,paste0(expid,"_CASY.png")),
    width=2*3.5, height=2*3.5, units="in", res=300)
par(mfcol=c(2,1),mai=c(.5,1,.1,.5),mgp=c(1.3,.4,0),tcl=-.25,xaxs="i",yaxs="i")
image(y=d2v(size[idx]),x=1:ncol(counts),z=t(counts.nrm), col=cols,breaks=brks,
      ylab=expression("cell volume, "*fL), xlab="",ylim=c(0,35),
      axes=FALSE)
axis(1, at=1:ncol(counts), label=sampleLabels,las=2)
axis(2)
box()
par(new=TRUE)
plot(1:ncol(counts), total/1e8, type="p",col=2,
     axes=FALSE,xlab=NA,ylab=NA,ylim=c(0,8),
     xlim=par("usr")[1:2],pch=19)
lines(1:ncol(counts), total/1e8,col=2)
axis(4,col=2,col.axis=2)
mtext("1e8 cells/mL",4, par("mgp")[1],col=2)
par(new=TRUE)
plot(1:ncol(counts), volume, type="p",col="white",
     axes=FALSE,xlab=NA,ylab=NA,ylim=c(0,4),
     xlim=par("usr")[1:2],pch=3,lwd=2)
lines(1:ncol(counts), volume,col="white")
axis(2,col="black",col.axis="black",line=par("mgp")[1]*2)
mtext(expression("total cell volume, "*mu*L/mL),2, 3*par("mgp")[1],col="black")
legend("topleft",c("cell count","total cell volume"),pch=c(19,3),
       col=c("red","white"),bty="n",text.col="white",pt.lwd=c(1,2))
#dev.off()

sample.cols <- rev(viridis::viridis(ncol(counts)))

#png(paste0(expid,"_raw.png"), width=400, height=200)
#par(mai=c(.5,.5,.1,.1),mgp=c(1.3,.4,0),tcl=-.25,xaxs="i",yaxs="i")
matplot(d2v(size), counts,type="l",lty=1,xlim=c(0,35),
        col=sample.cols,xlab=expression("cell volume, "*fL),ylim=c(0,1.5e7))
legend("topright", sampleLabels, col=sample.cols,
       lty=1,y.intersp=.6,cex=.6,bty="n",
       ncol=ifelse(length(sampleLabels)>12,2,1))
dev.off()



png(file.path(out.path,paste0(expid,"_CASY_diameter.png")),
    width=2*3.5, height=2*3.5, units="in", res=300)
par(mfcol=c(2,1),mai=c(.5,1,.1,.5),mgp=c(1.3,.4,0),tcl=-.25,xaxs="i",yaxs="i")
image(y=size,x=1:ncol(counts),z=t(counts.all), col=cols,breaks=brks,
      ylab=expression("cell diameter, "*mu*m), xlab="",ylim=c(0,5),
      axes=FALSE)
axis(1, at=1:ncol(counts), label=sampleLabels,las=2)
axis(2)
box()
par(new=TRUE)
plot(1:ncol(counts), total/1e8, type="p",col=2,
     axes=FALSE,xlab=NA,ylab=NA,ylim=c(0,8),
     xlim=par("usr")[1:2],pch=19)
lines(1:ncol(counts), total/1e8,col=2)
axis(4,col=2,col.axis=2)
mtext("1e8 cells/mL",4, par("mgp")[1],col=2)
par(new=TRUE)
plot(1:ncol(counts), volume, type="p",col="white",
     axes=FALSE,xlab=NA,ylab=NA,ylim=c(0,4),
     xlim=par("usr")[1:2],pch=3,lwd=2)
lines(1:ncol(counts), volume,col="white")
axis(2,col="black",col.axis="black",line=par("mgp")[1]*2)
mtext(expression("total cell volume, "*mu*L/mL),2, 3*par("mgp")[1],col="black")
legend("topleft",c("cell count","total cell volume"),pch=c(19,3),
       col=c("red","white"),bty="n",text.col="white",pt.lwd=c(1,2))
#dev.off()

sample.cols <- rev(viridis::viridis(ncol(counts)))

#png(paste0(expid,"_raw.png"), width=400, height=200)
#par(mai=c(.5,.5,.1,.1),mgp=c(1.3,.4,0),tcl=-.25,xaxs="i",yaxs="i")
matplot(size, counts,type="l",lty=1,xlim=c(0,5),
        col=sample.cols,xlab=expression("cell diameter, "*mu*m),ylim=c(0,1.5e7))
legend("topright", sampleLabels, col=sample.cols,
       lty=1,y.intersp=.6,cex=.6,bty="n",
       ncol=ifelse(length(sampleLabels)>12,2,1))
dev.off()

