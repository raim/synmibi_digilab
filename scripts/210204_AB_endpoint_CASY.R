options(stringsAsFactors=FALSE)

## moving average
ma <- function(x,n=10,circular=FALSE)
    stats::filter(x, rep(1/n, n), sides = 2, circular = circular)
## diameter -> volume
d2v <- function(x) (x/2)^3*pi *4/3


expid <- "210204_AB_endpoint"
PATH <- "/mnt/synmibi/Studierende/DATA/CASY/"
if ( Sys.info()["nodename"]=="intron" )
    PATH <- "/data/synmibi/CASY"

in.path <- file.path(PATH,expid)
out.path <- file.path("~/work/CoilHack/experiments/batch/",
                      expid,"analysis")
dat.path <- out.path

cat(paste("PARSING CASY DATA:", date(), "\n"))

## SAMPLE NOTES AND SETTINGS

## experiment parameters
dil <- 3000

## evc und topA KD war glaub ich immer 1:5000 verdünnt
## die anderen haben wir zunächst nur 1:1000 verdünnt,
## sobald die recovered sind dann auch 1:5000

## analysis&plot parameters
normalize <- TRUE
max.cnt <- 10e6 # max cell count for color scheme breaks
## size filter, min/max diameter to count
min.counts <- 1#.5
max.counts <- 5
min.norm <- min.counts
max.norm <- max.counts

## TODO
## *calculate median and peak sizes/volumes
## *calculate total cell volume

files <- list.files(path=in.path, pattern="coilhack_001_.*\\.TXT$")

sizes <- counts <- matrix(NA, nrow=1024, ncol=length(files))
sampleIDs <-  rep(NA, length(files))
times <- rep(list(NA), length(files))
## precalculated values
cvalues <- matrix(NA, nrow=length(files), ncol=c(6))
colnames(cvalues) <- c("Counts/ml","Volume/ml",
                       "Mean Volume (fl)", "Peak Volume (fl)",
                       "Mean Diameter (\xb5m)","Peak Diameter (\xb5m)")
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
    ## TODO: search correct dilution

    ## pre-caculated values
    cvalues[i,] <- as.numeric(trimws(data[match(colnames(cvalues),data[,1]),2]))
    cvalues[i,c("Counts/ml","Volume/ml")] <-
        cvalues[i,c("Counts/ml","Volume/ml")]*DIL

    ## time
    times[[i]] <- strptime(paste(data[7,2],data[8,2]),
                           format="%d.%m.%y %H:%M:%S")
    
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
filter <- grep("^EXP",sampleIDs)
#filter <- filter[grep("_[A-Z]", sampleIDs[filter], invert=TRUE)]
counts <- counts[,filter]
sizes <- sizes[,filter]
sampleIDs <- sampleIDs[filter]
times <-times[filter]
cvalues <- cvalues[filter,]

## TODO:
## * split into experiments (simple: set ID above and run script separately),
## * handle duplicates: average of all values

## fix known sample issues
## sample 08 is indexed as 20110116_18_15
sampleIDs <- sub("20110116_18_15", "08", sampleIDs)

sampleLabels <- sub("^EXP","",sampleIDs)
sampleLabels <- sub("^1_","",sampleLabels)

sampleRuns <- sub("_[0-9].*","", sampleLabels)

## split experiments!
RUNS <- unique(sampleRuns)

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

for ( i in 1:length(RUNS) ) {

    run <- RUNS[i]
    rid <- which(sampleRuns==run)

    ## TODO: loop through duplicates!
    cnts <- counts[,rid]
    cnts.nrm <- counts.nrm[,rid]
    cnts.all <- counts.all[,rid]

    ## sample times
    xtime <- times[rid]
    xtime <- unlist(lapply(xtime, function(x)
        difftime(x,xtime[[1]],units="days")))
    
    file.name <- file.path(out.path,paste0(expid,"_CASY_",run))

    png(paste0(file.name,"_volume.png"),
        width=2*3.5, height=3.5, units="in", res=300)
    par(mfcol=c(1,1),mai=c(.5,1,.1,.5),mgp=c(1.3,.4,0),tcl=-.25,
        xaxs="i",yaxs="i")
    image(y=d2v(size[idx]),x=xtime,z=t(cnts.nrm),
          col=cols,breaks=brks,
          ylab=expression("cell volume, "*fL), xlab="",ylim=c(0,35),
          axes=FALSE)
    axis(1)#, at=1:ncol(cnts), label=sampleLabels[rid],las=2,cex.axis=.7)
    axis(2)
    box()
    par(new=TRUE)
    plot(xtime, total[rid]/1e8, type="p",col=2,
         axes=FALSE,xlab=NA,ylab=NA,ylim=c(0,8),
         xlim=par("usr")[1:2],pch=1)
    ##lines(1:ncol(cnts), total/1e8,col=2)
    axis(4,col=2,col.axis=2)
    mtext("1e8 cells/mL",4, par("mgp")[1],col=2)
    par(new=TRUE)
    plot(xtime, volume[rid], type="p",col="white",
         axes=FALSE,xlab=NA,ylab=NA,ylim=c(0,4),
         xlim=par("usr")[1:2],pch=5,lwd=1)
    ##lines(1:ncol(cnts), volume,col="white")
    axis(2,col="black",col.axis="black",line=par("mgp")[1]*2)
    mtext(expression("total cell volume, "*mu*L/mL),2,
          3*par("mgp")[1],col="black")
    legend("topleft",c("cell count","total cell volume"),pch=c(1,5),
           col=c("red","white"),bty="n",text.col="white",pt.lwd=c(1,1))
    dev.off()
    
    sample.cols <- rev(viridis::viridis(ncol(cnts)))
    
    png(paste0(file.name,"_volume_raw.png"), 
        width=2*3.5, height=3.5, units="in", res=300)
    par(mai=c(.5,.5,.1,.1),mgp=c(1.3,.4,0),tcl=-.25,xaxs="i",yaxs="i")
    matplot(d2v(size), cnts,type="l",lty=1,xlim=c(0,35),
        col=sample.cols,xlab=expression("cell volume, "*fL),ylim=c(0,1.5e7))
    legend("topright", sampleLabels[rid], col=sample.cols,
           lty=1,y.intersp=.6,cex=.6,bty="n",
           ncol=max(c(1,round(length(sampleLabels[rid])/20))))
    dev.off()
    
    
    
    png(paste0(file.name,"_diameter.png"),
        width=2*3.5, height=3.5, units="in", res=300)
    par(mfcol=c(1,1),mai=c(.5,1,.1,.5),mgp=c(1.3,.4,0),tcl=-.25,
        xaxs="i",yaxs="i")
    image(y=size,x=xtime,z=t(cnts.all), col=cols,breaks=brks,
          ylab=expression("cell diameter, "*mu*m), xlab="",ylim=c(0,5),
          axes=FALSE)
    axis(1)#, at=1:ncol(cnts), label=sampleLabels[rid],las=2,cex.axis=.7)
    axis(2)
    box()
    par(new=TRUE)
    plot(xtime, total[rid]/1e8, type="p",col=2,
         axes=FALSE,xlab=NA,ylab=NA,ylim=c(0,8),
         xlim=par("usr")[1:2],pch=1)
#lines(1:ncol(cnts), total/1e8,col=2)
    axis(4,col=2,col.axis=2)
    mtext("1e8 cells/mL",4, par("mgp")[1],col=2)
    par(new=TRUE)
    plot(xtime, volume[rid], type="p",col="white",
         axes=FALSE,xlab=NA,ylab=NA,ylim=c(0,4),
         xlim=par("usr")[1:2],pch=5,lwd=1)
#lines(1:ncol(cnts), volume,col="white")
    axis(2,col="black",col.axis="black",line=par("mgp")[1]*2)
    mtext(expression("total cell volume, "*mu*L/mL),2, 3*par("mgp")[1],col="black")
    legend("topleft",c("cell count","total cell volume"),pch=c(1,5),
           col=c("red","white"),bty="n",text.col="white",pt.lwd=c(1,1))
    dev.off()
    
    sample.cols <- rev(viridis::viridis(ncol(cnts)))

    png(paste0(file.name,"_diameter_raw.png"), 
        width=2*3.5, height=3.5, units="in", res=300)
    par(mai=c(.5,.5,.1,.1),mgp=c(1.3,.4,0),tcl=-.25,xaxs="i",yaxs="i")
    matplot(size, cnts,type="l",lty=1,xlim=c(0,5),
            col=sample.cols,xlab=expression("cell diameter, "*mu*m),
            ylim=c(0,1.5e7))
    legend("topright", sampleLabels[rid], col=sample.cols,
           lty=1,y.intersp=.6,cex=.6,bty="n",
           ncol=max(c(1,round(length(sampleLabels[rid])/20))))
    dev.off()
    
    ## WRITE-OUT RESULTS
    ## TODO: summarize duplicates,
    ## TODO: write out count, total volume, peak volume
    colnames(cnts) <- sampleIDs[rid]
    
    allc <- cbind(diameter=size, volume=d2v(size), cnts)
    write.table(allc, file=paste0(file.name,".tsv"),
                quote=FALSE,sep="\t",row.names=FALSE)
    
    summary <- data.frame("sample"=sampleIDs[rid],
                          `cells/mL`=total[rid],
                          `volume,uL/mL`=volume[rid],
                          check.names=FALSE)
    summary <- cbind.data.frame(summary, cvalues)
    write.table(summary, file=paste0(file.name,"_summary.tsv"),
            quote=FALSE,sep="\t",row.names=FALSE)
}
