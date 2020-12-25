
options(stringsAsFactors=FALSE)

expid <- "20201204_RM_topA"
PATH <- "/mnt/synmibi/Studierende/DATA/Specord/"
if ( Sys.info()["nodename"]=="intron" )
    PATH <- "/data/synmibi/Specord"

in.path <- file.path(PATH,expid)
out.path <- file.path("~/work/CoilHack/experiments/reactor/pcc6803/",
                      expid,"analysis")


cat(paste("PARSING SPECORD DATA:", date(), "\n"))

## switch between normalization wavelengths
nrm <-443 # 684 # NULL # # reference wavelength for normalization
REF <- "SAMPLE_07" # as.character(7) ## reference sample
#REF2 <- as.character(-2) ## reference sample

files <- list.files(path=in.path, pattern="SAMPLE_.*\\.csv$")

smp <- sub("\\.csv","",files)
wl <- 400:750 # wavelengths
spc <- matrix(NA, ncol=length(wl), nrow=length(files))
rownames(spc) <- smp

## load spectra from Specord
for ( i in seq_along(files) ) {
    dat <- try(read.csv(file.path(in.path,files[i]), sep=";",skip=1))
    if ( inherits(dat, 'try-error')) next
    idx <- grep("SAMPLE",colnames(dat))
    if ( length(idx)==1 ) {
        dat <- dat[,idx]
        spc[i,] <- dat
    }
}

## remove NA samples
rm <- apply(spc,1,function(x) any(is.na(x)))
smp <- smp[!rm]
spc <- spc[!rm,]

## sample colors
cols <- sub("FF$","99",rev(viridis::viridis(nrow(spc))))
ltys <- rep(c(1,2), length(cols)/2)

nspc <- spc
    
## subtract 750 nm baseline
spc <- spc - spc[,wl==750]
    
nspc <- spc
## normalize to max
if ( !is.null(nrm) )
    nspc <- spc/spc[,wl==nrm] #apply(spc,1,max) #max at 443 nm

## calculate ratios
rspc <- t(t(as.matrix(nspc)) / unlist(nspc[smp==REF,]))
rspc2 <- t(t(as.matrix(spc))/unlist(spc[smp==REF,]))
    
    
### PLOTS
    
##
png(file.path(out.path,paste0(expid,"_Specord.png")),
    units="in", height=2*3.5, width=2*3.5, res=200)
par(mfcol=c(2,2),mai=c(.25,.35,.1,.1), mgp=c(1.1,.1,0), tcl=.25)
matplot(x=wl,t(spc),type="l",lwd=2,lty=ltys,col=cols,
        xlab=NA,ylab=expression("corrected absorbance"~A[corr]))
axis(3, labels=NA, tcl=.25)
matplot(x=wl,t(nspc),type="l",lwd=2,lty=ltys,col=cols,
        xlab=NA,ylab=expression("normalized absorbance"~A[norm]))
legend("bottomleft",smp, col=cols,lty=ltys,lwd=2, cex=.75,
           y.intersp=.7,bty="n",
       ncol=ifelse(length(smp)>12,2,1))
axis(3, labels=NA, tcl=.25)
mtext("wavelength, nm", 3, 0)
matplot(x=wl,t(rspc2),type="l",lwd=2,lty=ltys,col=cols,
        xlab=NA,ylab=bquote(ratio~A[coor]/A[.(REF)]),ylim=c(0,1.2))
axis(3, labels=NA, tcl=.25)
matplot(x=wl,t(rspc),type="l",lwd=2,lty=ltys,col=cols,
        xlab=NA,ylab=bquote(ratio~A[norm]/A[.(REF)]),ylim=c(.5,1.5))
axis(3, labels=NA, tcl=.25)
mtext("wavelength, nm", 3, 0)
dev.off()
