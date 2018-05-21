###############
# DEPENDENCIES
###############
library(extrafont)
library(fontcm)
loadfonts(quiet = TRUE)

#######################
# AUXILLIARY FUNCTIONS
#######################

toPdf <- function(expr, filename, ...) {
  toDev(expr, pdf, filename, ...)
}

figPath  <-  function(name) {
  file.path('output/figures', name)
}

toDev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf('Creating %s\n', filename))
  dev(filename, family='CM Roman', ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}



####################
# PLOTTING FUNCTIONS
####################

#' Plot text or points according to relative axis position.
#'
#' @title Plot text or points according to relative axis position
#' @param px Relative x-axis position (in proportion) where character is to be plotted.
#' @param py Relative y-axis position (in proportion) where character is to be plotted.
#' @param lab Plotted text. Works if argument \code{\link[graphics]{text}} is \code{TRUE}.
#' @param adj See argument of same name in R base function \code{\link[graphics]{par}}.
#' @param text Logical. Should text or points be plotted?
#' @param log Used if the original plot uses the argument log, e.g. \code{log='x'}, \code{log='y'} or \code{log='xy'}.
#' @param ... Additional arguments to R base function \code{\link[graphics]{text}}.
#' @export
proportionalLabel <- function(px, py, lab, adj=c(0, 1), text=TRUE, log=FALSE, ...) {
    usr  <-  par('usr')
    x.p  <-  usr[1] + px*(usr[2] - usr[1])
    y.p  <-  usr[3] + py*(usr[4] - usr[3])
    if(log=='x') {
        x.p<-10^(x.p)
    }
    if(log=='y') {
        y.p<-10^(y.p)
    }
    if(log=='xy') {
        x.p<-10^(x.p)
        y.p<-10^(y.p)
    }
    if(text){
        text(x.p, y.p, lab, adj=adj, ...)
    } else {
        points(x.p, y.p, ...)
    }
}

#' Draw equally-spaced white lines on plot window.
#'
#' @title Equally-spaced white lines on plot window
#' @param ... Additional arguments to internal function \code{\link{proportionalLabel}}.
#' @author Diego Barneche
#' @export
plotGrid  <-  function(lineCol='white',...) {
    proportionalLabel(rep(0.2, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.4, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.6, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.8, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.2, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.4, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.6, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.8, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
}


#' Internal. Create nice rounded numbers for plotting.
#'
#' @title Rounded numbers for plotting
#' @param value A numeric vector.
#' @param precision Number of rounding digits.
#' @return A character vector.
#' @author Diego Barneche.
rounded  <-  function(value, precision=1) {
  sprintf(paste0('%.', precision, 'f'), round(value, precision))
}


#' Creates transparent colours
#'
#' @title Creates transparent colours
#' @param col Colour.
#' @param opacity equivalent to alpha transparency parameter
#' @export
transparentColor <- function(col, opacity=0.5) {
    if (length(opacity) > 1 && any(is.na(opacity))) {
        n        <-  max(length(col), length(opacity))
        opacity  <-  rep(opacity, length.out=n)
        col      <-  rep(col, length.out=n)
        ok       <-  !is.na(opacity)
        ret      <-  rep(NA, length(col))
        ret[ok]  <-  Recall(col[ok], opacity[ok])
        ret
    } else {
        tmp  <-  col2rgb(col)/255
        rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
    }
}


##############################################################
##############################################################
##  Final figures for paper







##############################################################
##############################################################
##  Exploratory figures


# Simple figs for 2-patch Levene model for simultaneous hermaphrodites
proportionPolyMultiPatch  <-  function(h, delta, sMax=1) {
    
    # import data
    filename1  <-  paste("./output/data/simMultiPatchSgrad", "_C", 0,    "_delta", delta, "_hf", h, "_hm", h, "_sMax", sMax, ".csv", sep="")
    filename2  <-  paste("./output/data/simMultiPatchSgrad", "_C", 0.25, "_delta", delta, "_hf", h, "_hm", h, "_sMax", sMax, ".csv", sep="")
    filename3  <-  paste("./output/data/simMultiPatchSgrad", "_C", 0.5,  "_delta", delta, "_hf", h, "_hm", h, "_sMax", sMax, ".csv", sep="")
    filename4  <-  paste("./output/data/simMultiPatchSgrad", "_C", 0.75, "_delta", delta, "_hf", h, "_hm", h, "_sMax", sMax, ".csv", sep="")
    data1  <-  read.csv(filename1, header=TRUE)
    data2  <-  read.csv(filename2, header=TRUE)
    data3  <-  read.csv(filename3, header=TRUE)
    data4  <-  read.csv(filename4, header=TRUE)

    # Color Scheme
    COLS  <-  c("dodgerblue", "dodgerblue2", "dodgerblue3", "dodgerblue4", "black")

    # Set plot layout
    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A: C = 0
    par(omi=c(0.5, 0.5, 0.5, 0.5), mar = c(4,4,0.5,0.5), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        lines(poly5 ~ sMax, data=data1, lwd=3, col=COLS[1])
        lines(poly4 ~ sMax, data=data1, lwd=3, col=COLS[2])
        lines(poly3 ~ sMax, data=data1, lwd=3, col=COLS[3])
        lines(poly2 ~ sMax, data=data1, lwd=3, col=COLS[4])
        lines(poly1 ~ sMax, data=data1, lwd=3, col=COLS[5])
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(C)," = 0")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        if(h == 1/2) {
                legend(
                    x       =  usr[2]*0.425,
                    y       =  usr[4],
                    legend  =  c(
                                expression(paste("5 Patch")),
                                expression(paste("4 Patch")),
                                expression(paste("3 Patch")),
                                expression(paste("2 Patch")),
                                expression(paste("1 Patch"))),
                    col     =  COLS,
                    lty     =  c(1,1),
                    lwd     =  c(3,3),
                    cex     =  0.9,
                    xjust   =  1,
                    yjust   =  1,
                    bty     =  'n',
                    border  =  NA
            )
        }
        if(h == 1/4) {
                legend(
                    x       =  usr[2]*1,
                    y       =  usr[4]*0.4,
                    legend  =  c(
                                expression(paste("5 Patch")),
                                expression(paste("4 Patch")),
                                expression(paste("3 Patch")),
                                expression(paste("2 Patch")),
                                expression(paste("1 Patch"))),
                    col     =  COLS,
                    lty     =  c(1,1),
                    lwd     =  c(3,3),
                    cex     =  1,
                    xjust   =  1,
                    yjust   =  1,
                    bty     =  'n',
                    border  =  NA
            )
        }

## Panel B: C = 1/4
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        lines(poly5 ~ sMax, data=data2, lwd=3, col=COLS[1])
        lines(poly4 ~ sMax, data=data2, lwd=3, col=COLS[2])
        lines(poly3 ~ sMax, data=data2, lwd=3, col=COLS[3])
        lines(poly2 ~ sMax, data=data2, lwd=3, col=COLS[4])
        lines(poly1 ~ sMax, data=data2, lwd=3, col=COLS[5])
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(C)," = 0.25")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)

## Panel C: C = 1/2
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        lines(poly5 ~ sMax, data=data3, lwd=3, col=COLS[1])
        lines(poly4 ~ sMax, data=data3, lwd=3, col=COLS[2])
        lines(poly3 ~ sMax, data=data3, lwd=3, col=COLS[3])
        lines(poly2 ~ sMax, data=data3, lwd=3, col=COLS[4])
        lines(poly1 ~ sMax, data=data3, lwd=3, col=COLS[5])
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(C)," = 0.5")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  1.2,   expression(paste("Proportion parameter space")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel(1.2,  -0.3,   expression(paste("Maximum strength of selection: max(",italic(s[f]),", ", italic(s[m]),")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)        

## Panel D: C = 3/4
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        lines(poly5 ~ sMax, data=data4, lwd=3, col=COLS[1])
        lines(poly4 ~ sMax, data=data4, lwd=3, col=COLS[2])
        lines(poly3 ~ sMax, data=data4, lwd=3, col=COLS[3])
        lines(poly2 ~ sMax, data=data4, lwd=3, col=COLS[4])
        lines(poly1 ~ sMax, data=data4, lwd=3, col=COLS[5])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(D))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(C)," = 0.75")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.3,  0.5,   expression(paste("Proportion parameter space")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        
}






# Simple figs for 2-patch Levene model for simultaneous hermaphrodites
proportionPolyMultiPatchHDelta  <-  function(h, delta, sMax=1) {
    
    # import data
    filename1  <-  paste("./output/data/simMultiPatchSgrad", "_C", 0,    "_delta", delta[1], "_hf", h[1], "_hm", h[1], "_sMax", sMax, ".csv", sep="")
    filename2  <-  paste("./output/data/simMultiPatchSgrad", "_C", 0.25, "_delta", delta[1], "_hf", h[1], "_hm", h[1], "_sMax", sMax, ".csv", sep="")
    filename3  <-  paste("./output/data/simMultiPatchSgrad", "_C", 0.5,  "_delta", delta[1], "_hf", h[1], "_hm", h[1], "_sMax", sMax, ".csv", sep="")
    filename4  <-  paste("./output/data/simMultiPatchSgrad", "_C", 0.75, "_delta", delta[1], "_hf", h[1], "_hm", h[1], "_sMax", sMax, ".csv", sep="")
    filename5  <-  paste("./output/data/simMultiPatchSgrad", "_C", 0,    "_delta", delta[2], "_hf", h[2], "_hm", h[2], "_sMax", sMax, ".csv", sep="")
    filename6  <-  paste("./output/data/simMultiPatchSgrad", "_C", 0.25, "_delta", delta[2], "_hf", h[2], "_hm", h[2], "_sMax", sMax, ".csv", sep="")
    filename7  <-  paste("./output/data/simMultiPatchSgrad", "_C", 0.5,  "_delta", delta[2], "_hf", h[2], "_hm", h[2], "_sMax", sMax, ".csv", sep="")
    filename8  <-  paste("./output/data/simMultiPatchSgrad", "_C", 0.75, "_delta", delta[2], "_hf", h[2], "_hm", h[2], "_sMax", sMax, ".csv", sep="")
    data1      <-  read.csv(filename1, header=TRUE)
    data2      <-  read.csv(filename2, header=TRUE)
    data3      <-  read.csv(filename3, header=TRUE)
    data4      <-  read.csv(filename4, header=TRUE)
    data5      <-  read.csv(filename5, header=TRUE)
    data6      <-  read.csv(filename6, header=TRUE)
    data7      <-  read.csv(filename7, header=TRUE)
    data8      <-  read.csv(filename8, header=TRUE)

    # Color Scheme
    COLS  <-  c("dodgerblue", "dodgerblue2", "dodgerblue3", "dodgerblue4", "black")

    # Set plot layout
    layout.mat <- matrix(c(1:8), nrow=2, ncol=4, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

############
## Delta = 0
## Panel A: C = 0
    par(omi=c(0.75, 0.75, 0.75, 0.75), mar = c(2.5,2.5,0.5,0.5), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        lines(poly5 ~ sMax, data=data1, lwd=2, col=COLS[1])
        lines(poly4 ~ sMax, data=data1, lwd=2, col=COLS[2])
        lines(poly3 ~ sMax, data=data1, lwd=2, col=COLS[3])
        lines(poly2 ~ sMax, data=data1, lwd=2, col=COLS[4])
        lines(poly1 ~ sMax, data=data1, lwd=2, col=COLS[5])
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.1, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.3,   expression(paste(italic(C)," = 0")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( -0.7,  0.5,   expression(paste(italic(h[f])," = ",italic(h[m]), " = 1/2")), cex=1.5, adj=c(0.5, 0.5), xpd=NA,srt=90)


## Panel B: C = 1/4
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        lines(poly5 ~ sMax, data=data2, lwd=2, col=COLS[1])
        lines(poly4 ~ sMax, data=data2, lwd=2, col=COLS[2])
        lines(poly3 ~ sMax, data=data2, lwd=2, col=COLS[3])
        lines(poly2 ~ sMax, data=data2, lwd=2, col=COLS[4])
        lines(poly1 ~ sMax, data=data2, lwd=2, col=COLS[5])
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
#        proportionalLabel( 1.3,  1.3,   expression(paste(italic(h[f])," = ",italic(h[m]), " = 1/2")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.05,  1.1, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.3,   expression(paste(italic(C)," = 0.25")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)

## Panel C: C = 1/2
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        lines(poly5 ~ sMax, data=data3, lwd=2, col=COLS[1])
        lines(poly4 ~ sMax, data=data3, lwd=2, col=COLS[2])
        lines(poly3 ~ sMax, data=data3, lwd=2, col=COLS[3])
        lines(poly2 ~ sMax, data=data3, lwd=2, col=COLS[4])
        lines(poly1 ~ sMax, data=data3, lwd=2, col=COLS[5])
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.1, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.3,   expression(paste(italic(C)," = 0.5")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)

## Panel D: C = 3/4
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        lines(poly5 ~ sMax, data=data4, lwd=2, col=COLS[1])
        lines(poly4 ~ sMax, data=data4, lwd=2, col=COLS[2])
        lines(poly3 ~ sMax, data=data4, lwd=2, col=COLS[3])
        lines(poly2 ~ sMax, data=data4, lwd=2, col=COLS[4])
        lines(poly1 ~ sMax, data=data4, lwd=2, col=COLS[5])
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.1, expression(paste(bold(D))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.3,   expression(paste(italic(C)," = 0.75")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # Legend
                legend(
                    x       =  usr[2]*0.75,
                    y       =  usr[4],
                    legend  =  c(
                                expression(paste("5 Patch")),
                                expression(paste("4 Patch")),
                                expression(paste("3 Patch")),
                                expression(paste("2 Patch")),
                                expression(paste("1 Patch"))),
                    col     =  COLS,
                    lty     =  c(1,1),
                    lwd     =  2,
                    cex     =  0.9,
                    xjust   =  1,
                    yjust   =  1,
                    bty     =  'n',
                    border  =  NA
            )


## Delta = 1/2
## Panel E: C = 0
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        lines(poly5 ~ sMax, data=data5, lwd=2, col=COLS[1])
        lines(poly4 ~ sMax, data=data5, lwd=2, col=COLS[2])
        lines(poly3 ~ sMax, data=data5, lwd=2, col=COLS[3])
        lines(poly2 ~ sMax, data=data5, lwd=2, col=COLS[4])
        lines(poly1 ~ sMax, data=data5, lwd=2, col=COLS[5])
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.1, expression(paste(bold(E))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.4,  1.2,   expression(paste("Proportion parameter space")), cex=1.3, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( -0.7,  0.5,   expression(paste(italic(h[f])," = ",italic(h[m]), " = 1/4")), cex=1.5, adj=c(0.5, 0.5), xpd=NA,srt=90)

## Panel F: C = 1/4
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        lines(poly5 ~ sMax, data=data6, lwd=2, col=COLS[1])
        lines(poly4 ~ sMax, data=data6, lwd=2, col=COLS[2])
        lines(poly3 ~ sMax, data=data6, lwd=2, col=COLS[3])
        lines(poly2 ~ sMax, data=data6, lwd=2, col=COLS[4])
        lines(poly1 ~ sMax, data=data6, lwd=2, col=COLS[5])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
#        proportionalLabel( 1.3,  1.3,   expression(paste(italic(h[f])," = ",italic(h[m]), " = 1/2")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.05,  1.1, expression(paste(bold(F))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(1.25,  -0.4,   expression(paste("Maximum strength of selection: max(",italic(s[f]),", ", italic(s[m]),")")), cex=1.3, adj=c(0.5, 0.5), xpd=NA)        

## Panel G: C = 1/2
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        lines(poly5 ~ sMax, data=data7, lwd=2, col=COLS[1])
        lines(poly4 ~ sMax, data=data7, lwd=2, col=COLS[2])
        lines(poly3 ~ sMax, data=data7, lwd=2, col=COLS[3])
        lines(poly2 ~ sMax, data=data7, lwd=2, col=COLS[4])
        lines(poly1 ~ sMax, data=data7, lwd=2, col=COLS[5])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.1, expression(paste(bold(G))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)

## Panel H: C = 3/4
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        lines(poly5 ~ sMax, data=data8, lwd=2, col=COLS[1])
        lines(poly4 ~ sMax, data=data8, lwd=2, col=COLS[2])
        lines(poly3 ~ sMax, data=data8, lwd=2, col=COLS[3])
        lines(poly2 ~ sMax, data=data8, lwd=2, col=COLS[4])
        lines(poly1 ~ sMax, data=data8, lwd=2, col=COLS[5])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.1, expression(paste(bold(H))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
              
}



# Simple figs for 2-patch Levene model for simultaneous hermaphrodites
diffPolyMultiPatch  <-  function(h = 1/2, delta = 0, sMax = 1) {

    # import data
    filename1  <-  paste("./output/data/simMultiPatchSgrad", "_C", 0,    "_delta", delta, "_hf", h, "_hm", h, "_sMax", sMax, ".csv", sep="")
    filename2  <-  paste("./output/data/simMultiPatchSgrad", "_C", 0.25, "_delta", delta, "_hf", h, "_hm", h, "_sMax", sMax, ".csv", sep="")
    filename3  <-  paste("./output/data/simMultiPatchSgrad", "_C", 0.5,  "_delta", delta, "_hf", h, "_hm", h, "_sMax", sMax, ".csv", sep="")
    filename4  <-  paste("./output/data/simMultiPatchSgrad", "_C", 0.75, "_delta", delta, "_hf", h, "_hm", h, "_sMax", sMax, ".csv", sep="")
    data1  <-  read.csv(filename1, header=TRUE)
    data2  <-  read.csv(filename2, header=TRUE)
    data3  <-  read.csv(filename3, header=TRUE)
    data4  <-  read.csv(filename4, header=TRUE)

    # set ylims
    ymin  <-  min(c(data1$diffPoly12, data2$diffPoly12, data3$diffPoly12, data4$diffPoly12))
    ymax  <-  max(c(data1$diffPoly12, data2$diffPoly12, data3$diffPoly12, data4$diffPoly12), 0.275)

## Panel A: C = 0
    par(omi=c(0.5, 0.5, 0.5, 0.5), mar = c(4,4,0.5,0.5), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(ymin, ymax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        lines(diffPoly12 ~ sMax, lwd=3, col='grey70', data=data1)
        lines(diffPoly12 ~ sMax, lwd=3, col='grey60', data=data2)
        lines(diffPoly12 ~ sMax, lwd=3, col='grey50', data=data3)
        lines(diffPoly12 ~ sMax, lwd=3, col='grey40', data=data4)
        abline(h=0, lwd=2, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.4,  1.05,   expression(paste(italic(delta)," =")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.46, 1.05,  substitute(d,list(d=delta)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.56,  1.05,   expression(paste(", ", italic(h["m,f"])," =")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.69, 1.05,  substitute(dom,list(dom=h)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5,  -0.2,   expression(paste(italic(s[max]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
              x       =  usr[2]*0.25,
              y       =  usr[4],
              legend  =  c(
                          expression(paste(italic(C)," = ",0)),
                          expression(paste(italic(C)," = ",0.25)),
                          expression(paste(italic(C)," = ",0.50)),
                          expression(paste(italic(C)," = ",0.75))),
              col     =  c('grey70','grey60','grey50','grey40'),
              lty     =  1,
              lwd     =  3,
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )


}



invConditionsSA  <-  function() {
    # constants
    sm  <-  seq(0,1,by=0.0001)

    # Color Scheme
    COLS  <-  c("dodgerblue", "black")

    # Set plot layout
    layout.mat <- matrix(c(1:8), nrow=2, ncol=4, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

### Panels A -- D: Additive fitness effects
    h  <-  1/2
## Panel A: C = 0
    # Calculate invasion conditions
    UB   <-  InvB(C = 0, delta = 0, hf = h, hm = h, sm = sm)
    UB[UB > 1]  <-  1.00000001
    LB   <-  InvA(C = 0, delta = 0, hf = h, hm = h, sm = sm)
    UBd  <-  InvB(C = 0, delta = 0.5, hf = h, hm = h, sm = sm)
    UBd[UBd > 1]  <-  1.00000001
    LBd  <-  InvA(C = 0, delta = 0.5, hf = h, hm = h, sm = sm)

    par(omi=c(0.5, 0.5, 0.5, 0.5), mar = c(4,4,0.5,0.5), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, 1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        polygon(c(rev(sm),sm), c(rev(LB), UB), col=transparentColor('grey80', 0.6), border='grey70')
        polygon(c(rev(sm),sm), c(rev(LBd), UBd), col=transparentColor('dodgerblue', 0.15), border='grey70')
        lines(UB[UB<=1] ~ sm[UB<= 1], lwd=3, col=COLS[2])
        lines(LB ~ sm, lwd=3, col=COLS[2])
        lines(UBd[UBd<=1] ~ sm[UBd<=1], lwd=3, col=COLS[1])
        lines(LBd ~ sm, lwd=3, col=COLS[1])
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(C)," = 0")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste(italic(s[f]))), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( -0.7,  0.5,   expression(paste(italic(h[f])," = ",italic(h[m]), " = 1/2")), cex=1.5, adj=c(0.5, 0.5), xpd=NA,srt=90)

## Panel B: C = 1/4
    # Calculate invasion conditions
    UB   <-  InvB(C = 1/4, delta = 0, hf = h, hm = h, sm = sm)
    UB[UB > 1]  <-  1.00000001
    LB   <-  InvA(C = 1/4, delta = 0, hf = h, hm = h, sm = sm)
    UBd  <-  InvB(C = 1/4, delta = 0.5, hf = h, hm = h, sm = sm)
    UBd[UBd > 1]  <-  1.00000001
    LBd  <-  InvA(C = 1/4, delta = 0.5, hf = h, hm = h, sm = sm)

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, 1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        polygon(c(rev(sm),sm), c(rev(LB), UB), col=transparentColor('grey80', 0.3), border='grey70')
        polygon(c(rev(sm),sm), c(rev(LBd), UBd), col=transparentColor('dodgerblue', 0.15), border='grey70')
        lines(UB[UB<=1] ~ sm[UB<= 1], lwd=3, col=COLS[2])
        lines(LB ~ sm, lwd=3, col=COLS[2])
        lines(UBd[UBd<=1] ~ sm[UBd<=1], lwd=3, col=COLS[1])
        lines(LBd ~ sm, lwd=3, col=COLS[1])
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(C)," = 0.25")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)

## Panel C: C = 1/2
    # Calculate invasion conditions
    UB   <-  InvB(C = 1/2, delta = 0, hf = h, hm = h, sm = sm)
    UB[UB > 1]  <-  1.00000001
    LB   <-  InvA(C = 1/2, delta = 0, hf = h, hm = h, sm = sm)
    UBd  <-  InvB(C = 1/2, delta = 0.5, hf = h, hm = h, sm = sm)
    UBd[UBd > 1]  <-  1.00000001
    LBd  <-  InvA(C = 1/2, delta = 0.5, hf = h, hm = h, sm = sm)

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, 1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        polygon(c(rev(sm),sm), c(rev(LB), UB), col=transparentColor('grey80', 0.6), border='grey70')
        polygon(c(rev(sm),sm), c(rev(LBd), UBd), col=transparentColor('dodgerblue', 0.15), border='grey70')
        lines(UB[UB<=1] ~ sm[UB<= 1], lwd=3, col=COLS[2])
        lines(LB ~ sm, lwd=3, col=COLS[2])
        lines(UBd[UBd<=1] ~ sm[UBd<=1], lwd=3, col=COLS[1])
        lines(LBd ~ sm, lwd=3, col=COLS[1])
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(C)," = 0.5")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)

## Panel D: C = 3/4
    # Calculate invasion conditions
    UB   <-  InvB(C = 3/4, delta = 0, hf = h, hm = h, sm = sm)
    UB[UB > 1]  <-  1.00000001
    LB   <-  InvA(C = 3/4, delta = 0, hf = h, hm = h, sm = sm)
    UBd  <-  InvB(C = 3/4, delta = 0.5, hf = h, hm = h, sm = sm)
    UBd[UBd > 1]  <-  1.00000001
    LBd  <-  InvA(C = 3/4, delta = 0.5, hf = h, hm = h, sm = sm)

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, 1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        polygon(c(rev(sm),sm), c(rev(LBd), UBd), col=transparentColor('dodgerblue', 0.15), border='grey70')
        polygon(c(rev(sm),sm), c(rev(LB), UB), col=transparentColor('grey80', 0.6), border='grey70')
        lines(UBd[UBd<=1] ~ sm[UBd<=1], lwd=3, col=COLS[1])
        lines(LBd ~ sm, lwd=3, col=COLS[1])
        lines(UB[UB<=1] ~ sm[UB<= 1], lwd=3, col=COLS[2])
        lines(LB ~ sm, lwd=3, col=COLS[2])
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(D))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(C)," = 0.75")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # Legend
            legend(
                x       =  usr[2]*0.5,
                y       =  usr[4],
                legend  =  c(
                            expression(paste(delta, " = 1/2")),
                            expression(paste(delta, " = 0"))),
                col     =  COLS,
                lty     =  c(1,1),
                lwd     =  c(3,3),
                cex     =  0.9,
                xjust   =  1,
                yjust   =  1,
                bty     =  'n',
                border  =  NA
                    )

### Dominance reversal        
    h  <-  1/4
## Panel E: C = 0
    # Calculate invasion conditions
    UB   <-  InvB(C = 0, delta = 0, hf = h, hm = h, sm = sm)
    UB[UB > 1]  <-  1.00000001
    LB   <-  InvA(C = 0, delta = 0, hf = h, hm = h, sm = sm)
    UBd  <-  InvB(C = 0, delta = 0.5, hf = h, hm = h, sm = sm)
    UBd[UBd > 1]  <-  1.00000001
    LBd  <-  InvA(C = 0, delta = 0.5, hf = h, hm = h, sm = sm)

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, 1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        polygon(c(rev(sm),sm), c(rev(LB), UB), col=transparentColor('grey80', 0.6), border='grey70')
        polygon(c(rev(sm),sm), c(rev(LBd), UBd), col=transparentColor('dodgerblue', 0.15), border='grey70')
        lines(UB[UB<=1] ~ sm[UB<= 1], lwd=3, col=COLS[2])
        lines(LB ~ sm, lwd=3, col=COLS[2])
        lines(UBd[UBd<=1] ~ sm[UBd<=1], lwd=3, col=COLS[1])
        lines(LBd ~ sm, lwd=3, col=COLS[1])
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(E))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste(italic(s[f]))), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel(0.5,  -0.3,   expression(paste(italic(s[m]))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)        
        proportionalLabel( -0.7,  0.5,   expression(paste(italic(h[f])," = ",italic(h[m]), " = 1/4")), cex=1.5, adj=c(0.5, 0.5), xpd=NA,srt=90)

## Panel F: C = 1/4
    # Calculate invasion conditions
    UB   <-  InvB(C = 1/4, delta = 0, hf = h, hm = h, sm = sm)
    UB[UB > 1]  <-  1.00000001
    LB   <-  InvA(C = 1/4, delta = 0, hf = h, hm = h, sm = sm)
    UBd  <-  InvB(C = 1/4, delta = 0.5, hf = h, hm = h, sm = sm)
    UBd[UBd > 1]  <-  1.00000001
    LBd  <-  InvA(C = 1/4, delta = 0.5, hf = h, hm = h, sm = sm)

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, 1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        polygon(c(rev(sm),sm), c(rev(LB), UB), col=transparentColor('grey80', 0.3), border='grey70')
        polygon(c(rev(sm),sm), c(rev(LBd), UBd), col=transparentColor('dodgerblue', 0.15), border='grey70')
        lines(UB[UB<=1] ~ sm[UB<= 1], lwd=3, col=COLS[2])
        lines(LB ~ sm, lwd=3, col=COLS[2])
        lines(UBd[UBd<=1] ~ sm[UBd<=1], lwd=3, col=COLS[1])
        lines(LBd ~ sm, lwd=3, col=COLS[1])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(F))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5,  -0.3,   expression(paste(italic(s[m]))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)        

## Panel G: C = 1/2
    # Calculate invasion conditions
    UB   <-  InvB(C = 1/2, delta = 0, hf = h, hm = h, sm = sm)
    UB[UB > 1]  <-  1.00000001
    LB   <-  InvA(C = 1/2, delta = 0, hf = h, hm = h, sm = sm)
    UBd  <-  InvB(C = 1/2, delta = 0.5, hf = h, hm = h, sm = sm)
    UBd[UBd > 1]  <-  1.00000001
    LBd  <-  InvA(C = 1/2, delta = 0.5, hf = h, hm = h, sm = sm)

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, 1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        polygon(c(rev(sm),sm), c(rev(LB), UB), col=transparentColor('grey80', 0.6), border='grey70')
        polygon(c(rev(sm),sm), c(rev(LBd), UBd), col=transparentColor('dodgerblue', 0.15), border='grey70')
        lines(UB[UB<=1] ~ sm[UB<= 1], lwd=3, col=COLS[2])
        lines(LB ~ sm, lwd=3, col=COLS[2])
        lines(UBd[UBd<=1] ~ sm[UBd<=1], lwd=3, col=COLS[1])
        lines(LBd ~ sm, lwd=3, col=COLS[1])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(G))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5,  -0.3,   expression(paste(italic(s[m]))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)        

## Panel H: C = 3/4
    # Calculate invasion conditions
    UB   <-  InvB(C = 3/4, delta = 0, hf = h, hm = h, sm = sm)
    UB[UB > 1]  <-  1.00000001
    LB   <-  InvA(C = 3/4, delta = 0, hf = h, hm = h, sm = sm)
    UBd  <-  InvB(C = 3/4, delta = 0.5, hf = h, hm = h, sm = sm)
    UBd[UBd > 1]  <-  1.00000001
    LBd  <-  InvA(C = 3/4, delta = 0.5, hf = h, hm = h, sm = sm)

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, 1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot proportion of parameter space for 1 and 2 patches
        polygon(c(rev(sm),sm), c(rev(LBd), UBd), col=transparentColor('dodgerblue', 0.15), border='grey70')
        polygon(c(rev(sm),sm), c(rev(LB), UB), col=transparentColor('grey80', 0.6), border='grey70')
        lines(UBd[UBd<=1] ~ sm[UBd<=1], lwd=3, col=COLS[1])
        lines(LBd ~ sm, lwd=3, col=COLS[1])
        lines(UB[UB<=1] ~ sm[UB<= 1], lwd=3, col=COLS[2])
        lines(LB ~ sm, lwd=3, col=COLS[2])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(H))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5,  -0.3,   expression(paste(italic(s[m]))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)        
        

}










#' Fig.SX: Supplementary figure showing comparison between deterministic
#'         recursion simulations and invasion analysis based on eigenvalues
#'         for a single patch.
#' 
#'
#' @title Fig.S1: Supplementary figure showing comparison between deterministic
#'                recursion simulations and invasion analysis based on eigenvalues
#' @export
compareSimEig1PatchFig  <-  function(h = c(0.5, 0.25), delta = 0, sMax=1) {
    
    # import data
    filename1  <-  paste("./output/data/determFwdSimLoop", "_C", 0,    "_delta", delta, "_h", h[1], "_sMax", sMax, "_n", 10000, ".txt", sep="")
    filename2  <-  paste("./output/data/determFwdSimLoop", "_C", 0.25, "_delta", delta, "_h", h[1], "_sMax", sMax, "_n", 10000, ".txt", sep="")
    filename3  <-  paste("./output/data/determFwdSimLoop", "_C", 0.5,  "_delta", delta, "_h", h[1], "_sMax", sMax, "_n", 10000, ".txt", sep="")
    filename4  <-  paste("./output/data/determFwdSimLoop", "_C", 0.75, "_delta", delta, "_h", h[1], "_sMax", sMax, "_n", 10000, ".txt", sep="")
    filename5  <-  paste("./output/data/determFwdSimLoop", "_C", 0,    "_delta", delta, "_h", h[2], "_sMax", sMax, "_n", 10000, ".txt", sep="")
    filename6  <-  paste("./output/data/determFwdSimLoop", "_C", 0.25, "_delta", delta, "_h", h[2], "_sMax", sMax, "_n", 10000, ".txt", sep="")
    filename7  <-  paste("./output/data/determFwdSimLoop", "_C", 0.5,  "_delta", delta, "_h", h[2], "_sMax", sMax, "_n", 10000, ".txt", sep="")
    filename8  <-  paste("./output/data/determFwdSimLoop", "_C", 0.75, "_delta", delta, "_h", h[2], "_sMax", sMax, "_n", 10000, ".txt", sep="")
    data1      <-  read.table(filename1, header=TRUE)
    data2      <-  read.table(filename2, header=TRUE)
    data3      <-  read.table(filename3, header=TRUE)
    data4      <-  read.table(filename4, header=TRUE)
    data5      <-  read.table(filename5, header=TRUE)
    data6      <-  read.table(filename6, header=TRUE)
    data7      <-  read.table(filename7, header=TRUE)
    data8      <-  read.table(filename8, header=TRUE)

    # sm values for plotting invasion conditions
    sm  <-  seq(0,1,by=0.0001)

    # Color Scheme
    COLS  <-  c(transparentColor('seagreen3', opacity=0.2), 
                transparentColor('tomato2', opacity=0.2), 
                transparentColor('dodgerblue2', opacity=0.2),
                'black')

    # Set plot layout
    layout.mat <- matrix(c(1:8), nrow=2, ncol=4, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

############
## h = 1/2
## Panel A: C = 0

    # Calculate proportions of each outcome
    pAgree  <-  round(sum(data1$agree)/length(data1$agree), digits=3)
    pSim    <-  round(length(data1$sf[data1$simPoly == 1 & data1$agree == 0]) / 
                        length(data1$sf), digits=3)
    pEig    <-  round(length(data1$sf[data1$eigPoly == 1 & data1$agree == 0]) / 
                        length(data1$sf), digits=3)

    # Calcualte invasion conditions
    UB   <-  InvB(C = data1$C[1], delta = delta, hf = data1$hf[1], hm = data1$hm[1], sm = sm)
    UB[UB > 1]  <-  1.00000001
    LB   <-  InvA(C = data1$C[1], delta = delta, hf = data1$hf[1], hm = data1$hm[1], sm = sm)
    UBd  <-  InvB(C = data1$C[1], delta = delta, hf = data1$hf[1], hm = data1$hm[1], sm = sm)
    UBd[UBd > 1]  <-  1.00000001
    LBd  <-  InvA(C = data1$C[1], delta = delta, hf = data1$hf[1], hm = data1$hm[1], sm = sm)

    par(omi=c(0.75, 0.75, 0.75, 0.75), mar = c(2.5,2.5,0.5,0.5), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,sMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot 3 outcomes of comparisons between determ. simulations and eigenvalues
        points(sf[data1$agree == 1 & simPoly == 1 ] ~ sm[data1$agree == 1 & simPoly == 1 ], data=data1, lwd=3, col=COLS[1])
        points(sf[data1$simPoly == 1 & data1$agree == 0] ~ sm[data1$simPoly == 1 & data1$agree == 0], data=data1, lwd=3, col=COLS[2])
        points(sf[data1$eigPoly == 1 & data1$agree == 0] ~ sm[data1$eigPoly == 1 & data1$agree == 0], data=data1, lwd=3, col=COLS[3])
        # Plot SA invasion conditions 
        lines(UB[UB<=1] ~ sm[UB<= 1], lwd=3, col=COLS[4])
        lines(LB ~ sm, lwd=3, col=COLS[4])
        lines(UBd[UBd<=1] ~ sm[UBd<=1], lwd=3, col=COLS[4])
        lines(LBd ~ sm, lwd=3, col=COLS[4])
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.1, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.3,   expression(paste(italic(C)," = 0")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( -0.7,  0.5,   expression(paste(italic(h[f])," = ",italic(h[m]), " = 1/2")), cex=1.5, adj=c(0.5, 0.5), xpd=NA,srt=90)
        proportionalLabel( -0.4,  0.5,   expression(italic(s[f])), cex=1.3, adj=c(0.5, 0.5), xpd=NA,srt=90)
        points(0.02,0.99, pch=21, col=NA, cex=1, bg='seagreen3')
        points(0.02,0.91, pch=21, col=NA, cex=1, bg='tomato2')
        points(0.02,0.83, pch=21, col=NA, cex=1, bg='dodgerblue2')
        proportionalLabel(0.1, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.75, adj=c(0, 0.5), xpd=NA)


## Panel B: C = 1/4
    # Calculate proportions of each outcome
    pAgree  <-  round(sum(data2$agree)/length(data2$agree), digits=3)
    pSim    <-  round(length(data2$sf[data2$simPoly == 1 & data2$agree == 0]) / 
                        length(data2$sf), digits=3)
    pEig    <-  round(length(data2$sf[data2$eigPoly == 1 & data2$agree == 0]) / 
                        length(data2$sf), digits=3)

    # Calcualte invasion conditions
    UB   <-  InvB(C = data2$C[1], delta = delta, hf = data2$hf[1], hm = data2$hm[1], sm = sm)
    UB[UB > 1]  <-  1.00000001
    LB   <-  InvA(C = data2$C[1], delta = delta, hf = data2$hf[1], hm = data2$hm[1], sm = sm)
    UBd  <-  InvB(C = data2$C[1], delta = delta, hf = data2$hf[1], hm = data2$hm[1], sm = sm)
    UBd[UBd > 1]  <-  1.00000001
    LBd  <-  InvA(C = data2$C[1], delta = delta, hf = data2$hf[1], hm = data2$hm[1], sm = sm)

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,sMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot 3 outcomes of comparisons between determ. simulations and eigenvalues
        points(sf[data2$agree == 1 & simPoly == 1 ] ~ sm[data2$agree == 1 & simPoly == 1 ], data=data2, lwd=3, col=COLS[1])
        points(sf[data2$simPoly == 1 & data2$agree == 0] ~ sm[data2$simPoly == 1 & data2$agree == 0], data=data2, lwd=3, col=COLS[2])
        points(sf[data2$eigPoly == 1 & data2$agree == 0] ~ sm[data2$eigPoly == 1 & data2$agree == 0], data=data2, lwd=3, col=COLS[3])
        # Plot SA invasion conditions 
        lines(UB[UB<=1] ~ sm[UB<= 1], lwd=3, col=COLS[4])
        lines(LB ~ sm, lwd=3, col=COLS[4])
        lines(UBd[UBd<=1] ~ sm[UBd<=1], lwd=3, col=COLS[4])
        lines(LBd ~ sm, lwd=3, col=COLS[4])
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.05, 1.1, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5,  1.3,   expression(paste(italic(C)," = 0.25")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        points(0.02,0.99, pch=21, col=NA, cex=1, bg='seagreen3')
        points(0.02,0.91, pch=21, col=NA, cex=1, bg='tomato2')
        points(0.02,0.83, pch=21, col=NA, cex=1, bg='dodgerblue2')
        proportionalLabel(0.1, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.75, adj=c(0, 0.5), xpd=NA)

## Panel C: C = 1/2
    # Calculate proportions of each outcome
    pAgree  <-  round(sum(data3$agree)/length(data3$agree), digits=3)
    pSim    <-  round(length(data3$sf[data3$simPoly == 1 & data3$agree == 0]) / 
                        length(data3$sf), digits=3)
    pEig    <-  round(length(data3$sf[data3$eigPoly == 1 & data3$agree == 0]) / 
                        length(data3$sf), digits=3)

    # Calcualte invasion conditions
    UB   <-  InvB(C = data3$C[1], delta = delta, hf = data3$hf[1], hm = data3$hm[1], sm = sm)
    UB[UB > 1]  <-  1.00000001
    LB   <-  InvA(C = data3$C[1], delta = delta, hf = data3$hf[1], hm = data3$hm[1], sm = sm)
    UBd  <-  InvB(C = data3$C[1], delta = delta, hf = data3$hf[1], hm = data3$hm[1], sm = sm)
    UBd[UBd > 1]  <-  1.00000001
    LBd  <-  InvA(C = data3$C[1], delta = delta, hf = data3$hf[1], hm = data3$hm[1], sm = sm)

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,sMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot 3 outcomes of comparisons between determ. simulations and eigenvalues
        points(sf[data3$agree == 1 & simPoly == 1 ] ~ sm[data3$agree == 1 & simPoly == 1 ], data=data3, lwd=3, col=COLS[1])
        points(sf[data3$simPoly == 1 & data3$agree == 0] ~ sm[data3$simPoly == 1 & data3$agree == 0], data=data3, lwd=3, col=COLS[2])
        points(sf[data3$eigPoly == 1 & data3$agree == 0] ~ sm[data3$eigPoly == 1 & data3$agree == 0], data=data3, lwd=3, col=COLS[3])
        # Plot SA invasion conditions 
        lines(UB[UB<=1] ~ sm[UB<= 1], lwd=3, col=COLS[4])
        lines(LB ~ sm, lwd=3, col=COLS[4])
        lines(UBd[UBd<=1] ~ sm[UBd<=1], lwd=3, col=COLS[4])
        lines(LBd ~ sm, lwd=3, col=COLS[4])
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.1, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.3,   expression(paste(italic(C)," = 0.5")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        points(0.02,0.99, pch=21, col=NA, cex=1, bg='seagreen3')
        points(0.02,0.91, pch=21, col=NA, cex=1, bg='tomato2')
        points(0.02,0.83, pch=21, col=NA, cex=1, bg='dodgerblue2')
        proportionalLabel(0.1, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.75, adj=c(0, 0.5), xpd=NA)

## Panel D: C = 3/4
    # Calculate proportions of each outcome
    pAgree  <-  round(sum(data4$agree)/length(data4$agree), digits=3)
    pSim    <-  round(length(data4$sf[data4$simPoly == 1 & data4$agree == 0]) / 
                        length(data4$sf), digits=3)
    pEig    <-  round(length(data4$sf[data4$eigPoly == 1 & data4$agree == 0]) / 
                        length(data4$sf), digits=3)

    # Calcualte invasion conditions
    UB   <-  InvB(C = data4$C[1], delta = delta, hf = data4$hf[1], hm = data4$hm[1], sm = sm)
    UB[UB > 1]  <-  1.00000001
    LB   <-  InvA(C = data4$C[1], delta = delta, hf = data4$hf[1], hm = data4$hm[1], sm = sm)
    UBd  <-  InvB(C = data4$C[1], delta = delta, hf = data4$hf[1], hm = data4$hm[1], sm = sm)
    UBd[UBd > 1]  <-  1.00000001
    LBd  <-  InvA(C = data4$C[1], delta = delta, hf = data4$hf[1], hm = data4$hm[1], sm = sm)

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,sMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot 3 outcomes of comparisons between determ. simulations and eigenvalues
        points(sf[data4$agree == 1 & simPoly == 1 ] ~ sm[data4$agree == 1 & simPoly == 1 ], data=data4, lwd=3, col=COLS[1])
        points(sf[data4$simPoly == 1 & data4$agree == 0] ~ sm[data4$simPoly == 1 & data4$agree == 0], data=data4, lwd=3, col=COLS[2])
        points(sf[data4$eigPoly == 1 & data4$agree == 0] ~ sm[data4$eigPoly == 1 & data4$agree == 0], data=data4, lwd=3, col=COLS[3])
        # Plot SA invasion conditions 
        lines(UB[UB<=1] ~ sm[UB<= 1], lwd=3, col=COLS[4])
        lines(LB ~ sm, lwd=3, col=COLS[4])
        lines(UBd[UBd<=1] ~ sm[UBd<=1], lwd=3, col=COLS[4])
        lines(LBd ~ sm, lwd=3, col=COLS[4])
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.1, expression(paste(bold(D))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.3,   expression(paste(italic(C)," = 0.75")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        points(0.02,0.99, pch=21, col=NA, cex=1, bg='seagreen3')
        points(0.02,0.91, pch=21, col=NA, cex=1, bg='tomato2')
        points(0.02,0.83, pch=21, col=NA, cex=1, bg='dodgerblue2')
        proportionalLabel(0.1, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.75, adj=c(0, 0.5), xpd=NA)


## h = 1/4
## Panel E: C = 0
    # Calculate proportions of each outcome
    pAgree  <-  round(sum(data5$agree)/length(data5$agree), digits=3)
    pSim    <-  round(length(data5$sf[data5$simPoly == 1 & data5$agree == 0]) / 
                        length(data5$sf), digits=3)
    pEig    <-  round(length(data5$sf[data5$eigPoly == 1 & data5$agree == 0]) / 
                        length(data5$sf), digits=3)

    # Calcualte invasion conditions
    UB   <-  InvB(C = data5$C[1], delta = delta, hf = data5$hf[1], hm = data5$hm[1], sm = sm)
    UB[UB > 1]  <-  1.00000001
    LB   <-  InvA(C = data5$C[1], delta = delta, hf = data5$hf[1], hm = data5$hm[1], sm = sm)
    UBd  <-  InvB(C = data5$C[1], delta = delta, hf = data5$hf[1], hm = data5$hm[1], sm = sm)
    UBd[UBd > 1]  <-  1.00000001
    LBd  <-  InvA(C = data5$C[1], delta = delta, hf = data5$hf[1], hm = data5$hm[1], sm = sm)

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,sMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot 3 outcomes of comparisons between determ. simulations and eigenvalues
        points(sf[data5$agree == 1 & simPoly == 1 ] ~ sm[data5$agree == 1 & simPoly == 1 ], data=data5, lwd=3, col=COLS[1])
        points(sf[data5$simPoly == 1 & data5$agree == 0] ~ sm[data5$simPoly == 1 & data5$agree == 0], data=data5, lwd=3, col=COLS[2])
        points(sf[data5$eigPoly == 1 & data5$agree == 0] ~ sm[data5$eigPoly == 1 & data5$agree == 0], data=data5, lwd=3, col=COLS[3])
        # Plot SA invasion conditions 
        lines(UB[UB<=1] ~ sm[UB<= 1], lwd=3, col=COLS[4])
        lines(LB ~ sm, lwd=3, col=COLS[4])
        lines(UBd[UBd<=1] ~ sm[UBd<=1], lwd=3, col=COLS[4])
        lines(LBd ~ sm, lwd=3, col=COLS[4])
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.1, expression(paste(bold(E))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( -0.7,  0.5,   expression(paste(italic(h[f])," = ",italic(h[m]), " = 1/4")), cex=1.5, adj=c(0.5, 0.5), xpd=NA,srt=90)
        proportionalLabel( -0.4,  0.5,   expression(italic(s[f])), cex=1.3, adj=c(0.5, 0.5), xpd=NA,srt=90)
        proportionalLabel(0.5,  -0.4,   expression(italic(s[m])), cex=1.3, adj=c(0.5, 0.5), xpd=NA)        
        points(0.02,0.99, pch=21, col=NA, cex=1, bg='seagreen3')
        points(0.02,0.91, pch=21, col=NA, cex=1, bg='tomato2')
        points(0.02,0.83, pch=21, col=NA, cex=1, bg='dodgerblue2')
        proportionalLabel(0.1, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.75, adj=c(0, 0.5), xpd=NA)

## Panel F: C = 1/4    
    # Calculate proportions of each outcome
    pAgree  <-  round(sum(data6$agree)/length(data6$agree), digits=3)
    pSim    <-  round(length(data6$sf[data6$simPoly == 1 & data6$agree == 0]) / 
                        length(data6$sf), digits=3)
    pEig    <-  round(length(data6$sf[data6$eigPoly == 1 & data6$agree == 0]) / 
                        length(data6$sf), digits=3)

    # Calcualte invasion conditions
    UB   <-  InvB(C = data6$C[1], delta = delta, hf = data6$hf[1], hm = data6$hm[1], sm = sm)
    UB[UB > 1]  <-  1.00000001
    LB   <-  InvA(C = data6$C[1], delta = delta, hf = data6$hf[1], hm = data6$hm[1], sm = sm)
    UBd  <-  InvB(C = data6$C[1], delta = delta, hf = data6$hf[1], hm = data6$hm[1], sm = sm)
    UBd[UBd > 1]  <-  1.00000001
    LBd  <-  InvA(C = data6$C[1], delta = delta, hf = data6$hf[1], hm = data6$hm[1], sm = sm)

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,sMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot 3 outcomes of comparisons between determ. simulations and eigenvalues
        points(sf[data6$agree == 1 & simPoly == 1 ] ~ sm[data6$agree == 1 & simPoly == 1 ], data=data6, lwd=3, col=COLS[1])
        points(sf[data6$simPoly == 1 & data6$agree == 0] ~ sm[data6$simPoly == 1 & data6$agree == 0], data=data6, lwd=3, col=COLS[2])
        points(sf[data6$eigPoly == 1 & data6$agree == 0] ~ sm[data6$eigPoly == 1 & data6$agree == 0], data=data6, lwd=3, col=COLS[3])
        # Plot SA invasion conditions 
        lines(UB[UB<=1] ~ sm[UB<= 1], lwd=3, col=COLS[4])
        lines(LB ~ sm, lwd=3, col=COLS[4])
        lines(UBd[UBd<=1] ~ sm[UBd<=1], lwd=3, col=COLS[4])
        lines(LBd ~ sm, lwd=3, col=COLS[4])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.1, expression(paste(bold(F))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5,  -0.4,   expression(italic(s[m])), cex=1.3, adj=c(0.5, 0.5), xpd=NA)        
        points(0.02,0.99, pch=21, col=NA, cex=1, bg='seagreen3')
        points(0.02,0.91, pch=21, col=NA, cex=1, bg='tomato2')
        points(0.02,0.83, pch=21, col=NA, cex=1, bg='dodgerblue2')
        proportionalLabel(0.1, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.75, adj=c(0, 0.5), xpd=NA)

## Panel G: C = 1/2
    # Calculate proportions of each outcome
    pAgree  <-  round(sum(data7$agree)/length(data7$agree), digits=3)
    pSim    <-  round(length(data7$sf[data7$simPoly == 1 & data7$agree == 0]) / 
                        length(data7$sf), digits=3)
    pEig    <-  round(length(data7$sf[data7$eigPoly == 1 & data7$agree == 0]) / 
                        length(data7$sf), digits=3)

    # Calcualte invasion conditions
    UB   <-  InvB(C = data7$C[1], delta = delta, hf = data7$hf[1], hm = data7$hm[1], sm = sm)
    UB[UB > 1]  <-  1.00000001
    LB   <-  InvA(C = data7$C[1], delta = delta, hf = data7$hf[1], hm = data7$hm[1], sm = sm)
    UBd  <-  InvB(C = data7$C[1], delta = delta, hf = data7$hf[1], hm = data7$hm[1], sm = sm)
    UBd[UBd > 1]  <-  1.00000001
    LBd  <-  InvA(C = data7$C[1], delta = delta, hf = data7$hf[1], hm = data7$hm[1], sm = sm)

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,sMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot 3 outcomes of comparisons between determ. simulations and eigenvalues
        points(sf[data7$agree == 1 & simPoly == 1 ] ~ sm[data7$agree == 1 & simPoly == 1 ], data=data7, lwd=3, col=COLS[1])
        points(sf[data7$simPoly == 1 & data7$agree == 0] ~ sm[data7$simPoly == 1 & data7$agree == 0], data=data7, lwd=3, col=COLS[2])
        points(sf[data7$eigPoly == 1 & data7$agree == 0] ~ sm[data7$eigPoly == 1 & data7$agree == 0], data=data7, lwd=3, col=COLS[3])
        # Plot SA invasion conditions 
        lines(UB[UB<=1] ~ sm[UB<= 1], lwd=3, col=COLS[4])
        lines(LB ~ sm, lwd=3, col=COLS[4])
        lines(UBd[UBd<=1] ~ sm[UBd<=1], lwd=3, col=COLS[4])
        lines(LBd ~ sm, lwd=3, col=COLS[4])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5,  -0.4,   expression(italic(s[m])), cex=1.3, adj=c(0.5, 0.5), xpd=NA)        
        points(0.02,0.99, pch=21, col=NA, cex=1, bg='seagreen3')
        points(0.02,0.91, pch=21, col=NA, cex=1, bg='tomato2')
        points(0.02,0.83, pch=21, col=NA, cex=1, bg='dodgerblue2')
        proportionalLabel(0.1, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.75, adj=c(0, 0.5), xpd=NA)

## Panel H: C = 3/4
    # Calculate proportions of each outcome
    pAgree  <-  round(sum(data8$agree)/length(data8$agree), digits=3)
    pSim    <-  round(length(data8$sf[data8$simPoly == 1 & data8$agree == 0]) / 
                        length(data8$sf), digits=3)
    pEig    <-  round(length(data8$sf[data8$eigPoly == 1 & data8$agree == 0]) / 
                        length(data8$sf), digits=3)

    # Calcualte invasion conditions
    UB   <-  InvB(C = data8$C[1], delta = delta, hf = data8$hf[1], hm = data8$hm[1], sm = sm)
    UB[UB > 1]  <-  1.00000001
    LB   <-  InvA(C = data8$C[1], delta = delta, hf = data8$hf[1], hm = data8$hm[1], sm = sm)
    UBd  <-  InvB(C = data8$C[1], delta = delta, hf = data8$hf[1], hm = data8$hm[1], sm = sm)
    UBd[UBd > 1]  <-  1.00000001
    LBd  <-  InvA(C = data8$C[1], delta = delta, hf = data8$hf[1], hm = data8$hm[1], sm = sm)

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0, sMax), ylim = c(0,sMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot 3 outcomes of comparisons between determ. simulations and eigenvalues
        points(sf[data8$agree == 1 & simPoly == 1 ] ~ sm[data8$agree == 1 & simPoly == 1 ], data=data8, lwd=3, col=COLS[1])
        points(sf[data8$simPoly == 1 & data8$agree == 0] ~ sm[data8$simPoly == 1 & data8$agree == 0], data=data8, lwd=3, col=COLS[2])
        points(sf[data8$eigPoly == 1 & data8$agree == 0] ~ sm[data8$eigPoly == 1 & data8$agree == 0], data=data8, lwd=3, col=COLS[3])
        # Plot SA invasion conditions 
        lines(UB[UB<=1] ~ sm[UB<= 1], lwd=3, col=COLS[4])
        lines(LB ~ sm, lwd=3, col=COLS[4])
        lines(UBd[UBd<=1] ~ sm[UBd<=1], lwd=3, col=COLS[4])
        lines(LBd ~ sm, lwd=3, col=COLS[4])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5,  -0.4,   expression(italic(s[m])), cex=1.3, adj=c(0.5, 0.5), xpd=NA)        
        points(0.02,0.99, pch=21, col=NA, cex=1, bg='seagreen3')
        points(0.02,0.91, pch=21, col=NA, cex=1, bg='tomato2')
        points(0.02,0.83, pch=21, col=NA, cex=1, bg='dodgerblue2')
        proportionalLabel(0.1, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.75, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.75, adj=c(0, 0.5), xpd=NA)
              
}