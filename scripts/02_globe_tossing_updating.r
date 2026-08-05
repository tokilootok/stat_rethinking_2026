# animate globe tossing model

library(rethinking)
library(animation)
library(globe)
library(sf)
library(spData)

# globe dev code
if ( FALSE ) {

    library(globe)
    globeearth(eye=list(runif(1,-80,80),runif(1,-180,180)))

    # https://gis.stackexchange.com/questions/75033/given-lat-and-lon-identify-if-point-over-land-or-ocean-using-r
    # install.packages(c('sf','spData'))
    library(sf)
    library(spData) ## For `world`, an sf MULTIPOLYGON object

    ## Create an sf POINTS object
    set.seed(0)
    lat <- runif(10, -70, 70)
    lon <- runif(10, -180, 180)
    #points <- expand.grid(lon, lat)  # Note that I reversed OP's ordering of lat/long
    points <- data.frame(lon,lat)
    pts <- st_as_sf(points, coords=1:2, crs=4326)

    ## Find which points fall over land
    ii <- !is.na(as.numeric(st_intersects(pts, world)))

    # function to sample from globe and return whether water or land
    # returns TRUE for water
    globe_toss <- function(plot=FALSE) {
        lat <- runif(1,-70,70)
        lon <- runif(1,-180,180)
        pts <- st_as_sf( data.frame(lon,lat) , coords=1:2, crs=4326)
        if ( plot==TRUE ) globeearth(eye=list(lat,lon))
        return( c(lat,lon, is.na(as.numeric(st_intersects(pts, world))) ) )
    }

    sum( replicate(1e2,globe_toss()[3]) )

}

# routes
# c(lon,lat) format for each point
# https://www.r-graph-gallery.com/how-to-draw-connecting-routes-on-map-with-r-and-great-circles.html
make_route <- function( start , end , steps=50 ) {
    library(geosphere)
    route <- gcIntermediate(start, end, n=steps, addStartEnd=TRUE, breakAtDateLine=F)
    return(route)
}

myglobeearth <- function (gdata, runlen, eye, top, ..., do.plot = TRUE) 
{
    xylim <- 0.8
    if (missing(gdata)) {
        gdata <- globe::earth$coords
    }
    if (missing(runlen)) {
        runlen <- globe::earth$runlen
    }
    spos <- spatialpos(gdata[, 1], gdata[, 2])
    mpos <- orthogproj(eye = eye, top = top, spos)
    if (do.plot) 
        plot(c(-xylim, xylim), c(-xylim, xylim), type = "n", asp = 1, axes = FALSE, 
            xlab = "", ylab = "")
    x <- mpos[, 1]
    y <- mpos[, 2]
    ok <- (mpos[, 3] < 0)
    runlen <- runlen[runlen != 0]
    breaks <- cumsum(runlen)
    ok[breaks] <- FALSE
    s <- seq(x)[ok]
    par(xpd=NA)
    if (do.plot) {
        segments(x[s], y[s], x[s + 1], y[s + 1], ...)
        a <- seq(0, 2 * pi, length = 360)
        lines(cos(a), sin(a), lty = 2 , lwd=1 , col=1 )
    }
    result <- cbind(x[s], y[s], x[s + 1], y[s + 1])
    attr(result, "piece") <- as.integer(factor(cumsum(!ok)[ok]))
    return(invisible(result))
}

globe_toss <- function(plot=FALSE) {
    lat <- runif(1,-70,70)
    lon <- runif(1,-180,180)
    pts <- st_as_sf( data.frame(lon,lat) , coords=1:2, crs=4326)
    if ( plot==TRUE ) myglobeearth(eye=list(lon,lat))
    return( c(lon,lat, is.na(as.numeric(st_intersects(pts, world))) ) )
}

# globe_toss(TRUE); points(0,0,lwd=2,col=2)

blank(bty="n",h=0.8,w=2)
par(mfrow=c(1,2))

n <- 10
npts <- 100
mpts <- 10 # interpolation points

x <- seq( from=0 , to=1 , length.out=npts )

pf <- function(maxy) {
    plot(NULL,xlim=c(0,1),ylim=c(0,maxy),xlab="proportion water",ylab="density",xaxt="n")
    axis(1,at=c(0,0.5,1),labels=c(0,0.5,1))
}

f <- function(x,a,b) dbeta(x,a,b)

draw_samples <- function(samples) {
    mtext( paste(samples,collapse=" ") , adj=0 )
}

draw_posterior_panel <- function(maxy, samples, x, a_old, b_old, a_new = NULL, b_new = NULL, col = NULL) {
    pf(maxy)
    if ( length(samples) > 0 ) draw_samples(samples)
    lines( x , f(x,a_old,b_old) , lty=2 , lwd=2 , col=grau(0.5) )
    if ( !is.null(a_new) && !is.null(b_new) ) {
        lines( x , f(x,a_new,b_new) , lwd=5 , col=col )
    }
}

record_frame <- function(nf, do_save) {
    ani.record()
    if ( do_save ) {
        quartz.save( sprintf("frames/frame_%04d.png", 1000 + nf), type="png", device=dev.cur(), dpi=200 )
    }
    nf + 1
}

label_from_toss <- function(toss) if ( toss == 1 ) "W" else "L"
color_from_toss <- function(toss) if ( toss == 1 ) 4 else 2

kcols <- c( col.alpha(2,0.7) , col.alpha(4,0.7) )
gcol <- 1
glwd <- 0.5

# set.seed(1619)

# do all tosses first, so we can get max density for plot
tm <- t( replicate(n, globe_toss(FALSE)) )

# compute max density
a <- 1
b <- 1
samples <- c()
amax <- a + sum(tm[,3])
bmax <- b + sum(1-tm[,3])
dmax <- max( f(x,amax,bmax) )
lonlat <- c(0,0) # init globe position
do_save <- TRUE
nf <- 1

ani.record(reset = TRUE)  # clear history before recording
par(mfrow=c(1,2))

for ( i in seq_len(n) ) {

    # sample an observation
    # 1 water
    # 0 land
    toss_list <- tm[i,]
    toss <- toss_list[3]
    label <- label_from_toss(toss)
    toss_col <- color_from_toss(toss)
    toss_kcol <- ifelse(toss==1,kcols[2],kcols[1])
    samples_old <- samples
    samples <- c( samples , label )

    # update
    # new posterior is beta(a+water,b+land)
    a_old <- a
    b_old <- b
    a <- a + toss
    b <- b + (1-toss)

    # interpolate posterior
    aseq <- seq( from=a_old , to=a , length.out=mpts )
    bseq <- seq( from=b_old , to=b , length.out=mpts )

    # interpolate globe eye
    lonlat_old <- lonlat
    lonlat <- toss_list[1:2]
    route <- make_route( lonlat_old[1:2] , lonlat[1:2] , steps=mpts )

    # spin globe
    for ( j in seq_len(nrow(route)) ) {
        myglobeearth( eye=list( route[j,1] , route[j,2] ) , col=gcol , lwd=glwd )
        draw_posterior_panel(dmax+0.1, samples_old, x, a_old, b_old)
        nf <- record_frame(nf, do_save)
    }

    # update posterior
    print(toss_list)
    for ( j in seq_len(mpts) ) {

        # plot globe
        myglobeearth( eye=list( toss_list[1] , toss_list[2] ) , col=gcol , lwd=glwd )
        
        text(0,0,label,cex=0.8,col=toss_col,font=2)
        draw_posterior_panel(dmax+0.1, samples, x, a_old, b_old, aseq[j], bseq[j], toss_kcol)
        nf <- record_frame(nf, do_save)
    }
    
}#i

## now we can replay it, with an appropriate pause between
## frames
oopts = ani.options(interval = 0.05)
ani.replay()

# convert -alpha remove -background white -delay 10 -loop 0 frame*.png globe_tossing.gif

# convert -delay 0 animated.gif animated2.gif

#######
# globe posterior with some intervals on it

blank(h=0.8,bty="n")

z <- curve( dbeta(x,2,4) , lwd=4 , lty=1 , col="gray" , xlab="proportion water" , xaxt="n" , ylab="density" )
axis(1,at=c(0,0.5,1),labels=c(0,0.5,1))

p <- rbeta(1e4,2,4)
i <- PI(p,0.99)
lo <- z$x[ which.max(z$x > i[1]) ]
hi <- z$x[ which.min(z$x < i[2]) ]
polygon(c(z$x[z$x >= lo & z$x <= hi], hi, lo), c(z$y[z$x >= lo & z$x <= hi], 0, 0), col=2 , border=1 )
