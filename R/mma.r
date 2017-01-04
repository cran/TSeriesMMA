
#' Multiscale Multifractal Analysis of Time Series Data
#'
#' @param smin Minimal s scale used, when calculating Fq(s) functions family (default 10)
#' @param smax Maximal s scale used, when calculating Fq(s) functions family, has to be multiple of 5 (default 600; in general should be near to N/50, where N is a time series length)
#' @param qmin Minimal multifractal parameter q used (default -5)
#' @param qmax Maximal cmultifractal parameter q used (deafault 5)
#' @param data Time series data
#' @param col The color variation of the plot
#' @param theta Angle of view
#' @param phi Second angle of view
#' @examples \dontrun{
#' mma(smax=30, data=timeSeriesData)
#' }
#' @export
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics persp
#' @importFrom stats lm na.omit

mma<-function(smin=10,smax=600,qmin=-5,qmax=5,data,col='V1',theta=-45,phi=25)
{
  qlist <- seq(from=qmin, to=qmax, by=0.1)
  # We can't have q value = 0 because we take power of 1/q at one point
  for(i in 1:length(qlist)){
    if(qlist[i]==0){
      qlist[i]<-0.001
    }
  }
  # Read data from file - numbers alone, one per line in text file

  mean=mean(data[[col]])
  #mean
  # Find profile of the signal by calculating cumulative sum for the whole series
  prof = apply((data[1][1]-mean),2,cumsum)
  data[1][1]
  slength = length(prof)

  # fqs below is the matrix containing Fluctuations.
  # We store it in matrix form as each pair of values of q (multifractal parameter) and s (scale) has a particular value of fluctuation. These are all stored together.
  # We need to set the size of the fqs matrix beforehand to add rows, hence we calculate number of iterations and assign that
  iter <- length(qlist*(smax-smin))
  fqs <- matrix(NA, nrow=iter, ncol=3)


  i<-1

  for (s in smin:smax) {
    #Reshape the profiled signal into a matrix, so that we can extract segments
    ind <- c(1:length(prof))
    coordinates <- t(matrix(ind[1:(length(prof)-(length(prof)%%s))],s,(length(prof)-(length(prof)%%s))/s))
    segments<-matrix(prof, nrow(coordinates), ncol(coordinates), byrow=T)
    #Set of bases upto s
    xbase<-seq(1,s,1)
    f2nis<-c()
    for(ni in 1:nrow(segments)){
      # Obtain fits for all bases and calculate the mean variance for the same
      seg<-segments[ni,]
      fit<-lm(seg ~ xbase + I(xbase^2))
      sum<-0
      for(i in 1:length(xbase)){
        b<-(fit[[1]][[1]]*1+fit[[1]][[2]]*(xbase[i]^1)+fit[[1]][[3]]*(xbase[i]^2))
        sum<-sum+((seg[i]-b)^2)
      }
      variance<-(sum/length(xbase))
      f2nis<-c(f2nis, variance)
    }

    for(q in qlist){
      # print(i)
      # Calculating fluctuations and storing
      row<-c(q,s,mean((f2nis^(q/2)))^(1/q))
      fqs<-rbind(fqs, row)
      i<-i+1
      }
  }
  # Remove redundant rows
  fqs=na.omit(fqs)

  # Calculate log values of scales and fluctuations
  fqsll <- cbind(fqs[,1], fqs[,2], log(fqs[,2]), log(fqs[,3]));

  #We can increase the value of spacing and qlist to get more data points in the surface but the plot will be less smooth and calculations will be lengthy
  sspacing <- ((smax/5)-smin)/11
  qlist <- seq(from=qmin, to=qmax, by=1)
  for(i in 1:length(qlist)){
    if(qlist[i]==0){
      qlist[i]<-0.001
    }
  }

  # Calculating hurst exponent to be stored in hqs, for each unique pair of (q,s)
  # This will be used to make a 3D Hurst Surface
  hqRows<-(((smax/5)-smin)/sspacing)*length(qlist)
  hqs <- matrix(NA, nrow=hqRows, ncol=3)
  for(sit in seq(smin, (smax/5), sspacing)){
    for(qit in qlist){
      fittemp <- fqsll[fqsll[,1] == qit & fqsll[,2] >= sit & fqsll[,2] <= 5*sit,]
      htemp<-lm(fittemp[,4]~fittemp[,3])
      row<-c(qit, 3*sit, htemp[[1]][[2]])
      hqs<-rbind(hqs, row)
    }
  }

  # Remove redundant rows again
  hqs<-na.omit(hqs)
  # Reshape hqs values to a matrix so that we can plot it as a surface
  hplot<-t(matrix(hqs[,3],nrow=length(hqs[,3])/length(qlist),byrow=T))

  # Generate the desired number of colors from this palette
  jet.colors <- colorRampPalette( c( "#2c00cc", "#cc0000") )
  nbcol <- 100
  color <- jet.colors(nbcol)

  # Compute the z-value at the facet centres
  zfacet <- hplot[-1, -1] + hplot[-1, -ncol(hplot)] + hplot[-nrow(hplot), -1] + hplot[-nrow(hplot), -ncol(hplot)]
  # Recode facet z-values into color indices
  facetcol <- cut(zfacet, nbcol)
  persp(c(unique(hqs[,1])),c(unique(hqs[,2])), hplot, zlim=c(0,2.5),col = color[facetcol], xlab="Multifractal Parameters(Q)", ylab="Scales(S)", zlab = "Hurst Exponent", phi=phi, theta=theta, ticktype="detailed")
}
#signal <- read.table("data.txt")
#mma(smin=10,qmin=-5,qmax=5,data=signal)

