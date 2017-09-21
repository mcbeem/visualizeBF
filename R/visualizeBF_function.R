#' Function for exploring Bayes Factors.
#'
#' @param data A vector of data
#' @param pointsize Width of points for numerical integration, defaults to .001.
#' @param scale Scale parameter for Cauchy prior, defaults to .707
#' @param plot The type of plot to produce. 1=BF as weighted likelihood, 2=BF as comparison of densities for mu=0 under prior and posterior distributions.
#' @export
#' @return Returns a list with the following components:  
#' \itemize{
#'   \item \code{L.m0}    The likelihood evaluated at mu=0; the L of the null model (m0) \cr
#'   \item \code{L.m1}    The likelihood integrated over the whole space weighted by m1 prior; L of the alternative model (m1). \cr
#'   \item \code{prior.mu0} The density of the prior distribution evaluated at mu=0 (m1). \cr
#'   \item \code{posterior.m0} The density the posterior distribution evaluated at mu=0 (m0). \cr
#'   \item \code{BF10}    The Bayes Factor for alternative hypothesis relative to null.  \cr 
#'   \item \code{BF01}    The Bayes Factor for the null hypothesis relative to the alternative \cr 
#'   \item \code{figure}  Object containing the plot \cr
#'  }
#' @examples
#' set.seed(1)
#' data <- rnorm(n=50, mean=.3, sd=1)
#' visualizeBF(data, plot=1)
#' visualizeBF(data, plot=2)
#' 
#' a <- visualizeBF(data, plot=1)
#' 
#' # extract density of the prior distribution at mu=0
#' a$prior.mu0
#' 
#' # calculate BF10 as likelihood ratio
#' a$L.m1 / a$L.m0
#' 
#' # alternative method for calculating BF10 as density ratio
#' a$prior.mu0 / a$posterior.mu0
#' 
#' # view the plot
#' a$figure

visualizeBF <- function(data, pointsize=0.001, scale=.707, plot=1) {
    
  # do some checks
  if (!(plot %in% c(1, 2))) {stop("plot must be set to 1 or 2, see ?visualizeBF for details")}
  
  
  # standardize data so sample sd is exactly 1
  data <- data / sd(data)
  
  # calculate likelihood function for the data
  # there are two parameters
  
  pointsize <- .001
  
  b0s <- seq(-3,3, pointsize)
  
  # calculate LL for each parm value
  LL <- apply(apply(sapply(b0s, dnorm, x=data, sd=1), c(1,2), log), 2, sum)
  L <- apply(sapply(b0s, dnorm, x=data, sd=1), 2, prod)
  
  # prior for b0 based on Cauchy distribution
  prior.b0 <- dcauchy((b0s), scale=scale)
  
  # normalize it
  prior.b0 <- prior.b0 / (sum(prior.b0)*pointsize)
  log.prior.b0 <- log(prior.b0)
  
  # which b0 is zero?
  zeroloc <- which(b0s==0)
  
  #### This is calculating BF as the weighted average of the likelihood
  # numerator of Bayes theorem
  m0 <- exp(LL[zeroloc])
  m1 <- sum(exp(LL+log.prior.b0)) / sum(exp(log.prior.b0))
  m1/m0
  
  ### This is calculating BF as the ratio of prior and posterior at d=0 ###
  
  # posterior distribution
  posterior.b0 <- exp((LL+log.prior.b0)) / (sum(exp(LL+log.prior.b0))*pointsize)
  prior.b0[zeroloc] / posterior.b0[zeroloc]
  
  plot(b0s, posterior.b0, "l", col="blue")
  points(b0s, prior.b0, "l", col="red")
  abline(v=0)
  
  BF10 <- prior.b0[zeroloc] / posterior.b0[zeroloc]
  BF10
  BF10^-1
  
  m1/m0
  
  
  ### Make a fancy 2-panel plot
  plot.new()
  font.scale=.8
  par(mfrow=c(2,1), mar=rep(2,4), mgp=c(1.1,.3,0), cex=font.scale)

  
  if (plot==1) {
  
    # top panel: likelihood for full range of d
    plot(b0s, exp(LL), 'l', 
         main="", 
         xlab="", ylab="", cex.main=1)
    title(expression(bold(paste("The BF is the ratio of the weighted average likelihood to the likelihood at ", mu, " = 0", sep=""))), 
          line=1.3, cex=font.scale)
    mtext("Likelihood is black, weighted average likelihood is red", cex=font.scale)
    title(xlab=expression(mu), ylab="Likelihood")
    points(x=b0s[zeroloc], y=(m0), cex=1.5, col="red")
    points(x=b0s[zeroloc], y=(m1), cex=1.5, col="red")
    abline(h=(m1), col="red")
    if (m1 > m0) {
      points(x=c(0,0), y=c((m0), (m1)) , type='l', col="red", lwd=2)
      points(x=c(0,0), y=c(0, (m0)) , type='l', col="blue", lwd=2)
    }
  
    if (m0 > m1) {
      points(x=c(0, 0), y=c(0, (m0)) , type='l', col="blue", lwd=2)
      points(x=c(0, 0), y=c(0, (m1)) , type='l', col="red", lwd=2)
    }
    
    # bottom panel: zoomed in version of top panel
    plot(b0s, exp(LL), 'l', xlim=c(-.1, .1), ylim=c(0, max(m1,m0)*1.1), 
         main="Zoomed-in version of upper panel", 
         cex.main=1, cex.axis=font.scale, xlab="", ylab="")
    title(xlab=expression(mu), ylab="Likelihood")
    abline(h=0, col="darkgray")
    points(x=0, y=(m0), cex=1.5, col="red")
    points(x=0, y=(m1), cex=1.5, col="red")
    abline(h=(m1), col="red")
    
    if (m1 > m0) {
      points(x=c(0,0), y=c(0, (m1)) , type='l', col="red", lwd=2)
      points(x=c(0,0), y=c(0, (m0)) , type='l', col="blue", lwd=2)
    }
    
    if (m0 > m1) {
      points(x=c(0, 0), y=c(0, (m0)) , type='l', col="blue", lwd=2)
      points(x=c(0, 0), y=c(0, (m1)) , type='l', col="red", lwd=2)
    }
    
    figure <- recordPlot()
  }
  
  if (plot==2) {
  
    # top panel: prior for d
    ymax <- max(c(prior.b0, posterior.b0))
    
    plot(b0s, prior.b0, 'l', ylim=c(0,ymax), col="red", 
         main="", 
         cex.main=font.scale, cex.axis=font.scale, xlab="", ylab="")
    mtext("Prior is red, posterior is blue", cex=font.scale)
    title(expression(bold(paste("The BF is the ratio of the prior density at ", mu, " = 0 to the posterior density at ", mu, " = 0", sep=""))), 
          line=1.3, cex=font.scale)
    title(xlab=expression(mu), ylab="Density")
    abline(h=0, col="darkgray")
    points(b0s, posterior.b0, 'l', col="blue")
    
    points(x=0, y=prior.b0[zeroloc], cex=1.5, col="red")
    points(x=0, y=posterior.b0[zeroloc], cex=1.5, col="red")
    
    if (m1 > m0) {
      points(x=c(0,0), y=c(0, prior.b0[zeroloc]) , type='l', col="red", lwd=2)
      points(x=c(0,0), y=c(0, posterior.b0[zeroloc]) , type='l', col="blue", lwd=2)
    }
    
    if (m0 > m1) {
      points(x=c(0,0), y=c(0, posterior.b0[zeroloc]) , type='l', col="blue", lwd=2)
      points(x=c(0,0), y=c(0, prior.b0[zeroloc]) , type='l', col="red", lwd=2)
    }
  
    # bottom panel
    plot(b0s, prior.b0, 'l', xlim=c(-.1, .1), ylim=c(0,ymax), col="red", 
         main="(Zoomed-in version of upper panel)", 
         cex.axis=font.scale, xlab="", ylab="")
    title(xlab=expression(mu), ylab="Density")
    abline(h=0, col="darkgray")
    points(b0s, posterior.b0, 'l', col="blue")
    
    points(x=0, y=prior.b0[zeroloc], cex=1.5, col="red")
    points(x=0, y=posterior.b0[zeroloc], cex=1.5, col="red")
    
    if (m1 > m0) {
      points(x=c(0,0), y=c(0, prior.b0[zeroloc]) , type='l', col="red", lwd=2)
      points(x=c(0,0), y=c(0, posterior.b0[zeroloc]) , type='l', col="blue", lwd=2)
    }
    
    if (m0 > m1) {
      points(x=c(0,0), y=c(0, posterior.b0[zeroloc]) , type='l', col="blue", lwd=2)
      points(x=c(0,0), y=c(0, prior.b0[zeroloc]) , type='l', col="red", lwd=2)
    }
    
    figure <- recordPlot()
  }
  par(mfrow=c(1,1))
  
  BF10 <- m1/m0
  BF01 <- m0/m1
  
  figure
  
  return(list(figure=figure, L.m0=m0, L.m1=m1, prior.mu0=prior.b0[zeroloc], 
              posterior.mu0=posterior.b0[zeroloc],BF10=BF10, BF01=BF01))
}
