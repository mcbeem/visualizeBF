#' Function for exploring Bayes Factors.
#'
#' 
#' #' @param data A vector of data
#' @param pointsize Width of points for numerical integration, defaults to .001.
#' @param scale Scale parameter for Cauchy prior, defaults to .707
#' @param plot The type of plot to produce. 1=BF as weighted likelihood, 2=BF as comparison of likelihood for mu=0 under prior and posterior distributions.
#' @export
#' @return Returns a list with the following components:  
#' \itemize{
#'   \item \code{BF10}    The Bayes Factor for alternative hypothesis relative to null.  \cr 
#'   \item \code{BF01}    The Bayes Factor for the null hypothesis relative to the alternative \cr 
#'  }
#' @examples
#' set.seed(1)
#' data <- rnorm(n=50, mean=.3, sd=1)
#' visualizeBF(data, plot=1)
#' visualizeBF(data, plot=2)

visualizeBF <- function(data, pointsize=0.001, scale=.707, plot=1) {
    
  y <- data
  
  # standardize data so sample sd is exactly 1
  y <- y / sd(y)
  
  # calculate likelihood function for the data
  # there are two parameters
  
  pointsize <- .001
  
  b0s <- seq(-3,3, pointsize)
  
  # define a function to replace zero values with a tiny one to avoid infinities when
  #  calculating the log
  replacezeros <- function(x) {
    if(x==0) {x <- 1e-60}
    return(x)
  }
  
  # calculate LL for each parm value
  LL <- sapply(b0s, dnorm, x=y, sd=1) %>% apply(c(1,2), log) %>% apply(2, sum)
  L <- sapply(b0s, dnorm, x=y, sd=1) %>% apply(2, prod)
  
  
  # prior for b0 based on Cauchy distribution
  prior.b0 <- dcauchy((b0s), scale=scale)
  
  # normalize it
  prior.b0 <- prior.b0 / (sum(prior.b0)*pointsize)
  log.prior.b0 <- prior.b0 %>% log()
  
  
  # alternative normal prior
  #prior.b0 <- dnorm(x=b0s, mean=0, sd=.2) %>% log()
  
  # which b0 is zero?
  zeroloc <- which(b0s==0)
  
  #### This is calculating BF as the weighted average of the likelihood
  # numerator of Bayes theorem
  m0 <- LL[zeroloc] %>% exp()
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
  
  
  
  
  ### Make a fancy 3 panel plot
  par(mfrow=c(2,1))
  par(mar=rep(2.2,4))
  
  if (plot==1) {
  
    # top panel: likelihood for full range of d
    plot(b0s, exp(LL), 'l', main="Visualizing BF: Ratio of L at mu=0 to weighted average L")
    points(x=b0s[zeroloc]-.015, y=(m0), cex=1.5, col="red")
    points(x=b0s[zeroloc]+.015, y=(m1), cex=1.5, col="red")
    abline(h=(m1), col="red")
    if (m1 > m0) {
      points(x=c(0,0), y=c((m0), (m1)) , type='l', col="blue", lwd=2)
      points(x=c(0,0), y=c(0, (m0)) , type='l', col="red", lwd=2)
    }
  
    if (m0 > m1) {
      points(x=c(-.015,-.015), y=c(0, (m0)) , type='l', col="blue", lwd=2)
      points(x=c(.015,.015), y=c(0, (m1)) , type='l', col="red", lwd=2)
    }
    
    # bottom panel: zoomed in version of top panel
    plot(b0s, exp(LL), 'l', xlim=c(-.1, .1), ylim=c(0, max(m1,m0)*1.1), main="Zoomed-in version of upper panel")
    abline(h=0, col="darkgray")
    points(x=0, y=(m0), cex=1.5, col="red")
    points(x=.002, y=(m1), cex=1.5, col="red")
    abline(h=(m1), col="red")
    
    if (m1 > m0) {
      points(x=c(0,0), y=c(0, (m1)) , type='l', col="blue", lwd=2)
      points(x=c(.002,.002), y=c(0, (m0)) , type='l', col="red", lwd=2)
    }
    
    if (m0 > m1) {
      points(x=c(0,0), y=c(0, (m0)) , type='l', col="blue", lwd=2)
      points(x=c(.002,.002), y=c(0, (m1)) , type='l', col="red", lwd=2)
    }
  }
  
  if (plot==2) {
  
    # top panel: prior for d
    ymax <- max(c(prior.b0, posterior.b0))
    
    plot(b0s, prior.b0, 'l', ylim=c(0,ymax), col="red", main="Visualizing BF: Ratio of prior to posterior L at mu=0")
    abline(h=0, col="darkgray")
    points(b0s, posterior.b0, 'l', col="blue")
    
    points(x=.015, y=prior.b0[zeroloc], cex=1.5, col="red")
    points(x=-.015, y=posterior.b0[zeroloc], cex=1.5, col="red")
    
    if (m1 > m0) {
      points(x=c(-.015,-.015), y=c(0, posterior.b0[zeroloc]) , type='l', col="red", lwd=2)
      points(x=c(.015,.015), y=c(0, prior.b0[zeroloc]) , type='l', col="blue", lwd=2)
    }
    
    if (m0 > m1) {
      points(x=c(-.015,-.015), y=c(0, posterior.b0[zeroloc]) , type='l', col="blue", lwd=2)
      points(x=c(.015,.015), y=c(0, prior.b0[zeroloc]) , type='l', col="red", lwd=2)
    }
  
    # bottom panel
    plot(b0s, prior.b0, 'l', xlim=c(-.1, .1), ylim=c(0,ymax), col="red", main="Zoomed-in version of upper panel")
    abline(h=0, col="darkgray")
    points(b0s, posterior.b0, 'l', col="blue")
    
    points(x=.001, y=prior.b0[zeroloc], cex=1.5, col="red")
    points(x=-.001, y=posterior.b0[zeroloc], cex=1.5, col="red")
    
    if (m1 > m0) {
      points(x=c(-.001,-.001), y=c(0, posterior.b0[zeroloc]) , type='l', col="red", lwd=2)
      points(x=c(.001,.001), y=c(0, prior.b0[zeroloc]) , type='l', col="blue", lwd=2)
    }
    
    if (m0 > m1) {
      points(x=c(-.001,-.001), y=c(0, posterior.b0[zeroloc]) , type='l', col="blue", lwd=2)
      points(x=c(.001,.001), y=c(0, prior.b0[zeroloc]) , type='l', col="red", lwd=2)
    }
  }
  par(mfrow=c(1,1))
  
  BF10 <- m1/m0
  BF01 <- m0/m1
  
  return(list(BF10=BF10, BF01=BF01))
}