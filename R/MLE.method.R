#' Abundance Size Specturm with MLE estimate for b
#'
#' This function is modified from \url{https://github.com/andrew-edwards/fitting-size-spectra/blob/master/code/recommend/recommend.r} for re-producing Fig 6 in Edwards et al. (2017) with ggplot2
#'
#' @param x A vector of individual body masses
#' @author Chih-Lin Wei <chihlinwei@@gmail.com>, Jian-Xiang Liao <jianxiangliao@@gmail.com >
#' @return NBSS slope, confidence interval of slope; fitted values and upper CI and lower CI; individual rank and upper CI and lower CI
#' @references Edwards, A.M., Robinson, J.P.W., Plank, M.J., Baum, J.K., Blanchard, J.L., 2017. Testing and recommending methods for fitting size spectra to data. Methods Ecol Evol 8, 57–67. \url{https://doi.org/10.1111/2041-210X.12641}
#' @export
#' @examples
#' library(sizeSpectra)
#' MLE.method(mei$Size)

MLE.method<-
  function(x){
    log.x = log(x)                      # to avoid keep calculating
    sum.log.x = sum( log.x )
    xmin = min(x)
    xmax = max(x)

    # MLE (maximum likelihood method) calculations.

    # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
    #  as a starting point for nlm for MLE of b for PLB model.
    PL.bMLE = 1/( log(min(x)) - sum.log.x/length(x)) - 1

    PLB.minLL =  nlm(negLL.PLB, p=PL.bMLE, x=x, n=length(x),
                     xmin=xmin, xmax=xmax, sumlogx=sum.log.x) #, print.level=2 )

    PLB.bMLE = PLB.minLL$estimate

    # 95% confidence intervals for MLE method.

    PLB.minNegLL = PLB.minLL$minimum

    # Values of b to test to obtain confidence interval. For the real movement data
    #  sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
    #  symmetric interval here.

    bvec = seq(PLB.bMLE - 0.5, PLB.bMLE + 0.5, 0.00001)

    PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
    for(i in 1:length(bvec))
    {
      PLB.LLvals[i] = negLL.PLB(bvec[i], x=x, n=length(x), xmin=xmin,
                                xmax=xmax, sumlogx=sum.log.x)
    }
    critVal = PLB.minNegLL  + qchisq(0.95,1)/2
    # 1 degree of freedom, Hilborn and Mangel (1997) p162.
    bIn95 = bvec[ PLB.LLvals < critVal ]
    # b values in 95% confidence interval
    PLB.MLE.bConf = c(min(bIn95), max(bIn95))
    if(PLB.MLE.bConf[1] == min(bvec) | PLB.MLE.bConf[2] == max(bvec))
    { windows()
      plot(bvec, PLB.LLvals)
      abline(h = critVal, col="red")
      stop("Need to make bvec larger - see R window")   # Could automate
    }

    x.PLB = seq(min(x), max(x), length=1000)     # x values to plot PLB. Note
    # that these encompass the data, and are not based
    # on the binning (in MEE Figure 6 the line starts as
    # min(x), not the first bin.

    B.PLB = dPLB(x.PLB, b = PLB.bMLE, xmin=min(x.PLB),
                 xmax=max(x.PLB)) * length(x) * x.PLB
    B.PLB.UC <- dPLB(x.PLB, b = min(bIn95), xmin=min(x.PLB), xmax=max(x.PLB)) * length(x) * x.PLB
    B.PLB.LC <- dPLB(x.PLB, b = max(bIn95), xmin=min(x.PLB), xmax=max(x.PLB)) * length(x) * x.PLB
    # The biomass density, from equation (7), using the MLE for b.

    y.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE, xmin = min(x.PLB),
                      xmax = max(x.PLB))) * length(x)
    y.PLB.UC <- (1 - pPLB(x = x.PLB, b = min(bIn95), xmin = min(x.PLB), xmax = max(x.PLB))) * length(x)
    y.PLB.LC <- (1 - pPLB(x = x.PLB, b = max(bIn95), xmin = min(x.PLB),xmax = max(x.PLB))) * length(x)

    # NBSS slope, confidence interval of slope, fitted values and upper CI, lower CI, individual rank and upper CI, lower CI
    list(PLB.bMLE=PLB.bMLE, PLB.MLE.bConf=PLB.MLE.bConf, PLB=cbind(x.PLB, B.PLB, B.PLB.UC, B.PLB.LC, y.PLB, y.PLB.UC, y.PLB.LC))
  }
