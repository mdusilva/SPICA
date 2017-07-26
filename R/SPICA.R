#check overlap function

c <- 2.99792458e18 #speed of light in Angstron/sec
abzero <- -48.60    # magnitude zero points
stzero <- -21.10

cardelliLaw <- function(wavelength, A, Rv=3.1) {
  extinction <- wavelength
  X = wavelength * 1.e-4
  X = 1. / X
  for(i in 1:length(X)) {
    xi <- X[i]
    if (xi < 0.3) {

      stop("Outside Cardelli Law range")
    } else if (xi >= 0.3 && xi < 1.1) {

      a <- 0.574 * xi^1.61
      b <- -0.527 * xi^1.61

    } else if(xi >= 1.1 && xi < 3.3) {

      y <- xi - 1.82
      a <- 1. + 0.17699 * y - 0.50447 * y^2. -0.02427 * y^3. + 0.72085 * y^4. + 0.01979 * y^5. - 0.77530 * y^6. + 0.32999 * y^7.
      b <- 1.41338 * y + 2.28305 * y^2. + 1.07233 * y^3. - 5.38434 * y^4. - 0.62251 * y^5. + 5.30260 * y^6. - 2.09002 * y^7.
    } else if(xi >= 3.3 && xi < 8) {

      if (xi >= 5.9) {

        Fa <- -0.04473 * (xi - 5.9)^2. - 0.009779 * (xi - 5.9)^3.
        Fb <- 0.2130 * (xi - 5.9)^2. + 0.1207 * (xi - 5.9)^3.
      } else {

        Fa <- 0.
        Fb <- 0.
      }

      a <- 1.752 - 0.316 * xi - 0.104 / ((xi - 4.67)^2. + 0.341) + Fa
      b <- -3.090 + 1.825 * xi + 1.206 / ((xi - 4.62)^2. + 0.263) + Fb
    } else if (xi >= 8. && xi < 10.) {

      a <- -1.073 - 0.628 * (xi - 8.) + 0.137 * (xi - 8.)^2. - 0.070 * (xi - 8.)^3.
      b <- 13.670 + 4.257 * (xi - 8.) - 0.420 * (xi - 8.)^2. + 0.374 * (xi - 8.)^3.
    } else {

      stop("Outside Cardelli Law range")
    }
    extinction[i] <- A * (a + b / Rv)
  }

return (extinction)
}

checkOverlap <- function(rangefilt, rangespec){

  wbmax <- max(rangefilt)
  wbmin <- min(rangefilt)
  wsmax <- max(rangespec)
  wsmin <- min(rangespec)
  if (wbmin < wsmin || wbmax > wsmax) {
    #stop("must overlap")
    warning("Filter must be inside spectrum")
    return(FALSE)
  } else {
    return(TRUE)
  }

}

computeTotalMag <- function(filter, spectrum, A, magSys="Vega", magZero=0, rel.tol=.Machine$double.eps^0.25) {

apmag <- computeSyntheticMag(filter, spectrum, magSys=magSys, magZero=magZero, rel.tol=rel.tol)
extinction <- computeExtinction(filter, spectrum, A, rel.tol=rel.tol)

return(apmag+extinction)
}

computeExtinction <- function(filter, spectrum, A, rel.tol=.Machine$double.eps^0.25) {

  #Verify if filter and spectrum overlap
  if (checkOverlap(filter[,1], spectrum[,1]) == FALSE) {
    stop("must overlap")
  }
  intertpSpectrum <- approxfun(spectrum[,1], spectrum[,2], method="linear")
  interpFilter <- approxfun(filter[,1], filter[,2], method="linear")
  integratedSpectrum <- integrate(function(x) {intertpSpectrum(x) * interpFilter(x) * x * 10.^(-0.4 * cardelliLaw(x,A))}, lower = min(filter[,1]), upper = max(filter[,1]), rel.tol=rel.tol)
  integratedFilter <- integrate(function(x) {intertpSpectrum(x) * interpFilter(x) * x}, lower = min(filter[,1]), upper = max(filter[,1]), rel.tol=rel.tol)
  result <- -2.5 * log10(integratedSpectrum$value) + 2.5 * log10(integratedFilter$value)
  return(result)
}


computeSyntheticMag <- function(filter, spectrum, magSys="Vega", magZero=0, rel.tol=.Machine$double.eps^0.25) {

  #Verify if filter and spectrum overlap
  if (checkOverlap(filter[,1], spectrum[,1]) == FALSE) {
    stop("must overlap")
  }

  #waveVec <- seq(from=wsmin, to=wsmax, length.out=1000)
  intertpSpectrum <- approxfun(spectrum[,1], spectrum[,2], method="linear")
  interpFilter <- approxfun(filter[,1], filter[,2], method="linear")

  #Compute magnitude
  if (magSys == "AB") {
   integratedSpectrum <- integrate(function(x) {intertpSpectrum(x) * interpFilter(x) * x}, lower = min(filter[,1]), upper = max(filter[,1]), rel.tol=rel.tol)
   integratedFilter <- integrate(function(x) {interpFilter(x) * c / x}, lower=min(filter[,1]), upper=max(filter[,1]))
   result <- -2.5 * log10(integratedSpectrum$value) + 2.5 * log10(integratedFilter$value) + abzero
     } else if (magSys == "ST") {
       integratedSpectrum <- integrate(function(x) {intertpSpectrum(x) * interpFilter(x) * x}, lower = min(filter[,1]), upper = max(filter[,1]), rel.tol=rel.tol)
       integratedFilter <- integrate(function(x) {interpFilter(x) * x}, lower=min(filter[,1]), upper=max(filter[,1]))
       result <- -2.5 * log10(integratedSpectrum$value) + 2.5 * log10(integratedFilter$value) + stzero
  } else if (magSys == "Vega") {
    if (checkOverlap(filter[,1], vegaSpectrum[,1]) == FALSE) {
      stop("must overlap")
    }
    intertpVegaSpectrum <- approxfun(vegaSpectrum[,1], vegaSpectrum[,2], method="linear")
    integratedSpectrum <- integrate(function(x) {intertpSpectrum(x) * interpFilter(x) * x}, lower = min(filter[,1]), upper = max(filter[,1]), rel.tol=rel.tol)
    #set correct Vea interation limits
    integratedVegaSpectrum <- integrate(function(x) {intertpVegaSpectrum(x) * interpFilter(x) * x}, lower = min(filter[,1]), upper = max(filter[,1]), rel.tol=rel.tol)
    result <- -2.5 * log10(integratedSpectrum$value) + 2.5 * log10(integratedVegaSpectrum$value) + magZero
  }

return(result)
}
