#check overlap function

c <- 2.99792458e18 #speed of light in Angstron/sec
abzero <- -48.60    # magnitude zero points
stzero <- -21.10

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


computeSyntheticMag <- function(filter, spectrum, magSys="Vega", magZero=0) {

  #Verify if filter and spectrum overlap
  if (checkOverlap(filter[,1], spectrum[,1]) == FALSE) {
    stop("must overlap")
  }

  #waveVec <- seq(from=wsmin, to=wsmax, length.out=1000)
  intertpSpectrum <- approxfun(spectrum[,1], spectrum[,2], method="linear")
  interpFilter <- approxfun(filter[,1], filter[,2], method="linear")

  #Compute magnitude
  if (magSys == "AB") {
   integratedSpectrum <- integrate(function(x) {intertpSpectrum(x) * interpFilter(x) * x}, lower = wsmin, upper = wsmax) / integrate(function(x) {interpFilter(x) * c / x}, lower=wsmin, upper=wsmax)
    result <- -2.5 * log10(integratedSpectrum) + abzero
     } else if (magSys == "ST") {
       integratedSpectrum <- integrate(function(x) {intertpSpectrum(x) * interpFilter(x) * x}, lower = wsmin, upper = wsmax) / integrate(function(x) {interpFilter(x) * x}, lower=wsmin, upper=wsmax)
       resukt <- -2.5 * log10(integratedSpectrum) + stzero
  } else if (magSys == "Vega") {
    if (checkOverlap(filter[,1], vegaSpectrum[,1]) == FALSE) {
      stop("must overlap")
    }
    intertpVegaSpectrum <- approxfun(vegaSpectrum[,1], vegaSpectrum[,2], method="linear")
    integratedSpectrum <- integrate(function(x) {intertpSpectrum(x) * interpFilter(x) * x}, lower = wsmin, upper = wsmax)
    #set correct Vea interation limits
    integratedVegaSpectrum <- integrate(function(x) {intertpVegaSpectrum(x) * interpFilter(x) * x}, lower = , upper = )
    result <- -2.5 * log10(integratedSpectrum) + 2.5 * log10(integratedVegaSpectrum) + magZero
  }

return(result)
}
