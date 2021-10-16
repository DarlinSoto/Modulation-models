##########################
# AUTHOR: DARLIN SOTO    #
# DATE: OCTOBER 16, 2021 #
# EMAIL: dmsoto1@uc.cl   #
##########################


bspline <- function(x,xl,xr,ndx,bdeg){
  dx <- (xr-xl)/ndx
  knots <- seq(xl-bdeg*dx,xr+bdeg*dx,by=dx)
  B <- spline.des(knots,x,bdeg+1,0*x,outer.ok=TRUE)$design
  output <- list(knots=knots,B=B)
  return(output)
}