subsampleX <- function(X,ims) {
 # assumes X$val has a matrix format and select subsamples in the ensemble according to input ensemble indices
 Xs <- X
 for (iv in 1:length(Xs)) {
   for (it in 1:length(Xs[[iv]]$val)) {
     Xs[[iv]]$val[[it]] <-  X[[iv]]$val[[it]][,ims]
   }
 }
 return(Xs)
}
