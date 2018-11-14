nnGpos <- function(G, pos) {

# G   :: gmeta6 object
# pos :: [x,y] coordinate matrix
# return: vector indexes within G [byrow=TRUE; from upper-left corner], for each coordinate in pos

 nnid <- function(x,v) {as.integer(which.min(abs(x-v)))}
 xids <- sapply(pos[,1],nnid,v=G$xseq)
 yids <- sapply(pos[,2],nnid,v=G$ryseq)
 ids  <- as.integer(G$cols*(yids-1)+xids)
 return(ids)
} # end function nnGpos()
