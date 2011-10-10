.gcLoess <- function(counts,gc) {
  ff <- function(y,x) {
    xx <- x[(y>0)&(y<=quantile(y,probs=0.99))]
    yy <- log(y[(y>0)&(y<=quantile(y,probs=0.99))])
    
    l <- loess(yy~xx)
    y.fit <- predict(l,newdata=x)
    names(y.fit) <- names(y)
    y.fit[is.na(y.fit)] <- 0
  
    retval <- y/exp(y.fit-median(yy))
    return(retval)
  }
  apply(counts,2,ff,x=gc)
}

.gcQuant <- function(counts,gc,num.bins=10,which=c("full","median","upper")) {
  which <- match.arg(which)
  bins <- cut(gc,breaks=quantile(gc,probs=seq(0,1,length.out=num.bins+1)))
  bins[is.na(bins)] <- levels(bins)[1]
  names(bins) <- names(gc)
  f <- function(y) {
    if(is.null(names(y))) {
      names(y) <- 1:length(y)
    }
    tmp <- tapply(y,bins,function(x) x)
    switch(which,
           full = {y.norm <- normalizeQuantileRank(tmp)},
           median = {y.norm <- lapply(tmp,function(x) exp(log(x) - median(log(x)) + median(log(unlist(tmp)))))},
           upper = {y.norm <- lapply(tmp,function(x) exp(log(x) - quantile(log(x),probs=.75) + quantile(log(unlist(tmp)),prob=.75)))}
           )
    names <- unlist(sapply(y.norm,names))
    y.norm <- unlist(y.norm)
    names(y.norm) <- names
    y.norm[names(y)]
  }
  apply(counts,2,f)
}
