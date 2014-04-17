## initialize
setMethod(
          f = "initialize",
          signature = "SeqExpressionSet",
          definition = function(.Object, ..., assayData=assayDataNew(
                                   counts = matrix(0L, 0, 0),
                                   normalizedCounts =  matrix(0L, 0, 0),
                                   offset = matrix(0L, 0, 0)))
          {
            callNextMethod(.Object, ..., assayData=assayData)
          })

setValidity(
            Class= "SeqExpressionSet",
            method = function(object) {
              is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
                abs(x - round(x)) < tol
              }
              msg <- validMsg(NULL, assayDataValidMembers(assayData(object), c("counts", "normalizedCounts", "offset")))
              if(!is.null(assayDataElement(object, "counts"))) {
                if(class(assayDataElement(object, "counts"))!="matrix") {
                  msg <- validMsg(msg, "'counts' must be an integer or numeric matrix")
                }
                if(!all(is.wholenumber(assayDataElement(object, "counts")), na.rm=TRUE)) {
                  warning("'counts' contains non-integer numbers")
                } 
              }
              if (is.null(msg))
                TRUE
              else msg
            }
            )

setMethod("updateObject",
          signature = signature(object = "SeqExpressionSet"),
          definition = function(object, ..., verbose = FALSE) {
            if (verbose) {
              message("updateObject(object = 'SeqExpressionSet')")
            }
            object <- callNextMethod()
            if (isCurrent(object)["SeqExpressionSet"]) {
              return(object)
            } else {
              classVersion(object)["SeqExpressionSet"] <- classVersion("SeqExpressionSet")["SeqExpressionSet"]
              object
            }
          })

## constructor
newSeqExpressionSet <- function(counts,
                                normalizedCounts = matrix(data=NA, nrow=nrow(counts), ncol=ncol(counts),
                                  dimnames=dimnames(counts)),
                                offset = matrix(data=0, nrow=nrow(counts), ncol=ncol(counts),
                                  dimnames=dimnames(counts)),
                                phenoData = annotatedDataFrameFrom(counts, FALSE),
                                featureData = annotatedDataFrameFrom(counts, TRUE),
                                ...) {
  if(class(phenoData) == "data.frame") {
    phenoData <- AnnotatedDataFrame(phenoData)
  }
  if(class(featureData) == "data.frame") {
    featureData <- AnnotatedDataFrame(featureData)
  }
  
  new("SeqExpressionSet",
      assayData = assayDataNew(counts=counts, normalizedCounts=normalizedCounts, offset=offset),
      phenoData = phenoData,
      featureData = featureData, ...)
}


## exprs DEPRECATED
setMethod(
          f = "exprs",
          signature = "SeqExpressionSet",
          definition = function(object) {
            .Deprecated("counts")
            
            if(all(is.na(normCounts(object)))) {
              counts <- counts(object)
            } else {
              counts <- normCounts(object)
            }
            return(counts)
          }
        )

setReplaceMethod(
                 f = "exprs",
                 signature = "SeqExpressionSet",
                 definition = function(object, value) {
                   .Deprecated("counts<-")
                   assayDataElement(object, "counts") <- as.matrix(value)
                   validObject(object)
                   object
                 }
                 )

## counts
setMethod(
          f = "counts",
          signature = "SeqExpressionSet",
          definition = function(object) {
            assayDataElement(object, "counts")
          }
        )

setReplaceMethod(
                 f = "counts",
                 signature = "SeqExpressionSet",
                 definition = function(object,value) {
                   assayDataElement(object, "counts") <- as.matrix(value)
                   validObject(object)
                   object
                 }
                 )

## normCounts
setMethod(
          f = "normCounts",
          signature = "SeqExpressionSet",
          definition = function(object) {
            assayDataElement(object, "normalizedCounts")
          }
        )

setReplaceMethod(
                 f = "normCounts",
                 signature = "SeqExpressionSet",
                 definition = function(object,value) {
                   assayDataElement(object, "normCounts") <- as.matrix(value)
                   validObject(object)
                   object
                 }
                 )

## offst
setMethod(
          f = "offst",
          signature = "SeqExpressionSet",
          definition = function(object) {
            assayDataElement(object, "offset")
          }
        )


setReplaceMethod(
                 f = "offst",
                 signature = "SeqExpressionSet",
                 definition = function(object, value) {
                   assayDataElement(object, "offset") <- as.matrix(value)
                   validObject(object)
                   object
                 }
                 )


## some graphics
setMethod(
          f = "boxplot",
          signature = "SeqExpressionSet",
          definition = function(x, ...) {
            if(all(is.na(normCounts(x)))) {
              boxplot(as.data.frame(log(counts(x) + 0.1)),...)
            } else {
              boxplot(as.data.frame(log(normCounts(x) + 0.1)),...)
            }
          }
          )


setMethod(
          f = "meanVarPlot",
          signature = "SeqExpressionSet",
          definition = function(x,log=FALSE,...) {
            if(ncol(counts(x))<=1) {
              stop("At least two samples are needed to compute variance.")
            }
            if(all(is.na(normCounts(x)))) {
              counts <- counts(x)
            } else {
              counts <- normCounts(x)
            }
            m <- apply(counts, 1, mean)
            v <- apply(counts, 1, var)
            if(log) {
              mm <- pmax(0, log(m))
              vv <- pmax(0, log(v))
            } else {
              mm <- m[m<=quantile(m,probs=.9)]
              vv <- v[m<=quantile(m,probs=.9)]
            }
            smoothScatter(mm,vv,xlab="mean",ylab="variance",...)
            lines(abline(0,1))
            lines(lowess(mm,vv),col=2)
          }
          )


setMethod(
          f = "biasPlot",
          signature = signature(x="matrix", y="numeric"),
          definition = function(x, y, cutoff=1000, log=FALSE, col=NULL, ...) {
            if(log) {
              x <- log(x + 0.1)
            }
            if(is.null(col)) {
              col <- 1:ncol(x)
            } else if(length(col) < ncol(x)) {
              col = rep(col,length.out=ncol(x))
            }
            plot(lowess(y[x[,1]<=cutoff],x[x[,1]<=cutoff,1]),type='l',col=col[1],...)
            if(ncol(x)>1) {
              for(i in 2:ncol(x)) {
                lines(lowess(y[x[,i]<=cutoff],x[x[,i]<=cutoff,i]),col=col[i],type='l',...)
              }
            }
          }
          )

setMethod(
          f = "biasPlot",
          signature = signature(x="SeqExpressionSet", y="character"),
          definition = function(x, y, cutoff=1000, log=FALSE, color_code=NULL, legend=TRUE, col=NULL, ..., xlab = y, ylab = "gene counts") {
            flag <- FALSE
            if(is.null(col)) {
              if(is.null(color_code)) {
                color_code <- 1
              }
              col <- as.numeric(as.factor(pData(x)[,color_code]))
              flag <- TRUE
            } else if(!is.null(color_code)) {
              warning("If both col and color_code are specified, col overrides color_code")
            }
            if(all(is.na(normCounts(x)))) {
              counts <- counts(x)
            } else {
              counts <- normCounts(x)
            }
            biasPlot(counts, fData(x)[,y], log=log, col=col, ..., xlab=xlab, ylab=ylab)
            if(ncol(counts)>1 & legend & flag) {
            legend("topleft",unique(as.character(pData(x)[,color_code])),fill=unique(pData(x)[,color_code]))
          }
          }
          )

setMethod(
          f = "biasBoxplot",
          signature = signature(x="numeric",y="numeric",num.bins="ANY"),
          definition = function(x, y, num.bins,...) {
            if(missing(num.bins)) {
              num.bins <- 10
            }
            bins <- cut(y,breaks=quantile(y,probs=seq(0,1,length.out=num.bins+1)))
            bins[is.na(bins)] <- levels(bins)[1]
            names(bins) <- names(y)
            boxplot(x~bins,...)
            abline(h=0,col=2)
          }
          )

setMethod(
          f = "MDPlot",
          signature = signature(x="matrix",y="numeric"),
          definition = function(x, y, controls=NULL, ...) {
            if(ncol(x)<=1) {
              stop("At least a two-column matrix needed for the mean-difference plot.")
            } else {
              mean <- rowMeans(log(x + 0.1))
              difference <- log(x[,y[2]] + 0.1)-log(x[,y[1]] + 0.1)
              smoothScatter(mean,difference,...)
              lines(lowess(mean,difference),col=1)
              abline(h=0,lty=2)
              if(!is.null(controls)) {
                points(mean[controls], difference[controls], pch=20, col=2)
                lines(lowess(mean[controls], difference[controls]), col=2)
              }
            }
          }
          )

setMethod(
          f = "MDPlot",
          signature = signature(x="SeqExpressionSet",y="numeric"),
          definition = function(x, y, controls=NULL, ...) {
            if(ncol(counts(x))<=1) {
              stop("At least two samples needed for the mean-difference plot.")
            } else {
              if(all(is.na(normCounts(x)))) {
                m <- counts(x)[,y]
              } else {
                m <- normCounts(x)[,y]
              }
              mean <- rowMeans(log(m + 0.1))
              difference <- log(m[,2] + 0.1) - log(m[,1] + 0.1)
              smoothScatter(mean,difference,...)
              lines(lowess(mean,difference),col=1)
              abline(h=0,lty=2)
              if(!is.null(controls)) {
                points(mean[controls], difference[controls], pch=20, col=2)
                lines(lowess(mean[controls], difference[controls]), col=2)
              }
            }
          }
          )


setMethod(
          f = "withinLaneNormalization",
          signature = signature(x="matrix",y="numeric"),
          definition = function(x, y, which=c("loess", "median", "upper", "full"), offset=FALSE, num.bins=10, round=TRUE) {
            which <- match.arg(which)
            if(which=="loess") {
              retval <- .gcLoess(x,y)
            } else {
              retval <- .gcQuant(x,y,num.bins,which)
            }
            if(!offset) {
              if(round) {
                retval <- round(retval)
              }
              return(retval)
            } else {
              ret <- log(retval + 0.1)-log(x + 0.1)
              return(ret)
            }
          }
          )

setMethod(
          f = "withinLaneNormalization",
          signature = signature(x="SeqExpressionSet", y="character"),
          definition = function(x, y, which=c("loess", "median", "upper", "full"),
                                offset=FALSE, num.bins=10, round=TRUE) {
            if(offset) {
              o <- withinLaneNormalization(counts(x),
                                           fData(x)[,y],
                                           which, offset, num.bins, round)
            } else {
              o <- offst(x)
            }
            
            newSeqExpressionSet(counts=counts(x),
                                normalizedCounts=withinLaneNormalization(counts(x),
                                  fData(x)[,y], which, offset=FALSE, num.bins, round),
                                offset=o,
                                  phenoData=phenoData(x), featureData=featureData(x))
          }
          )

#between-lane

setMethod(
          f = "betweenLaneNormalization",
          signature = signature(x="matrix"),
          definition = function(x, which=c("median", "upper", "full"), offset=FALSE, round=TRUE) {
            which <- match.arg(which)
            if(which=="full") {
              retval <- normalizeQuantileRank(as.matrix(x), robust=TRUE)
            } else {
              if(which=="upper") {
                sum <- apply(x, 2, quantile, 0.75)
              } else {
                sum <- apply(x, 2, median)
              }
              sum <- sum/mean(sum)
              retval <- scale(x, center=FALSE, scale=sum)
            }
            if(!offset) {
              if(round) {
                retval <- round(retval)
              }
              return(retval)
            } else {
              ret <- log(retval + 0.1) - log(x + 0.1)
              return(ret)
            }
          }
          )

setMethod(
          f = "betweenLaneNormalization",
          signature = signature(x="SeqExpressionSet"),
          definition = function(x, which=c("median","upper","full"), offset=FALSE, round=TRUE) {
            if(all(is.na(normCounts(x)))) {
              counts <- counts(x)
            } else {
              counts <- normCounts(x)
            }
            if(offset) {
              o = offst(x) + betweenLaneNormalization(counts, which, offset=TRUE, round)
            }
            else {
              o = offst(x)
            }
            newSeqExpressionSet(counts=counts(x),
                                  normalizedCounts=betweenLaneNormalization(counts,
                                    which, offset=FALSE, round),
                                  offset=o,
                                  phenoData=phenoData(x), featureData=featureData(x))
          }
          )

setAs("SeqExpressionSet",
      "CountDataSet",
      function(from) {
        if(!("conditions" %in% colnames(pData(from)))) {
          stop("phenoData must contain a column named 'conditions'")
        }
        if(all(is.na(normCounts(from)))) {
          counts <- counts(from)
        } else {
          counts <- normCounts(from)
        }
        if(NCOL(pData(from))==1 & length(levels(pData(from)$conditions))==2) {
          newCountDataSet(round(counts),pData(from)[,1])
        } else {
          newCountDataSet(round(counts), pData(from), sizeFactors = NULL, featureData = featureData(from))
        }
      }
      )

setMethod(
          f = "plotRLE",
          signature = signature(x="matrix"),
          definition = function(x,...) {
              y <- log(x+1)
              median <- apply(y, 1, median)
              rle <- apply(y, 2, function(x) x - median)

              boxplot(rle, ...)
              abline(h=0, lty=2)
          }
          )

setMethod(
          f = "plotRLE",
          signature = signature(x="SeqExpressionSet"),
          definition = function(x, ...) {
            if(ncol(counts(x))<=1) {
              stop("At least two samples needed for the mean-difference plot.")
            } else {
              if(all(is.na(normCounts(x)))) {
                counts <- counts(x)
              } else {
                counts <- normCounts(x)
              }
              plotRLE(counts, ...)
            }
          }
          )


setMethod(
          f = "plotPCA",
          signature = signature(x="matrix"),
          definition = function(x, k=2, ...) {
              Y <- apply(log(x+1), 1, function(y) scale(y, center=TRUE, scale=FALSE))
              s <- svd(Y)

              if(k>ncol(x)) {
                stop("The number of PCs must be less than the number of samples.")
              }
              
              if(k<2) {
                stop("The number of PCs must be at least 2.")
              } else if (k==2) {
                plot(s$u[,1], s$u[,2], type='n', ..., xlab="PC1", ylab="PC2")
                text(s$u[,1], s$u[,2], labels=colnames(x), ...)
              } else {
                colnames(s$u) <- paste("PC", 1:ncol(s$u), sep="")
                pairs(s$u[,1:k], ...)
              }
            }
          )

setMethod(
          f = "plotPCA",
          signature = signature(x="SeqExpressionSet"),
          definition = function(x, k=2, ...) {
            if(ncol(counts(x))<=1) {
              stop("At least two samples needed for the PCA plot.")
            } else {
              if(all(is.na(normCounts(x)))) {
                counts <- counts(x)
              } else {
                counts <- normCounts(x)
              }
              plotPCA(counts, ...)
            }
          }
          )

setMethod(
          f = "plotRLE",
          signature = signature(x="matrix"),
          definition = function(x,...) {
              y <- log(x+1)
              median <- apply(y, 1, median)
              rle <- apply(y, 2, function(x) x - median)

              boxplot(rle, ...)
              abline(h=0, lty=2)
          }
          )

setMethod(
          f = "plotRLE",
          signature = signature(x="SeqExpressionSet"),
          definition = function(x, ...) {
            if(ncol(counts(x))<=1) {
              stop("At least two samples needed for the mean-difference plot.")
            } else {
              if(all(is.na(normCounts(x)))) {
                counts <- counts(x)
              } else {
                counts <- normCounts(x)
              }             
              plotRLE(counts, ...)
            }
          }
          )
