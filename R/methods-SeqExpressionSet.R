#initialize
setMethod(
          f = "initialize",
          signature = "SeqExpressionSet",
          definition = function(.Object, ..., assayData=assayDataNew(
                                   exprs=matrix(0L, 0, 0),
                                   offset=matrix(0L, 0, 0)))
          {
            callNextMethod(.Object, ..., assayData=assayData)
          })

setValidity(
            Class= "SeqExpressionSet",
            method = function(object) {
              is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
                abs(x - round(x)) < tol
              }
              msg <- validMsg(NULL, assayDataValidMembers(assayData(object),c("exprs","offset")))
              if(!is.null(assayDataElement(object, "exprs"))) {
                if(class(assayDataElement(object, "exprs"))!="matrix") {
                  msg <- validMsg(msg,"exprs must be an integer or numeric matrix")
                }
                if(!all(is.wholenumber(assayDataElement(object, "exprs")),na.rm=T)) {
                  warning("exprs contains non-integer numbers")
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

#constructor
newSeqExpressionSet <- function(exprs,
                                offset=matrix(data=0,nrow=nrow(exprs),ncol=ncol(exprs),
                                  dimnames=dimnames(exprs)),
                                phenoData=annotatedDataFrameFrom(exprs, FALSE),
                                featureData=annotatedDataFrameFrom(exprs, TRUE),
                                ...) {
  if(class(phenoData)=="data.frame") {
    phenoData <- AnnotatedDataFrame(phenoData)
  }
  if(class(featureData)=="data.frame") {
    featureData <- AnnotatedDataFrame(featureData)
  }
  
  new("SeqExpressionSet",
      assayData=assayDataNew(exprs=exprs,offset=offset),
      phenoData=phenoData,featureData=featureData, ...)
}


#exprs and offset
setMethod(
          f = "exprs",
          signature = "SeqExpressionSet",
          definition = function(object) {
            assayDataElement(object, "exprs")
          }
        )

setReplaceMethod(
                 f = "exprs",
                 signature = "SeqExpressionSet",
                 definition = function(object,value) {
                   assayDataElement(object, "exprs") <- as.matrix(value)
                   validObject(object)
                   object
                 }
                 )

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
                 definition = function(object,value) {
                   assayDataElement(object, "offset") <- as.matrix(value)
                   validObject(object)
                   object
                 }
                 )


#some graphics
setMethod(
          f = "boxplot",
          signature = "SeqExpressionSet",
          definition = function(x, ...) {
            boxplot(as.data.frame(log(exprs(x) + 0.1)),...)
          }
          )


setMethod(
          f = "meanVarPlot",
          signature = "SeqExpressionSet",
          definition = function(x,log=FALSE,...) {
            if(ncol(exprs(x))<=1) {
              stop("At least two lanes are needed to compute variance.")
            }
            m <- apply(exprs(x),1,mean)
            v <- apply(exprs(x),1,var)
            if(log) {
              mm <- pmax(0,log(m))
              vv <- pmax(0,log(v))
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
          signature = signature(x="matrix",y="numeric"),
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
          signature = signature(x="SeqExpressionSet",y="character"),
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
            biasPlot(exprs(x), fData(x)[,y], log=log, col=col, ..., xlab=xlab, ylab=ylab)
            if(ncol(exprs(x))>1 & legend & flag) {
            legend("topleft",unique(as.character(pData(x)[,color_code])),fill=unique(pData(x)[,color_code]))
          }
          }
          )

setMethod(
          f = "biasBoxplot",
          signature = signature(x="numeric",y="numeric",num.bins="ANY"),
          definition = function(x,y,num.bins,...) {
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
          definition = function(x,y,...) {
            if(ncol(x)<=1) {
              stop("At least a two-column matrix needed for the mean-difference plot.")
            } else {
              mean <- rowMeans(log(x + 0.1))
              difference <- log(x[,y[2]] + 0.1)-log(x[,y[1]] + 0.1)
              smoothScatter(mean,difference,...)
              lines(lowess(mean,difference),col=2)
              abline(h=0,lty=2)
            }
          }
          )

setMethod(
          f = "MDPlot",
          signature = signature(x="SeqExpressionSet",y="numeric"),
          definition = function(x,y,...) {
            if(ncol(exprs(x))<=1) {
              stop("At least two lanes needed for the mean-difference plot.")
            } else {
              m <- exprs(x)[,y]
              mean <- rowMeans(log(m + 0.1))
              difference <- log(m[,2] + 0.1) - log(m[,1] + 0.1)
              smoothScatter(mean,difference,...)
              lines(lowess(mean,difference),col=2)
              abline(h=0,lty=2)
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
          signature = signature(x="SeqExpressionSet",y="character"),
          definition = function(x,y,which=c("loess","median","upper","full"),offset=FALSE,num.bins=10, round=TRUE) {
            if(offset) {
              newSeqExpressionSet(exprs=exprs(x), phenoData=phenoData(x), featureData=featureData(x), offset=withinLaneNormalization(exprs(x), fData(x)[,y], which, offset, num.bins, round))
            } else {
              newSeqExpressionSet(exprs=withinLaneNormalization(exprs(x), fData(x)[,y], which, offset, num.bins, round), phenoData=phenoData(x), featureData=featureData(x))
            }
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
          definition = function(x,which=c("median","upper","full"), offset=FALSE, round=TRUE) {
            if(offset) {
              if(all(offst(x)==0)) {
                newSeqExpressionSet(exprs=exprs(x), phenoData=phenoData(x), featureData=featureData(x), offset=betweenLaneNormalization(exprs(x), which, offset, round))
              } else {
                counts <- exp(log(exprs(x) + 0.1) + offst(x)) - 0.1
                newSeqExpressionSet(exprs=exprs(x), phenoData=phenoData(x), featureData=featureData(x), offset=offst(x) + betweenLaneNormalization(counts,which,offset,round))
              }
            }
            else {
              newSeqExpressionSet(exprs=betweenLaneNormalization(exprs(x), which, offset, round), phenoData=phenoData(x), featureData=featureData(x), offset=offst(x))
            }
          }
          )

setAs("SeqExpressionSet",
      "CountDataSet",
      function(from) {
        if(!("conditions" %in% colnames(pData(data)))) {
          stop("phenoData must contain a column named 'conditions'")
        }
        if(NCOL(pData(data))==1 & length(levels(pData(data)$conditions))==2) {
          newCountDataSet(round(exprs(from)),pData(from)[,1])
        } else {
          newCountDataSet(round(exprs(from)), pData(from), sizeFactors = NULL, featureData = featureData(from))
        }
      }
      )

