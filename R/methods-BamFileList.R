#plot methods
setMethod(
          f = "barplot",
          signature = signature(height="BamFileList"),
          definition = function(height,unique=FALSE,...) {
              bars <- countBam(height)$records
              names(bars) <- names(height)
              barplot(bars,...)
            }
          )

setMethod(
          f = "plotQuality",
          signature = "BamFileList", 
          definition = function(x, ..., col=rainbow(length(x)))
          {
            param <- ScanBamParam(what=c("strand","qual"))
            quals <- lapply(x, function(x) {
              bam <- scanBam(x, param=param)
              strand <- bam[[1]]$strand
              fq <- bam[[1]]$qual
              wd <- unique(width(fq))
              if (1L != length(wd)) {
                minwd <- min(wd)
                message("reducing width to trailing ", minwd,
                        "\n  path: ", x)
                fq <- narrow(fq, end(fq) - minwd + 1L, end(fq))
              }
              quality <- as(fq[strand=="+" | is.na(strand)], "matrix")
              quality <- rbind(quality,as(fq[strand=="-" & !is.na(strand)], "matrix")[,NCOL(quality):1])
              colMeans(quality)
            })
            wd <- max(sapply(quals, length))
            min <- min(sapply(quals, min))
            max <- max(sapply(quals, max))
            plot(quals[[1]], xlim=c(1, wd), ylim=c(min, max), type="l",
                 ..., xlab="Cycle", ylab="Quality", col=col[1])
            for (i in 1 + seq_along(quals[-1]))
              lines(quals[[i]], type="l", col=col[i])
            invisible()
          })

setMethod(
          f = "plot",
          signature = signature(x="BamFileList",y="FastqFileList"),
          definition = function(x,y,unique=FALSE,...) {
            if(!is.null(elementMetadata(x)) & !is.null(elementMetadata(y))) {
              if(!identical(elementMetadata(x),elementMetadata(y)))
                stop("x and y must have the same elementMetadata")
            }
            countBam(x)$records
            param <- ScanBamParam(what="seq")
            fq <- sapply(y,function(y) length(readFastq(y)))
            all <-countBam(x)$records
            if(unique) {
              unique <-  sapply(x,function(x) length(unique(unlist(scanBam(x,param=param))$seq)))
              nonunique <- all - unique
              bars <- t(cbind(unique,nonunique))
              barplot(t(t(bars)/fq),...)
            } else {
              barplot(all/fq,...)
            }
          }
          )

setMethod(
          f = "plotQuality",
          signature = "BamFile", 
          definition = function(x, ..., col=rainbow(length(x)))
          {
            param <- ScanBamParam(what=c("strand","qual"))
            bam <- scanBam(x, param=param)
            fq <- bam[[1]]$qual
            strand <- bam[[1]]$strand
            wd <- unique(width(fq))
            if (1L != length(wd)) {
              minwd <- min(wd)
              message("reducing width to trailing ", minwd,
                      "\n  path: ", x)
              fq <- narrow(fq, end(fq) - minwd + 1L, end(fq))
            }
            quality <- as(fq[strand=="+" | is.na(strand)], "matrix")
            quality <- rbind(quality,as(fq[strand=="-" & !is.na(strand)], "matrix")[,NCOL(quality):1])
            boxplot(quality,...)
          }
          )

setMethod(
          f = "barplot",
          signature = signature(height="BamFile"),
          definition = function(height,strata=c("rname","strand"),...) {
            st <- match.arg(strata)
            param <- ScanBamParam(what=st)
            bam <- scanBam(height, param=param)
            res <- bam[[1]][[1]]
            barplot(table(res),...)
          }
          )

setMethod(
          f="plotNtFrequency",
          signature= signature(x="BamFile"),
          definition = function(x,...) {
            param <- ScanBamParam(what=c("strand","seq"))
            bam <- scanBam(x, param=param)
            res <- bam[[1]]$seq
            strand <- bam[[1]]$strand
            nt <- t(alphabetByCycle(res[strand=="+" | is.na(strand)]))
            nt <- nt + t(alphabetByCycle(res[strand=="-" & !is.na(strand)]))[NROW(nt):1,]
            nt <- nt[,c("A","C","G","T","N")]
            nt <- t(scale(t(nt), center=FALSE, scale=rowSums(nt)))
            matplot(nt, xlab = "Cycle", ylab = "Nt frequency", type = "l", lty=1, col=1:5, ...)
            abline(h=1/4, col="gray")
            legend("bottomright",colnames(nt), fill=1:5,bg="white")
          }
          )
