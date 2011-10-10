setMethod(
          f = "plotQuality",
          signature = "FastqFileList", 
          definition = function(x, ..., col=rainbow(length(x)))
          {
            quals <- lapply(x, function(x) {
              fq <- quality(readFastq(x))
              wd <- unique(width(fq))
              if (1L != length(wd)) {
                minwd <- min(wd)
                message("reducing width to trailing ", minwd,
                        "\n  path: ", x)
                fq <- narrow(fq, end(fq) - minwd + 1L, end(fq))
              }
              colMeans(as(fq, "matrix"))
            })
            wd <- max(sapply(quals, length))
            min <- min(sapply(quals, min))
            max <- max(sapply(quals, max))
            plot(quals[[1]], xlim=c(1, wd), ylim=c(min, max), type="b",
                 ..., xlab="Cycle", ylab="Quality", col=col[1])
            for (i in 1 + seq_along(quals[-1]))
              lines(quals[[i]], type="l", col=col[i])
            invisible()
          })


setMethod(
          f = "barplot",
          signature = signature(height="FastqFileList"),
          definition = function(height,...) {
            x <- lapply(height,function(x) length(readFastq(x)))
            barplot(unlist(x),...)
          }
          )
