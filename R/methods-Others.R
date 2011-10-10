setMethod(
          f="plotNtFrequency",
          signature= signature(x="ShortRead"),
          definition = function(x,...) {
            nt <- t(alphabetByCycle(sread(x)))
            nt <- nt[,c("A","C","G","T","N")]
            nt <- t(scale(t(nt), center=FALSE, scale=rowSums(nt)))
            matplot(nt, xlab = "Cycle", ylab = "Nt frequency", type = "l", lty=1, col=1:5, ...)
            abline(h=1/4, col="gray")
            legend("topright",colnames(nt), fill=1:5,bg="white")
          }
          )

setMethod(
          f = "boxplot",
          signature = signature(x="FastqQuality"),
          definition = function(x,...) {
            bp <- as.data.frame(as(x,"matrix"))
            colnames(bp) <- as.character(1:NCOL(bp))
            boxplot(bp,...)
          }
          )




