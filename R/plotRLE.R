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
            if(ncol(exprs(x))<=1) {
              stop("At least two samples needed for the mean-difference plot.")
            } else {
              plotRLE(exprs(x), ...)
            }
          }
          )
