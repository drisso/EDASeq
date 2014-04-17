setGeneric(
           name = "plotQuality",
           def = function(x,...) {standardGeneric("plotQuality")}
           )

setGeneric(
           name = "biasPlot",
           def = function(x,y,...) {standardGeneric("biasPlot")}
           )

setGeneric(
           name = "meanVarPlot",
           def = function(x,...) {standardGeneric("meanVarPlot")}
           )

setGeneric(
           name = "MDPlot",
           def = function(x,y,...) {standardGeneric("MDPlot")}
           )

setGeneric(
           name = "biasBoxplot",
           def = function(x,y,num.bins,...) {standardGeneric("biasBoxplot")}
           )

setGeneric(
           name = "offst",
           def = function(object) {standardGeneric("offst")}
           )

setGeneric(
           name = "offst<-",
           def = function(object,value) {standardGeneric("offst<-")}
           )

setGeneric(
           name = "betweenLaneNormalization",
           def = function(x, which=c("median","upper","full"), offset=FALSE, round=TRUE) {standardGeneric("betweenLaneNormalization")}
           )

setGeneric(
           name = "withinLaneNormalization",
           def = function(x, y, which=c("loess","median","upper","full"), offset=FALSE, num.bins=10, round=TRUE) {standardGeneric("withinLaneNormalization")}
           )

setGeneric(
           name = "plotNtFrequency",
           def = function(x,...) {standardGeneric("plotNtFrequency")}
           )

setGeneric(
           name = "plotRLE",
           def = function(x,...) {standardGeneric("plotRLE")}
           )

setGeneric(
           name = "normCounts",
           def = function(object) {standardGeneric("normCounts")}
           )

setGeneric(
           name = "normCounts<-",
           def = function(object, value) {standardGeneric("normCounts<-")}
           )

setGeneric(
           name = "plotPCA",
           def = function(x, k=2, ...) {standardGeneric("plotPCA")}
           )
