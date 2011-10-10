setClass(
         Class = "SeqExpressionSet",
         contains = "eSet",
         prototype = prototype(new("VersionedBiobase",
           versions = c(classVersion("eSet"), SeqExpressionSet = "1.0.0")))
         )
