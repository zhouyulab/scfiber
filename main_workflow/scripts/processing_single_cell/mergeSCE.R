library(batchelor)
load("countSCE.RData")
no_correction2 <- function(...){
  return(correctExperiments(..., PARAM=NoCorrectParam()))
}

mergeSCE <- do.call(no_correction2, countSCE)
save(mergeSCE, file=file.path("mergeSCE.RData"))