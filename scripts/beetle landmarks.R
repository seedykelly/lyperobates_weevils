library(geomorph)
rm(list=ls())
setwd("/Users/clintkelly1/Documents/Projects/Fungus weevils/nts")
digitize2d("03F-Pair2.jpg", 11, 5)
digitize2d("01F-Pair1.jpg", 11, 5)

test.list<-list.files("/Users/clintkelly1/Documents/Projects/Fungus weevils/nts",pattern=".nts")

readland.nts("01F-Pair12dcoords.nts")

readmulti.nts(filelist="test.list")

readmulti.nts(test.list)
readmulti.nts()