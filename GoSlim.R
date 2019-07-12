library("GSEABase")
dir <- setwd("~/directory")

go <- read.csv("../file", sep = "\t", header = FALSE ) # Reading in big annotation file
goint <- subset(go, go$V1 %in% candidate) # Filtering by  a list of candidates
vect <- as.character(goint$V2) # Creating a list of your GO terms
myCollection <- GOCollection(vect)
fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
slim <- getOBOCollection(fl)
goSlim(myCollection, slim, "MF") # MF - molecular function, BP - biological process and CC - cellular component

#-----GO IDs to abstracts-----
library(GO.db)
goterms <- as.data.frame(Term(GOTERM))

files <- c("axsGO","exsGO","pasGO")
for (i in files) 
{
output <- as.data.frame(matrix(0, ncol = 2))
golist <- readLines(paste("../../19.04.16_clusters/newclusters/styles/", i, sep=""))
myCollection <- GOCollection(golist)
shady <- goSlim(myCollection, slim, "MF")
for (j in golist) 
  {
output[nrow(output)+1,1] <- j
tryCatch({
output[nrow(output),2] <- paste(goterms[grep(j, row.names(goterms)),])
}, error=function(e){})
  }
output <- output[-1,]
write.table(output, file = paste("../../", i,  sep=""), sep="\t", col.names=FALSE, row.names = FALSE, quote = FALSE)
write.table(shady, file = paste("../../", i, "_slim",  sep=""), sep="\t", col.names=FALSE, row.names = FALSE, quote = FALSE)
}