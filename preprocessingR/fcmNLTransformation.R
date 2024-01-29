# ************************************
# FCM Non-Linear Transformation SCRIPT
# ************************************

# **********************
# Libraries importations
# **********************
library(flowCore)
library(lattice)
library(flowViz)
library(flowClust)
library(flowTrans)

# ***************************************
# Transformation of one specific FCS file
# ***************************************

# File name of the file to be transformed (the file must be in the working directory)
# Note that the non-linear is applied after compensation and cleaning
filenameIn <- "file1.fcs"

# Reads the FCS file
fcsRaw <- read.FCS(filenameIn, truncate_max_range = FALSE) 

# Plot the principal metadata
pData(parameters(fcsRaw))

# Define the transformation fitted for the data
trans <-estimateLogicle(fcsRaw,colnames(fcsRaw)) 

# The transformation is applied
fcsTrans <- transform(fcsRaw,trans) 

# File name of the transformed file (it will appear in the working directory)
filenameOut <- sprintf("%s%s",substr(filenameIn,1,nchar(filenameIn)-4),"_Trans.fcs")

# Writing the output transformed file (it will appear in the working directory)
write.FCS(fcsTrans,filenameOut ,delimiter = "$",endian="big")
# *************************************
# End of transformation of one FCS file
# *************************************


# *******************************************************
# Transformation of all FCS file in the working directory
# *******************************************************
listFCS <- list.files(pattern = ".fcs")
for(i in 1:length(listFCS)){
   print(sprintf("File #%d over %d",i,length(listFCS)))
   filenameIn <- listFCS[i]
   fcsRaw<- read.FCS(filenameIn, truncate_max_range = FALSE)
   trans <-estimateLogicle(fcsRaw,colnames(fcsRaw)) 
   fcsTrans<- transform(fcsRaw,trans)
   filenameOut <- sprintf("%s%s",substr(filenameIn,1,nchar(filenameIn)-4),"_Trans.fcs")
   write.FCS(fcsTrans,filenameOut ,delimiter = "$",endian="big")
}
# *******************************
# End of automatic transformation
# *******************************
