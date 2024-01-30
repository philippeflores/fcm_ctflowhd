# *************************
# FCM pre-processing SCRIPT
# *************************

# **********************
# Libraries importations
# **********************
library(flowCore)
library(flowAI)
library(lattice)
library(flowViz)
library(flowClust)
library(flowTrans)

# *************************************
# Pre-process one specific FCS file
# *************************************

# File name of the file to be pre-processed (the file must be in the working directory)
# Note that compensation is performed first then cleaning and finally transformation
filenameIn <- "file1.fcs"

# Reads the FCS file
fcsRaw <- read.FCS(filenameIn, truncate_max_range = FALSE) 

# Plot the principal metadata
pData(parameters(fcsRaw))

# Find the SPILLOVER matrix (while NULL, try next line)
if (is.null(fcsRaw@description[["SPILL"]])) {
   if (is.null(fcsRaw@description[["SPILLOVER"]])) {
      if (is.null(fcsRaw@description[["$SPILL"]])) {
         if (is.null(fcsRaw@description[["$SPILLOVER"]])) {
         } else { strSpill <- "$SPILLOVER" }
      } else { strSpill <- "$SPILL" }
   } else { strSpill <- "SPILLOVER" }
} else { strSpill <- "SPILL" }

# Plot the compensation matrix
fcsRaw@description[[strSpill]]

# Creation of the compensated FCS file
fcsComp<- compensate(fcsRaw,fcsRaw@description[[strSpill]]) 

# Creation of the clean FCS file from the compensated file
fcsCompClean<- flow_auto_qc(fcsComp) 

# Define the transformation fitted for the data
trans <-estimateLogicle(fcsCompClean,colnames(fcsRaw)) 

# The transformation is applied to the compensated and clean file
fcsCompCleanTrans <- transform(fcsCompClean,trans) 

# File name of the pre-processed file (it will appear in the working directory)
filenameOut <- sprintf("%s%s",substr(filenameIn,1,nchar(filenameIn)-4),"_CCT.fcs")

# Writing the output pre-processed file (it will appear in the working directory)
write.FCS(fcsCompCleanTrans,filenameOut ,delimiter = "$",endian="big")
# *************************************
# End of pre-processing of one FCS file
# *************************************


# *****************************************************
# Cleaning of all FCS file in the working directory
# *****************************************************
listFCS <- list.files(pattern = ".fcs")
for(i in 1:length(listFCS)){
   print(sprintf("File #%d over %d",i,length(listFCS)))
   filenameIn <- listFCS[i]
   fcsRaw<- read.FCS(filenameIn, truncate_max_range = FALSE)
   if (is.null(fcsRaw@description[["SPILL"]])) {
     if (is.null(fcsRaw@description[["SPILLOVER"]])) {
       if (is.null(fcsRaw@description[["$SPILL"]])) {
         if (is.null(fcsRaw@description[["$SPILLOVER"]])) {
         } else { strSpill <- "$SPILLOVER" }
       } else { strSpill <- "$SPILL" }
     } else { strSpill <- "SPILLOVER" }
   } else { strSpill <- "SPILL" }
   fcsComp <- compensate(fcsRaw,fcsRaw@description[[strSpill]])
   fcsCompClean <- flow_auto_qc(fcsComp)
   trans <-estimateLogicle(fcsCompClean,colnames(fcsRaw)) 
   fcsCompCleanTrans<- transform(fcsCompClean,trans)
   filenameOut <- sprintf("%s%s",substr(filenameIn,1,nchar(filenameIn)-4),"_CCT.fcs")
   write.FCS(fcsCompCleanTrans,filenameOut ,delimiter = "$",endian="big")
}
# *****************************
# End of automatic Cleaning
# *****************************
