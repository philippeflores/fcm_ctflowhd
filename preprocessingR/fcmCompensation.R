# ***********************
# FCM Compensation SCRIPT
# ***********************

# **********************
# Libraries importations
# **********************
library(flowCore)

# *************************************
# Compensation of one specific FCS file
# *************************************

# File name of the file to be compensated (the file must be in the working directory)
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

# File name of the compensated file (it will appear in the working directory)
filenameOut <- sprintf("%s%s",substr(filenameIn,1,nchar(filenameIn)-4),"_Comp.fcs")

# Writing the output compensated file (it will appear in the working directory)
write.FCS(fcsComp,filenameOut ,delimiter = "$",endian="big")
# ***********************************
# End of compensation of one FCS file
# ***********************************

# *****************************************************
# Compensation of all FCS file in the working directory
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
   fcsComp<- compensate(fcsRaw,fcsRaw@description[[strSpill]])
   filenameOut <- sprintf("%s%s",substr(filenameIn,1,nchar(filenameIn)-4),"_Comp.fcs")
   write.FCS(fcsComp,filenameOut ,delimiter = "$",endian="big")
}
# *****************************
# End of automatic compensation
# *****************************
