# ***********************
# FCM Cleaning SCRIPT
# ***********************

# **********************
# Libraries importations
# **********************
library(flowCore)
library(flowAI)

# *************************************
# Cleaning one specific FCS file
# *************************************

# File name of the file to be cleaned (the file must be in the working directory)
# Note that cleaning is performed after compensation
filenameIn <- "file1.fcs"

# Reads the FCS file
fcsRaw <- read.FCS(filenameIn, truncate_max_range = FALSE) 

# Plot the principal metadata
pData(parameters(fcsRaw))

# Creation of the clean FCS file
fcsClean<- flow_auto_qc(fcsRaw) 

# File name of the clean file (it will appear in the working directory)
filenameOut <- sprintf("%s%s",substr(filenameIn,1,nchar(filenameIn)-4),"_Clean.fcs")

# Writing the output clean file (it will appear in the working directory)
write.FCS(fcsClean,filenameOut ,delimiter = "$",endian="big")
# ***********************************
# End of cleaning of one FCS file
# ***********************************


# *****************************************************
# Cleaning of all FCS file in the working directory
# *****************************************************
listFCS <- list.files(pattern = ".fcs")
for(i in 1:length(listFCS)){
  print(sprintf("File #%d over %d",i,length(listFCS)))
  filenameIn <- listFCS[i]
  fcsRaw<- read.FCS(filenameIn, truncate_max_range = FALSE)
  fcsClean<- flow_auto_qc(fcsRaw)
  filenameOut <- sprintf("%s%s",substr(filenameIn,1,nchar(filenameIn)-4),"_Clean.fcs")
  write.FCS(fcsClean,filenameOut ,delimiter = "$",endian="big")
}
# *****************************
# End of automatic Cleaning
# *****************************
