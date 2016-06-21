source("data_to_list.R")

## ask user if they have SMPS data, if 'no' plot MAC, if 'yes' plot both BrCcorr & MAC
SMPS_check <- readline(prompt="Do you have SMPS file corresponding to Daily Spectra? Enter yes/no:")

#ExptDay <- readline(prompt="Enter the expt day as CESAM_YYMMDD: ")  # removed pattern requirement from list files function

require(tcltk)
library(tcltk)
SpectraFolder <- tk_choose.dir(default = "", caption = "Select directory for absorbance spectra files")

# Get the names of the data files including the full path
files <- list.files(path=SpectraFolder, full.names=TRUE) # see full.names=TRUE
numFiles <- length(files)

#Before this lines runs, we need to get the nrow for our spectra rather than hard code it
#then we can use the wavelength range cutoffs. To do this, we need only read in the first column 
#wavelength, not spectra.
#require("R.utils")
#nheaders <- 17 # 17 lines of header
#numRows <- countLines(files[1]) - nheaders
numRows <- 3648
data_matrixAll <- matrix(NA,nrow=numRows,ncol=numFiles+1) # a place holder (an empty data frame obj), 
TimeSeries <- vector(length=numFiles+1)

i <- 1 # index for data_list in the following for loop
for(file in files) 
{
  # we are trying to pull, file by file, absorbances into a frame, which are contained in the first element of newlist from function above
  tempList <- data_to_list(f=file) # use wl_low and wl_high
  
  data_matrixAll[,i+1] <- tempList$absorbance
  data_matrixAll[,1] <- tempList$wavelength
  
  TimeSeries[i+1] <- tempList$timeStamp
  # data_list_uv_mean[[i]] <- lapply(data_list_uv[[i]],mean)
  i <- i + 1
}
class(TimeSeries) <- c('POSIXt','POSIXct') #reestablishes as time series

