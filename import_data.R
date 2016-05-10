data_to_dframe <- function(file='no_input', wl_low=170, wl_high=900) {
  # Process raw data file in tab delimited format 
  # Args:
  #   file: data file name 
  #   wl_low: lower bound of wavelength
  #   wl_high: upper bound of wavelength
  #
  # Returns: 
  #   A data Frame containing three coloumns, wavelength, 
  #    absorbance, and timestamp
  # 
  #######################
  ## Opening the individual text files, getting wavelength and absorbance but not date
  ###########################
  if(file=='no_input') # in case there's no file name input
    file = file.choose()
  # read the datafile line by line
  txt <- readLines(file) 
  # find all the lines starting with a number
  ind <- grepl("^[0-9]",txt)
  # create a new array only with the lines starting with a number
  data_str <- txt[ind]
  # split each line into pieces
  fieldList <- strsplit(data_str, split = "\t")
  # create a matrix 
  data_mat <- matrix(
    unlist(fieldList), 
    nrow=length(fieldList),
    byrow=TRUE)
 
  ########################################
  ## Date and time import
  ## result: date_time (date-time object)
  ########################################
  # extract the date line at line 3 and save it in the POSIXct (date/time) class format
  datetime_line <- strsplit(txt[3], split=" ")
  datetime_str <- paste( # Date format %b%d%Y 
    datetime_line[[1]][3], # month, abbreviated %b
    datetime_line[[1]][4], # day, [01-31] %d
    datetime_line[[1]][7], # year, four digits %Y
    sep="")
  # first convert the strings into date objects.
  datetime_obj <- as.Date(datetime_str, "%b%d%Y"); 
  # now merge the date and time into a single POSIXct object
  date_time <- as.POSIXct(paste(datetime_obj, datetime_line[[1]][5]), format="%Y-%m-%d %H:%M:%S")
  
  ########################################
  ## Adding the timestamp to list to return in function
  ## Result: data_frm (a data frame object)
  ########################################
 
  # As a result, list has three items: wavelength, absorbance, Timestamp

  
  returnlist <- list("wavelength"=as.numeric(data_mat[,1]),"absorbance"=as.numeric(data_mat[,2]), "timeStamp"= date_time)
  return(returnlist) 

}



ExptDay <- readline(prompt="Enter the expt day as CESAM_YYMMDD: ")
SpectraFolder <- tk_choose.dir(default = "", caption = "Select directory for absorbance spectra files")
files <- list.files(path=SpectraFolder, pattern=ExptDay, full.names=TRUE)

numFiles <- length(files)

#Before this lines runs, we need to get the nrow for our spectra rather than hard code it
#then we can use the wavelength range cutoffs. To do this, we need only read in the first column 
#wavelength, not spectra.

data_matrixAll <- matrix(NA,nrow=3648,ncol=numFiles+1) # a place holder (an empty data frame obj), 
TimeSeries <- vector(length=numFiles+1)


i = 1 # index for data_list in the following for loop
for(file in files) 
{
  # we are trying to pull, file by file, absorbances into a frame, which are contained in the first element of newlist from function above
  tempList <- data_to_dframe(f=file) 

  data_matrixAll[,i+1] <- tempList$absorbance
  data_matrixAll[,1] <- tempList$wavelength
  
  TimeSeries[i+1] <- tempList$timeStamp
 # data_list_uv_mean[[i]] <- lapply(data_list_uv[[i]],mean)
  i = i + 1
}
class(TimeSeries) <- c('POSIXt','POSIXct') #reestablishes as time series

#Now we will take the average of the absorption around the brown carbon region, keeping the time series intact.
BrC365 <- colMeans(subset(data_matrixAll,data_matrixAll[,1] > 360 & data_matrixAll[,1] < 370))
#Now we take the average absorption in some long wavelength reference region, usually 700 nm but recently at CESAM
#we found that 600 or 550 is better due to fluctuating signal 700
BrCref <- colMeans(subset(data_matrixAll,data_matrixAll[,1] > 695 & data_matrixAll[,1] < 705))

#We must now subtract the absorbance at the reference wavelength from the BrC wavelength
BrCcorr <- BrC365-BrCref #closer to actual signal we want

#############################
## Now we begin using SMPS data to normalize absorbance. We read the file in, format time
#########################

#read in particle concentration data, by allowing user to select the file
SMPSfile <- tk_choose.files(default="",caption="Select a tab-delimited SMPS file with mm/dd/yyyy format")
SMPS <- read.table(SMPSfile, sep="\t", header=TRUE,stringsAsFactors = FALSE)

SMPS$smpstimeFormatted <- as.POSIXct(SMPS$smpstime, format="%m/%d/%Y %H:%M")

#use the approx function to interpolate
InterSMPS <- approx(SMPS$smpstimeFormatted, SMPS$smpsconc, TimeSeries, method = "linear", rule = 1, f = 0, ties = mean)

#####################
#Finally, we normalize the signal of brown carbon coming from the difference of Abs at 365 nm and 
#some reference wavelength (typically 700 nm but in Paris 2015 we found that 550 was more suitable
#due to fluctuations at 700 for unknown reasons)
########################

MAC <- BrCcorr/InterSMPS$y
plot(TimeSeries[2:numFiles+1],BrCcorr[2:numFiles+1], col="red", type = "l")
plot(TimeSeries[2:numFiles+1], MAC[2:numFiles+1], col="green", type = "l",ylim = c(0,0.005))




