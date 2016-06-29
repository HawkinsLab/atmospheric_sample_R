##########################
# In order to decide which plots you would like to see , set these variables
# to either yes or no before runnning the file.
##########################
PTR_plot <- "no"
rainbow_plot <- "yes"
corrected_rainbow_plot <- "no"
log_log_plot <- "no"

##########################
# These variables can be changed. Please change the values with the variables here
# and do not change the variables lower in the script.
##########################
ref_wave <- 550   # reference wavelength subtracted to correct baseline drift
sd_limit <- 1000   # standard deviation used to remove noisy spectra







##########################
# Import the spectral data. Function returns a list with wavelengths, absorbance data
# for all files, and timestamps
##########################

data_to_dframe <- function(file, wl_low=170, wl_high=900) {
  
  if(missing(file)) # in case there's no file name input
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
    unlist(fieldList), # change a list to a vector
    nrow=length(fieldList),
    byrow=TRUE)

  ## Date and time import
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

  # As a result, returnlist has three items: wavelength, absorbance, Timestamp
  returnlist <- list("wavelength"=as.numeric(data_mat[,1]),"absorbance"=as.numeric(data_mat[,2]), "timeStamp"= date_time)
  return(returnlist) 
  
}







##########################
# Choose the data to import and format it correctly.
##########################

require(tcltk)
library(tcltk)

# Make window appear that allows user to choose which absorbance spectra to analyze
SpectraFolder <- tk_choose.dir(default = "", caption = "Select directory for absorbance spectra files")

# Get the names of the data files (including the full path) and how many files there are
files <- list.files(path=SpectraFolder, full.names=TRUE) # see full.names=TRUE
numFiles <- length(files)
# number of rows hard coded in
numRows <- 3648

# Make a matrix (currently an empty placeholder) that will store all the spectra and wavelengths
data_matrixAll <- matrix(NA,nrow=numRows,ncol=numFiles+1)
# Make a vector (currently an empty placeholder) that will store all the times
TimeSeries <- vector(length=numFiles+1)

i <- 1
for(file in files) {
  # Put the wavelengths, absorbance data, and time from one file into a temporary 2-column data frame
  tempDframe <- data_to_dframe(f=file) 
  # Put the absorbance data from this file into the full data-containing matrix
  data_matrixAll[,i+1] <- tempDframe$absorbance
  # Put the timestamp from the current file into the time-containing vector
  TimeSeries[i+1] <- tempDframe$timeStamp
  # Change the index so we can move to the next position in data_matrix to store new file info
  i <- i + 1
}

# Put the wavelength values in the data-containing matrix
data_matrixAll[,1] <- tempList$wavelength

# Make the data in the TimeSeries vector into POSIXct time class (instead of original string)
class(TimeSeries) <- c('POSIXt','POSIXct')







##########################
# Choose the type of baseline correction to apply, and apply it.
##########################

# Average the absorbance values from 360 to 370 to get an absorbance value for 365 nm
BrC365 <- colMeans(subset(data_matrixAll,data_matrixAll[,1] > 360 & data_matrixAll[,1] < 370))
# Average the absorbance values around a reference point (that is set above)
BrCref <- colMeans(subset(data_matrixAll,data_matrixAll[,1] > (ref_wave - 5) & data_matrixAll[,1] < (ref_wave + 5)))
# Remove effect of baseline drift by subtracting baseline absorbance from 365 nm absorbance
BrCcorr <- BrC365-BrCref #closer to actual signal we want






##########################
# Ask for SMPS data. If it exists, read it in and format it correctly in a data
# frame, then interpolate so that we can get data points at the times we sampled at.
# Calculate MAC at 365 nm.
##########################

SMPS_check <- readline(prompt="Do you have SMPS file corresponding to Daily Spectra? y or n: ")

if (SMPS_check == "y"){
  
  # Ask what type of file it is (raw or processed)
  SMPS_analysis_type <- readline(prompt = "Is the SMPS data file raw data (r) or processed data (p)? ")

  }

if (SMPS_check == "y" && SMPS_analysis_type == "p") {
  
  # read in SMPS processed data by allowing user to select the file
  SMPSfile <- tk_choose.files(default="",caption="Select a tab-delimited SMPS file with mm/dd/yyyy format")
  
  # format SMPS data
  SMPS <- read.table(SMPSfile, sep="\t", header=TRUE,stringsAsFactors = FALSE)
  SMPS$smpstimeFormatted <- as.POSIXct(SMPS$smpstime, format="%m/%d/%Y %H:%M")
  
  # use the approx function to interpolate
  InterSMPS <- approx(SMPS$smpstimeFormatted, SMPS$smpsconc, TimeSeries, method = "linear", rule = 1, f = 0, ties = mean)
  
}

if (SMPS_check == "y" && SMPS_analysis_type == "r") {
  
  # read in SMPS raw data by allowing user to select the file
  SMPS_testFile <- tk_choose.files(default="",caption="Select a raw SMPS file")
  
  # format SMPS data
  #### following note seems to be irrelevant as written, with colnames removed
  #### NOTE: This assumes length of file is 136, thus including all the columns because otherwise it would think there were only 2 columns
  SMPS <- read.csv(SMPS_testFile, skip = 17, header = FALSE, sep="\t", fill=TRUE)   #col.names = paste0(seq_len(136))
  SMPS_conc <- as.numeric(as.character(SMPS$V136 ))    # total concentration, need to convert to class form 
  SMPS_datetime<- as.POSIXct(paste(SMPS$V2, SMPS$V3), format="%m/%d/%y %H:%M")
  SMPS.df <- data.frame(SMPS_datetime, SMPS_conc)
  
  # check if line 15 column 2 is "dw/dlogDp"
  if (SMPS$V2[15]=="dw/dlogDp"){
    num <- grepl("^[0-9]",SMPS$V1)    # if this is true, find all the lines starting with a number
    SMPS.df <- SMPS.df[num,]        # update date frame to just include lines starting with a number
    } 
  
  # use the approx function to interpolate
  InterSMPS <- approx(SMPS.df$SMPS_datetime, SMPS.df$SMPS_conc, TimeSeries, method = "linear", rule = 1, f = 0, ties = mean)
  
}

if (SMPS_check == "y") {

  # calculate MAC at 365 nm using baseline-corrected absorbance  
  MAC <- (BrCcorr*1329787)/InterSMPS$y  # obscure number comes from unit and dilution correction (page 39) in HGW lab notebook
  
}







############
## We create a POSIXct date/time with the experiment date and entered desired start time. The difference 
## between actual start time and desired start time is subtracted from every date/time.
############

# Ask user to enter an experiment reference start time in HH:MM:SS format. 
# Ex: to start at 2 pm , enter 14:00:00. 
time_plot <- readline(prompt="Reference time? For 00:00:00, press enter. For other time, enter as HH:MM:SS. ")

if (time_plot == "\n") {
  time_plot <- "00:00:00"
}

Time <- TimeSeries[2:numFiles+1]    # rename time series vector
getDate <- as.Date(head(Time, n=1))    # get date of first datetime stamp
date_time <- paste(getDate, time_plot)    # create a string with the experiment date and user entered start time
refTime <- as.POSIXct(date_time)    # convert to POSIXct date / time 
difference <- as.numeric(difftime(head(Time, n=1), refTime), units="secs")    # difference between actual and desired start; seconds can be subtracted directly from POSIX class
CorrectedTime_Ref <- Time - difference    # new referenced times; original time vector minus the difference between actual and desired start      







#######################
## plot for BrCcorr and MAC at local Paris time
########################
require(ggplot2)
require(grid)
library(ggplot2)
library(grid)

# create a grob (grid graphical object) with the current date to add to plots 
date_grob = grobTree(textGrob(getDate, x=0.1,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=15, fontface="italic")))

# time series vs BrC correction 
qplot1 <- qplot(Time, BrCcorr[2:numFiles+1], colour="red", geom = "line", 
                xlab="Local Time (Paris)", ylab="Corrected Absorbance at 365 nm", show.legend=FALSE) + annotation_custom(date_grob) + theme_bw()
print(qplot1) 

# time series vs MAC, only if you have access to SMPS file
if (SMPS_check == "y") {
  qplot2 <- qplot(Time, MAC[2:numFiles+1], colour="green", geom = "line",
                  xlab="Local Time (Paris)", ylim=c(-100,3500), ylab="MAC", show.legend=FALSE) + annotation_custom(date_grob) + theme_bw()
  print(qplot2)
} 

#######################
## plots for BrCcorr and MAC at reference time 
########################

qplot3 <- qplot(CorrectedTime_Ref, BrCcorr[2:numFiles+1], colour ="red", geom="line",
                xlab=paste("Time since ", time_plot), ylab="Corrected Absorbance at 365 nm", show.legend=FALSE)+ annotation_custom(date_grob) + theme_bw()
print(qplot3) 

if (SMPS_check == "y"){
  qplot4 <- qplot(CorrectedTime_Ref, MAC[2:numFiles+1], colour ="red", geom="line",
                  xlab=paste("Time since", time_plot), ylab="MAC", show.legend=FALSE)+ annotation_custom(date_grob) + theme_bw()
  print(qplot4) 
}









