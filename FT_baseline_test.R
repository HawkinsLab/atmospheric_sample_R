## Elyse Pennington
## new baselining routine: Fourier transform to remove light brightening
## 1) perform a Fourier transform on each wavelength over time
## 2) remove low frequency values (caused by slow changes in the lamp brightness
##    over time)
## 3) perform a reverse Fourier transform to convert the frequency data for each
##    wavelength back to time data



# Hawkins Lab Data Analysis

##########################
# In order to decide which plots you would like to see , set these variables
# to either yes or no before runnning the file.
##########################
PTR_plot <- "no"
rainbow_plot <- "yes"
corrected_rainbow_plot <- "yes"    # cannot markt this as yes if rainbow_plot is no
log_log_plot <- "no"

one_wave_bl <- "yes"
FT_bl <- "yes"

##########################
# These variables can be changed. Please change the values with the variables here
# and do not change the variables lower in the script.
##########################
ref_wave <- 550   # reference wavelength subtracted to correct baseline drift
sd_limit <- 0.005   # standard deviation used to remove noisy spectra







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

########################## old
# Average the absorbance values from 360 to 370 to get an absorbance value for 365 nm
BrC365 <- colMeans(subset(data_matrixAll,data_matrixAll[,1] > 360 & data_matrixAll[,1] < 370))
# Average the absorbance values around a reference point (that is set above)
BrCref <- colMeans(subset(data_matrixAll,data_matrixAll[,1] > (ref_wave - 5) & data_matrixAll[,1] < (ref_wave + 5)))
# Remove effect of baseline drift by subtracting baseline absorbance from 365 nm absorbance
BrCcorr <- BrC365-BrCref #closer to actual signal we want
########################## old


for (i in 1:numRows) {
  freq_data <- fft(data_matrixAll[i,])
  mag0 <- sqrt( Re(freq_data)^2 + Im(freq_data)^2 ) * 2 / length(freq_data)
  mag <- fftshift(mag0)
  #shifted_time <- fftshift(TimeSeries)
  FTdata_matrixAll[i] <- fft(freq_data, reverse = TRUE)
}










##########################
## We create a POSIXct date/time with the experiment date and entered desired start time. The difference 
## between actual start time and desired start time is subtracted from every date/time.
##########################

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






##########################
# Rainbow plot -- Absorbance versus wavelength for each time
# Plots only if user has SMPS file
##########################

library("fields")
library("maps")
library("spam")

if (rainbow_plot == "yes" && one_wave_bl == "yes") {
  # create new time vector of just times, but keeping the first time stamp as a placeholder
  justTime <- strftime(TimeSeries, format = "%H:%M:%S")
  
  # create new matrix with same dimensions as data matrix
  spectra_corr <- matrix(NA, nrow=numRows, ncol=numFiles+1)   
  
  # rename columns using time series of just times 
  colnames(spectra_corr) <- justTime
  
  # Fill spectra_corr with wavelength vector and corrected spectra
  spectra_corr[,1] <- data_matrixAll[,1]    # put wavelengths in first column
  for (i in 2:(length(BrCref))){
    # put corrected absorbances in other columns
    spectra_corr[,i] <- data_matrixAll[,i] - BrCref[i]
  }
  
  # function to convert HH:MM:SS to hours in order to plot as legend
  Time_asHours <- sapply(strsplit(justTime,":"),
                         function(x) {
                           x <- as.numeric(x)
                           ((x[1]+x[2]/60) + (x[3]/3600)) })
  
  # plot
  grid.newpage()    # create new page for plot
  layout(t(1:2), widths=c(10,2))    # set up the layout of the plot
  my.palette <- rainbow(length(justTime), start=0, end=4/6)      # create rainbow colors with length of time vector, from red (start=0) to blue (end=4/6)
  matplot(spectra_corr[,1], spectra_corr[,-1], type="l", xlim=c(300,700), ylim=c(-0.15,.45), 
          xlab="wavelength", ylab="absorbance", col = my.palette) 
  image.plot(smallplot= c(.99,1,0.1,.9), zlim=c(Time_asHours[2],Time_asHours[length(Time_asHours)]), 
             legend.only=TRUE, horizontal = FALSE, col=my.palette, legend.lab="Local Time")   # add color bar to plot 
  mtext(getDate, side=3)    # add date to the plot
}







if (rainbow_plot == "yes" && FT_bl == "yes") {
  # create new time vector of just times, but keeping the first time stamp as a placeholder
  justTime <- strftime(TimeSeries, format = "%H:%M:%S")
  
  # create new matrix with same dimensions as data matrix
  spectra_corr <- matrix(NA, nrow=numRows, ncol=numFiles+1)   
  
  # rename columns using time series of just times 
  colnames(spectra_corr) <- justTime
  
  # Fill spectra_corr with wavelength vector and corrected spectra
  spectra_corr[,1] <- data_matrixAll[,1]    # put wavelengths in first column
  for (i in 2:(length(BrCref))){
    # put corrected absorbances in other columns
    spectra_corr[,i] <- data_matrixAll[,i] - BrCref[i]
  }
  
  # function to convert HH:MM:SS to hours in order to plot as legend
  Time_asHours <- sapply(strsplit(justTime,":"),
                         function(x) {
                           x <- as.numeric(x)
                           ((x[1]+x[2]/60) + (x[3]/3600)) })
  
  # plot
  grid.newpage()    # create new page for plot
  layout(t(1:2), widths=c(10,2))    # set up the layout of the plot
  my.palette <- rainbow(length(justTime), start=0, end=4/6)      # create rainbow colors with length of time vector, from red (start=0) to blue (end=4/6)
  matplot(spectra_corr[,1], spectra_corr[,-1], type="l", xlim=c(300,700), ylim=c(-0.15,.45), 
          xlab="wavelength", ylab="absorbance", col = my.palette) 
  image.plot(smallplot= c(.99,1,0.1,.9), zlim=c(Time_asHours[2],Time_asHours[length(Time_asHours)]), 
             legend.only=TRUE, horizontal = FALSE, col=my.palette, legend.lab="Local Time")   # add color bar to plot 
  mtext(getDate, side=3)    # add date to the plot
}






# plot the original spectra
grid.newpage()    # create new page for plot
layout(t(1:2), widths=c(10,2))    # set up the layout of the plot
my.palette <- rainbow(length(justTime), start=0, end=4/6)      # create rainbow colors with length of time vector, from red (start=0) to blue (end=4/6)
matplot(data_matrixAll[,1], data_matrixAll[,-1], type="l", xlim=c(300,700), ylim=c(-0.15,.45), 
        xlab="wavelength", ylab="absorbance", col = my.palette) 
image.plot(smallplot= c(.99,1,0.1,.9), zlim=c(Time_asHours[2],Time_asHours[length(Time_asHours)]), 
           legend.only=TRUE, horizontal = FALSE, col=my.palette, legend.lab="Local Time")   # add color bar to plot 
mtext(getDate, side=3)    # add date to the plot

