# Hawkins Lab Data Analysis

##########################
# In order to decide which plots you would like to see , set these variables
# to either yes or no before runnning the file.
##########################

PTR_plot <- "no"
rainbow_plot <- "yes"
corrected_rainbow_plot <- "no"    # cannot mark this as yes if rainbow_plot is no
log_log_plot <- "no"
save_graphs <- "yes"

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
  SMPSfile <- tk_choose.files(default="", caption="Select a tab-delimited SMPS file with mm/dd/yyyy format")
  
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
## plot BrCcorr and MAC at local Paris time
##########################
require(ggplot2)
require(grid)
library(ggplot2)
library(grid)

# create a PDF file to save plot to
if (save_graphs == "yes" ) {
  pdf("corrected abs at 365 vs real time.pdf")
}

# create a grob (grid graphical object) with the current date to add to plots 
date_grob = grobTree(textGrob(getDate, x=0.1,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=15, fontface="italic")))
# corrected absorbance vs time series
qplot1 <- qplot(Time, BrCcorr[2:numFiles+1], colour="red", geom = "line", 
                xlab="Local Time (Paris)", ylab="Corrected Absorbance at 365 nm", show.legend=FALSE) + annotation_custom(date_grob) + theme_bw()
print(qplot1) 

# turn off PDF save
if (save_graphs == "yes" ) {
  dev.off()
}

# create a PDF file to save plot to
if (save_graphs == "yes" ) {
  pdf("MAC vs real time.pdf")
}

# MAC vs time series, only if you have access to SMPS file
if (SMPS_check == "y") {
  qplot2 <- qplot(Time, MAC[2:numFiles+1], colour="green", geom = "line",
                  xlab="Local Time (Paris)", ylim = c(0,2000), ylab="MAC", show.legend=FALSE) + annotation_custom(date_grob) + theme_bw()
  print(qplot2)
} 

# turn off PDF save
if (save_graphs == "yes" ) {
  dev.off()
}






##########################
## plots for BrCcorr and MAC at reference time 
##########################

# create a PDF file to save plot to
if (save_graphs == "yes" ) {
  pdf("corrected absorbance at 365 vs referenced time.pdf")
}

# corrected absorbance vs time series
qplot3 <- qplot(CorrectedTime_Ref, BrCcorr[2:numFiles+1], colour ="red", geom="line",
                xlab=paste("Time since ", time_plot), ylab="Corrected Absorbance at 365 nm", show.legend=FALSE)+ annotation_custom(date_grob) + theme_bw()
print(qplot3) 

# turn off PDF save
if (save_graphs == "yes" ) {
  dev.off()
}

# create a PDF file to save plot to
if (save_graphs == "yes" ) {
  pdf("MAC vs referenced time.pdf")
}

# MAC vs time series, only if you have access to SMPS file
if (SMPS_check == "y"){
  qplot4 <- qplot(CorrectedTime_Ref, MAC[2:numFiles+1], colour ="red", geom="line",
                  xlab=paste("Time since", time_plot), ylab="MAC", show.legend=FALSE)+ annotation_custom(date_grob) + theme_bw()
  print(qplot4) 
}

# turn off PDF save
if (save_graphs == "yes" ) {
  dev.off()
}







##########################
# Plot of MAC and ion counts (from PTR Corr Series file) over time
# Plots only if user has SMPS file
##########################

if (PTR_plot=="yes") {
  if (SMPS_check == "y"){
    
    require("reshape2")
    library("reshape2")
    library("fields")
    
    # read in PTR Corr Series file, must be tab-demilited
    PTR_Series <- tk_choose.files(default="",caption="Select a tab-delimited PTR Corr Series file")
    PTR <- read.table(PTR_Series, sep="\t", header=TRUE, stringsAsFactors = FALSE)
    
    # format time to POSIXct
    # Note: PTR files have different formatting  
    #PTR_time <- as.POSIXct(PTR$Time, format="%m/%d/%y %H:%M")  # if date/time is date & time 
    PTR_time <- as.POSIXct(PTR$Time, format="%H:%M:%S")   # if date/time is just time 
    
    # create data frame with time, MA, MG, Imine
    PTR.df <- data.frame(PTR_time, PTR$MA, PTR$MG, PTR$Imine)
    
    # create data frame of MAC data with Time Series
    MAC.df <- data.frame(Time, MAC[2:numFiles+1])
    colnames(MAC.df) <- c("Time", "MAC") # rename columns 
    
    # create a PDF file to save plot to
    if (save_graphs == "yes" ) {
      pdf("PTR.pdf")
    }
    
    par(new=F)
    layout(mat=1)
    # set margins to fit both y-axis
    par(mar = c(4,6,4,6))
    
    # get min and max counts of PTR counts for plotting 
    min_counts <- c(min(PTR.df$PTR.MA, na.rm = TRUE), max(PTR.df$PTR.MG, na.rm = TRUE), max(PTR.df$PTR.Imine, na.rm = TRUE))
    max_counts <- c(max(PTR.df$PTR.MA, na.rm = TRUE), max(PTR.df$PTR.MG, na.rm = TRUE), max(PTR.df$PTR.Imine, na.rm = TRUE))
    
    # Plot of PTR counts versus time
    # Each line is plotted 
    with(PTR.df, plot(PTR_time, PTR.MA, col="blue", type="l", axes=F, xlab=NA, ylab=NA, ylim=c(min(min_counts), max(max_counts))))
    lines(PTR_time, PTR.df$PTR.MG, col="red")
    lines(PTR_time, PTR.df$PTR.Imine, col="green")
    
    # Add axis with label on right side of plot 
    axis(side = 4)
    mtext("Ion count", side = 4, line = 2, las=2)
    
    par(new = T) # indicates that next plotting command should not clear the plot
    
    # plot MAC versus time, add in y-axis label after
    with(MAC.df, plot(Time, MAC, type="l", col="black", xlab="Time", ylab=NA))
    mtext("MAC", line=3, las=2, side=2)
    # add to an existing plot using new=T
    
    # Add in legend on bottomright of plot
    legend("bottomright",
           col = c("blue", "red", "green", "black"),
           legend=c("MA", "MG", "Imine", "MAC"),
           lty=c(1))
    # add date on top of plot
    mtext(getDate, side=3, line=2) 
    
    par(new=F)
    
    # turn off PDF save
    if (save_graphs == "yes" ) {
      dev.off()
    }
    
  }} # end of check for PTR plots








##########################
# Rainbow plot -- Absorbance versus wavelength for each time
# Plots only if user has SMPS file
##########################

library("fields")
library("maps")
library("spam")

if (rainbow_plot == "yes") {
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
  
  # create a PDF file to save plot to
  if (save_graphs == "yes" ) {
    pdf("rainbow_plot_all.pdf")
  }

  # plot
  grid.newpage()    # create new page for plot
  layout(t(1:2), widths=c(10,2))    # set up the layout of the plot
  my.palette <- rainbow(length(justTime), start=0, end=4/6)      # create rainbow colors with length of time vector, from red (start=0) to blue (end=4/6)
  matplot(spectra_corr[,1], spectra_corr[,-1], type="l", xlim=c(300,700), ylim=c(-0.15,.45), 
          xlab="wavelength", ylab="absorbance", col = my.palette) 
  image.plot(smallplot= c(.99,1,0.1,.9), zlim=c(Time_asHours[2],Time_asHours[length(Time_asHours)]), 
             legend.only=TRUE, horizontal = FALSE, col=my.palette, legend.lab="Local Time")   # add color bar to plot 
  mtext(getDate, side=3)    # add date to the plot
  
  # turn off PDF save
  if (save_graphs == "yes" ) {
    dev.off()
  }
}






##########################
# Rainbow plot with noisy spectra removed
# Currently removing spectra with a standard deviation over sd_limit (assigned
# above) over the region 480 nm to 530 nm
# Also makes a plot to show which spectra were removed
# Plots only if user has SMPS file
##########################

if (corrected_rainbow_plot == "yes") {
  
  num_cols <- ncol(spectra_corr) - 1    # Number of columns to loop through
  times <- vector(length=num_cols)    # empty vector for times
  removed <- vector(length=num_cols)     # assigns 0 or 1 to each time depending on whether or not it is removed
  matrix_corr <- spectra_corr    # make new matrix filled with all the corrected spectra
  
  for (i in 2:num_cols) {
    
    times[i-1] <- colnames(spectra_corr)[i]    # fill in matrix with times
    
    if (sd(spectra_corr[1450:1700,i]) > sd_limit) {   # choose spectra with standard deviation over sd_limit
      matrix_corr <- matrix_corr[,-i]    # remove that spectrum from the matrix of all spectra
      removed[i-1] <- 1    # assign 1 to removed time
    } else {
      removed[i-1] <- 0    # assign 0 to all other times
    }
    
  }
  
  # create a PDF file to save plot to
  if (save_graphs == "yes" ) {
    pdf("rainbow_plot_corrected.pdf")
  }
  
  # plot
  grid.newpage()
  par(mar= c(5, 4, 4, 2))
  plot(removed, xlab="Time", yaxt="n")    # no y-axis
  axis(2, at = seq(0, 1, by = 1), las=2)    # add y-axis
  
  # add wavelength vector to corrected matrix
  matrix_corr <- cbind(spectra_corr[,1], matrix_corr)
  
  #### Plots new matrix with same parameters as above 
  # set new plot
  grid.newpage()
  
  # layout to fit both plot & color bar legend
  layout(t(1:2), widths=c(10,2))
  
  # create rainbow colors with length of time vector, from red (start=0) to blue (end=4/6)
  my.palette <- rainbow(length(justTime), start=0, end=4/6)
  
  # plot the matrix
  matplot(matrix_corr[,1], matrix_corr[,-1], type="l", xlim=c(300,700), ylim=c(-0.15,.45), 
          xlab="wavelength", ylab="absorbance", col = my.palette) 
  
  # add color bar to plot 
  image.plot(smallplot= c(.99,1,0.1,.9), zlim=c(Time_asHours[2],Time_asHours[length(Time_asHours)]), 
             legend.only=TRUE, horizontal = FALSE, col=my.palette, legend.lab="Local Time")
  # add date to the plot
  mtext(getDate, side=1) 
  
  # create a PDF file to save plot to
  if (save_graphs == "yes" ) {
    dev.off()
  }
  
} # end of check for corrected_rainbow_plot






##########################
# log-log plot of absorptivity vs wavelength
# Plots only if user has SMPS file
##########################

if (log_log_plot == "yes"){
  
  matrix_log <- spectra_corr    # create new matrix by copying the corrected spectra 
  matrix_log[matrix_log < 0] <- 0    # set any value less than 0 to 0
  matrix_log <- log(matrix_log)    # take natural log of entire matrix
  matrix_log[is.infinite(matrix_log)] <- 0    # set any infinite (negative) value to 0 in order to take mean of each row
  row_Means <- rowMeans(matrix_log[,2:ncol(matrix_log)], na.rm = TRUE)    # take mean along rows (at each wavelength)
  
  # create a PDF file to save plot to
  if (save_graphs == "yes" ) {
    pdf("log-log.pdf")
  }
  
  # plot
  grid.newpage()    # make new page to plot on
  layout(t(1:2), widths=c(10,2))    # set up the layout of the graph
  
  # plot this log/log matrix
  matplot(matrix_log[,1], matrix_log[,-1], type="l", xlim = c(log(300), log(500)),
          xlab="log(wavelength)", ylab="log(absorptivity)", col = my.palette) 
  
  # add line of mean absorbance vs log wavelength
  lines(matrix_log[,1], row_Means) # line of wavelength versus mean
  
  # add color bar to plot 
  image.plot(smallplot= c(.99,1,0.1,.9), zlim=c(Time_asHours[2],Time_asHours[length(Time_asHours)]), 
             legend.only=TRUE, horizontal = FALSE, col=my.palette, legend.lab="Local Time")
  
  # create a PDF file to save plot to
  if (save_graphs == "yes" ) {
    dev.off()
  }
  
} # end of check for log/log plot





