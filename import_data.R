# Hawkins Lab Data Analysis

##########################
# In order to decide which plots you would like to see, set these variables
# to either yes or no before runnning the file.
##########################

PTR_plot <- "no"
rainbow_plot <- "yes"
corrected_rainbow_plot <- "yes"    # cannot mark this as yes if rainbow_plot is no
log_log_plot <- "yes"    # cannot mark this as yes if rainbow_plot & corrected_rainbow_plot are no
save_graphs <- "yes"

##########################
# These variables can be changed. Please change the values with the variables here
# and do not change the variables lower in the script.
##########################
ref_wave <- 550   # reference wavelength subtracted to correct baseline drift
sd_limit <- 0.005   # standard deviation used to remove noisy spectra

##########################
# Graphing paramters
##########################
MAC_ymin <- 0
MAC_ymax <- 800
corr_abs_ymin <- 0
corr_abs_ymax <- 0.022
rainbow_ymin <- -0.02
rainbow_ymax <- 0.09
log_ymin <- 0.0005   # this number must be larger than 0
log_ymax <- 0.1

# import -----------
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

# import and format -----------
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
data_matrixAll[,1] <- tempDframe$wavelength

# Make the data in the TimeSeries vector into POSIXct time class (instead of original string)
class(TimeSeries) <- c('POSIXt','POSIXct')

# baseline calcs -----------
##########################
# Choose the type of baseline correction to apply, and apply it.
##########################

# Average the absorbance values from 360 to 370 to get an absorbance value for 365 nm
BrC365 <- colMeans(subset(data_matrixAll,data_matrixAll[,1] > 360 & data_matrixAll[,1] < 370))
# Average the absorbance values around a reference point (that is set above)
BrCref <- colMeans(subset(data_matrixAll,data_matrixAll[,1] > (ref_wave - 5) & data_matrixAll[,1] < (ref_wave + 5)))
# Remove effect of baseline drift by subtracting baseline absorbance from 365 nm absorbance
BrCcorr <- BrC365-BrCref #closer to actual signal we want

# time calcs -----------
##########################
## We create a POSIXct date/time with the experiment date and entered desired start time. The difference 
## between actual start time and desired start time is subtracted from every date/time.
##########################

# Ask user to enter an experiment reference start time in HH:MM:SS format. 
# Ex: to start at 2 pm , enter 14:00:00. 
#time_plot <- readline(prompt="Reference time? For 00:00:00, press enter. For other time, enter as HH:MM:SS. ")

#if (time_plot == "\n") {
#  time_plot <- "00:00:00"
#}

time_plot <- "00:00:00"

Time <- TimeSeries[2:numFiles+1]    # rename time series vector
getDate <- as.Date(head(Time, n=1))    # get date of first datetime stamp
date_time <- paste(getDate, time_plot)    # create a string with the experiment date and user entered start time
refTime <- as.POSIXct(date_time)    # convert to POSIXct date / time 
difference <- as.numeric(difftime(head(Time, n=1), refTime), units="secs")    # difference between actual and desired start; seconds can be subtracted directly from POSIX class
CorrectedTime_Ref <- Time - difference    # new referenced times; original time vector minus the difference between actual and desired start      

# gets the date of the files being analyzed from the reference date above
date_str <- as.character(getDate)   # make the date a string
split <- strsplit(date_str, "-")   # split the string at the -
expt_date <- paste("16", split[[1]][2], split[[1]][3], sep = "")   # puts the date in YYMMDD format
split_path <- strsplit(getwd(), split = .Platform$file.sep)   # gets the path - split into pieces - of the computer that this code is being run on
path_prelim <- paste(.Platform$file.sep, file.path(split_path[[1]][2], split_path[[1]][3], "Dropbox (Hawkins Research Lab)", "Hawkins Research Lab Team Folder", "Paris CESAM study 2016", "Preliminary Graphs"), sep = "")   # concatenate the computer's user and the Dropbox folders to save to

# SMPS calcs -----------
##########################
# Ask for SMPS data. If it exists, read it in and format it correctly in a data
# frame, then interpolate so that we can get data points at the times we sampled at.
# Calculate MAC at 365 nm.
##########################

path_prelim_SMPS <- paste(.Platform$file.sep, file.path(split_path[[1]][2], split_path[[1]][3], "Dropbox (Hawkins Research Lab)", "Hawkins Research Lab Team Folder", "Paris CESAM study 2016", "SMPS"), sep = "")   # concatenate the computer's user and the Dropbox folders to save to
SMPS_check <- "y"

#SMPS_check <- readline(prompt="Do you have SMPS file corresponding to Daily Spectra? y or n: ")

#if (SMPS_check == "y"){
  
  # Ask what type of file it is (raw or processed)
#  SMPS_analysis_type <- readline(prompt = "Is the SMPS data file raw data (r) or processed data (p)? ")

#  }

#if (SMPS_check == "y" && SMPS_analysis_type == "p") {
  
  # read in SMPS processed data by allowing user to select the file
#  SMPSfile <- tk_choose.files(default="", caption="Select a tab-delimited SMPS file with mm/dd/yyyy format")
  
  # format SMPS data
#  SMPS <- read.table(SMPSfile, sep="\t", header=TRUE,stringsAsFactors = FALSE)
#  SMPS$smpstimeFormatted <- as.POSIXct(SMPS$smpstime, format="%m/%d/%Y %H:%M")
  
  # use the approx function to interpolate
#  InterSMPS <- approx(SMPS$smpstimeFormatted, SMPS$smpsconc, TimeSeries, method = "linear", rule = 1, f = 0, ties = mean)
  
#}

#&& SMPS_analysis_type == "r"
if (SMPS_check == "y") {
  
  # read in SMPS raw data by allowing user to select the file
  SMPS_testFile <- file.path(path_prelim_SMPS, paste("CESAM_", expt_date, "_SMPS.txt", sep=""))
  #SMPS_testFile <- tk_choose.files(default="",caption="Select a raw SMPS file")
  
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
  MAC <- (BrCcorr*2014286)/InterSMPS$y  # obscure number comes from unit and dilution correction (page 39) in HGW lab notebook
  
}

# create new theme -----------

require(ggplot2)
require(grid)
library(ggplot2)
library(grid)

theme_best <- function(base_size = 12, base_family = "Helvetica") {
  theme(
  line = element_line(colour = "black", size = 0.5, linetype = 1, lineend = "butt"), 
  rect = element_rect(fill = "white", colour = "black", size = 0.5, linetype = 1), 
  #text = element_text(family = base_family, face = "plain", colour = "black", size = base_size, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9), 
  
  axis.line = element_blank(),   # line parallel to axis; don't need this
  axis.text = element_text(size = rel(0.8), colour = "black"),   # tick labels
  axis.text.x = element_text(vjust = 0.5, margin=margin(9,5,11,5,"pt")),   # x-axis tick labels
  axis.text.y = element_text(hjust = 1, margin=margin(5,8,10,11,"pt")),   # y-axis tick labels
  axis.ticks = element_line(colour = "black", size = 0.2),   # axis tick marks
  axis.title = element_text(colour = "black"),   # axis titles
  axis.title.x = element_text(vjust = 0.1),   # x-axis title
  axis.title.y = element_text(angle = 90),   # y-axis title
  axis.ticks.length = unit(-0.3, "lines"),   # length of tick marks

  
  legend.background = element_rect(colour = NA), 
  legend.margin = unit(0.2, "cm"), 
  legend.key = element_rect(fill = "black", colour = "white"), 
  legend.key.size = unit(1.2, "lines"), 
  legend.key.height = NULL, 
  legend.key.width = NULL, 
  legend.text = element_text(size = rel(0.8), colour = "white"), 
  legend.text.align = NULL, 
  legend.title = element_text(size = rel(0.8), face = "bold", hjust = 0, colour = "white"), 
  legend.title.align = NULL, 
  legend.position = "right", 
  legend.direction = "vertical", 
  legend.justification = "center", 
  legend.box = NULL, 
  
  panel.background = element_rect(fill = "white", colour = NA),
  panel.border = element_rect(fill = NA, colour = "black", size = 0.7),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.margin = unit(0.25, "lines"), 
  
  strip.background = element_blank(), 
  strip.text.x = element_blank(), 
  strip.text.y = element_blank(), 
  strip.text = element_blank(),
  
  plot.background = element_rect(colour = "white", fill = "white"), 
  plot.title = element_text(size = rel(1.3), margin=margin(5,8,13,11,"pt")), 
  plot.margin = unit(c(1,1,1,1), "lines"),
  
  complete = TRUE
  )
} # Check that it is a complete theme attr(theme_black(), "complete")

# rainbow -----------
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
  
  # plot
  #grid.newpage()    # create new page for plot
  layout(t(1:2), widths=c(10,2))    # set up the layout of the plot
  my.palette <- rainbow(length(justTime), start=0, end=4/6)      # create rainbow colors with length of time vector, from red (start=0) to blue (end=4/6)
  
  # create a PDF file to save plot to
  if (save_graphs == "yes" ) {
    pdf(file.path(path_prelim, paste(expt_date, "_CESAM", sep = ""), "rainbow plot.pdf"))
  }
  
  matplot(spectra_corr[,1], spectra_corr[,-1], type="l", xlim=c(300,700), ylim=c(rainbow_ymin, rainbow_ymax), 
          xlab="Wavelength (nm)", ylab="Corrected Absorbance",
          main = paste("Corrected Absorbance on ", split[[1]][2], "/", split[[1]][3], "/", split[[1]][1], sep = ""),
          col = my.palette) 
  image.plot(smallplot= c(.99,1,0.1,.9), zlim=c(Time_asHours[2],Time_asHours[length(Time_asHours)]), 
             legend.only=TRUE, horizontal = FALSE, col=my.palette, legend.lab="Local Time")   # add color bar to plot 
  mtext(getDate, side=3)    # add date to the plot
  
  # turn off PDF save
  if (save_graphs == "yes" ) {
    dev.off()
  }
}

# corrected rainbow -----------
##########################
# Rainbow plot with noisy spectra removed
# Currently removing spectra with a standard deviation over sd_limit (assigned
# above) over the region 390 nm to 410 nm
# Also makes a plot to show which spectra were removed
# Plots only if user has SMPS file
##########################

if (corrected_rainbow_plot == "yes") {
  
  num_cols <- ncol(spectra_corr) - 1    # number of columns to loop through
  times <- vector(length=num_cols)    # empty vector for times
  removed <- vector(length=num_cols)     # assigns 0 or 1 to each time depending on whether or not it is removed
  matrix_corr <- spectra_corr    # copy all of the spectra into a new matrix; wavelengths in first column 
  index_log <- vector(length=0)   # empty matrix that will be used to store the indices of bad spectra
  
  for (i in 2:num_cols) {
    
    times[i-1] <- colnames(spectra_corr)[i]    # fill in matrix with times
    
    if (abs(mean(spectra_corr[920:1020,i])) > abs(5 * mean(spectra_corr[920:1020,2:num_cols]))) {   # choose spectra with absorbances (averaged over 390 and 410 nm) with magnitudes 3+ times greater than the average of all spectra over that range; this is because bubbles have much greater magnitude "signal"
      index_log <- append(index_log, i)   # add the index of the bad spectrum to our log of indices
      removed[i-1] <- 1    # assign 1 to the removed time log
    } else {
      removed[i-1] <- 0    # assign 0 to all other times
    }
  }
  
  index_log_rev <- rev(index_log)   # reverse the order of the indices in the vector
  # should add explanation of why these need to be reversed
  
  for (i in index_log_rev) {
    matrix_corr <- matrix_corr[,-i]   # remove the bad spectra
  }
  
  # plot
  grid.newpage()
  par(mar= c(5, 2, 4, 2) + 0.1)
  plot(removed, xlab="Time", yaxt="n")    # no y-axis
  axis(2, at = seq(0, 1, by = 1), las=2)    # add y-axis
  
  #### Plots new matrix with same parameters as above 
  
  # layout to fit both plot & color bar legend
  layout(t(1:2), widths=c(10,2))
  
  # create rainbow colors with length of time vector, from red (start=0) to blue (end=4/6)
  my.palette <- rainbow(length(justTime), start=0, end=4/6)
  
  # create a PDF file to save plot to
  if (save_graphs == "yes" ) {
    pdf(file.path(path_prelim, paste(expt_date, "_CESAM", sep = ""), "rainbow plot corrected.pdf"))
  }
  
  # plot the matrix
  matplot(matrix_corr[,1], matrix_corr[,-1], type="l", xlim=c(300,700), ylim=c(rainbow_ymin, rainbow_ymax), 
          main = paste("Corrected Absorbance on ", split[[1]][2], "/", split[[1]][3], "/", split[[1]][1], sep = ""),
          xlab="Wavelength", ylab="Absorbance", col = my.palette) 
  
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

# PTR -----------
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
      pdf(file.path(path_prelim, paste(expt_date, "_CESAM", sep = ""), "PTR.pdf"))
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

# ref time plots -----------
##########################
## plots for BrCcorr and MAC at reference time 
##########################

# create a PDF file to save plot to
if (save_graphs == "yes" ) {
  pdf(file.path(path_prelim, paste(expt_date, "_CESAM", sep = ""), "corr abs vs ref time.pdf"))
}

# corrected absorbance vs time series
qplot3 <- qplot(CorrectedTime_Ref, BrCcorr[2:numFiles+1], colour ="red", geom="line",
                xlab=paste("Time since ", format(Time[1], format =  "%H:%M:%S")), ylab="Corrected Absorbance at 365 nm",
                main = paste("Corrected Absorbance at 365 nm on ", split[[1]][2], "/", split[[1]][3], "/", split[[1]][1], sep = ""),
                show.legend=FALSE) + theme_best() + 
                scale_x_datetime(limits = c(CorrectedTime_Ref[1], CorrectedTime_Ref[numFiles]), expand = c(0, 0)) +
                scale_y_continuous(limits = c(corr_abs_ymin, corr_abs_ymax), expand = c(0.001,0.001))
# scale_x_datetime and scale_y_continuous force the axes to intersect at the limits (aka it removes the offset from 0)
print(qplot3) 

# turn off PDF save
if (save_graphs == "yes" ) {
  dev.off()
}

# create a PDF file to save plot to
if (save_graphs == "yes" ) {
  pdf(file.path(path_prelim, paste(expt_date, "_CESAM", sep = ""), "MAC vs ref time.pdf"))
}

# MAC vs time series, only if you have access to SMPS file
if (SMPS_check == "y") {
  qplot4 <- qplot(CorrectedTime_Ref, MAC[2:numFiles+1], colour ="red", geom="line",
                  xlab=paste("Time since", format(Time[1], format =  "%H:%M:%S")), ylab=expression(paste(MAC[365], " (c", m^{2}, "/g OM)")), 
                  main = paste("MAC at 365 nm on ", split[[1]][2], "/", split[[1]][3], "/", split[[1]][1], sep = ""),
                  show.legend=FALSE) + theme_best() + 
                  scale_x_datetime(limits = c(CorrectedTime_Ref[1], CorrectedTime_Ref[numFiles]), expand = c(0, 0)) +
                  scale_y_continuous(limits = c(MAC_ymin, MAC_ymax), expand = c(0.001,0.001))
  # scale_x_datetime and scale_y_continuous force the axes to intersect at the limits (aka it removes the offset from 0)
  print(qplot4) 
}

# turn off PDF save
if (save_graphs == "yes" ) {
  dev.off()
}

# local time plots -----------
##########################
## plot BrCcorr and MAC at local Paris time
##########################

# create a PDF file to save plot to
if (save_graphs == "yes" ) {
  pdf(file.path(path_prelim, paste(expt_date, "_CESAM", sep = ""), "corr abs vs real time.pdf"))
}

# create a grob (grid graphical object) with the current date to add as an annotation to plots 
#date_grob = grobTree(textGrob(getDate, gp=gpar(col="black", fontsize=15, fontface="italic")))
#x=0.1,  y=0.95, hjust=0

# corrected absorbance vs time series
qplot1 <- qplot(Time, BrCcorr[2:numFiles+1], colour="red", geom = "line",
                xlab="Local Time (Paris)", ylab="Corrected Absorbance",
                main = paste("Corrected Absorbance at 365 nm on ", split[[1]][2], "/", split[[1]][3], "/", split[[1]][1], sep = ""),
                show.legend=FALSE) + theme_best() + 
                scale_x_datetime(limits = c(Time[1], Time[numFiles]), expand = c(0, 0)) +
                scale_y_continuous(limits = c(corr_abs_ymin, corr_abs_ymax), expand = c(0.001,0.001))
# scale_x_datetime and scale_y_continuous force the axes to intersect at the limits (aka it removes the offset from 0)
print(qplot1) 

# turn off PDF save
if (save_graphs == "yes" ) {
  dev.off()
}

# create a PDF file to save plot to
if (save_graphs == "yes" ) {
  pdf(file.path(path_prelim, paste(expt_date, "_CESAM", sep = ""), "MAC vs real time.pdf"))
}

# MAC vs time series, only if you have access to SMPS file
if (SMPS_check == "y") {
  qplot2 <- qplot(Time, MAC[2:numFiles+1], colour="green", geom = "line",
                  xlab="Local Time (Paris)", ylab=expression(paste(MAC[365], " (c", m^{2}, "/g OM)")),
                  main = paste("MAC at 365 nm on ", split[[1]][2], "/", split[[1]][3], "/", split[[1]][1], sep = ""),
                  show.legend=FALSE) + theme_best() + 
                  scale_x_datetime(limits = c(Time[1], Time[numFiles]), expand = c(0, 0)) + 
                  scale_y_continuous(limits = c(MAC_ymin, MAC_ymax), expand = c(0.001,0.001))
  # scale_x_datetime and scale_y_continuous force the axes to intersect at the limits (aka it removes the offset from 0)
  print(qplot2)
}



# turn off PDF save
if (save_graphs == "yes" ) {
  dev.off()
}

# log log -----------
##########################
# log-log plot of absorptivity vs wavelength
# Plots only if user has SMPS file
##########################

if (log_log_plot == "yes"){

  matrix_ll <- matrix_corr   # copy the spectra (with noisy ones removed) into a new matrix
  row_Means_ll <- rowMeans(matrix_ll[,2:ncol(matrix_ll)], na.rm = TRUE)    # take the mean of all the absorbances at each wavelength 
    
  # plot
  layout(t(1:2), widths=c(10,2))    # set up the layout of the graph
  
  # create a PDF file to save plot to
  if (save_graphs == "yes" ) {
    pdf(file.path(path_prelim, paste(expt_date, "_CESAM", sep = ""), "log log.pdf"))
  }
 
  # plot on a log-log axis
  matplot(matrix_ll[,1], matrix_ll[,-1], type="l", xlim = c(300, 500), ylim = c(log_ymin, log_ymax),
          xlab="Wavelength (nm)", ylab="Absorbance", col = my.palette, log = "xy")
  title(paste("Absorbance vs. Wavelength on ", split[[1]][2], "/", split[[1]][3], "/", split[[1]][1], sep = ""))
  
  # add line of mean absorbance vs wavelength
  lines(matrix_ll[,1], row_Means_ll)
  
  # add color bar to plot 
  image.plot(smallplot= c(.99,1,0.1,.9), zlim=c(Time_asHours[2],Time_asHours[length(Time_asHours)]), 
             legend.only=TRUE, horizontal = FALSE, col=my.palette, legend.lab="Local Time")
  
  # create a PDF file to save plot to
  if (save_graphs == "yes" ) {
    dev.off()
  }
  
} # end of check for log/log plot







