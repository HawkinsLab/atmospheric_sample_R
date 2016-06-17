
data_to_dframe <- function(file, wl_low=170, wl_high=900) {
  # Process raw data file in tab delimited format 
  # Args:
  #   file: data file name 
  #   wl_low: lower bound of wavelength
  #   wl_high: upper bound of wavelength
  #
  # Returns: 
  #   A list containing three items, wavelength, absorbance, and timestamp
  # 
  #Lelia test commit 6/1/2016
  #######################
  ## Opening the individual text files, getting wavelength and absorbance but not date
  ###########################
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
  tempList <- data_to_dframe(f=file) 

  data_matrixAll[,i+1] <- tempList$absorbance
  data_matrixAll[,1] <- tempList$wavelength
  
  TimeSeries[i+1] <- tempList$timeStamp
 # data_list_uv_mean[[i]] <- lapply(data_list_uv[[i]],mean)
  i <- i + 1
}
class(TimeSeries) <- c('POSIXt','POSIXct') #reestablishes as time series

#Now we will take the average of the absorption around the brown carbon region, keeping the time series intact.
BrC365 <- colMeans(subset(data_matrixAll,data_matrixAll[,1] > 360 & data_matrixAll[,1] < 370))


# Now we take the average absorption in some long wavelength reference region, usually 700 nm but recently at CESAM
# We found that 600 or 550 is better due to fluctuating signal 700

#BrCref <- colMeans(subset(data_matrixAll,data_matrixAll[,1] > 695 & data_matrixAll[,1] < 705))
BrCref <- colMeans(subset(data_matrixAll,data_matrixAll[,1] > 545 & data_matrixAll[,1] < 555))
 

#We must now subtract the absorbance at the reference wavelength from the BrC wavelength
BrCcorr <- BrC365-BrCref #closer to actual signal we want

#############################
## Now we begin using SMPS data to normalize absorbance. First we ask user if they have SMPS data. 
## If 'no' just plot MAC, if 'yes' plot MAC and calculate/plot BrC correction. 
## We read the file in and format time.
#########################

SMPS_check <- readline(prompt="Do you have SMPS file corresponding to Daily Spectra? Enter yes/no:")

if (SMPS_check == "yes"){
  
#read in particle concentration data, by allowing user to select the file
SMPSfile <- tk_choose.files(default="",caption="Select a tab-delimited SMPS file with mm/dd/yyyy format")
SMPS <- read.table(SMPSfile, sep="\t", header=TRUE,stringsAsFactors = FALSE)

SMPS$smpstimeFormatted <- as.POSIXct(SMPS$smpstime, format="%m/%d/%Y %H:%M")

#use the approx function to interpolate
InterSMPS <- approx(SMPS$smpstimeFormatted, SMPS$smpsconc, TimeSeries, method = "linear", rule = 1, f = 0, ties = mean)

####################################
### This commented out section can be used to read in raw SMPS csv file

# Allow user to select csv file with raw SMPS data
#SMPS_testFile <- tk_choose.files(default="",caption="Select a raw SMPS file")

## NOTE: This assumes length of file is 136, thus including all the columns because otherwise it would think there were only 2 columns
#SMPS <- read.csv(SMPS_testFile, head=FALSE, sep=",", col.names = paste0(seq_len(136)), fill=TRUE)

# check if line 15 column 2 is "dw/dlogDp"
#if (SMPS$X2[15]=="dw/dlogDp"){
# num <- grepl("^[0-9]",SMPS$X1)    # if this is true, find all the lines starting with a number
  #SMPS.df <- SMPS.df[num,]        # update date frame to just include lines starting with a number
#  SMPS <- SMPS[num,]  
#} 

## create a date frame with date, time, and total smps concentration
#SMPS_conc <- SMPS$X136    # total concentration
#SMPS_datetime<- as.POSIXct(paste(SMPS$X2, SMPS$X3), format="%m/%d/%y %H:%M")
#SMPS.df <- data.frame(SMPS_datetime, SMPS_conc)

#InterSMPS <- approx(SMPS_datetime, SMPS_conc, TimeSeries, method = "linear", rule = 1, f = 0, ties = mean)

##################################
# NOTE
# To workaround the "Error in plot.new() : figure margins too large
# par(mar=c(1,1,1,1))

#####################
#Finally, we normalize the signal of brown carbon coming from the difference of Abs at 365 nm and 
#some reference wavelength (typically 700 nm but in Paris 2015 we found that 550 was more suitable
#due to fluctuations at 700 for unknown reasons)
########################

MAC <- BrCcorr/InterSMPS$y
} # end of SMPS calculations 

############
## We create a POSIXct date/time with the experiment date and entered desired start time. The difference 
## between actual start time and desired start time is subtracted from every date/time.
############

# Ask user to enter an experiment start time in HH:MM:SS format. Ex: to start at 9 am , enter 09:00:00. 
time_plot <- readline(prompt="What is the reference time for start of experiment? (enter as HH:MM:SS):")

# rename time series vector
Time <- TimeSeries[2:numFiles+1]   

# get date of first date/time stamp
getDate <- as.Date(head(Time, n=1))  
  
# create a string with the experiment date and user entered start time
date_time <- paste(getDate, time_plot)   

# convert to POSIXct date / time
refTime <- as.POSIXct(date_time)         

## In order to adjust the time to start at reference time, we find the difference between the actual start and reference start time.
## That difference is subtracted from the entire time series.

# difference between actual and desired start
difference <- as.numeric(difftime(head(Time, n=1), refTime), units="secs") # seconds can be subtracted directly from POSIX class

# time vector minus difference between actual and desired start 
CorrectedTime_Ref <- Time-difference           


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
      xlab="Local Time (Paris)", ylab="Brown Carbon correction", show.legend=FALSE) + annotation_custom(date_grob) + theme_bw()
print(qplot1) 

# time series vs MAC, only if you have access to SMPS file 
if (SMPS_check == "yes"){
  qplot2 <- qplot(Time, MAC[2:numFiles+1], colour="green", geom = "line",
        xlab="Local Time (Paris)", ylab="MAC", show.legend=FALSE) + annotation_custom(date_grob) + theme_bw()
  print(qplot2)
} 

#######################
## plots for BrCcorr and MAC at reference time 
########################

qplot3 <- qplot(CorrectedTime_Ref, BrCcorr[2:numFiles+1], colour ="red", geom="line",
      xlab=paste("Time since", time_plot), ylab="Brown Carbon correction", show.legend=FALSE)+ annotation_custom(date_grob) + theme_bw()
print(qplot3) 
  
if (SMPS_check == "yes"){
  qplot4 <- qplot(CorrectedTime_Ref, MAC[2:numFiles+1], colour ="red", geom="line",
        xlab=paste("Time since", time_plot), ylab="MAC", show.legend=FALSE)+ annotation_custom(date_grob) + theme_bw()
  print(qplot4) 
  }


#################
# Original plots (without ggplot)  
#################

#plot(TimeSeries[2:numFiles+1], BrCcorr[2:numFiles+1], col="red", type = "l")
#plot(TimeSeries[2:numFiles+1], MAC[2:numFiles+1], col="green", type = "l", ylim = c(0,0.005))

#################
## Plot of MAC over time with ion counts included, only if user has SMPS file
## File is PTR Corr Series file 
#################

if (SMPS_check == "yes"){
  
require("reshape2")
library("reshape2")
  
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

} # end of SMPS check for PTR plots

#################
## Absorbance versus wavelength for each time 
## 'Rainbow plot'
## We create a matrix with the same dimensions as data_matrixAll 
## For each time, we subtract the BrCref value from all absorbace measurements 
#################

# create new time vector of just times, but keeping the first time stamp
# first time stamp is place holder 
justTime <- strftime(TimeSeries, format = "%H:%M:%S")

# create new matrix with same dimensions as data matrix
matrix_noBrC <- matrix(NA, nrow=numRows, ncol=numFiles+1)   

# rename columns using time series of just times 
colnames(matrix_noBrC) <- justTime

# Fill matrix_noBrC with wavelength vector and with BrCref subtracted from every absorbance in data_matrix
i<-1
for (i in 2:(length(BrCref))){
  # Add wavelength vector from data_matrix to first column
  matrix_noBrC[,1] <- data_matrixAll[,1]           
  # Subtract BrCref from all absorbance measurements from data_matrix
  matrix_noBrC[,i] <- data_matrixAll[,i]-BrCref[i] 
  i<-i+1
}

# function to convert HH:MM:SS to hours in order to plot as legend
Time_asHours <- sapply(strsplit(justTime,":"),
                       function(x) {
                         x <- as.numeric(x)
                         ((x[1]+x[2]/60) + (x[3]/3600)) })

# layout to fit both plot & color bar legend
layout(t(1:2), widths=c(8,1))
# set new plot
grid.newpage()

# create rainbow colors with length of time vector, from red (start=0) to blue (end=4/6)
my.palette <- rainbow(length(justTime), start=0, end=4/6)

# plot the matrix
  # matrix_noBrC[,1] is the wavelength vector
  # matrix_noBrC[,-1] is all of the absorbance vectors
matplot(matrix_noBrC[,1], matrix_noBrC[,-1], type="l", xlim=c(300,700), ylim=c(-0.15,.45), 
        xlab="wavelength", ylab="absorbance", col = my.palette) 

# add color bar to plot 
image.plot(smallplot= c(.99,1,0.1,.9), zlim=c(Time_asHours[2],Time_asHours[length(Time_asHours)]), 
           legend.only=TRUE, horizontal = FALSE, col=my.palette, legend.lab="Local Time")
# add date to the plot
mtext(getDate, side=3, line=2) 

##########
# In order to eliminate some of the noise, we remove any measurement with a 
# standard deviation above 0.05 between wavelengths of approximately 500 and 550 
# Stan dev is measured from approx. 480 at matrix_noBrC[1450,] to 530 at matrix_noBrC[1700,]
# We assign each time stamp a 0/1. A 1 is given to times vectors that are removed.
#########

# Number of columns to loop through
num_cols <- ncol(matrix_noBrC)

# Create two vectors to fill with length of number of absorbance measurements. 
# First vector is times
# Second vector assigns 0 or 1 to each time depending on whether or not it is removed
times <- vector(length=num_cols)
removed <- vector(length=num_cols)

# Fill new matrix with values of matrix_noBrC 
matrix_corr <- matrix_noBrC

i <- 1
for (i in 2:num_cols){
  # vector with every time stamp
  times[i] <- colnames(matrix_noBrC)[i]  
  
  if (sd(matrix_noBrC[1450:1700,i]) > 0.05){
    # remove column with value i from matrix
    matrix_corr <- matrix_corr[,-i]    
    # assign 1 to removed time
    removed[i] <- 1     
  }
  else {
    # assign 0 to all other times
    removed[i] <- 0
  }
  i < i+1
}

# Create a data frame with time vector and 0/1 vector
binTimes <- data.frame(times, removed)
# Plot to determine if there is a pattern to which times are removed
plot(binTimes$times, binTimes$removed)

# add wavelength vector to corrected matrix
matrix_corr <- cbind(matrix_noBrC[,1], matrix_corr)


#### Plots new matrix with same parameters as above 
# layout to fit both plot & color bar legend
layout(t(1:2), widths=c(8,1))
# set new plot
grid.newpage()

# create rainbow colors with length of time vector, from red (start=0) to blue (end=4/6)
my.palette <- rainbow(length(justTime), start=0, end=4/6)

# plot the matrix
matplot(matrix_corr[,1], matrix_corr[,-1], type="l", xlim=c(300,700), ylim=c(-0.15,.45), 
        xlab="wavelength", ylab="absorbance", col = my.palette) 

# add color bar to plot 
image.plot(smallplot= c(.99,1,0.1,.9), zlim=c(Time_asHours[2],Time_asHours[length(Time_asHours)]), 
           legend.only=TRUE, horizontal = FALSE, col=my.palette, legend.lab="Local Time")
# add date to the plot
mtext(getDate, side=3, line=2) 






