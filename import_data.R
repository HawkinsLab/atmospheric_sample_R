
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
  #browser()
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
 
  # As a result, returnlist has three items: wavelength, absorbance, Timestamp

  
  returnlist <- list("wavelength"=as.numeric(data_mat[,1]),"absorbance"=as.numeric(data_mat[,2]), "timeStamp"= date_time)
  return(returnlist) 

}

## ask user if they have SMPS data, if 'no' plot MAC, if 'yes' plot both BrCcorr & MAC
SMPS_check <- readline(prompt="Do you have SMPS file corresponding to Daily Spectra? Enter yes/no:")

#ExptDay <- readline(prompt="Enter the expt day as CESAM_YYMMDD: ")  # removed pattern requirement from list files function

require(tcltk)
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
#Now we take the average absorption in some long wavelength reference region, usually 700 nm but recently at CESAM
#we found that 600 or 550 is better due to fluctuating signal 700
BrCref <- colMeans(subset(data_matrixAll,data_matrixAll[,1] > 695 & data_matrixAll[,1] < 705))

#We must now subtract the absorbance at the reference wavelength from the BrC wavelength
BrCcorr <- BrC365-BrCref #closer to actual signal we want

#############################
## Now we begin using SMPS data to normalize absorbance. We read the file in, format time
#########################

if (SMPS_check == "yes"){
#read in particle concentration data, by allowing user to select the file
SMPSfile <- tk_choose.files(default="",caption="Select a tab-delimited SMPS file with mm/dd/yyyy format")
SMPS <- read.table(SMPSfile, sep="\t", header=TRUE,stringsAsFactors = FALSE)

SMPS$smpstimeFormatted <- as.POSIXct(SMPS$smpstime, format="%m/%d/%Y %H:%M")

#use the approx function to interpolate
InterSMPS <- approx(SMPS$smpstimeFormatted, SMPS$smpsconc, TimeSeries, method = "linear", rule = 1, f = 0, ties = mean)

####################################
## this is for reading in raw SMPS csv file
## NOTE: this assumes length of file is 136, thus including all the columns because otherwise it would think there were only 2 columns
#SMPS <- read.csv(file="June_24_SMPS.csv",head=FALSE,sep=",", col.names = paste0(seq_len(136)), fill=TRUE)

## create a date frame with date, time, and total smps concentration
#SMPS_conc <- SMPS$X136    # total concentration
#SMPS_datetime<- as.POSIXct(paste(SMPS$X2, SMPS$X3), format="%m/%d/%y %H:%M:%S")
#SMPS.df <- data.frame(SMPS_datetime, SMPS_conc)

# check if line 15 column 2 is "dw/dlogDp"
#if (SMPS$X2[15]=="dw/dlogDp"){
#  num <- grepl("^[0-9]",SMPS$X1)    # if this is true, find all the lines starting with a number
#  SMPS.df <- SMPS.df[num,]        # update date frame to just include lines starting with a number
#} 
#InterSMPS <- approx(SMPS.df$SMPS_datetime, SMPS.df$SMPS_conc, TimeSeries, method = "linear", rule = 1, f = 0, ties = mean)
##################################


# To workaround the "Error in plot.new() : figure margins too large
# par(mar=c(1,1,1,1))

#####################
#Finally, we normalize the signal of brown carbon coming from the difference of Abs at 365 nm and 
#some reference wavelength (typically 700 nm but in Paris 2015 we found that 550 was more suitable
#due to fluctuations at 700 for unknown reasons)
########################

MAC <- BrCcorr/InterSMPS$y
}    

## ask user to enter an experiment start time with 0 for midnight 
time_plot <- readline(prompt="What is the reference time for start of experiment? (enter as HH:MM:SS):")

Time <- TimeSeries[2:numFiles+1]  # rename time series vector 

#########
## This function converts hours into seconds, because seconds can be subtracted directly from POSIXct date/time format
#########

sec <- function(u) {
  x <- u * 3600
  return(x)
}

############
## We create a POSIXct date/time with the experiment date and entered desired start time. The difference 
## between actual start time and desired start time is subtracted from every date/time.
############

  # get date of first date/time stamp
getDate <- as.Date(head(TimeSeries[2:numFiles+1], n=1))        
  # create a string with the experiment date and user entered start time
date_time <- paste(getDate, time_plot)   
  #  POSIXct date/time
refTime <- as.POSIXct(date_time)                                
  # difference between actual and desired start
start_difference <- as.numeric(difftime(head(TimeSeries[2:numFiles+1], n=1), refTime), units="hours")  
  # time vector minus difference between actual and desired start 
CorrectedTime_Ref <- Time-sec(start_difference)           


#######################
## plot for BrCcorr and MAC at local Paris time
########################
require(ggplot2)
library(ggplot2)
library(grid)
require(grid)

my_grob = grobTree(textGrob(getDate, x=0.1,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=15, fontface="italic")))

   # time series vs BrC correction 
qplot1 <- qplot(TimeSeries[2:numFiles+1], BrCcorr[2:numFiles+1], colour="red", geom="line", 
      xlab="Local Time (Paris)", ylab="Brown Carbon correction", show.legend=FALSE) + annotation_custom(my_grob) + theme_bw()
print(qplot1) 

    # time series vs MAC, only if you have access to SMPS file 
if (SMPS_check == "yes"){
  qplot2 <- qplot(TimeSeries[2:numFiles+1], MAC[2:numFiles+1], colour="green", 
        xlab="Local Time (Paris)", ylab="MAC", show.legend=FALSE, geom = "line") + annotation_custom(my_grob) + theme_bw()
  print(qplot2)
} 


#######################
## plots for BrCcorr and MAC at reference time 
########################

  qplot3 <- qplot(CorrectedTime_Ref, BrCcorr[2:numFiles+1], colour ="red", geom="line",
        xlab=paste("Time since", time_plot), ylab="Brown Carbon correction", show.legend=FALSE)+ annotation_custom(my_grob) + theme_bw()
  print(qplot3) 
  
  if (SMPS_check == "yes"){
  qplot4 <- qplot(CorrectedTime_Ref, MAC[2:numFiles+1], colour ="red", geom="line",
        xlab=paste("Time since", time_plot), ylab="MAC", show.legend=FALSE)+ annotation_custom(my_grob) + theme_bw()
  print(qplot4) 
  }


#################
# Original plots (without ggplot)  
#################

#plot(TimeSeries[2:numFiles+1],BrCcorr[2:numFiles+1], col="red", type = "l")
#plot(TimeSeries[2:numFiles+1], MAC[2:numFiles+1], col="green", type = "l", ylim = c(0,0.005))

#################
## Plot of MAC over time with ion counts included
## File is PTR Corr Series file 
#################
#require("reshape2")
library("reshape2")
  
  
  ## read in PTR Corr Series file, must be tab-demilited
  PTR_Series <- tk_choose.files(default="",caption="Select a tab-delimited PTR Corr Series file")
  PTR <- read.table(PTR_Series, sep="\t", header=TRUE, stringsAsFactors = FALSE)
  
  # format time to POSIXct
  PTR_timeFormatted <- as.POSIXct(PTR$Time, format="%m/%d/%y %H:%M")
  #PTR_timeFormatted <- as.POSIXct(PTR$Time, format="%m/%d/%Y %H:%M")
  
  # create data frame with time, MA, MG, Imine
  PTR.df <- data.frame(PTR_timeFormatted, PTR$MA, PTR$MG, PTR$Imine)
  
  # create vector of MAC data, along with data frame of MAC data with Time Series
  MAC_data <- MAC[2:numFiles+1]
  MAC.df <- data.frame(Time, MAC_data)
  
  library(ggplot2)
  library(gtable)
  library(grid)
  
  grid.newpage()
  
  ## create ggplot for MAC versus time and for PTR data versus time 
  p1 <- ggplot(MAC.df, aes(Time, MAC_data)) + geom_line(colour = "black") + theme_bw() + ylab("MAC")
  p2 <- ggplot() +
    geom_line(data = PTR.df, aes(x = PTR_timeFormatted, y = PTR.MA, colour="red")) +   
    geom_line(data = PTR.df, aes(x = PTR_timeFormatted, y = PTR.MG, colour="blue")) + 
    geom_line(data = PTR.df, aes(x = PTR_timeFormatted, y = PTR.Imine, colour = "green")) +
    theme_bw()  +  
    theme(panel.background = element_rect(fill = NA)) +
    ylab("Ion count")
  

  # extract gtable  
  g1 <- ggplot_gtable(ggplot_build(p1))
  g2 <- ggplot_gtable(ggplot_build(p2))
  
  # overlap the panel of 2nd plot on that of 1st plot
  pp <- c(subset(g1$layout, name == "panel", se = t:r))
  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)
  
  # formatting for right axis
  ia <- which(g2$layout$name == "axis-l")
  ga <- g2$grobs[[ia]]
  ax <- ga$children[[2]]
  ax$widths <- rev(ax$widths)  
  ax$grobs <- rev(ax$grobs)     
  
  g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1) # add columns to right axis
  g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
  g <- gtable_add_grob(g, g2$grob[[7]], pp$t, length(g$widths), pp$b)
  
  # draw it
  grid.draw(g) 

  
  # set margins to fit both y-axis
par(mar = c(4,6,4,6))
  # plot MAC versus time, add in y-axis label after
with(MAC.df, plot(Time, MAC_data, type="l", col="black", xlab="Time", ylab=NA))
mtext("MAC", line=3, las=2, side=2)
  # add to an existing plot using new=T
par(new = T)
  # create second plot, with lines added in for all three ion types
with(PTR.df, plot(PTR_timeFormatted, PTR.MA, col="blue", type="l", axes=F, xlab=NA, ylab=NA))
lines(PTR_timeFormatted, PTR.df$PTR.MG, col="red")
lines(PTR_timeFormatted, PTR.df$PTR.Imine, col="green")
  # add axis on right side of plot and axis label 
axis(side = 4)
mtext("Ion count", side = 4, line = 2, las=2)
  # add in legend 
legend("bottomright",
       col = c("blue", "red", "green", "black"),
       legend=c("MA", "MG", "Imine", "MAC"),
       lty=c(1))


