#I am writing a comment on 6/1/16 -Jason asdjfadsfa;slkf
# test comment by Moira 6/1/16

file = file.choose()
# read the datafile line by line
txt <- readLines(file) 
# find all the lines starting with a number
ind <- grepl("^[0-9]",txt)
# create a new array only with the lines starting with a number
#  because there are several header lines (we could skip a fixed number of lines)
data_str <- txt[ind]
# split each line into pieces
fieldList <- strsplit(data_str, split = "\t")
# create a matrix 
data_mat <- matrix(
  unlist(fieldList), 
  nrow=length(fieldList),
  byrow=TRUE)
# name the columns
colnames(data_mat) <- c("wavelength","absorbance")
## 
##NOW don't covert the matrix to a data frame
## convert a matrix (with colnames) into a frame
# data_frm <- as.data.frame(data_mat, stringsAsFactors=FALSE)

# type casting to numeric values
#data_frm$wavelength <- as.numeric(data_frm$wavelength)
#data_frm$absorbance <- as.numeric(data_frm$absorbance)
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
## Adding the timestamp to corresponding obs as a new column
## Result: data_frm (a data frame object)
########################################
data_frm$timestamp <- rep(date_time, nrow(data_frm)) # repeat date_time value
# As a result, data_frm has three coloumns: wavelength, absorbance, Timestamp

data_frm <- subset(data_frm, data_frm$wavelength > wl_low & data_frm$wavelength < wl_high)
