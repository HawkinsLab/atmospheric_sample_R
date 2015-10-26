####################################
## Data import
## Result: data_frm (frame object)
####################################
# read the datafile line by line
# txt <- readLines("CESAM_150615_00000.txt") 
txt <- readLines(file.choose()) # brings up a file browser
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
# name the columns
colnames(data_mat) <- c("wavelength","absorbance")
# convert a matrix (with colnames) into a frame
data_frm <- as.data.frame(data_mat, stringsAsFactors=FALSE)
# type casting to numeric values
data_frm$wavelength <- as.numeric(data_frm$wavelength)
data_frm$absorbance <- as.numeric(data_frm$absorbance)
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
## Result: data_frm (frame object)
########################################
data_frm$timestamp <- rep(date_time, nrow(data_frm)) # repeat date_time value
# As a result, data_frm has three coloumns: wavelength, absorbance, Timestamp

########################################
## Plot
########################################
# A simple absorbance plot (for now)
plot(data_frm$wavelength,data_frm$absorbance,xlab="Wavelength",ylab="Absorbance",type="l")

########################################
## Clear intermediate variables
########################################
rm(data_mat,data_str,date_time,datetime_line,
   datetime_obj,datetime_str,fieldList,ind,txt)
