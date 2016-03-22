data_to_dframe <- function(file='no_input', wl_low=170, wl_high=900) {
  # Process raw data file in CSV format 
  # Args:
  #   file: data file name 
  #   wl_low: lower bound of wavelength
  #   wl_high: upper bound of wavelength
  #
  # Returns: 
  #   A data Frame containing three coloumns, wavelength, 
  #    absorbance, and timestamp
  # 
  # Lelia changed this Feb 17 2016
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
 # data_mat <- subset(data_mat, data_mat[,1] > wl_low & data_mat[,1] < wl_high)
  # name the columns
  #colnames(data_mat) <- c("wavelength","absorbance")
  # convert a matrix (with colnames) into a frame
  #data_frm <- as.data.frame(data_mat, stringsAsFactors=FALSE)
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
  #print(date_time)
  ########################################
  ## Adding the timestamp to corresponding obs as a new column
  ## Result: data_frm (a data frame object)
  ########################################
 # data_frm$timestamp <- rep(date_time, nrow(data_frm)) # repeat date_time value
  # As a result, data_frm has three coloumns: wavelength, absorbance, Timestamp

  #data_frm <- subset(data_frm, data_frm$waclasvelength > wl_low & data_frm$wavelength < wl_high)
  returnlist <- list("wavelength"=as.numeric(data_mat[,1]),"absorbance"=as.numeric(data_mat[,2]), "timeStamp"= date_time)
  return(returnlist) #I am trying to deliver a frame with wavelength and absorbance, and also a single time stamp variable per file

}

files <- list.files(pattern="CESAM_150619*")
numFiles <- length(files)
#a_dframe <- data_to_frame()
data_matrixAll <- matrix(NA,nrow=3648,ncol=numFiles+1) # a place holder (an empty data frame obj), 
TimeSeries <- vector(length=numFiles+1)

#data_list_uv_mean <- list() #am empty list
#data_list_uv <- lapply(data_to_dframe(f=file, 360, 370), get)
i = 1 # index for data_list in the following for loop
for(file in files) 
{
  # we are trying to pull, file by file, absorbances into a frame, which are contained in the first element of newlist from function above
  tempList <- data_to_dframe(f=file) 
 # print(tempList)
  data_matrixAll[,i+1] <- tempList$absorbance
  data_matrixAll[,1] <- tempList$wavelength
  #print(tempList$timestamp)
  TimeSeries[i+1] <- tempList$timeStamp
 # data_list_uv_mean[[i]] <- lapply(data_list_uv[[i]],mean)
  i = i + 1
}
class(TimeSeries) <- c('POSIXt','POSIXct')
#data_matrixAll <- as.numeric(data_matrixAll) #turns the data matrix from characters into numeric data so we can do maths
#TimeSeries[1]=0
#print(TimeSeries)
#data_list_ir <- list() # a place holder (an empty list obj)
#data_list_ir_mean <- list()
#i = 1 # index for data_list in the following for loop
# for(file in files)
# {
#   # this might not be the best way...
#   data_list_ir[[i]] <- data_to_dframe(f=file, 695, 705)
#   data_list_ir_mean[[i]] <- lapply(data_list_ir[[i]],mean)
#   i = i + 1
# }

