data_to_list <- function(file, wl_low=195, wl_high=907) {
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
  
  # Change the format of data_mat elements to numeric and subset
  # according to wl_low and wl_high
  data_mat[,1] <- as.numeric(data_mat[,1])
  data_mat[,2] <- as.numeric(data_mat[,2])
  data_mat <- subset(data_mat, data_mat[,1] > wl_low & data_mat[,1] < wl_high)
  
#  returnlist <- list("wavelength"=as.numeric(data_mat[,1]),"absorbance"=as.numeric(data_mat[,2]), "timeStamp"= date_time)
  returnlist <- list("wavelength"=data_mat[,1],"absorbance"=data_mat[,2], "timeStamp"= date_time)
  return(returnlist) 
}