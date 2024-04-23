### Take in ESF data from a single channel (as a CSV file),
### binarise the data above and below a threshold (the default is 90th centile).
### We assume that the spectrogramData has the frequency components in the first row, and
### the first column is the time data.
### ARGUMENTS
### spectrogramData :- the data from the csv file, as a numeric matrix
### threshold       :- the threshold value to binarise the power (default at 90th centile of the entire dataset)
### freqSelect      :- vector of indicies for the frequencies to keep
### freqSelect      :- a min and max value of the frequency components to keep
binariseESF <- function(spectrogramData, freqSelect = NULL, threshold = NULL){
 
  
  ### parse the frequencies and limit if required
  if(!is.null(freqSelect)){
    freqBands = round(as.numeric(gsub("X","",colnames(spectrogramData))))
    spectrogramData = spectrogramData[,
                                      which(freqBands >= freqSelect[1] & freqBands < freqSelect[2])
                                        ]
    
  }

  ## work out the threshold if required
  if(is.null(threshold)){
    threshold = quantile(spectrogramData[,seq(2,ncol(spectrogramData))], 0.9)
  }
  ### do the thresholding
  spectrogramData[which(spectrogramData < threshold)] = 0
  spectrogramData[which(spectrogramData >= threshold)] = 1
  
  return(spectrogramData)
}


# 
# ### test the code
# out = list()
# for(i in seq(6)){
#   tic("read in data and convert to matrix")
#   test = read.csv(file.path(paste0("~/Documents/Circadian/BeeSpy/Exemplar Data_30minInclVideo/ExemplarDataChannel",i-1,".csv")))  
#   test = as.matrix(test)  
#   mode(test) = "numeric"
#   toc()
#   
#   tic("binarise")
#   out[[i]] = binariseESF(test, c(0, 750))
#   toc()
# }
# 
# 
# ## plot each channel
# par(mfcol=c(3,1))
# for(i in seq(6)){
#   image(
#     x = seq(dim(out[[i]])[1]/100)*0.2,
#     y = seq(5,750,by=5),
#     out[[i]][], 
#     main=paste0("channel", i-1),
#     xlab="time",
#     ylab="freq"
#     )
# }
# ## plot the difference between all 'good' signals
# par(mfcol=c(1,1))
# image(
#     x = seq(dim(out[[i]])[1])*0.2,
#     y = seq(5,750,by=5),
#     xor(out[[6]], xor(out[[4]], xor(out[[3]], xor(out[[1]], out[[2]]) ))), 
#     main=paste0("channel", i-1),
#     xlab="time",
#     ylab="freq"
# )

