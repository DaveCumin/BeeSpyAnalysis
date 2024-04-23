### Take in ESF data from a single channel (as a CSV file),
### binarise the data above and below a threshold (the default is 90th centile).
### We assume that the spectrogramData has the frequency components in the first row, and
### the first column is the time data.
### ARGUMENTS
### spectrogramData :- the data from the csv file, as a matrix
### threshold       :- the threshold value to binarise the power (default at 90th centile of the entire dataset)
### freqSelect      :- vector of indicies for the frequencies to keep
### freqSelect      :- a min and max value of the frequency components to keep
binariseESF <- function(spectrogramData, threshold = NULL, freqSelect = NULL){
 
  
  ### parse the frequencies and limit if required
  if(!is.null(freqSelect)){
    freqBands = round(as.numeric(gsub("X","",colnames(spectrogramData))))
    spectrogramData = spectrogramData[,
                                      which(freqBands >= freqSelect[1] & freqBands < freqSelect[2])
                                        ]
    
  }
  
  ## convert to numeric
  mode(spectrogramData) = "numeric"
  
  ## work out the threshold if required
  if(is.null(threshold)){
    threshold = quantile(spectrogramData[,seq(2,ncol(spectrogramData))], 0.9)
  }
  ### do the thresholding
  spectrogramData[which(spectrogramData < threshold)] = 0
  spectrogramData[which(spectrogramData >= threshold)] = 1
  
  return(spectrogramData)
}



### test the code
test = read.csv(file.path("~/Documents/Circadian/BeeSpy/Exemplar Data_30minInclVideo/ExemplarDataChannel0.csv"))
test = as.matrix(test)
out = binariseESF(test)
