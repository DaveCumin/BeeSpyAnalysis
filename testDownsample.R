rm(list=ls()) #remove all objects
dev.off()
library(tictoc)

tic("read data in")
data = read.csv(file.path("~/Documents/Circadian/BeeSpy/2024 04 Germany/20240417 hive test with colony/(array-1d denoised, chn=0) 17042400.000...005 1.csv"))
toc()
head(data)
plot(data$values[seq(5000)])

spec = read.csv(file.path("~/Documents/Circadian/BeeSpy/2024 04 Germany/20240417 hive test with colony/(array-2d spec) Unbenannt 1.csv"))
head(spec)


## Show that the spectrogram from Manu's program is similar to the spectrogram created from the raw data
par(mfcol=c(4,1))

## Manu's program
image(seq(nrow(spec))/5, ##divide by 5 to make the time numbers correct
      seq(ncol(spec))*5, ##divide by 5 to make the freq numbers correct
      as.matrix(spec), col = hcl.colors(100, "temps"),
      main = "Data from Manu's program")

## Make a spectrogram from the data
require(signal)
##first do a linear detrend of the data
tic("detrend")
data_detrended = data$values - predict(lm(data$values~seq(length(data$values))))
toc()

tic("Make spectrogram")
ssd = specgram(data_detrended, 
         n = 1000, 
         Fs = 5000, 
         overlap = 0
         )
toc()

plot(ssd, 
     col = hcl.colors(100, "temps"), 
     main="created spectrogram in R")



## Now show that subsampling makes little difference to the spectrogram

plot(ssd, 
     ylim=c(0,500), 
     col = hcl.colors(100, "temps"), 
     main="Only up to 500Hz")


downsampled = data$values[seq(1, length(data$values), by=5)]
downsampled_detrended = downsampled - predict(lm(downsampled~seq(length(downsampled))))
ssd_1 = specgram(downsampled_detrended, 
                 n = 200, 
                 Fs = 1000, 
                 overlap = 0
)


plot(ssd_1, 
     ylim=c(0,500), 
     col = hcl.colors(100, "temps"), 
     main="downsampled to 1k")
