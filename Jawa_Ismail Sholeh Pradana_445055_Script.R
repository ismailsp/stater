#load library
library(raster)
library(zoo)
rasterOptions(progress = "text", tmpdir = paste0(getwd(), "/tmp"))

#stack raster data
ch <- stack("precip_chirps_jawa_8215.tif")
ch <- ch[[1:215]]
ndvi <- stack("ndvi_gimms_res_jawa_8215.tif")
ndvi.res <- raster::resample(ndvi, ch, "bilinear")
writeRaster(ndvi.res, "ndvi_resample_crop_1.tif", format = "GTiff", datatype = "FLT4S", overwrite=TRUE)

#membuat vektor rerata dari semua bulan
ndvi.vc <- cellStats(ndvi.res, stat = mean, na.rm = T)
ch.vc <- cellStats(ch, stat = mean, na.rm = T)

#konversi menjadi TS
#month.idx <- seq(as.Date("2002/7/1"), as.Date("2020/6/1"), by = "month") 
ndvi.ts <- ts(ndvi.vc, start = c(2002, 7), frequency = 12)
ch.ts <- ts(ch.vc, start = c(2002, 7), frequency = 12)

#plot data ts
plot(ndvi.ts, type = "l")
plot(ch.ts, type = "l")

#lakukan perhitungan acf dan pacf
par(mfrow = c(1,2))
acf(c(ndvi.ts),main='ACF ndvi')
pacf(c(ndvi.ts),main='PACF ndvi')

par(mfrow = c(1,2))
acf(c(ch.ts),main='ACF CH')
pacf(c(ch.ts),main='PACF CH')

#check apakah ada seasonality di situ
#buat data menjadi stationary dan cek acf dan pacf
ndvi.st <- diff(log10(ndvi.ts))
ch.st <- diff(log10(ch.ts))

par(mfrow = c(1,2))
acf(c(ch.st), main='ACF CH Diff Log')
pacf(c(ch.st), main='PACF CH Diff Log')

par(mfrow = c(1,2))
acf(c(ndvi.st), main='ACF ndvi Diff Log')
pacf(c(ndvi.st), main='PACF ndvi Diff Log')

#hitung cross-correlation data
par(mfrow = c(1,1))
ccf.res <- ccf(c(ndvi.st), c(ch.st), lag.max = 12)
print(ccf.res)

#cek maksimum lag pada bulan berapa 
max(ccf.res$acf)
ccf.res$lag[which.max(ccf.res$acf)]

#buat peta lag dan nindvi korelasi maksimal dari curah hujan dan ndvi
fun.ccfmax <- function(x){
  if (isTRUE(any(is.na(x)))){ 
    return(c(NA, NA))
  } else if (isTRUE(any(x <= 0))){
    return(c(NA, NA))
  } else {
    chr2 <- x[1:215]
    ndvi2 <- x[216:430]
    ndvi.st <- diff(log10(chr2))
    ch.st <- diff(log10(ndvi2))
    ccf <- ccf(c(ndvi.st), c(ch.st), lag.max = 12, plot = F)
    mx <- max(ccf$acf)
    lag <- ccf$lag[which.max(ccf$acf)]
    return(c(mx, lag))
  }                             
}                               

ccf.raster <- calc(stack(ch, ndvi.res), fun.ccfmax)
plot(ccf.raster)
writeRaster(ccf.raster, "ccf_max_lag.tif", format = "GTiff", datatype = "FLT4S", overwrite=TRUE)

