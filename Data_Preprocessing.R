library(MALDIquant)
library(MALDIquantForeign)
setwd('./R/Sample_Data_Folder_for_xcms/')
CJ1 <- importMzMl('./CJ1.mzML')
CJ2 <- importMzMl('./CJ2.mzML')
CJ3 <- importMzMl('./CJ3.mzML')
CS1 <- importMzMl('./CS1.mzML')
CS2 <- importMzMl('./CS2.mzML')
CS4 <- importMzMl('./CS4.mzML')
GM1 <- importMzMl('./GM_1.mzML')
GM2 <- importMzMl('./GM_2.mzML')
GM3 <- importMzMl('./GM_3.mzML')
GP1 <- importMzMl('./GP_1.mzML')
GP2 <- importMzMl('./GP_2.mzML')
GP3 <- importMzMl('./GP_3.mzML')

avgCJ1 <- averageMassSpectra(CJ1)
avgCJ2 <- averageMassSpectra(CJ2)
avgCJ3 <- averageMassSpectra(CJ3)
avgCS1 <- averageMassSpectra(CS1)
avgCS2 <- averageMassSpectra(CS2)
avgCS4 <- averageMassSpectra(CS4)
avgGM1 <- averageMassSpectra(GM1)
avgGM2 <- averageMassSpectra(GM2)
avgGM3 <- averageMassSpectra(GM3)
avgGP1 <- averageMassSpectra(GP1)
avgGP2 <- averageMassSpectra(GP2)
avgGP3 <- averageMassSpectra(GP3)

LAESIspectra <- c(avgCJ1, avgCJ2, avgCJ3, avgCS1, avgCS2, avgCS4, avgGM1, avgGM2, avgGM3, avgGP1, avgGP2, avgGP3)
LAESIspectra.info <- read.table("spectra_info.csv", sep = ",", header = TRUE)

par(mfrow = c(3,4))
plot(LAESIspectra[[1]], main = "CJ1")
plot(LAESIspectra[[4]], main = "CS1")
plot(LAESIspectra[[7]], main = "GM1")
plot(LAESIspectra[[10]], main = "GP1")
plot(LAESIspectra[[2]], main = "CJ2")
plot(LAESIspectra[[5]], main = "CS2")
plot(LAESIspectra[[8]], main = "GM2")
plot(LAESIspectra[[11]], main = "GP2")
plot(LAESIspectra[[3]], main = "CJ3")
plot(LAESIspectra[[6]], main = "CS4")
plot(LAESIspectra[[9]], main = "GM3")
plot(LAESIspectra[[12]], main = "GP3")
dev.off()

LAESIspectraSmoothed <- smoothIntensity(LAESIspectra, method="SavitzkyGolay", halfWindowSize=10)
baseline <- estimateBaseline(LAESIspectraSmoothed[[1]], method="SNIP",iterations=150)
plot(LAESIspectraSmoothed[[1]])
lines(baseline, col="red", lwd=2)
LAESIspectraSmoothedBaselineCorrected <- removeBaseline(LAESIspectraSmoothed, method="SNIP",iterations=150)
LAESIspectraSmoothedBaselineCorrectedNormalized <- calibrateIntensity(LAESIspectraSmoothedBaselineCorrected, method="TIC")
LAESIspectraSmoothedBaselineCorrectedNormalizedAligned <- alignSpectra(LAESIspectraSmoothedBaselineCorrectedNormalized)

avgSpectraReplicates <- averageMassSpectra(LAESIspectraSmoothedBaselineCorrectedNormalizedAligned, labels=LAESIspectra.info$species)
avgSpectraReplicates.info <-
LAESIspectra.info[!duplicated(LAESIspectra.info$species), ]

plot(avgSpectraReplicates[[1]])
lines(noise, col="red")
lines(noise[, 1], 5*noise[, 2], col="blue")
LAESIPeaks <- detectPeaks(avgSpectraReplicates, SNR=5, halfWindowSize=10)
plot(avgSpectraReplicates[[1]], xlim = c(400, 450))
points(LAESIPeaks[[1]], col="red", pch=4)
LAESIPeaksBinned <- binPeaks(LAESIPeaks)
LAESIPeaksBinnedFiltered <- filterPeaks(LAESIPeaksBinned, minFrequency=c(0.5, 0.5),labels=avgSpectraReplicates.info$type, mergeWhitelists=TRUE)

LAESIfeatureMatrix <- intensityMatrix(LAESIPeaksBinnedFiltered, avgSpectraReplicates)
rownames(LAESIfeatureMatrix) <- avgSpectraReplicates.info$species

LAESIdistanceMatrix <- dist(LAESIfeatureMatrix, method="euclidean")
LAESIhClust <- hclust(LAESIdistanceMatrix, method="complete")
plot(LAESIhClust, hang = -1)





