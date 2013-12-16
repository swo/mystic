library(RColorBrewer)
col <- rev(colorRampPalette(brewer.pal(9, 'RdBu'))(50))

files <- c('ma_1.csv','ma_2.csv','ma_3.csv','ma_4.csv','ma_5.csv', 'ma_6.csv', 'tea_1.csv','tea_2.csv','tea_3.csv','tea_4.csv')
rates <- c('iron_oxidation_(oxygen)', 'ammonia_oxidation', 'sulfur_oxidation', 'iron_oxidation_(nitrate)', 'methanotrophy_(oxygen)', 'methanotrophy_(sulfate)', 'aerobic_heterotrophy', 'nitrate_reduction', 'iron_reduction', 'sulfur_reduction')
imax <- length(files)

par(mfrow=c(2, 5), mar=c(5,1,1,1), new=F)
for (i in 1: imax) {
    file <- files[i]
    rate <- rates[i]
    
    cat(c(file, rate, "\n"))
    raw <- read.csv(file, header=F)
    mat <- as.matrix(raw)
    
    # need some fancy footwork with matrix so that time=depth=0 appears at top left
    image(t(apply(mat, 2, rev)), col=col, bty='n', xaxt='n', yaxt='n', xlab=rate, cex.lab=0.5)
}
