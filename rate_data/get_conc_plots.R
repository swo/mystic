library(RColorBrewer)
col <- rev(colorRampPalette(brewer.pal(9, 'RdBu'))(50))

files <- c('conc_1.csv','conc_2.csv','conc_3.csv','conc_4.csv','conc_5.csv','conc_6.csv','conc_7.csv','conc_8.csv', 'conc_9.csv', 'conc_10.csv')
names <- c('O', 'C', 'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-', 'CO2', 'CH4')
imax <- length(files)

par(mfrow=c(2, 5), mar=c(4,1,1,1), new=F)
for (i in 1: imax) {
    file <- files[i]
    name <- names[i]
    
    cat(c(file, rate, "\n"))
    raw <- read.csv(file, header=F)
    mat <- as.matrix(raw)
    
    # need some fancy footwork with matrix so that time=depth=0 appears at top left
    image(t(apply(mat, 2, rev)), col=col, bty='n', xaxt='n', yaxt='n', xlab=name, cex.lab=1.0)
}
