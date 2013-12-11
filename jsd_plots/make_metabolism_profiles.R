# produces the curves for each metabolism

names <- c('nitrification', 'sulfur oxidation', 'nit on Fe', 'aerobic heterotrophy', 'nitrate reduction', 'iron reduction', 'sulfate reduction', 'methanog');
rates <- read.csv('rates.csv', header=F)

prop.rates <- prop.table(as.matrix(rates), margin=1)

#pdf('figs/rates.pdf')
par(mfrow=c(3, 3))
for (i in 1:8) {
    dat <- prop.rates[, i]
    name <- names[i]
    plot(dat, type='l', main=name, xlab='depth below thermocline', ylab='fraction of total metabolism')
}
#dev.off()