library(gplots)
library(RColorBrewer)

rates <- read.csv('rates.csv', header=F, col.names=c('ammox', 'HS- ox', 'nit on Fe', 'O2 resp', 'NO3- resp', 'Fe3+ resp', 'SO42- resp', 'methanog'))
d <- dist(rates)
dm <- as.matrix(d)
m <- as.matrix(rates)

# check for proportional rates (margin=1 means rowwise)
prop.rates <- as.matrix(prop.table(m, margin=1))
dpm <- as.matrix(dist(prop.rates))

# get log proportional
log.prop.rates <- log10(prop.rates + 1e-12)
ldpm <- as.matrix(dist(log.prop.rates))

# grab the JSD matrix
jsd <- read.csv('rates_jsd.csv', header=F)
jd <- as.dist(jsd)
jdm <- as.matrix(jd)
jsd.hc <- hclust(jd)

col <- brewer.pal(11, 'RdBu')
#heatmap.2(prop.rates, trace='none', col=col)

# the proportional rates heatmap
#heatmap.2(log.prop.rates, trace='none', col=col)

# the jsd distances
library(fields)
#pdf('jsd_plot.pdf')
image.plot(jdm, asp=1, legend.lab='JSD')
#dev.off()

# the jsd clustering
#plot(jsd.hc)