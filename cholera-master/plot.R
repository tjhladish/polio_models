d = read.table('output', header=F, sep=' ')
y_max = max(unlist(d))

png('cholera.png', width=2000, height=1200, res=150)
plot(d[,1], ylim=c(0,y_max), type = 'l', col=1, xlab='Time (3 mo periods)', ylab='Compartment size')
lines(d[,2], col=2)
lines(d[,3], col=3)
lines(d[,4], col=4)
dev.off()
