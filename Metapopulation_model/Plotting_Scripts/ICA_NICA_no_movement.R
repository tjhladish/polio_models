dir = '/Users/Celeste/Desktop/multipatch_model/SC_statistic/'
dir1 = '/Users/Celeste/Desktop/polio-small-pop/polio-small-pop/'
dir2 = '/Users/Celeste/Desktop/multipatch_model/absRisk/'
dir3 = '/Users/Celeste/Desktop/multipatch_model/Figures/'


library(RColorBrewer)

marker = list(color = brewer.pal(8,"YlOrRd"))

#one village of size 16,000
movement_comp_1x16000_0pcase = read.table(paste0(dir,'sc_0pcase_N_16000_numVil_1_migRate_0_det_1_beta_135_fast_multipatch.out'), col.names=c('time', 'E&D'))
one_pcase_1x16000 = read.table(paste0(dir,'sc_1pcase_N_16000_det_1_beta_135_fast_multipatch.out'),col.names=c('time','SC_stat'))

plot(movement_comp_1x16000_0pcase,type='l',col='black')
lines(one_pcase_1x16000,type='l',col='red')

#8 villages of size 2,000
movement_comp_8x2000_move_0_0pcase = read.table(paste0(dir,'sc_0pcase_N_2000_numVil_8_migRate_0_det_1_beta_135_fast_multipatch.out'), col.names=c('time', 'E&D'))
one_pcase_8x2000 = read.table(paste0(dir,'sc_1pcase_N_2000_numVil_8_det_1_beta_135_fast_multipatch.out'),col.names=c('time','SC_stat'))

plot(movement_comp_8x2000_move_0_0pcase,type='l',col='black')
lines(one_pcase_8x2000,type='l',col='red')

movement_comp_3x10000_move_0_0pcase = read.table(paste0(dir,'sc_0pcase_N_10000_numVil_3_migRate_0_det_1_beta_135_fast_multipatch.out'), col.names=c('time', 'E&D'))
one_pcase_3x10000 = read.table(paste0(dir,'sc_1pcase_N_10000_numVil_3_migRate_0_det_1_beta_135_fast_multipatch.out'), col.names=c('time', 'E&D'))

plot(movement_comp_3x10000_move_0_0pcase,type='l',col='black')
lines(one_pcase_3x10000,col='red')

movement_comp_2x15000_move_0_0pcase = read.table(paste0(dir,'sc_0pcase_N_15000_numVil_2_migRate_0_det_1_beta_135_fast_multipatch.out'), col.names=c('time', 'E&D'))
one_pcase_2x15000 = read.table(paste0(dir,'sc_1pcase_N_15000_numVil_2_migRate_0_det_1_beta_135_fast_multipatch.out'), col.names=c('time', 'E&D'))

plot(movement_comp_2x15000_move_0_0pcase,type='l',col='black')
lines(one_pcase_2x15000,col='red')

movement_comp_3x10000_move_0.1_0pcase = read.table(paste0(dir,'sc_0pcase_N_10000_numVil_3_migRate_0.1_det_1_beta_135_fast_multipatch.out'), col.names=c('time', 'E&D'))
one_pcase_3x10000_move_0.1 = read.table(paste0(dir,'sc_1pcase_N_10000_numVil_3_migRate_0.1_det_1_beta_135_fast_multipatch.out'), col.names=c('time', 'E&D'))

plot(movement_comp_3x10000_move_0.1_0pcase,type='l',col='black')
lines(one_pcase_3x10000_move_0.1,col='red')

movement_comp_8x2000_move_0.1_0pcase = read.table(paste0(dir,'sc_0pcase_N_2000_numVil_8_migRate_0.1_det_1_beta_135_fast_multipatch.out'), col.names=c('time', 'E&D'))
one_pcase_8x2000_move_0.1 = read.table(paste0(dir,'sc_1pcase_N_2000_numVil_8_migRate_0.1_det_1_beta_135_fast_multipatch.out'),col.names=c('time','SC_stat'))

plot(movement_comp_8x2000_move_0.1_0pcase,type='l',col='black')
lines(one_pcase_8x2000_move_0.1,type='l',col='red')

movement_comp_8x2000_move_0.5_0pcase = read.table(paste0(dir,'sc_0pcase_N_2000_numVil_8_migRate_0.5_det_1_beta_135_fast_multipatch.out'), col.names=c('time', 'E&D'))
one_pcase_8x2000_move_0.5 = read.table(paste0(dir,'sc_1pcase_N_2000_numVil_8_migRate_0.5_det_1_beta_135_fast_multipatch.out'),col.names=c('time','SC_stat'))

plot(movement_comp_8x2000_move_0.5_0pcase,type='l',col='black')
lines(one_pcase_8x2000_move_0.5,type='l',col='red')

absRisk_10000_moveRate_0 = read.csv(file=paste0(dir2,'absRisk_N_10000_numVil_3_migRate_0_det_1_beta_135_fast_multipatch.csv'),header=FALSE)
absRisk_15000_moveRate_0 = read.csv(file=paste0(dir2,'absRisk_N_15000_numVil_3_migRate_0_det_1_beta_135_fast_multipatch.csv'),header=FALSE)
absRisk_10000_moveRate_0.1 = read.csv(file=paste0(dir2,'absRisk_N_10000_numVil_3_migRate_0.1_det_1_beta_135_fast_multipatch.csv'),header=FALSE)

png(paste0(dir3,'ICA_NICA_no_movement_N_10000.png'), width=1400, height=800, res=150)
par(mfrow=c(2,1),mar=c(0,4,0,0),oma=c(3,3,1,1),mgp=c(2,1,0))
plot(movement_comp_3x10000_move_0_0pcase, type='l',lty=1, xlab = '',ylab = 'Silent circulation statistic',xlim=c(0,5),xaxt='n',col=marker$color[1])
lines(one_pcase_3x10000,type='l', col=marker$color[2])
legend('topright', legend = c(expression(paste("ICA: No movement (3 villages, ","N"["sub"],'=10000)')),expression(paste("NICA: No movement (3 villages, ","N"["sub"],'=10000)'))),col=marker$color[1:2],lwd=2, bty='n' )
plot(absRisk_10000_moveRate_0, type='l',lty=1, xlab="Time since last paralytic case (years)",ylab="Differential probability\n of circulation",xlim=c(0,5),col=marker$color[1])
abline(0,0,lty=3)
legend('bottomright', legend = c(expression(paste("NICA - ICA: No movement (3 villages, ","N"["sub"],'=10000)'))),col=marker$color[1],lwd=2, bty='n' )
mtext("Time since last detected paralytic case (years)",side=1,line=2)
dev.off()

png(paste0(dir3,'ICA_NICA_no_movement_N_15000.png'), width=1400, height=800, res=150)
par(mfrow=c(2,1),mar=c(0,4,0,0),oma=c(3,3,1,1),mgp=c(2,1,0))
plot(movement_comp_2x15000_move_0_0pcase, type='l',lty=1, xlab = '',ylab = 'Silent circulation statistic',xlim=c(0,5),xaxt='n',col=marker$color[1])
lines(one_pcase_2x15000,type='l', col=marker$color[2])
legend('topright', legend = c(expression(paste("ICA: No movement (2 villages, ","N"["sub"],'=15000)')),expression(paste("NICA: No movement (2 villages, ","N"["sub"],'=15000)'))),col=marker$color[1:2],lwd=2, bty='n' )
plot(absRisk_15000_moveRate_0, type='l',lty=1, xlab="Time since last paralytic case (years)",ylab="Differential probability\n of circulation",xlim=c(0,5),col=marker$color[1],yaxt='n')
ticks<-c(-.01,-.005,0)
axis(2,at=ticks,labels=ticks)
abline(0,0,lty=3)
legend('bottomright', legend = c(expression(paste("NICA - ICA: No movement (2 villages, ","N"["sub"],'=15000)'))),col=marker$color[1],lwd=2, bty='n' )
mtext("Time since last detected paralytic case (years)",side=1,line=2)
dev.off()

png(paste0(dir3,'ICA_NICA_moveRate_1_per_10_N_10000.png'), width=1400, height=800, res=150)
par(mfrow=c(2,1),mar=c(0,4,0,0),oma=c(3,3,1,1),mgp=c(2,1,0))
plot(movement_comp_3x10000_move_0.1_0pcase, type='l',lty=1, xlab = '',ylab = 'Silent circulation statistic',xlim=c(0,5),xaxt='n',col=marker$color[1])
lines(one_pcase_3x10000_move_0.1,type='l', col=marker$color[2])
legend('topright', legend = c(expression(paste("ICA: Move rate 0.1 (3 villages, ","N"["sub"],'=10000)')),expression(paste("NICA: Move rate 0.1 (3 villages, ","N"["sub"],'=10000)'))),col=marker$color[1:2],lwd=2, bty='n' )
plot(absRisk_10000_moveRate_0.1, type='l',lty=1, xlab="Time since last paralytic case (years)",ylab="Differential probability\n of circulation",xlim=c(0,5),col=marker$color[1])
abline(0,0,lty=3)
legend('bottomright', legend = c(expression(paste("NICA - ICA: Move rate 0.1 (3 villages, ","N"["sub"],'=10000)'))),col=marker$color[1],lwd=2, bty='n' )
mtext("Time since last detected paralytic case (years)",side=1,line=2)
dev.off()


