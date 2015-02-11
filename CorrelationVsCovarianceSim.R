library('ggplot2')
library(sapa)
require(plyr)

#Simulations use convention:
#s=shared signal, ns=non-shared signal, nsN=non-shared noise
#example: s1_ns1_nsNPt25_X = shared signal with amplitude 1, non-shared signal with amplitude 1, non-shared noise with amplitude 1, time series "X"

sink(file="CorrelationVsCovarianceSim_output.txt")

#25 subjects modeled, 200 time points each
numSubjs=25;
numTimePoints=200

##Simulating 
#Simulating two time series (RegionX & RegionY) consisting of equal-parts shared signal and independent noise...
shared=rnorm(numTimePoints*numSubjs)
nonSharedActivityX=rnorm(numTimePoints*numSubjs)
nonSharedActivityY=rnorm(numTimePoints*numSubjs)
nonSharedNoiseX=rnorm(numTimePoints*numSubjs)
nonSharedNoiseY=rnorm(numTimePoints*numSubjs)
s1_ns1_nsNPt25_X=shared+nonSharedActivityX+.25*nonSharedNoiseX
s1_ns1_nsNPt25_Y=shared+nonSharedActivityY+.25*nonSharedNoiseY
simActivity=data.frame(s1_ns1_nsNPt25_X=s1_ns1_nsNPt25_X, s1_ns1_nsNPt25_Y=s1_ns1_nsNPt25_Y, subject=factor(rep(1:numSubjs,each=numTimePoints)))
simActivity_orig=simActivity

print('Original data:')
print('Pearson correlation, Fz:')
correlationsBySubj_orig=ddply(simActivity, .(subject), function(x) {atanh(cor(x$s1_ns1_nsNPt25_X, x$s1_ns1_nsNPt25_Y))})
print(mean(correlationsBySubj_orig$V1))
print('Covariance:')
covBySubj_orig=ddply(simActivity, .(subject), function(x) {cov(x$s1_ns1_nsNPt25_X, x$s1_ns1_nsNPt25_Y)})
print(mean(covBySubj_orig$V1))




print('==Increased shared variance amplitude (2x), both regions:')
shared=rnorm(numTimePoints*numSubjs)
nonSharedActivityX=rnorm(numTimePoints*numSubjs)
nonSharedActivityY=rnorm(numTimePoints*numSubjs)
nonSharedNoiseX=rnorm(numTimePoints*numSubjs)
nonSharedNoiseY=rnorm(numTimePoints*numSubjs)
s2_ns1_nsNPt25_X=(2*shared)+nonSharedActivityX+.25*nonSharedNoiseX
s2_ns1_nsNPt25_Y=(2*shared)+nonSharedActivityY+.25*nonSharedNoiseY
simActivity$s2_ns1_nsNPt25_X=s2_ns1_nsNPt25_X
simActivity$s2_ns1_nsNPt25_Y=s2_ns1_nsNPt25_Y
#Tests
print("Pearson corr diff, Fz:")
correlationsBySubj=ddply(simActivity, .(subject), function(x) {atanh(cor(x$s2_ns1_nsNPt25_X, x$s2_ns1_nsNPt25_Y))})
print(mean(correlationsBySubj$V1)-mean(correlationsBySubj_orig$V1))
tvals=t.test(correlationsBySubj$V1,correlationsBySubj_orig$V1)
cat("T-value: ",tvals$statistic,", p-value:",tvals$p.value,"\n")
print("Covariance diff:")
covBySubj=ddply(simActivity, .(subject), function(x) {cov(x$s2_ns1_nsNPt25_X, x$s2_ns1_nsNPt25_Y)})
print(mean(covBySubj$V1)-mean(covBySubj_orig$V1))
tvals=t.test(covBySubj$V1,covBySubj_orig$V1)
cat("T-value: ",tvals$statistic,", p-value:",tvals$p.value,"\n")


#Graph variables
graphVar_X=c(s1_ns1_nsNPt25_X[1:numTimePoints],s2_ns1_nsNPt25_X[1:numTimePoints])
graphVar_Y=c(s1_ns1_nsNPt25_Y[1:numTimePoints],s2_ns1_nsNPt25_Y[1:numTimePoints])
graphVar_groups=c(rep('A',numTimePoints),rep('B',numTimePoints))
graphVar_manipulation=c(rep("Increased shared signal",2*numTimePoints))
graphVar_numRegions=c(rep("Both regions",2*numTimePoints))






print('==Increased unshared variance amplitude (2x) in both regions:')
shared=rnorm(numTimePoints*numSubjs)
nonSharedActivityX=rnorm(numTimePoints*numSubjs)
nonSharedActivityY=rnorm(numTimePoints*numSubjs)
nonSharedNoiseX=rnorm(numTimePoints*numSubjs)
nonSharedNoiseY=rnorm(numTimePoints*numSubjs)
s1_ns2_nsNPt25_X=shared+2*nonSharedActivityX+.25*nonSharedNoiseX
s1_ns2_nsNPt25_Y=shared+2*nonSharedActivityY+.25*nonSharedNoiseY
simActivity$s1_ns2_nsNPt25_X=s1_ns2_nsNPt25_X
simActivity$s1_ns2_nsNPt25_Y=s1_ns2_nsNPt25_Y
#Tests
print("Pearson corr diff, Fz:")
correlationsBySubj=ddply(simActivity, .(subject), function(x) {atanh(cor(x$s1_ns2_nsNPt25_X, x$s1_ns2_nsNPt25_Y))})
print(mean(correlationsBySubj$V1)-mean(correlationsBySubj_orig$V1))
tvals=t.test(correlationsBySubj$V1,correlationsBySubj_orig$V1)
cat("T-value: ",tvals$statistic,", p-value:",tvals$p.value,"\n")
print("Covariance diff:")
covBySubj=ddply(simActivity, .(subject), function(x) {cov(x$s1_ns2_nsNPt25_X, x$s1_ns2_nsNPt25_Y)})
print(mean(covBySubj$V1)-mean(covBySubj_orig$V1))
tvals=t.test(covBySubj$V1,covBySubj_orig$V1)
cat("T-value: ",tvals$statistic,", p-value:",tvals$p.value,"\n")

#Graph variables
graphVar_X=c(graphVar_X,s1_ns1_nsNPt25_X[1:numTimePoints],s1_ns2_nsNPt25_X[1:numTimePoints])
graphVar_Y=c(graphVar_Y,s1_ns1_nsNPt25_Y[1:numTimePoints],s1_ns2_nsNPt25_Y[1:numTimePoints])
graphVar_groups=c(graphVar_groups,rep('A',numTimePoints),rep('B',numTimePoints))
graphVar_manipulation=c(graphVar_manipulation,rep("Increased unshared signal",2*numTimePoints))
graphVar_numRegions=c(graphVar_numRegions,rep("Both regions",2*numTimePoints))





print('==Equal increase in shared & unshared variance amplitude (2x) in both regions (2x shared signal, 2x unshared signal, 1x noise):')
shared=rnorm(numTimePoints*numSubjs)
nonSharedActivityX=rnorm(numTimePoints*numSubjs)
nonSharedActivityY=rnorm(numTimePoints*numSubjs)
nonSharedNoiseX=rnorm(numTimePoints*numSubjs)
nonSharedNoiseY=rnorm(numTimePoints*numSubjs)
s2_ns2_nsNPt25_X=(2*shared)+(2*nonSharedActivityX)+.25*nonSharedNoiseX
s2_ns2_nsNPt25_Y=(2*shared)+(2*nonSharedActivityY)+.25*nonSharedNoiseY
simActivity$s2_ns2_nsNPt25_X=s2_ns2_nsNPt25_X
simActivity$s2_ns2_nsNPt25_Y=s2_ns2_nsNPt25_Y
#Tests
print("Pearson corr diff, Fz:")
correlationsBySubj=ddply(simActivity, .(subject), function(x) {atanh(cor(x$s2_ns2_nsNPt25_X, x$s2_ns2_nsNPt25_Y))})
print(mean(correlationsBySubj$V1)-mean(correlationsBySubj_orig$V1))
tvals=t.test(correlationsBySubj$V1,correlationsBySubj_orig$V1)
cat("T-value: ",tvals$statistic,", p-value:",tvals$p.value,"\n")
print("Covariance diff:")
covBySubj=ddply(simActivity, .(subject), function(x) {cov(x$s2_ns2_nsNPt25_X, x$s2_ns2_nsNPt25_Y)})
print(mean(covBySubj$V1)-mean(covBySubj_orig$V1))
tvals=t.test(covBySubj$V1,covBySubj_orig$V1)
cat("T-value: ",tvals$statistic,", p-value:",tvals$p.value,"\n")


#Graph variables
graphVar_X=c(graphVar_X,s1_ns1_nsNPt25_X[1:numTimePoints],s2_ns2_nsNPt25_X[1:numTimePoints])
graphVar_Y=c(graphVar_Y,s1_ns1_nsNPt25_Y[1:numTimePoints],s2_ns2_nsNPt25_Y[1:numTimePoints])
graphVar_groups=c(graphVar_groups,rep('A',numTimePoints),rep('B',numTimePoints))
graphVar_manipulation=c(graphVar_manipulation,rep("Increased shared & unshared signals",2*numTimePoints))
graphVar_numRegions=c(graphVar_numRegions,rep("Both regions",2*numTimePoints))



#Scatter plots
customColors=c('darkblue','red')
gtbl=data.frame(X=graphVar_X,Y=graphVar_Y,groups=graphVar_groups,manipulation=factor(graphVar_manipulation,ordered=T,levels=c("Increased shared signal","Increased unshared signal","Increased shared & unshared signals")),numRegions=factor(graphVar_numRegions))
ggplot(gtbl, aes(x=X,y=Y,color=groups)) + geom_point() + theme_bw() + scale_colour_manual(values = customColors) + xlab("Region X activity") + ylab("Region Y activity") + facet_grid( ~ manipulation) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(file="Rplot_scatter.pdf", width=8, height=4)


sink()
