#Make sure these libraries are installed prior to running this code
library('ggplot2')
library(sapa)
require(plyr)
library('entropy')
library(grid)
#Included functions
source('cohFunc.R')
source('specCov.R')


#Simulations from Cole et al. "Reconceptualizing brain network change as shared signal dynamics".
#Contact: Michael W. Cole, mwcole@mwcole.net
#Eight simulations manipulating shared and unshared variance between two time series, simulating 25 subjects each. Seven functional connectivity measures are used: covariance, spectral covariance, Pearson correlation, coherence, mutual information, regression/PPI, and Spearman correlation

#Saving output to a text file
sink(file="CorrelationVsCovarianceSim_WithOtherMeasures_output.txt")


#25 subjects modeled, 200 time points each
numSubjs=25;
numTimePoints=200

#Manipulations to time series
sharedSignalLevelsToTest=c(.5,1,2)
unsharedSignalLevelsToTest=c(.5,1,2)

#Data structures to save to
simRAWData=data.frame(matrix(ncol=5,nrow=numTimePoints*numSubjs))
names(simRAWData)=c('signalLevel','unsharedSignalLevel','subjectNum','regionXdata','regionYdata')

##Simulating
matRowNum=1
timeIndex=1
for(signalLevel in sharedSignalLevelsToTest) {
  for(unsharedSignalLevel in unsharedSignalLevelsToTest) {
    
    #Simulating two time series (X & Y), consisting of shared signal, unshared signal, and noise
    shared=vector(length=numTimePoints*numSubjs)
    nonSharedActivityX=vector(length=numTimePoints*numSubjs)
    nonSharedActivityY=vector(length=numTimePoints*numSubjs)
    nonSharedNoiseX=vector(length=numTimePoints*numSubjs)
    nonSharedNoiseY=vector(length=numTimePoints*numSubjs)
    indexStart=1
    for(subjNum in 1:numSubjs) {
      #Generating Gaussian random time points, separately for each component (shared, unshared, noise)
      shared[indexStart:(indexStart+numTimePoints-1)]=rnorm(numTimePoints)
      nonSharedActivityX[indexStart:(indexStart+numTimePoints-1)]=rnorm(numTimePoints)
      nonSharedActivityY[indexStart:(indexStart+numTimePoints-1)]=rnorm(numTimePoints)
      nonSharedNoiseX[indexStart:(indexStart+numTimePoints-1)]=rnorm(numTimePoints)
      nonSharedNoiseY[indexStart:(indexStart+numTimePoints-1)]=rnorm(numTimePoints)
      indexStart=indexStart+numTimePoints
    }
    #Combining randomly generated components to create time series
    timeseriesX=(signalLevel*shared)+(unsharedSignalLevel*nonSharedActivityX)+(.25*nonSharedNoiseX)
    timeseriesY=(signalLevel*shared)+(unsharedSignalLevel*nonSharedActivityY)+(.25*nonSharedNoiseY)
    
    #Saving simulated data
    simRAWData[timeIndex:(timeIndex+numTimePoints*numSubjs-1),]$regionXdata=timeseriesX
    simRAWData[timeIndex:(timeIndex+numTimePoints*numSubjs-1),]$regionYdata=timeseriesY
    simRAWData[timeIndex:(timeIndex+numTimePoints*numSubjs-1),]$subjectNum=rep(1:numSubjs,each=numTimePoints)
    simRAWData[timeIndex:(timeIndex+numTimePoints*numSubjs-1),]$signalLevel=rep(signalLevel,each=numTimePoints*numSubjs)
    simRAWData[timeIndex:(timeIndex+numTimePoints*numSubjs-1),]$unsharedSignalLevel=rep(unsharedSignalLevel,each=numTimePoints*numSubjs)
    
    matRowNum=matRowNum+1
    timeIndex=timeIndex+numTimePoints*numSubjs
  }
}



##Testing for functional connectivity effects

#Function to test all of the functional connectivity measures
fctests <- function(originalInputDFrame, changedInputDFrame) {
  
  outputText=''
  
  #Covariance
  outputText=paste(outputText,'==Covariance difference:==\n')
  #Original signal and noise
  fcBySubj_orig=ddply(originalInputDFrame, .(subjectNum), function(x) {cov(x$regionXdata, x$regionYdata)})
  #Changed data
  fcBySubj=ddply(changedInputDFrame, .(subjectNum), function(x) {cov(x$regionXdata, x$regionYdata)})
  meanDiff_cov=mean(fcBySubj$V1)-mean(fcBySubj_orig$V1)
  tvals_cov=t.test(fcBySubj$V1,fcBySubj_orig$V1)
  outputText=paste(outputText,'Mean difference:',meanDiff_cov,'\n')
  outputText=paste(outputText,'T-value:',tvals_cov$statistic,'\n')
  outputText=paste(outputText,'P-value:',tvals_cov$p.value,'\n')
  if(tvals_cov$p.value<0.05){
    outputText=paste(outputText,'**Significant (P<0.05)**','\n')
  }
  
  #Spectral covariance
  outputText=paste(outputText,'==Spectral covariance difference:==\n')
  #Original signal and noise
  fcBySubj_orig=ddply(originalInputDFrame, .(subjectNum), function(x) {specCov(x$regionXdata, x$regionYdata)})
  #Changed data
  fcBySubj=ddply(changedInputDFrame, .(subjectNum), function(x) {specCov(x$regionXdata, x$regionYdata)})
  meanDiff_scov=mean(fcBySubj$V1)-mean(fcBySubj_orig$V1)
  tvals_scov=t.test(fcBySubj$V1,fcBySubj_orig$V1)
  outputText=paste(outputText,'Mean difference:',meanDiff_scov,'\n')
  outputText=paste(outputText,'T-value:',tvals_scov$statistic,'\n')
  outputText=paste(outputText,'P-value:',tvals_scov$p.value,'\n')
  if(tvals_scov$p.value<0.05){
    outputText=paste(outputText,'**Significant (P<0.05)**','\n')
  }
  
  #Pearson correlation
  outputText=paste(outputText,'==Pearson correlation difference:==\n')
  #Original signal and noise
  fcBySubj_orig=ddply(originalInputDFrame, .(subjectNum), function(x) {atanh(cor(x$regionXdata, x$regionYdata))})
  #Changed data
  fcBySubj=ddply(changedInputDFrame, .(subjectNum), function(x) {atanh(cor(x$regionXdata, x$regionYdata))})
  meanDiff_corr=mean(fcBySubj$V1)-mean(fcBySubj_orig$V1)
  tvals_corr=t.test(fcBySubj$V1,fcBySubj_orig$V1)
  outputText=paste(outputText,'Mean difference:',meanDiff_corr,'\n')
  outputText=paste(outputText,'T-value:',tvals_corr$statistic,'\n')
  outputText=paste(outputText,'P-value:',tvals_corr$p.value,'\n')
  if(tvals_corr$p.value<0.05){
    outputText=paste(outputText,'**Significant (P<0.05)**','\n')
  }

  #Coherence
  outputText=paste(outputText,'==Coherence difference:==\n')
  #Original signal and noise
  fcBySubj_orig=ddply(originalInputDFrame, .(subjectNum), function(x) {cohFunc(x$regionXdata, x$regionYdata)})
  #Changed data
  fcBySubj=ddply(changedInputDFrame, .(subjectNum), function(x) {cohFunc(x$regionXdata, x$regionYdata)})
  meanDiff_coh=mean(fcBySubj$V1)-mean(fcBySubj_orig$V1)
  tvals_coh=t.test(fcBySubj$V1,fcBySubj_orig$V1)  
  outputText=paste(outputText,'Mean difference:',meanDiff_coh,'\n')
  outputText=paste(outputText,'T-value:',tvals_coh$statistic,'\n')
  outputText=paste(outputText,'P-value:',tvals_coh$p.value,'\n')
  if(tvals_coh$p.value<0.05){
    outputText=paste(outputText,'**Significant (P<0.05)**','\n')
  }
  
  #Mutual information
  outputText=paste(outputText,'==Mutual information difference:==\n')
  #Original signal and noise
  fcBySubj_orig=ddply(originalInputDFrame, .(subjectNum), function(x) {mi.empirical(discretize2d(x$regionXdata, x$regionYdata,16,16))})
  #Changed data
  fcBySubj=ddply(changedInputDFrame, .(subjectNum), function(x) {mi.empirical(discretize2d(x$regionXdata, x$regionYdata,16,16))})
  meanDiff_mi=mean(fcBySubj$V1)-mean(fcBySubj_orig$V1)
  tvals_mi=t.test(fcBySubj$V1,fcBySubj_orig$V1)
  outputText=paste(outputText,'Mean difference:',meanDiff_mi,'\n')
  outputText=paste(outputText,'T-value:',tvals_mi$statistic,'\n')
  outputText=paste(outputText,'P-value:',tvals_mi$p.value,'\n')
  if(tvals_mi$p.value<0.05){
    outputText=paste(outputText,'**Significant (P<0.05)**','\n')
  }
  
  #Regression/PPI, Beta X on Y
  outputText=paste(outputText,'==Regression/PPI, Beta X on Y:==\n')
  #Original signal and noise
  fcBySubj_orig=ddply(originalInputDFrame, .(subjectNum), function(x) {lm(x$regionXdata ~ x$regionYdata, na.action = NULL)$coefficients[2]})
  #Changed data
  fcBySubj=ddply(changedInputDFrame, .(subjectNum), function(x) {lm(x$regionXdata ~ x$regionYdata, na.action = NULL)$coefficients[2]})
  meanDiff_regXonY=mean(fcBySubj[,2])-mean(fcBySubj_orig[,2])
  tvals_regXonY=t.test(fcBySubj[,2],fcBySubj_orig[,2])
  outputText=paste(outputText,'Mean difference:',meanDiff_regXonY,'\n')
  outputText=paste(outputText,'T-value:',tvals_regXonY$statistic,'\n')
  outputText=paste(outputText,'P-value:',tvals_regXonY$p.value,'\n')
  if(tvals_regXonY$p.value<0.05){
    outputText=paste(outputText,'**Significant (P<0.05)**','\n')
  }
  
  #Regression/PPI, Beta Y on X
  outputText=paste(outputText,'==Regression/PPI, Beta Y on X:==\n')
  #Original signal and noise
  fcBySubj_orig=ddply(originalInputDFrame, .(subjectNum), function(x) {lm(x$regionYdata ~ x$regionXdata)$coefficients[2]})
  #Changed data
  fcBySubj=ddply(changedInputDFrame, .(subjectNum), function(x) {lm(x$regionYdata ~ x$regionXdata)$coefficients[2]})
  meanDiff_regYonX=mean(fcBySubj[,2])-mean(fcBySubj_orig[,2])
  tvals_regYonX=t.test(fcBySubj[,2],fcBySubj_orig[,2])
  outputText=paste(outputText,'Mean difference:',meanDiff_regYonX,'\n')
  outputText=paste(outputText,'T-value:',tvals_regYonX$statistic,'\n')
  outputText=paste(outputText,'P-value:',tvals_regYonX$p.value,'\n')
  if(tvals_regYonX$p.value<0.05){
    outputText=paste(outputText,'**Significant (P<0.05)**','\n')
  }

  #Spearman correlation
  outputText=paste(outputText,'==Spearman correlation:==\n')
  #Original signal and noise
  fcBySubj_orig=ddply(originalInputDFrame, .(subjectNum), function(x) {cor(x$regionXdata, x$regionYdata,method="spearman")})
  #Changed data
  fcBySubj=ddply(changedInputDFrame, .(subjectNum), function(x) {cor(x$regionXdata, x$regionYdata,method="spearman")})
  meanDiff_spcorr=mean(fcBySubj$V1)-mean(fcBySubj_orig$V1)
  tvals_spcorr=t.test(fcBySubj$V1,fcBySubj_orig$V1)  
  outputText=paste(outputText,'Mean difference:',meanDiff_spcorr,'\n')
  outputText=paste(outputText,'T-value:',tvals_spcorr$statistic,'\n')
  outputText=paste(outputText,'P-value:',tvals_spcorr$p.value,'\n')
  if(tvals_spcorr$p.value<0.05){
    outputText=paste(outputText,'**Significant (P<0.05)**','\n')
  }
  
  outputText=paste(outputText,'\n')
  
  result = list(outputText=outputText, meanDiff_cov=meanDiff_cov, tvals_cov=tvals_cov, meanDiff_corr=meanDiff_corr, tvals_corr=tvals_corr, meanDiff_coh=meanDiff_coh, tvals_coh=tvals_coh, meanDiff_scov=meanDiff_scov, tvals_scov=tvals_scov, meanDiff_mi=meanDiff_mi, tvals_mi=tvals_mi, meanDiff_regXonY=meanDiff_regXonY, tvals_regXonY=tvals_regXonY, meanDiff_regYonX=meanDiff_regYonX, tvals_regYonX=tvals_regYonX, meanDiff_spcorr=meanDiff_spcorr, tvals_spcorr=tvals_spcorr)
  return(result)
}
  
##Going through each scenario, showing key points in the shared vs. unshared 2D parameter space (corresponding to 8 scenarios in Figure S3)

print('==Equal increase in shared & unshared variance amplitude (2x) in both regions (2x shared signal, 2x unshared signal, 1x noise):')
fcResults=fctests(simRAWData[simRAWData$signalLevel==1&simRAWData$unsharedSignalLevel==1,], simRAWData[simRAWData$signalLevel==2&simRAWData$unsharedSignalLevel==2,])
cat(fcResults$outputText)

print('==Increased shared variance amplitude (2x), both regions:')
fcResults=fctests(simRAWData[simRAWData$signalLevel==1&simRAWData$unsharedSignalLevel==1,], simRAWData[simRAWData$signalLevel==2&simRAWData$unsharedSignalLevel==1,])
cat(fcResults$outputText)

print('==Increase in shared variance amplitude (2x), decrease in unshared variance amplitude (0.5x) in both regions (2x shared signal, 0.5x unshared signal, 1x noise):')
fcResults=fctests(simRAWData[simRAWData$signalLevel==1&simRAWData$unsharedSignalLevel==1,], simRAWData[simRAWData$signalLevel==2&simRAWData$unsharedSignalLevel==.5,])
cat(fcResults$outputText)

print('==Increased unshared variance amplitude (2x) in both regions:')
fcResults=fctests(simRAWData[simRAWData$signalLevel==1&simRAWData$unsharedSignalLevel==1,], simRAWData[simRAWData$signalLevel==1&simRAWData$unsharedSignalLevel==2,])
cat(fcResults$outputText)

print('==Decreased unshared variance amplitude (0.5x) in both regions:')
fcResults=fctests(simRAWData[simRAWData$signalLevel==1&simRAWData$unsharedSignalLevel==1,], simRAWData[simRAWData$signalLevel==1&simRAWData$unsharedSignalLevel==.5,])
cat(fcResults$outputText)

print('==Decrease in shared variance amplitude (0.5x), increase in unshared variance amplitude (2x) in both regions (0.5x shared signal, 2x unshared signal, 1x noise):')
fcResults=fctests(simRAWData[simRAWData$signalLevel==1&simRAWData$unsharedSignalLevel==1,], simRAWData[simRAWData$signalLevel==.5&simRAWData$unsharedSignalLevel==2,])
cat(fcResults$outputText)

print('==Decreased shared variance amplitude (0.5x), both regions:')
fcResults=fctests(simRAWData[simRAWData$signalLevel==1&simRAWData$unsharedSignalLevel==1,], simRAWData[simRAWData$signalLevel==.5&simRAWData$unsharedSignalLevel==1,])
cat(fcResults$outputText)

print('==Equal decrease in shared & unshared variance amplitude (0.5x) in both regions (0.5x shared signal, 0.5x unshared signal, 1x noise):')
fcResults=fctests(simRAWData[simRAWData$signalLevel==1&simRAWData$unsharedSignalLevel==1,], simRAWData[simRAWData$signalLevel==.5&simRAWData$unsharedSignalLevel==.5,])
cat(fcResults$outputText)


sink()
