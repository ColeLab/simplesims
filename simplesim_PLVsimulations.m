%This script requires PLV code from:
%https://praneethnamburi.wordpress.com/2011/08/10/plv/


%=======
numTrials=50;
numSubjs=25;

%Create sources of non-shared noise
nsNoiseA=.25*randn(200,numTrials,numSubjs);
nsNoiseB=.25*randn(200,numTrials,numSubjs);
%Create sources of non-shared activity
nsActivityA=randn(200,numTrials,numSubjs);
nsActivityB=randn(200,numTrials,numSubjs);
%Create source of shared activity
sharedActivity=randn(200,numTrials,numSubjs);
%Create time series for RegionA
sAct1_nsAct1_nsN1_A=sharedActivity+nsActivityA+nsNoiseA;
sAct1_nsAct1_nsN1_B=sharedActivity+nsActivityB+nsNoiseB;





disp('==Increased shared activity amplitude (2x) in both regions==')
sAct2_nsAct1_nsN1_A=(2*sharedActivity)+nsActivityA+nsNoiseA;
sAct2_nsAct1_nsN1_B=(2*sharedActivity)+nsActivityB+nsNoiseB;
plvBefore=zeros(numSubjs,1);
plvAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    srate = 100; %Hz
    filtSpec.order = 50;
    filtSpec.range = [10 20]; %Hz
    dat=zeros(2,200,numTrials);
    dat(1,:,:)=sAct1_nsAct1_nsN1_A(:,:,subjNum);
    dat(2,:,:)=sAct1_nsAct1_nsN1_B(:,:,subjNum);
    plvB = pn_eegPLV(dat, srate, filtSpec);
    plvBefore(subjNum)=mean(mean(plvB(25:end, 1, 2)));
    dat=zeros(2,200,numTrials);
    dat(1,:,:)=sAct2_nsAct1_nsN1_A(:,:,subjNum);
    dat(2,:,:)=sAct2_nsAct1_nsN1_B(:,:,subjNum);
    plvA = pn_eegPLV(dat, srate, filtSpec);
    plvAfter(subjNum)=mean(mean(plvA(25:end, 1, 2)));
end
[h,pval,ci,stats]=ttest2(plvAfter, plvBefore);
disp(['PLV before: ' num2str(mean(plvBefore)) ', After: ' num2str(mean(plvAfter)) ';  T-test t-value: ' num2str(stats.tstat) ', p-value: ' num2str(pval)])


disp('==Increased non-shared activity amplitude (2x) in both regions==')
sAct1_nsAct2_nsN1_A=sharedActivity+(2*nsActivityA)+nsNoiseA;
sAct1_nsAct2_nsN1_B=sharedActivity+(2*nsActivityB)+nsNoiseB;
plvBefore=zeros(numSubjs,1);
plvAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    srate = 100; %Hz
    filtSpec.order = 50;
    filtSpec.range = [10 20]; %Hz
    dat=zeros(2,200,numTrials);
    dat(1,:,:)=sAct1_nsAct1_nsN1_A(:,:,subjNum);
    dat(2,:,:)=sAct1_nsAct1_nsN1_B(:,:,subjNum);
    plvB = pn_eegPLV(dat, srate, filtSpec);
    plvBefore(subjNum)=mean(mean(plvB(25:end, 1, 2)));
    dat=zeros(2,200,numTrials);
    dat(1,:,:)=sAct1_nsAct2_nsN1_A(:,:,subjNum);
    dat(2,:,:)=sAct1_nsAct2_nsN1_B(:,:,subjNum);
    plvA = pn_eegPLV(dat, srate, filtSpec);
    plvAfter(subjNum)=mean(mean(plvA(25:end, 1, 2)));
end
[h,pval,ci,stats]=ttest2(plvAfter, plvBefore);
disp(['PLV before: ' num2str(mean(plvBefore)) ', After: ' num2str(mean(plvAfter)) ' T-test t-value: ' num2str(stats.tstat) ', p-value: ' num2str(pval)])



disp('==Increased shared (2x) & non-shared (2x) activity amplitude in both regions==')
sAct2_nsAct2_nsN1_A=(2*sharedActivity)+(2*nsActivityA)+nsNoiseA;
sAct2_nsAct2_nsN1_B=(2*sharedActivity)+(2*nsActivityB)+nsNoiseB;
plvBefore=zeros(numSubjs,1);
plvAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    srate = 100; %Hz
    filtSpec.order = 50;
    filtSpec.range = [10 20]; %Hz
    dat=zeros(2,200,numTrials);
    dat(1,:,:)=sAct1_nsAct1_nsN1_A(:,:,subjNum);
    dat(2,:,:)=sAct1_nsAct1_nsN1_B(:,:,subjNum);
    plvB = pn_eegPLV(dat, srate, filtSpec);
    plvBefore(subjNum)=mean(mean(plvB(25:end, 1, 2)));
    dat=zeros(2,200,numTrials);
    dat(1,:,:)=sAct2_nsAct2_nsN1_A(:,:,subjNum);
    dat(2,:,:)=sAct2_nsAct2_nsN1_B(:,:,subjNum);
    plvA = pn_eegPLV(dat, srate, filtSpec);
    plvAfter(subjNum)=mean(mean(plvA(25:end, 1, 2)));
end
[h,pval,ci,stats]=ttest2(plvAfter, plvBefore);
disp(['PLV before: ' num2str(mean(plvBefore)) ', After: ' num2str(mean(plvAfter)) ' T-test t-value: ' num2str(stats.tstat) ', p-value: ' num2str(pval)])

