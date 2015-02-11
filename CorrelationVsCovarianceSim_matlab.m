
%Simulations use convention:
%s=shared signal, ns=non-shared signal, nsN=non-shared noise
%example: sAct1_nsAct1_nsN1_A = shared signal with amplitude 1, non-shared signal with amplitude 1, non-shared noise with amplitude 1, time series "A"

%25 subjects modeled, 200 time points each
numSubjs=25;

disp('====Simulating data for 2 brain regions from 25 subjects by creating random normally distributed time series for each===')
disp('==Each region''s time series consists of: Shared (across regions) activity, non-shared activity, non-shared noise (MR & physiological)==')


%Create sources of non-shared noise
nsNoiseA=0.25*randn(200,numSubjs);
nsNoiseB=0.25*randn(200,numSubjs);
%Create sources of non-shared activity
nsActivityA=randn(200,numSubjs);
nsActivityB=randn(200,numSubjs);
%Create source of shared activity
sharedActivity=randn(200,numSubjs);
%Create time series for RegionA
sAct1_nsAct1_nsN1_A=sharedActivity+nsActivityA+nsNoiseA;
sAct1_nsAct1_nsN1_B=sharedActivity+nsActivityB+nsNoiseB;
disp(' ')

disp('==Increased shared activity amplitude (2x) in both regions==')
sAct2_nsAct1_nsN1_A=(2*sharedActivity)+nsActivityA+nsNoiseA;
sAct2_nsAct1_nsN1_B=(2*sharedActivity)+nsActivityB+nsNoiseB;
corrBefore=zeros(numSubjs,1);
corrAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    corrBefore(subjNum)=corr(sAct1_nsAct1_nsN1_A(:,subjNum),sAct1_nsAct1_nsN1_B(:,subjNum));
    corrAfter(subjNum)=corr(sAct2_nsAct1_nsN1_A(:,subjNum),sAct2_nsAct1_nsN1_B(:,subjNum));
end
[h,pval]=ttest2(corrAfter, corrBefore);
disp(['Pearson correlation before: ' num2str(mean(corrBefore)) ', After: ' num2str(mean(corrAfter)) '; T-test p-value: ' num2str(pval)])
covBefore=zeros(numSubjs,1);
covAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    cB=cov(sAct1_nsAct1_nsN1_A(:,subjNum),sAct1_nsAct1_nsN1_B(:,subjNum));
    covBefore(subjNum)=cB(1,2);
    cA=cov(sAct2_nsAct1_nsN1_A(:,subjNum),sAct2_nsAct1_nsN1_B(:,subjNum));
    covAfter(subjNum)=cA(1,2);
end
[h,pval]=ttest2(covAfter, covBefore);
disp(['Covariance before: ' num2str(mean(covBefore)) ', After: ' num2str(mean(covAfter)) '; T-test p-value: ' num2str(pval)])
disp(' ')


disp('==Increased unshared activity amplitude (2x) in one region==')
sAct1_nsAct2_nsN1_A=sharedActivity+(2*nsActivityA)+nsNoiseA;
corrBefore=zeros(numSubjs,1);
corrAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    corrBefore(subjNum)=corr(sAct1_nsAct1_nsN1_A(:,subjNum),sAct1_nsAct1_nsN1_B(:,subjNum));
    corrAfter(subjNum)=corr(sAct1_nsAct2_nsN1_A(:,subjNum),sAct1_nsAct1_nsN1_B(:,subjNum));
end
[h,pval]=ttest2(corrAfter, corrBefore);
disp(['Pearson correlation before: ' num2str(mean(corrBefore)) ', After: ' num2str(mean(corrAfter)) '; T-test p-value: ' num2str(pval)])
covBefore=zeros(numSubjs,1);
covAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    cB=cov(sAct1_nsAct1_nsN1_A(:,subjNum),sAct1_nsAct1_nsN1_B(:,subjNum));
    covBefore(subjNum)=cB(1,2);
    cA=cov(sAct1_nsAct2_nsN1_A(:,subjNum),sAct1_nsAct1_nsN1_B(:,subjNum));
    covAfter(subjNum)=cA(1,2);
end
[h,pval]=ttest2(covAfter, covBefore);
disp(['Covariance before: ' num2str(mean(covBefore)) ', After: ' num2str(mean(covAfter)) '; T-test p-value: ' num2str(pval)])
disp(' ')

disp('==Increased non-shared activity amplitude (2x) in both regions==')
sAct1_nsAct2_nsN1_A=sharedActivity+(2*nsActivityA)+nsNoiseA;
sAct1_nsAct2_nsN1_B=sharedActivity+(2*nsActivityB)+nsNoiseB;
corrBefore=zeros(numSubjs,1);
corrAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    corrBefore(subjNum)=corr(sAct1_nsAct1_nsN1_A(:,subjNum),sAct1_nsAct1_nsN1_B(:,subjNum));
    corrAfter(subjNum)=corr(sAct1_nsAct2_nsN1_A(:,subjNum),sAct1_nsAct2_nsN1_B(:,subjNum));
end
[h,pval]=ttest2(corrAfter, corrBefore);
disp(['Pearson correlation before: ' num2str(mean(corrBefore)) ', After: ' num2str(mean(corrAfter)) '; T-test p-value: ' num2str(pval)])
covBefore=zeros(numSubjs,1);
covAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    cB=cov(sAct1_nsAct1_nsN1_A(:,subjNum),sAct1_nsAct1_nsN1_B(:,subjNum));
    covBefore(subjNum)=cB(1,2);
    cA=cov(sAct1_nsAct2_nsN1_A(:,subjNum),sAct1_nsAct2_nsN1_B(:,subjNum));
    covAfter(subjNum)=cA(1,2);
end
[h,pval]=ttest2(covAfter, covBefore);
disp(['Covariance before: ' num2str(mean(covBefore)) ', After: ' num2str(mean(covAfter)) '; T-test p-value: ' num2str(pval)])
% figure;scatter(sAct1_nsAct1_nsN1_A(:,1),sAct1_nsAct1_nsN1_B(:,1));
% hold on;scatter(sAct1_nsAct2_nsN1_A(:,1),sAct1_nsAct2_nsN1_B(:,1),4); title('Example: Increased non-shared activity amplitude (3x) in both regions'); hold off;
disp(' ')


disp('==Increased shared (2x) & non-shared (2x) activity amplitude in both regions==')
sAct2_nsAct2_nsN1_A=(2*sharedActivity)+(2*nsActivityA)+nsNoiseA;
sAct2_nsAct2_nsN1_B=(2*sharedActivity)+(2*nsActivityB)+nsNoiseB;
corrBefore=zeros(numSubjs,1);
corrAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    corrBefore(subjNum)=corr(sAct1_nsAct1_nsN1_A(:,subjNum),sAct1_nsAct1_nsN1_B(:,subjNum));
    corrAfter(subjNum)=corr(sAct2_nsAct2_nsN1_A(:,subjNum),sAct2_nsAct2_nsN1_B(:,subjNum));
end
[h,pval]=ttest2(corrAfter, corrBefore);
disp(['Pearson correlation before: ' num2str(mean(corrBefore)) ', After: ' num2str(mean(corrAfter)) '; T-test p-value: ' num2str(pval)])
covBefore=zeros(numSubjs,1);
covAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    cB=cov(sAct1_nsAct1_nsN1_A(:,subjNum),sAct1_nsAct1_nsN1_B(:,subjNum));
    covBefore(subjNum)=cB(1,2);
    cA=cov(sAct2_nsAct2_nsN1_A(:,subjNum),sAct2_nsAct2_nsN1_B(:,subjNum));
    covAfter(subjNum)=cA(1,2);
end
[h,pval]=ttest2(covAfter, covBefore);
disp(['Covariance before: ' num2str(mean(covBefore)) ', After: ' num2str(mean(covAfter)) '; T-test p-value: ' num2str(pval)])
disp(' ')


