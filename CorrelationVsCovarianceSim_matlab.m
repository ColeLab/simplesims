%MATLAB version of the simple simulations
%
%Note: Rather than amplifying each time series component, new time series
%are generated and added to the baseline time series. This is somewhat more
%realistic than amplifying the baseline time series components. The results
%are the same using this approach.


disp('====Simulating data for 2 brain regions from 25 subjects by creating random normally distributed time series for each===')
disp('==Each region''s time series consists of: Shared (across regions) activity, unshared activity, unshared noise (MR & physiological)==')
disp(' ')


%% The basic components of the simulations

%25 subjects modeled, 200 time points each
numSubjs=25;
numTimepoints=200;

%%Baselines
%The baseline level of communication between the regions
shared_baseline=randn(numTimepoints,numSubjs);
%The baseline level of internal activity within region X
unsharedX_baseline=randn(numTimepoints,numSubjs);
%The baseline level of internal activity within region Y
unsharedY_baseline=randn(numTimepoints,numSubjs);
%The noise in X
noiseX=.25*randn(numTimepoints,numSubjs);
%The noise in Y
noiseY=.25*randn(numTimepoints,numSubjs);

%%Increases
%The increased signal that is communicated between X and Y
shared_inc=randn(numTimepoints,numSubjs);
%The increased amount of internal activity within region X
unsharedX_inc=randn(numTimepoints,numSubjs);
%The increased amount of internal activity within region Y
unsharedY_inc=randn(numTimepoints,numSubjs);

%% Constructing baseline time series
X_baseline = shared_baseline + unsharedX_baseline + noiseX;
Y_baseline = shared_baseline + unsharedY_baseline + noiseY;
disp('--Constructing baseline time series--')
disp('X_baseline = shared_baseline + unsharedX_baseline + noiseX')
disp('Y_baseline = shared_baseline + unsharedY_baseline + noiseY')
disp(' ')

%% Increased shared activity amplitude (2x) in both regions
disp('==Increased shared activity amplitude (2x) in both regions==')

%Increasing coupling between X and Y, resulting in increased communication in both directions
X_sim = X_baseline + shared_inc;
Y_sim = Y_baseline + shared_inc;
disp('X_sim = X_baseline + shared_inc')
disp('Y_sim = Y_baseline + shared_inc')

%Group level statistical tests

%Pearson correlation
corrBefore=zeros(numSubjs,1);
corrAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    corrBefore(subjNum)=corr(X_baseline(:,subjNum),Y_baseline(:,subjNum));
    corrAfter(subjNum)=corr(X_sim(:,subjNum),Y_sim(:,subjNum));
end
[h,pval]=ttest2(corrAfter, corrBefore);
disp(['Pearson correlation before: ' num2str(mean(corrBefore)) ', After: ' num2str(mean(corrAfter)) '; T-test p-value: ' num2str(pval)])

%Covariance
covBefore=zeros(numSubjs,1);
covAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    cB=cov(X_baseline(:,subjNum),Y_baseline(:,subjNum));
    covBefore(subjNum)=cB(1,2);
    cA=cov(X_sim(:,subjNum),Y_sim(:,subjNum));
    covAfter(subjNum)=cA(1,2);
end
[h,pval]=ttest2(covAfter, covBefore);
disp(['Covariance before: ' num2str(mean(covBefore)) ', After: ' num2str(mean(covAfter)) '; T-test p-value: ' num2str(pval)])
disp(' ')



%% Increased unshared activity amplitude (2x) in one region
disp('==Increased unshared activity amplitude (2x) in one region==')

%Increasing coupling between X and Y, resulting in increased communication in both directions
X_sim = X_baseline + unsharedX_inc;
Y_sim = Y_baseline;
disp('X_sim = X_baseline + unsharedX_inc')
disp('Y_sim = Y_baseline')

%Group level statistical tests

%Pearson correlation
corrBefore=zeros(numSubjs,1);
corrAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    corrBefore(subjNum)=corr(X_baseline(:,subjNum),Y_baseline(:,subjNum));
    corrAfter(subjNum)=corr(X_sim(:,subjNum),Y_sim(:,subjNum));
end
[h,pval]=ttest2(corrAfter, corrBefore);
disp(['Pearson correlation before: ' num2str(mean(corrBefore)) ', After: ' num2str(mean(corrAfter)) '; T-test p-value: ' num2str(pval)])

%Covariance
covBefore=zeros(numSubjs,1);
covAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    cB=cov(X_baseline(:,subjNum),Y_baseline(:,subjNum));
    covBefore(subjNum)=cB(1,2);
    cA=cov(X_sim(:,subjNum),Y_sim(:,subjNum));
    covAfter(subjNum)=cA(1,2);
end
[h,pval]=ttest2(covAfter, covBefore);
disp(['Covariance before: ' num2str(mean(covBefore)) ', After: ' num2str(mean(covAfter)) '; T-test p-value: ' num2str(pval)])
disp(' ')




%% Increased unshared activity amplitude (2x) in both regions
disp('==Increased non-shared activity amplitude (2x) in both regions==')

%Increasing coupling between X and Y, resulting in increased communication in both directions
X_sim = X_baseline + unsharedX_inc;
Y_sim = Y_baseline + unsharedY_inc;
disp('X_sim = X_baseline + unsharedX_inc')
disp('Y_sim = Y_baseline + unsharedY_inc')

%Group level statistical tests

%Pearson correlation
corrBefore=zeros(numSubjs,1);
corrAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    corrBefore(subjNum)=corr(X_baseline(:,subjNum),Y_baseline(:,subjNum));
    corrAfter(subjNum)=corr(X_sim(:,subjNum),Y_sim(:,subjNum));
end
[h,pval]=ttest2(corrAfter, corrBefore);
disp(['Pearson correlation before: ' num2str(mean(corrBefore)) ', After: ' num2str(mean(corrAfter)) '; T-test p-value: ' num2str(pval)])

%Covariance
covBefore=zeros(numSubjs,1);
covAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    cB=cov(X_baseline(:,subjNum),Y_baseline(:,subjNum));
    covBefore(subjNum)=cB(1,2);
    cA=cov(X_sim(:,subjNum),Y_sim(:,subjNum));
    covAfter(subjNum)=cA(1,2);
end
[h,pval]=ttest2(covAfter, covBefore);
disp(['Covariance before: ' num2str(mean(covBefore)) ', After: ' num2str(mean(covAfter)) '; T-test p-value: ' num2str(pval)])
disp(' ')




%% Increased shared (2x) & unshared (2x) activity amplitude in both regions
disp('==Increased shared (2x) & non-shared (2x) activity amplitude in both regions==')

%Increasing coupling between X and Y, resulting in increased communication in both directions
X_sim = X_baseline + unsharedX_inc + shared_inc;
Y_sim = Y_baseline + unsharedY_inc + shared_inc;
disp('X_sim = X_baseline + unsharedX_inc + shared_inc')
disp('Y_sim = Y_baseline + unsharedY_inc + shared_inc')

%Group level statistical tests

%Pearson correlation
corrBefore=zeros(numSubjs,1);
corrAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    corrBefore(subjNum)=corr(X_baseline(:,subjNum),Y_baseline(:,subjNum));
    corrAfter(subjNum)=corr(X_sim(:,subjNum),Y_sim(:,subjNum));
end
[h,pval]=ttest2(corrAfter, corrBefore);
disp(['Pearson correlation before: ' num2str(mean(corrBefore)) ', After: ' num2str(mean(corrAfter)) '; T-test p-value: ' num2str(pval)])

%Covariance
covBefore=zeros(numSubjs,1);
covAfter=zeros(numSubjs,1);
for subjNum=1:numSubjs
    cB=cov(X_baseline(:,subjNum),Y_baseline(:,subjNum));
    covBefore(subjNum)=cB(1,2);
    cA=cov(X_sim(:,subjNum),Y_sim(:,subjNum));
    covAfter(subjNum)=cA(1,2);
end
[h,pval]=ttest2(covAfter, covBefore);
disp(['Covariance before: ' num2str(mean(covBefore)) ', After: ' num2str(mean(covAfter)) '; T-test p-value: ' num2str(pval)])
disp(' ')



