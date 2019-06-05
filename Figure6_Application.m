clc
clear
close all

for shuffle=0:2
    for MonkeyNum=1:2
        [ccTheo2_cross_F{shuffle+1,MonkeyNum},ccSim2_cross_F{shuffle+1,MonkeyNum},ccTheo2_square_F{shuffle+1,MonkeyNum},ccSim2_square_F{shuffle+1,MonkeyNum},ccTheo1_F{shuffle+1,MonkeyNum},ccSim1_F{shuffle+1,MonkeyNum}] = ccComputeCombined (MonkeyNum,shuffle);
    end
end


%% Filtered fitting plot
figure
jj=0;
for MonkeyNum=1:2 
    for shuffle=0:2
        jj=jj+1;
        Pointsize=5;
        subplot(2,3,jj)
        plot(ccTheo2_cross_F{shuffle+1,MonkeyNum},ccSim2_cross_F{shuffle+1,MonkeyNum},'r.','markersize', 1);hold on;
        plot(ccTheo1_F{shuffle+1,MonkeyNum},ccSim1_F{shuffle+1,MonkeyNum},'b.','markersize', 1);hold on;
        plot(ccTheo2_square_F{shuffle+1,MonkeyNum},ccSim2_square_F{shuffle+1,MonkeyNum},'g.','markersize', 1);hold on;
        [slope_fit1_F{shuffle+1,MonkeyNum}] = Fitting_cc (ccTheo1_F{shuffle+1,MonkeyNum},ccSim1_F{shuffle+1,MonkeyNum});
        [slope_fit2_cross_F{shuffle+1,MonkeyNum}] = Fitting_cc (ccTheo2_cross_F{shuffle+1,MonkeyNum},ccSim2_cross_F{shuffle+1,MonkeyNum});
        [slope_fit2_square_F{shuffle+1,MonkeyNum}] = Fitting_cc (ccTheo2_square_F{shuffle+1,MonkeyNum},ccSim2_square_F{shuffle+1,MonkeyNum});
        plot([-1,1],slope_fit1_F{shuffle+1,MonkeyNum}.*[-1,1],'b-','markersize', 6);hold on;
        plot([-1,1],slope_fit2_cross_F{shuffle+1,MonkeyNum}.*[-1,1],'r-','markersize', 6);hold on;
        plot([-1,1],slope_fit2_square_F{shuffle+1,MonkeyNum}.*[-1,1],'g-','markersize', 6);hold on;
        plot([-1,1],[-1,1] ,'k-','markersize', 6);hold on;
        axis square;
        set(gca,'XTick',[-1,1]);
        set(gca,'YTick',[-1,1]);
        set(gca,'Yticklabel',[-1,1]);
        set(gca,'Xticklabel',[-1,1]);
        axis([-1 1 -1 1]);
        if jj==1
            title('Original','FontSize',10);
            legend('cross','linear','square');
        elseif jj==2
            title('Shuffle internal noise','FontSize',10);
        elseif jj==3
            title('Shuffle external noise','FontSize',10);    
        elseif jj==4
            xlabel('Optimal NACCC','FontSize',10);
            ylabel('Measured NACCC','FontSize',10);
        end
    end
end




function [ccTheo2_cross_F,ccSim2_cross_F,ccTheo2_square_F,ccSim2_square_F,ccTheo1_F,ccSim1_F] = ccComputeCombined (MonkeyNum,shuffle)
if MonkeyNum==1
    SessionTotalNum=59;
    load('monkey1.mat');
else
    SessionTotalNum=71;
    load('monkey2.mat');
end
MonkeySessionList=[1:SessionTotalNum];
contrastAll=0; % flag to denote if we want to pick the all the contrasts or just the highest contrast.
%% Compute Theo and experimental CCs for each session
for SelectedSession=1:SessionTotalNum   
    if MonkeyNum==1
        SessionData=monkey1{SelectedSession};
    else
        SessionData=monkey2{SelectedSession};
    end
    %% Find same contrast  
    Contrast=SessionData.contrast;
    Union_Contrast=unique(Contrast);
    if contrastAll==1
        Index_SelectContrast=Contrast==Contrast;
    else
        Index_SelectContrast=Contrast==Union_Contrast(end);
    end
    %% Pick data under same contrast.
    Stimulus=double(SessionData.stimulus_class=='B');
    Choice=double(SessionData.selected_class=='B');
    Orientation=double(SessionData.orientation);
    Orientation=Orientation-mean(Orientation);
    response=double(SessionData.Counts_matrix);
    n=size(response,2);
    Stimulus1=Stimulus(Index_SelectContrast);
    Choice1=Choice(Index_SelectContrast);    
    Stimulus1=Stimulus1.*2-1;
    Choice1=Choice1.*2-1;
    Orientation1=Orientation(Index_SelectContrast);
    response1=response(Index_SelectContrast,:);
    NumOfNA_r(SelectedSession)=0;
    for ij=1:n
        if sum(response1(:,ij)==0)==n
            NumOfNA_r(SelectedSession)=ij;
        end
    end
    m=size(response1,1);
    Orientation_p=Orientation1(Stimulus1==1);
    Orientation_n=Orientation1(Stimulus1==-1);
    sigma_p=15;sigma_n=3;
    sp=sigma_p^2;
    sn=sigma_n^2;
    thetahat_DB=sqrt((sp+sn)/2);
     sbar=thetahat_DB^2;
     clearvars sT
     for j=1:m
        if Stimulus1(j)==1
            sT(j)=sp;
        else
            sT(j)=sn;
        end
     end
     
[ccTheo,ccSim,CCPredRatio(SelectedSession),zeta(SelectedSession),delta(SelectedSession),n_half_fit(SelectedSession),us_S1(SelectedSession),us_S0(SelectedSession),OrientationMidMeaningful{SelectedSession},Choice_A_Ratio{SelectedSession},predRatio{SelectedSession}] = v1_ccAnalysis (Stimulus1,Choice1,Orientation1,response1,shuffle,sT,sbar,sp,sn);
sp_all(SelectedSession)=sp;
sn_all(SelectedSession)=sn;
ccTheo1(SelectedSession,:)=ccTheo(1:n);
ccTheo2_square(SelectedSession,:)=ccTheo(n+1:2*n);
ccTheo2_cross(SelectedSession,:)=ccTheo(2*n+1:end);
ccSim1(SelectedSession,:)=ccSim(1:n);
ccSim2_square(SelectedSession,:)=ccSim(n+1:2*n);
ccSim2_cross(SelectedSession,:)=ccSim(2*n+1:end);

%% Compute Accuracy
thetahat_0=sqrt(2*log(sigma_p/sigma_n)/(1/sigma_n^2-1/sigma_p^2));
Choice_opt=double((Orientation1-0).^2>thetahat_0^2);
Choice_opt=Choice_opt.*2-1;
Accuracy(SelectedSession)=sum(Choice1==Stimulus1)./numel(Choice1);
Accuracy_opt(SelectedSession)=0.5+0.25*I0compute(thetahat_0,sqrt(sp))-0.25*I0compute(thetahat_0,sqrt(sn));
end

%% linear fit for every session
for SelectedSession=1:SessionTotalNum     
    [slope_fit1(SelectedSession)] = Fitting_cc (ccTheo1(SelectedSession,:),ccSim1(SelectedSession,:))
    [slope_fit2_cross(SelectedSession)] = Fitting_cc (ccTheo2_cross(SelectedSession,:),ccSim2_cross(SelectedSession,:))
    [slope_fit2_square(SelectedSession)] = Fitting_cc (ccTheo2_square(SelectedSession,:),ccSim2_square(SelectedSession,:))
end

%% Filtering sessions according to Accuracy
if MonkeyNum==1
    Accu_thre=0.7
else
    Accu_thre=0.75
end
MonkeySessionList_Filtered=MonkeySessionList(Accuracy>=Accu_thre);
ccSim2_cross_F=reshape(ccSim2_cross(MonkeySessionList_Filtered,:),1,[]);
ccTheo2_cross_F=reshape(ccTheo2_cross(MonkeySessionList_Filtered,:),1,[]);
ccSim2_square_F=reshape(ccSim2_square(MonkeySessionList_Filtered,:),1,[]);
ccTheo2_square_F=reshape(ccTheo2_square(MonkeySessionList_Filtered,:),1,[]);
ccSim1_F=reshape(ccSim1(MonkeySessionList_Filtered,:),1,[]);
ccTheo1_F=reshape(ccTheo1(MonkeySessionList_Filtered,:),1,[]);

end






%% Defined function
function [ccTheo,ccSim,CCpredRatio,zeta,delta,n_half_pred,us_S1,us_S0,OrientationMidMeaningful,Choice_A_Ratio,predRatio] = v1_ccAnalysis (Stimulus1,Choice1,Orientation1,response1,shuffle,sT,sbar,sp,sn)
Orientation_Sp=Orientation1(Stimulus1==1);
Orientation_Sn=Orientation1(Stimulus1==-1);
Choice_Sp=Choice1(Stimulus1==1);
Choice_Sn=Choice1(Stimulus1==-1);
Stim_Sp=Stimulus1(Stimulus1==1);
Stim_Sn=Stimulus1(Stimulus1==-1);
sT_Sp=sT(Stimulus1==1);
sT_Sn=sT(Stimulus1==-1);
response1_Sp= response1(Stimulus1==1,:);
response1_Sn= response1(Stimulus1==-1,:);

%% Combine
Stimulus1=[Stim_Sp;Stim_Sn];
Choice1=[Choice_Sp;Choice_Sn];
Orientation1=[Orientation_Sp,Orientation_Sn];
sT=[sT_Sp,sT_Sn];
response1=[response1_Sp;response1_Sn];
n=size(response1,2);
m=size(response1,1);

%% PsychometricThreshold and estimate CCsim/CCTheo
BoundarySize=3;
OrientationMid=-100:BoundarySize:100;
MinTrialNum=3;
k=0;
AllChoice=double.empty;
AllFixedRb=double.empty;
AllStimulus=double.empty;
for j=1:numel(OrientationMid)
    SelectedIndex=abs(Orientation1-OrientationMid(j))<BoundarySize/2;
    if sum(double(SelectedIndex))>=MinTrialNum
        k=k+1;
        OrientationMidMeaningful(k)=OrientationMid(j);
        Choice_Range=Choice1(SelectedIndex);
        Choice_A_Ratio(k)=sum(Choice_Range(Choice_Range==1))./numel(Choice_Range);
    end
end
theta_fit = glmfit(abs(OrientationMidMeaningful'), [Choice_A_Ratio' ones(numel(Choice_A_Ratio),1)], 'binomial', 'link', 'logit');
tmp=theta_fit(1)+theta_fit(2)*abs(OrientationMidMeaningful);        %theta'x
predRatio=1./(1+exp(-tmp));  %logistic function
us_S1=(sum(Choice_Sp))./numel(Choice_Sp);
us_S0=(sum(Choice_Sn))./numel(Choice_Sn);
n_half_pred=(-log(1)-theta_fit(1))/theta_fit(2);
us=(us_S1^2+us_S0^2)/2;
zeta=2/sqrt(pi)*n_half_pred/sqrt(sp)*1/sqrt(1-us);
delta=sqrt(5/4)/(1-sn/sp)*sqrt((us_S1-us_S0)^2/(1-0.25*(us_S1+us_S0)^2));
CCpredRatio=zeta./delta;
%% Shuffle
 Index_Sp_Cp=find(Choice_Sp==1);
 Index_Sp_Cn=find(Choice_Sp==-1);
 Index_Sn_Cp=find(Choice_Sn==1);
 Index_Sn_Cn=find(Choice_Sn==-1);

 Orientation_Sp_Cp=Orientation_Sp(Index_Sp_Cp);
 Orientation_Sp_Cn=Orientation_Sp(Index_Sp_Cn);
 Orientation_Sn_Cp=Orientation_Sn(Index_Sn_Cp);
 Orientation_Sn_Cn=Orientation_Sn(Index_Sn_Cn);

 response1_Sp_Cp=response1_Sp(Index_Sp_Cp,:);
 response1_Sp_Cn=response1_Sp(Index_Sp_Cn,:);
 response1_Sn_Cp=response1_Sn(Index_Sn_Cp,:);
 response1_Sn_Cn=response1_Sn(Index_Sn_Cn,:);

 if shuffle==1 % shuffle r given orientation, s and choice
       
     [response_shuff_Sp_Cp] = Shuffle_func (Orientation_Sp_Cp,response1_Sp_Cp,shuffle);
     [response_shuff_Sp_Cn] = Shuffle_func (Orientation_Sp_Cn,response1_Sp_Cn,shuffle);
     [response_shuff_Sn_Cp] = Shuffle_func (Orientation_Sn_Cp,response1_Sn_Cp,shuffle);
     [response_shuff_Sn_Cn] = Shuffle_func (Orientation_Sn_Cn,response1_Sn_Cn,shuffle);   
     response1_Sp_shuff=response1_Sp;
     response1_Sp_shuff(Index_Sp_Cp,:)=response_shuff_Sp_Cp;
     response1_Sp_shuff(Index_Sp_Cn,:)=response_shuff_Sp_Cn;
     response1_Sn_shuff=response1_Sn;
     response1_Sn_shuff(Index_Sn_Cp,:)=response_shuff_Sn_Cp;
     response1_Sn_shuff(Index_Sn_Cn,:)=response_shuff_Sn_Cn;
    response1_shuff=[response1_Sp_shuff;response1_Sn_shuff];    
    [R2_square,R2_cross] = rMomentsGenerate_v1_exp (response1,response1_shuff);
    [R1_square_Sp,R1_cross_Sp] = rMomentsGenerate_v1_exp (response1_Sp,response1_Sp_shuff);
    [R1_square_Sn,R1_cross_Sn] = rMomentsGenerate_v1_exp (response1_Sn,response1_Sn_shuff);  
 
 elseif shuffle==2 % shuffle r given s and choice.
     [response_shuff_Sp_Cp] = Shuffle_func (Orientation_Sp_Cp,response1_Sp_Cp,shuffle);
     [response_shuff_Sp_Cn] = Shuffle_func (Orientation_Sp_Cn,response1_Sp_Cn,shuffle);
     [response_shuff_Sn_Cp] = Shuffle_func (Orientation_Sn_Cp,response1_Sn_Cp,shuffle);
     [response_shuff_Sn_Cn] = Shuffle_func (Orientation_Sn_Cn,response1_Sn_Cn,shuffle);   
     response1_Sp_shuff=response1_Sp;
     response1_Sp_shuff(Index_Sp_Cp,:)=response_shuff_Sp_Cp;
     response1_Sp_shuff(Index_Sp_Cn,:)=response_shuff_Sp_Cn;
     response1_Sn_shuff=response1_Sn;
     response1_Sn_shuff(Index_Sn_Cp,:)=response_shuff_Sn_Cp;
     response1_Sn_shuff(Index_Sn_Cn,:)=response_shuff_Sn_Cn;
    response1_shuff=[response1_Sp_shuff;response1_Sn_shuff];    
    [R2_square,R2_cross] = rMomentsGenerate_v1_exp (response1,response1_shuff);
    [R1_square_Sp,R1_cross_Sp] = rMomentsGenerate_v1_exp (response1_Sp,response1_Sp_shuff);
    [R1_square_Sn,R1_cross_Sn] = rMomentsGenerate_v1_exp (response1_Sn,response1_Sn_shuff);
 
 elseif shuffle==0 % No shuffling
    [R2_square,R2_cross] = rMomentsGenerate_v1_exp (response1,response1);
    [R1_square_Sp,R1_cross_Sp] = rMomentsGenerate_v1_exp (response1_Sp,response1_Sp);
    [R1_square_Sn,R1_cross_Sn] = rMomentsGenerate_v1_exp (response1_Sn,response1_Sn);
 end

%% Compute cctheo and ccsim
Corr_s_R=corr(sT',[response1,R2_square,R2_cross]);
std_shat_R_sDB_approx=sqrt(2)*sbar;
Corr_s_shatR=corr(sT',Choice1);
ccTheo=CCpredRatio.*Corr_s_R./Corr_s_shatR;

%% Compute coarse-discrimination CC
trialNumRatio_Sp=numel(Choice_Sp)./m;
trialNumRatio_Sn=numel(Choice_Sn)./m;
R_all_Sp=[response1_Sp,R1_square_Sp,R1_cross_Sp];
R_all_Sn=[response1_Sn,R1_square_Sn,R1_cross_Sn];
varShat_ave=trialNumRatio_Sp.*var(Choice_Sp)+trialNumRatio_Sn.*var(Choice_Sn);
for j=1:size(R_all_Sn,2)
    cov_mat_Sp=cov(Choice_Sp,R_all_Sp(:,j));
    cov_mat_Sn=cov(Choice_Sn,R_all_Sn(:,j));
    cov_ave=trialNumRatio_Sp.*cov_mat_Sp(1,2)+trialNumRatio_Sn.*cov_mat_Sn(1,2);
    varR_ave=trialNumRatio_Sp.*var(R_all_Sp(:,j))+trialNumRatio_Sn.*var(R_all_Sn(:,j));    
    ccSim(j)=(cov_ave)./sqrt(varR_ave.*varShat_ave);
    if isnan(ccSim(j))
        ccSim(j)=0;
    end
    if isnan(ccTheo(j))
        ccTheo(j)=0;
    end
end
end


function [R2_square,R2_cross] = rMomentsGenerate_v1_exp (r1,r2)
m = size(r1,1);
n = size(r1,2);
rmean1=repmat(mean(r1,1),m,1);
rmean2=repmat(mean(r2,1),m,1);
z1=r1-rmean1;
z2=r2-rmean2;
k=0;
for i=1:n
        k=k+1;
        R2_square(:,k)=z1(:,i).*z2(:,i);
end
k=0;
for i=1:n
    for j=i+1:n
        k=k+1;
        R2_cross(:,k)=z1(:,i).*z2(:,j);
    end
end
end


function I=I0compute(theta,sigma)
    I=2*erfc(theta/sqrt(2)/sigma)-1;
end

function [slopePCA] = Fitting_cc (ccTheo,ccSim)
    ccdat=[ccTheo;ccSim]';
    coeff = pca(ccdat);
    slopePCA = abs(coeff(1,2) / coeff(1,1));
end
        
function [response_shuff] = Shuffle_func (Orientation,response,shuffle)
    if shuffle==1
        [~,index]=sort(Orientation);
    elseif shuffle==2
         shufflingtimes=100000;
         index=randi([1 numel(Orientation)],1,shufflingtimes);
    end
    response_shuff=response;
    for i=1:numel(Orientation)/2
        temp=response_shuff(index(2*i),:);
        response_shuff(index(2*i),:)=response_shuff(index(2*i-1),:);
        response_shuff(index(2*i-1),:)=temp;
    end
end