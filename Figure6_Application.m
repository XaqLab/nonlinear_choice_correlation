clc
clear
close all

for shuffle=0:2
    for MonkeyNum=1:2
        [ccTheo2_cross_F{shuffle+1,MonkeyNum},ccSim2_cross_F{shuffle+1,MonkeyNum},ccTheo2_square_F{shuffle+1,MonkeyNum},ccSim2_square_F{shuffle+1,MonkeyNum},ccTheo1_F{shuffle+1,MonkeyNum},ccSim1_F{shuffle+1,MonkeyNum},pSignificant_Theo1_F{shuffle+1,MonkeyNum},pSignificant_Theo2sq_F{shuffle+1,MonkeyNum},pSignificant_Theo2cr_F{shuffle+1,MonkeyNum},pSignificant_Sim1_F{shuffle+1,MonkeyNum},pSignificant_Sim2sq_F{shuffle+1,MonkeyNum},pSignificant_Sim2cr_F{shuffle+1,MonkeyNum}] = ccComputeCombined (MonkeyNum,shuffle);
        pSignificant_1_Comb{shuffle+1,MonkeyNum}=pValueSepToComb(pSignificant_Theo1_F{shuffle+1,MonkeyNum},pSignificant_Sim1_F{shuffle+1,MonkeyNum});
        pSignificant_sq_Comb{shuffle+1,MonkeyNum}=pValueSepToComb(pSignificant_Theo2sq_F{shuffle+1,MonkeyNum},pSignificant_Sim2sq_F{shuffle+1,MonkeyNum});
        pSignificant_cross_Comb{shuffle+1,MonkeyNum}=pValueSepToComb(pSignificant_Theo2cr_F{shuffle+1,MonkeyNum},pSignificant_Sim2cr_F{shuffle+1,MonkeyNum});
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

%%

fname = sprintf('CC_test_AllData.mat');
save(fname,'pSignificant_1_Comb','pSignificant_sq_Comb','pSignificant_cross_Comb','pSignificant_Theo1_F','pSignificant_Theo2sq_F','pSignificant_Theo2cr_F','pSignificant_Sim1_F','pSignificant_Sim2sq_F','pSignificant_Sim2cr_F','ccSim1_F','ccSim2_cross_F','ccSim2_square_F','ccTheo1_F','ccTheo2_cross_F','ccTheo2_square_F','slope_fit1_F','slope_fit2_cross_F','slope_fit2_square_F');

save('CC_test_PlotsAndData.mat');


function [ccTheo2_cross_F,ccSim2_cross_F,ccTheo2_square_F,ccSim2_square_F,ccTheo1_F,ccSim1_F,pSignificant_Theo1_F,pSignificant_Theo2sq_F,pSignificant_Theo2cr_F,pSignificant_Sim1_F,pSignificant_Sim2sq_F,pSignificant_Sim2cr_F] = ccComputeCombined (MonkeyNum,shuffle)
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
    SelectedSession
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
%%
[ccTheo,ccSim,CCPredRatio(SelectedSession),zeta(SelectedSession),delta(SelectedSession),n_half_fit(SelectedSession),us_S1(SelectedSession),us_S0(SelectedSession),OrientationMidMeaningful{SelectedSession},Choice_A_Ratio{SelectedSession},predRatio{SelectedSession}] = v1_ccAnalysis (Stimulus1,Choice1,Orientation1,response1,shuffle,sT,sbar,sp,sn);
sp_all(SelectedSession)=sp;
sn_all(SelectedSession)=sn;
ccTheo1(SelectedSession,:)=ccTheo(1:n);
ccTheo2_square(SelectedSession,:)=ccTheo(n+1:2*n);
ccTheo2_cross(SelectedSession,:)=ccTheo(2*n+1:end);
ccSim1(SelectedSession,:)=ccSim(1:n);
ccSim2_square(SelectedSession,:)=ccSim(n+1:2*n);
ccSim2_cross(SelectedSession,:)=ccSim(2*n+1:end);

%%     
[ccTheo1BS,ccTheo2_squareBS,ccTheo2_crossBS,ccSim1BS,ccSim2_squareBS,ccSim2_crossBS] = ccSignificance (Stimulus1,Choice1,Orientation1,response1,shuffle,sT,sp,sn);

% figure
% histogram(ccTheo1BS(:,10),'Normalization','probability','BinWidth',0.03); hold on;
% histogram(ccTheo2_squareBS(:,10),'Normalization','probability','BinWidth',0.03); hold on;
% histogram(ccTheo2_crossBS(:,10),'Normalization','probability','BinWidth',0.03); 
% legend('linear','square','cross');
% title('CC theo')
% xlabel('value')
% ylabel('probability')
% 
% figure
% histogram(ccSim1BS(:,10),'Normalization','probability','BinWidth',0.03); hold on;
% histogram(ccSim2_squareBS(:,10),'Normalization','probability','BinWidth',0.03); hold on;
% histogram(ccSim2_crossBS(:,10),'Normalization','probability','BinWidth',0.03); 
% legend('linear','square','cross');
% title('CC Sim')
% xlabel('value')
% ylabel('probability')


mu_Theo1(SelectedSession,:)=mean(ccTheo1BS);
std_Theo_1(SelectedSession,:)=std(ccTheo1BS);
mu_Theo2square(SelectedSession,:)=mean(ccTheo2_squareBS);
std_Theo2square(SelectedSession,:)=std(ccTheo2_squareBS);
mu_Theo2cross(SelectedSession,:)=mean(ccTheo2_crossBS);
std_Theo2cross(SelectedSession,:)=std(ccTheo2_crossBS);
mu_Sim1(SelectedSession,:)=mean(ccSim1BS);
std_Sim_1(SelectedSession,:)=std(ccSim1BS);
mu_Sim2square(SelectedSession,:)=mean(ccSim2_squareBS);
std_Sim2square(SelectedSession,:)=std(ccSim2_squareBS);
mu_Sim2cross(SelectedSession,:)=mean(ccSim2_crossBS);
std_Sim2cross(SelectedSession,:)=std(ccSim2_crossBS);

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

%% Filtered CCs and Significance
ccSim2_cross_F=reshape(ccSim2_cross(MonkeySessionList_Filtered,:),1,[]);
ccTheo2_cross_F=reshape(ccTheo2_cross(MonkeySessionList_Filtered,:),1,[]);
ccSim2_square_F=reshape(ccSim2_square(MonkeySessionList_Filtered,:),1,[]);
ccTheo2_square_F=reshape(ccTheo2_square(MonkeySessionList_Filtered,:),1,[]);
ccSim1_F=reshape(ccSim1(MonkeySessionList_Filtered,:),1,[]);
ccTheo1_F=reshape(ccTheo1(MonkeySessionList_Filtered,:),1,[]);

%% 
for jj=1:numel(MonkeySessionList_Filtered)
    FilteredSessionNum=MonkeySessionList_Filtered(jj);
    pSignificant_Theo1(jj,:) = 1-normcdf(ccTheo1(FilteredSessionNum,:),mu_Theo1(FilteredSessionNum,:),std_Theo_1(FilteredSessionNum,:));
    pSignificant_Theo2sq(jj,:)  = 1-normcdf(ccTheo2_square(FilteredSessionNum,:),mu_Theo2square(FilteredSessionNum,:),std_Theo2square(FilteredSessionNum,:));
    pSignificant_Theo2cr(jj,:)  = 1-normcdf(ccTheo2_cross(FilteredSessionNum,:),mu_Theo2cross(FilteredSessionNum,:),std_Theo2cross(FilteredSessionNum,:));
    pSignificant_Sim1(jj,:)  = 1-normcdf(ccSim1(FilteredSessionNum,:),mu_Sim1(FilteredSessionNum,:),std_Sim_1(FilteredSessionNum,:));
    pSignificant_Sim2sq(jj,:)  = 1-normcdf(ccSim2_square(FilteredSessionNum,:),mu_Sim2square(FilteredSessionNum,:),std_Sim2square(FilteredSessionNum,:));
    pSignificant_Sim2cr(jj,:)  = 1-normcdf(ccSim2_cross(FilteredSessionNum,:),mu_Sim2cross(FilteredSessionNum,:),std_Sim2cross(FilteredSessionNum,:));
end

pSignificant_Theo1_F=reshape(pSignificant_Theo1,1,[]);
pSignificant_Theo2sq_F=reshape(pSignificant_Theo2sq,1,[]);
pSignificant_Theo2cr_F=reshape(pSignificant_Theo2cr,1,[]);
pSignificant_Sim1_F=reshape(pSignificant_Sim1,1,[]);
pSignificant_Sim2sq_F=reshape(pSignificant_Sim2sq,1,[]);
pSignificant_Sim2cr_F=reshape(pSignificant_Sim2cr,1,[]);


end


function [ccTheo1BS,ccTheo2_squareBS,ccTheo2_crossBS,ccSim1BS,ccSim2_squareBS,ccSim2_crossBS] = ccSignificance (Stimulus1,Choice1,Orientation1,response1,shuffle,sT,sp,sn)
     %% Compute the null distribution
     n=size(response1,2);
     Bootstrap_Times=100;
    for kk=1 : Bootstrap_Times
        % Compute ccTheo shuffled to estimate the shuffled distribution
        % parameters
        Sample_Index_1 = randi([1,numel(Stimulus1)],1,numel(Stimulus1));% randomly chose from data with replacement. 
        sT_shuffled=sT(Sample_Index_1);
        [ccTheoBS(1,kk,:)] = v1_ccTheo_Compute (Stimulus1,Choice1,Orientation1,response1,shuffle,sT_shuffled,sT,sp,sn);                  
        ccTheo1BS(kk,:)=ccTheoBS(1,kk,1:n);
        ccTheo2_squareBS(kk,:)=ccTheoBS(1,kk,n+1:2*n);
        ccTheo2_crossBS(kk,:)=ccTheoBS(1,kk,2*n+1:end);              
        % Compute ccSim shuffled to estimate the shuffled distribution
        % parameters
        Sample_Index_2 = randi([1,numel(Stimulus1)],1,numel(Stimulus1));    
        Choice1_shuffled=Choice1(Sample_Index_2);
        [ccSimBS(1,kk,:)] = v1_ccSim_Compute (Stimulus1,Choice1_shuffled,Orientation1,response1,shuffle);            
        ccSim1BS(kk,:)=ccSimBS(1,kk,1:n);
        ccSim2_squareBS(kk,:)=ccSimBS(1,kk,n+1:2*n);
        ccSim2_crossBS(kk,:)=ccSimBS(1,kk,2*n+1:end);   
    end
%     figure
%     for jj=1:96
%         histogram(ccTheo2_crossBS(:,jj)); hold on;
%     end
% 
%     figure
%     for jj=1:96
%         histogram(ccSim1BS(:,jj)); hold on;
%     end

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
 if shuffle==1 % shuffle r given orientation and s
    [~,indexp]=sort(Orientation_Sp);
    [~,indexn]=sort(Orientation_Sn);
    response1_Sn_shuff=response1_Sn;
    response1_Sp_shuff=response1_Sp;
    for i=1:numel(Orientation_Sp)/2
        temp=response1_Sp_shuff(indexp(2*i),:);
        response1_Sp_shuff(indexp(2*i),:)=response1_Sp_shuff(indexp(2*i-1),:);
        response1_Sp_shuff(indexp(2*i-1),:)=temp;
    end
    
    for i=1:numel(Orientation_Sn)/2
        temp=response1_Sn_shuff(indexn(2*i),:);
        response1_Sn_shuff(indexn(2*i),:)=response1_Sn_shuff(indexn(2*i-1),:);
        response1_Sn_shuff(indexn(2*i-1),:)=temp;
    end
    response1_shuff=[response1_Sp_shuff;response1_Sn_shuff];    
    [rmean_Shatn_shuff,rmean_Shatp_shuff] = R_Centralization_decoded(response1_Sn_shuff, response1_Sp_shuff, Stimulus1);
    [rmean_Shatn,rmean_Shatp] = R_Centralization_decoded(response1_Sn, response1_Sp, Stimulus1);
    [R2_square,R2_cross] = rMomentsGenerate_v1_exp (response1,response1_shuff);
    [R1_square_Sp,R1_cross_Sp] = rMomentsGenerate_v1_centra (response1_Sp,response1_Sp_shuff,rmean_Shatp,rmean_Shatp_shuff);
    [R1_square_Sn,R1_cross_Sn] = rMomentsGenerate_v1_centra (response1_Sn,response1_Sn_shuff,rmean_Shatn,rmean_Shatn_shuff);
 elseif shuffle==2 % shuffle r given s and n.
     shufflingtimes=100000;
     indexp=randi([1 numel(Choice_Sp)],1,shufflingtimes);
     indexn=randi([1 numel(Choice_Sn)],1,shufflingtimes);
    response1_Sn_shuff=response1_Sn;
    response1_Sp_shuff=response1_Sp;
    for i=1:numel(indexp)/2
        temp=response1_Sp_shuff(indexp(2*i),:);
        response1_Sp_shuff(indexp(2*i),:)=response1_Sp_shuff(indexp(2*i-1),:);
        response1_Sp_shuff(indexp(2*i-1),:)=temp;
    end    
    for i=1:numel(indexn)/2
        temp=response1_Sn_shuff(indexn(2*i),:);
        response1_Sn_shuff(indexn(2*i),:)=response1_Sn_shuff(indexn(2*i-1),:);
        response1_Sn_shuff(indexn(2*i-1),:)=temp;
    end
    response1_shuff=[response1_Sp_shuff;response1_Sn_shuff];
    [rmean_Shatn_shuff,rmean_Shatp_shuff] = R_Centralization_decoded(response1_Sn_shuff, response1_Sp_shuff, Stimulus1);
    [rmean_Shatn,rmean_Shatp] = R_Centralization_decoded(response1_Sn, response1_Sp, Stimulus1);
    [R2_square,R2_cross] = rMomentsGenerate_v1_exp (response1,response1_shuff);
    [R1_square_Sp,R1_cross_Sp] = rMomentsGenerate_v1_centra (response1_Sp,response1_Sp_shuff,rmean_Shatp,rmean_Shatp_shuff);
    [R1_square_Sn,R1_cross_Sn] = rMomentsGenerate_v1_centra (response1_Sn,response1_Sn_shuff,rmean_Shatn,rmean_Shatn_shuff);
 elseif shuffle==0 % No shuffling
    [rmean_Shatn,rmean_Shatp] = R_Centralization_decoded(response1_Sn, response1_Sp, Stimulus1);
    [R2_square,R2_cross] = rMomentsGenerate_v1_exp (response1,response1);
    [R1_square_Sp,R1_cross_Sp] = rMomentsGenerate_v1_centra (response1_Sp,response1_Sp, rmean_Shatp,rmean_Shatp);
    [R1_square_Sn,R1_cross_Sn] = rMomentsGenerate_v1_centra (response1_Sn,response1_Sn, rmean_Shatn,rmean_Shatn);
    
    
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


function [ccTheo] = v1_ccTheo_Compute (Stimulus1,Choice1,Orientation1,response1,shuffle,sT_shuffled,sT,sp,sn)
    
    
    
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
 if shuffle==1 % shuffle r given orientation and s
    [~,indexp]=sort(Orientation_Sp);
    [~,indexn]=sort(Orientation_Sn);
    response1_Sn_shuff=response1_Sn;
    response1_Sp_shuff=response1_Sp;
    for i=1:numel(Orientation_Sp)/2
        temp=response1_Sp_shuff(indexp(2*i),:);
        response1_Sp_shuff(indexp(2*i),:)=response1_Sp_shuff(indexp(2*i-1),:);
        response1_Sp_shuff(indexp(2*i-1),:)=temp;
    end
    
    for i=1:numel(Orientation_Sn)/2
        temp=response1_Sn_shuff(indexn(2*i),:);
        response1_Sn_shuff(indexn(2*i),:)=response1_Sn_shuff(indexn(2*i-1),:);
        response1_Sn_shuff(indexn(2*i-1),:)=temp;
    end
    response1_shuff=[response1_Sp_shuff;response1_Sn_shuff];    
    [R2_square,R2_cross] = rMomentsGenerate_v1_exp (response1,response1_shuff);
 elseif shuffle==2 % shuffle r given s and n.
     shufflingtimes=100000;
     indexp=randi([1 numel(Choice_Sp)],1,shufflingtimes);
     indexn=randi([1 numel(Choice_Sn)],1,shufflingtimes);
    response1_Sn_shuff=response1_Sn;
    response1_Sp_shuff=response1_Sp;
    for i=1:numel(indexp)/2
        temp=response1_Sp_shuff(indexp(2*i),:);
        response1_Sp_shuff(indexp(2*i),:)=response1_Sp_shuff(indexp(2*i-1),:);
        response1_Sp_shuff(indexp(2*i-1),:)=temp;
    end    
    for i=1:numel(indexn)/2
        temp=response1_Sn_shuff(indexn(2*i),:);
        response1_Sn_shuff(indexn(2*i),:)=response1_Sn_shuff(indexn(2*i-1),:);
        response1_Sn_shuff(indexn(2*i-1),:)=temp;
    end
    response1_shuff=[response1_Sp_shuff;response1_Sn_shuff];    
    [R2_square,R2_cross] = rMomentsGenerate_v1_exp (response1,response1_shuff);
 elseif shuffle==0 % No shuffling
    [R2_square,R2_cross] = rMomentsGenerate_v1_exp (response1,response1);
 end

%% Compute cctheo and ccsim

Corr_s_R=corr(sT_shuffled',[response1,R2_square,R2_cross]);
Corr_s_R_true=corr(sT',[response1,R2_square,R2_cross]);
Corr_s_shatR=corr(sT',Choice1);
ccTheo=CCpredRatio.*Corr_s_R./Corr_s_shatR;


end

%%
function [ccSim] = v1_ccSim_Compute (Stimulus1,Choice1,Orientation1,response1,shuffle)
Orientation_Sp=Orientation1(Stimulus1==1);
Orientation_Sn=Orientation1(Stimulus1==-1);
Choice_Sp=Choice1(Stimulus1==1);
Choice_Sn=Choice1(Stimulus1==-1);
Stim_Sp=Stimulus1(Stimulus1==1);
Stim_Sn=Stimulus1(Stimulus1==-1);
response1_Sp= response1(Stimulus1==1,:);
response1_Sn= response1(Stimulus1==-1,:);

%% Combine
Stimulus1=[Stim_Sp;Stim_Sn];
Choice1=[Choice_Sp;Choice_Sn];
Orientation1=[Orientation_Sp,Orientation_Sn];
response1=[response1_Sp;response1_Sn];
n=size(response1,2);
m=size(response1,1);

%% Shuffle
 if shuffle==1 % shuffle r given orientation and s
    [~,indexp]=sort(Orientation_Sp);
    [~,indexn]=sort(Orientation_Sn);
    response1_Sn_shuff=response1_Sn;
    response1_Sp_shuff=response1_Sp;
    for i=1:numel(Orientation_Sp)/2
        temp=response1_Sp_shuff(indexp(2*i),:);
        response1_Sp_shuff(indexp(2*i),:)=response1_Sp_shuff(indexp(2*i-1),:);
        response1_Sp_shuff(indexp(2*i-1),:)=temp;
    end
    
    for i=1:numel(Orientation_Sn)/2
        temp=response1_Sn_shuff(indexn(2*i),:);
        response1_Sn_shuff(indexn(2*i),:)=response1_Sn_shuff(indexn(2*i-1),:);
        response1_Sn_shuff(indexn(2*i-1),:)=temp;
    end
    [rmean_Shatn_shuff,rmean_Shatp_shuff] = R_Centralization_decoded(response1_Sn_shuff, response1_Sp_shuff, Stimulus1);
    [rmean_Shatn,rmean_Shatp] = R_Centralization_decoded(response1_Sn, response1_Sp, Stimulus1);
    [R1_square_Sp,R1_cross_Sp] = rMomentsGenerate_v1_centra (response1_Sp,response1_Sp_shuff,rmean_Shatp,rmean_Shatp_shuff);
    [R1_square_Sn,R1_cross_Sn] = rMomentsGenerate_v1_centra (response1_Sn,response1_Sn_shuff,rmean_Shatn,rmean_Shatn_shuff);
 elseif shuffle==2 % shuffle r given s and n.
     shufflingtimes=100000;
     indexp=randi([1 numel(Choice_Sp)],1,shufflingtimes);
     indexn=randi([1 numel(Choice_Sn)],1,shufflingtimes);
    response1_Sn_shuff=response1_Sn;
    response1_Sp_shuff=response1_Sp;
    for i=1:numel(indexp)/2
        temp=response1_Sp_shuff(indexp(2*i),:);
        response1_Sp_shuff(indexp(2*i),:)=response1_Sp_shuff(indexp(2*i-1),:);
        response1_Sp_shuff(indexp(2*i-1),:)=temp;
    end    
    for i=1:numel(indexn)/2
        temp=response1_Sn_shuff(indexn(2*i),:);
        response1_Sn_shuff(indexn(2*i),:)=response1_Sn_shuff(indexn(2*i-1),:);
        response1_Sn_shuff(indexn(2*i-1),:)=temp;
    end
    [rmean_Shatn_shuff,rmean_Shatp_shuff] = R_Centralization_decoded(response1_Sn_shuff, response1_Sp_shuff, Stimulus1);
    [rmean_Shatn,rmean_Shatp] = R_Centralization_decoded(response1_Sn, response1_Sp, Stimulus1);
    [R1_square_Sp,R1_cross_Sp] = rMomentsGenerate_v1_centra (response1_Sp,response1_Sp_shuff,rmean_Shatp,rmean_Shatp_shuff);
    [R1_square_Sn,R1_cross_Sn] = rMomentsGenerate_v1_centra (response1_Sn,response1_Sn_shuff,rmean_Shatn,rmean_Shatn_shuff);

    
 elseif shuffle==0 % No shuffling
    [rmean_Shatn,rmean_Shatp] = R_Centralization_decoded(response1_Sn, response1_Sp, Stimulus1);
    [R1_square_Sp,R1_cross_Sp] = rMomentsGenerate_v1_centra (response1_Sp,response1_Sp, rmean_Shatp,rmean_Shatp);
    [R1_square_Sn,R1_cross_Sn] = rMomentsGenerate_v1_centra (response1_Sn,response1_Sn, rmean_Shatn,rmean_Shatn);

 end

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
end
end
%%
function [R2_square,R2_cross] = rMomentsGenerate_v1_centra (r1,r2,rmean1,rmean2)
    n = size(r1,2);
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


function pComb=pValueSepToComb(px,py)
    pComb=0.5.*erfc(-sqrt((erfcinv(2.*px)).^2+(erfcinv(2.*py)).^2));
end


%% Estimate Stimulus from linear response, use that to centralize the r to get deltar
function [rmean_Shatn,rmean_Shatp] = R_Centralization_decoded(response1_Sn, response1_Sp, Stimulus1)
    mp=sum(Stimulus1(:) == 1)   ;
    mn=sum(Stimulus1(:) == -1)   ;
    m=mp+mn;
    Fb_Sn=mean(response1_Sn,1);
    Fb_Sp=mean(response1_Sp,1);
    response1=[response1_Sp;response1_Sn];
    srange=2;
    fp=(Fb_Sp-Fb_Sn)./srange;
    F0=mean(response1,1);
    Rref=repmat(F0,[m,1]);
    wq=pinv(response1-Rref)*Stimulus1;
    wq=wq./(fp*wq);
    shat_1=(response1-Rref)*wq;
    Stimulus_est_linear=((shat_1>0)-0.5).*2;

    response1_Shatn=response1(Stimulus_est_linear==-1,:);
    response1_Shatp=response1(Stimulus_est_linear==1,:);

    rmean_Shatn=repmat(mean(response1_Shatn,1),mn,1);
    rmean_Shatp=repmat(mean(response1_Shatp,1),mp,1);
end


