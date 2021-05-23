%% This code generates Figure S5C: Histograms of choice correlations show that purely internal noise is not significantly correlated with choices 
%% (p<0.01, two-sample Kolmogorov-Smirnov test for a choice-shuffled null distribution), unlike the nuisance-generated fluctuations seen in Figure 6D. 
%% To isolate the correlation of internal noise on choice, we compute the Normalized Average Conditional Choice Correlation (NACCC) where we condition on, 
%% and then average over, the complete stimulus (s,?) rather than just on the task-relevant stimulus s as in Eq.17 Individual choice correlations 
%% within the histograms are each colored by their significance according to their own null distribution (Methods 4.4).

clc
clear
close all

for shuffle=0:0
    for MonkeyNum=1:2
        [ccSim1BS_F{shuffle+1,MonkeyNum},ccSim2_squareBS_F{shuffle+1,MonkeyNum},ccSim2_crossBS_F{shuffle+1,MonkeyNum},ccTheo2_cross_F{shuffle+1,MonkeyNum},ccSim2_cross_F{shuffle+1,MonkeyNum},ccTheo2_square_F{shuffle+1,MonkeyNum},ccSim2_square_F{shuffle+1,MonkeyNum},ccTheo1_F{shuffle+1,MonkeyNum},ccSim1_F{shuffle+1,MonkeyNum},pSignificant_Theo1_F{shuffle+1,MonkeyNum},pSignificant_Theo2sq_F{shuffle+1,MonkeyNum},pSignificant_Theo2cr_F{shuffle+1,MonkeyNum},pSignificant_Sim1_F{shuffle+1,MonkeyNum},pSignificant_Sim2sq_F{shuffle+1,MonkeyNum},pSignificant_Sim2cr_F{shuffle+1,MonkeyNum}] = ccComputeCombined (MonkeyNum,shuffle);
        corcor_cross(shuffle+1,MonkeyNum)=corr(ccTheo2_cross_F{shuffle+1,MonkeyNum}',ccSim2_cross_F{shuffle+1,MonkeyNum}');
        corcor_square(shuffle+1,MonkeyNum)=corr(ccTheo2_square_F{shuffle+1,MonkeyNum}',ccSim2_square_F{shuffle+1,MonkeyNum}');
        corcor_1(shuffle+1,MonkeyNum)=corr(ccTheo1_F{shuffle+1,MonkeyNum}',ccSim1_F{shuffle+1,MonkeyNum}');
    end
end


%% Figures
figure
jj=0;
for MonkeyNum=1:2
        jj=jj+1;
        Pointsize=5;
        subplot(1,2,jj)
        plot(ccTheo2_cross_F{shuffle+1,MonkeyNum},ccSim2_cross_F{shuffle+1,MonkeyNum},'r.','markersize', 1);hold on;
        plot(ccTheo1_F{shuffle+1,MonkeyNum},ccSim1_F{shuffle+1,MonkeyNum},'b.','markersize', 1);hold on;
        plot(ccTheo2_square_F{shuffle+1,MonkeyNum},ccSim2_square_F{shuffle+1,MonkeyNum},'g.','markersize', 1);hold on;
        plot([-1,1],[-1,1] ,'k-','markersize', 6);hold on;
        axis square;
        set(gca,'XTick',[-1,1]);
        set(gca,'YTick',[-1,1]);
        set(gca,'Yticklabel',[-1,1]);
        set(gca,'Xticklabel',[-1,1]);
        axis([-1 1 -1 1]);
        if jj==1
            title('M1','FontSize',10);
            legend('cross','linear','square');
            xlabel('Optimal NACCC','FontSize',10);
            ylabel('Measured NACCC','FontSize',10);
        elseif jj==2
            title('M2','FontSize',10);
        end
end





        


%% Generate Null distributions of CCS
RandIndexM1=datasample(1:numel(ccSim2_cross_F{1,1}),numel(ccSim1_F{1,1}),'Replace',false)
RandIndexM2=datasample(1:numel(ccSim2_cross_F{1,2}),numel(ccSim1_F{1,2}),'Replace',false)

% RandIndexM1=[1:numel(ccSim2_cross_F{1,1})];
% RandIndexM2=[1:numel(ccSim2_cross_F{1,2})];


ccSim1BSall=ccSim1BS_F{1,1};
ccSim1BS1{1,1}=ccSim1BSall(1,:);

ccSim2_squareBSall=ccSim2_squareBS_F{1,1};
ccSim2_squareBS1{1,1}=ccSim2_squareBSall(1,:);

ccSim2_crossBSall=ccSim2_crossBS_F{1,1};
ccSim2_crossBS1{1,1}=ccSim2_crossBSall(1,RandIndexM1);

ccSim1BSall=ccSim1BS_F{1,2};
ccSim1BS1{1,2}=ccSim1BSall(1,:);

ccSim2_squareBSall=ccSim2_squareBS_F{1,2};
ccSim2_squareBS1{1,2}=ccSim2_squareBSall(1,:);

ccSim2_crossBSall=ccSim2_crossBS_F{1,2};
ccSim2_crossBS1{1,2}=ccSim2_crossBSall(1,RandIndexM2);

ccSim2_cross_M1=ccSim2_cross_F{1,1};
ccSim2_cross_M2=ccSim2_cross_F{1,2};

%% KS TEST
[h1_M1_all,p1_M1_all] = kstest2(reshape(ccSim1BS_F{1,1},1,[]),ccSim1_F{1,1},'Alpha',0.01)
[h2s_M1_all,p2s_M1_all] = kstest2(reshape(ccSim2_squareBS_F{1,1},1,[]),ccSim2_square_F{1,1},'Alpha',0.01)
[h2c_M1_all,p2c_M1_all] = kstest2(reshape(ccSim2_crossBS_F{1,1},1,[]),ccSim2_cross_F{1,1},'Alpha',0.01)

[h1_M2_all,p1_M2_all] = kstest2(reshape(ccSim1BS_F{1,2},1,[]),ccSim1_F{1,2},'Alpha',0.01)
[h2s_M2_all,p2s_M2_all] = kstest2(reshape(ccSim2_squareBS_F{1,2},1,[]),ccSim2_square_F{1,2},'Alpha',0.01)
[h2c_M2_all,p2c_M2_all] = kstest2(reshape(ccSim2_crossBS_F{1,2},1,[]),ccSim2_cross_F{1,2},'Alpha',0.01)


[h1_M1,p1_M1,Z_M1] = kstest2(ccSim1BS1{1,1},ccSim1_F{1,1},'Alpha',0.01)
[h2s_M1,p2s_M1,Z2s_M1] = kstest2(ccSim2_squareBS1{1,1},ccSim2_square_F{1,1},'Alpha',0.01)
[h2c_M1,p2c_M1,Z2c_M1] = kstest2(ccSim2_crossBS1{1,1},ccSim2_cross_M1(RandIndexM1),'Alpha',0.01)

[h1_M2,p1_M2,Z_M2] = kstest2(ccSim1BS1{1,2},ccSim1_F{1,2},'Alpha',0.01)
[h2s_M2,p2s_M2,Z2s_M2] = kstest2(ccSim2_squareBS1{1,2},ccSim2_square_F{1,2},'Alpha',0.01)
[h2c_M2,p2c_M2,Z2c_M2] = kstest2(ccSim2_crossBS1{1,2},ccSim2_cross_M2(RandIndexM2),'Alpha',0.01)

%% Plot histograms of both null and original CC distributions.
ccSim1BSallM1=ccSim1BS_F{1,1};
ccSim2_squareBSallM1=ccSim2_squareBS_F{1,1};
ccSim2_crossBSallM1=ccSim2_crossBS_F{1,1};

ccSim1BSallM2=ccSim1BS_F{1,2};
ccSim2_squareBSallM2=ccSim2_squareBS_F{1,2};
ccSim2_crossBSallM2=ccSim2_crossBS_F{1,2};


ccSim2_cross_F{1,1}= ccSim2_cross_M1(RandIndexM1);
ccSim2_cross_F{1,2}= ccSim2_cross_M2(RandIndexM2);

p_sim2cr1=pSignificant_Sim2cr_F{1,1};
p_sim2cr2=pSignificant_Sim2cr_F{1,2};
pSignificant_Sim2cr_F{1,1}=p_sim2cr1(RandIndexM1);
pSignificant_Sim2cr_F{1,2}=p_sim2cr2(RandIndexM2);

ccSim2_crossNull=ccSim2_crossBS1;
ccSim2_squareNull=ccSim2_squareBS1;
ccSim1Null=ccSim1BS1;

BinW=0.05;
shuffleIndex=1;

MonkeyNum=1;
figure
Pointsize=5;
subplot(3,1,1)
histogram(ccSim1BSallM1(shuffleIndex,:),'Normalization','probability','BinWidth',BinW,'BinLimits',[-0.5,0.5]);hold on;
histogram(ccSim1_F{shuffle+1,MonkeyNum},'Normalization','probability','BinWidth',BinW,'BinLimits',[-0.5,0.5]);
axis([-0.5,0.5,0,inf])
legend('Null','True')
subplot(3,1,2)
histogram(ccSim2_squareBSallM1(shuffleIndex,:),'Normalization','probability','BinWidth',BinW,'BinLimits',[-0.5,0.5]);hold on;
histogram(ccSim2_square_F{shuffle+1,MonkeyNum},'Normalization','probability','BinWidth',BinW,'BinLimits',[-0.5,0.5]);
axis([-0.5,0.5,0,inf])
subplot(3,1,3)
histogram(ccSim2_crossBS1{1,1},'Normalization','probability','BinWidth',BinW,'BinLimits',[-0.5,0.5]);hold on;
histogram(ccSim2_cross_M1(RandIndexM1),'Normalization','probability','BinWidth',BinW,'BinLimits',[-0.5,0.5]);
axis([-0.5,0.5,0,inf])

MonkeyNum=2;
figure
Pointsize=5;
subplot(3,1,1)
histogram(ccSim1BSallM2(shuffleIndex,:),'Normalization','probability','BinWidth',BinW,'BinLimits',[-0.5,0.5]);hold on;
histogram(ccSim1_F{shuffle+1,MonkeyNum},'Normalization','probability','BinWidth',BinW,'BinLimits',[-0.5,0.5]);
axis([-0.5,0.5,0,inf])
legend('Null','True')
subplot(3,1,2)
histogram(ccSim2_squareBSallM2(shuffleIndex,:),'Normalization','probability','BinWidth',BinW,'BinLimits',[-0.5,0.5]);hold on;
histogram(ccSim2_square_F{shuffle+1,MonkeyNum},'Normalization','probability','BinWidth',BinW,'BinLimits',[-0.5,0.5]);
axis([-0.5,0.5,0,inf])
subplot(3,1,3)
histogram(ccSim2_crossBS1{1,2},'Normalization','probability','BinWidth',BinW,'BinLimits',[-0.5,0.5]);hold on;
histogram(ccSim2_cross_M2(RandIndexM2),'Normalization','probability','BinWidth',BinW,'BinLimits',[-0.5,0.5]);
axis([-0.5,0.5,0,inf])



%% Save data for Pvalue plots

fname = sprintf('CC_Pvalue_Conditioned_Data.mat');
save(fname,'pSignificant_Sim1_F','pSignificant_Sim2sq_F','pSignificant_Sim2cr_F','ccSim1_F','ccSim2_cross_F','ccSim2_square_F','ccSim1Null','ccSim2_crossNull','ccSim2_squareNull');

function [ccSim1BS_F,ccSim2_squareBS_F,ccSim2_crossBS_F,ccTheo2_cross_F,ccSim2_cross_F,ccTheo2_square_F,ccSim2_square_F,ccTheo1_F,ccSim1_F,pSignificant_Theo1_F,pSignificant_Theo2sq_F,pSignificant_Theo2cr_F,pSignificant_Sim1_F,pSignificant_Sim2sq_F,pSignificant_Sim2cr_F] = ccComputeCombined (MonkeyNum,shuffle)

if MonkeyNum==1
    load('monkey1_conditioned_nuisance.mat');
else
    load('monkey2_conditioned_nuisance.mat');
end
SessionTotalNum=numel(data);
MonkeySessionList=[1:SessionTotalNum];
%% Compute Theo and experimental CCs for each session
count=0;
for SelectedSession=1:SessionTotalNum   
    SelectedSession
    SessionData=data{SelectedSession};
    %% Find same contrast  
    Contrast=SessionData.contrast;
    Index_SelectContrast=Contrast==Contrast;
    %% Pick data under same contrast.
    Stimulus=double(SessionData.stimulus_class=='B');
    Choice=double(SessionData.selected_class=='B');
    Orientation=double(SessionData.orientation);
    Orientation=Orientation-mean(Orientation);
    response=double(SessionData.Counts_matrix);
    [m,n]=size(response);
    
    if n<2
        continue;
    end
    
    Stimulus1=Stimulus(Index_SelectContrast);
    Choice1=Choice(Index_SelectContrast);    
    Stimulus1=Stimulus1.*2-1;
    Choice1=Choice1.*2-1;
    S1trial=sum(Stimulus1==1);
    S0trial=sum(Stimulus1==-1);
    C1trial=sum(Choice1==1);
    C0trial=sum(Choice1==-1);
    Accuracy=sum(Choice1==Stimulus1)./numel(Choice1);
    Corr_s_shatR=corr(Stimulus1,Choice1);

    
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
    if min([C1trial,C0trial,S1trial,S0trial])<5||m<50||Accuracy<0.60||Accuracy>0.99||Corr_s_shatR<0.25
        continue
    end
%     if abs(mean(Orientation_p)-mean(Orientation_n))>0.5||min([C1trial,C0trial,S1trial,S0trial])<3||m<50||Accuracy<0.40||Accuracy>0.99||Corr_s_shatR<0.25
%         continue
%     end

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
count=count+1;
[ccTheo,ccSim,CCPredRatio(SelectedSession),zeta(SelectedSession),delta(SelectedSession),n_half_fit(SelectedSession),us_S1(SelectedSession),us_S0(SelectedSession),OrientationMidMeaningful{SelectedSession},Choice_A_Ratio{SelectedSession},predRatio{SelectedSession}] = v1_ccAnalysis (Stimulus1,Choice1,Orientation1,response1,shuffle,sT,sbar,sp,sn);

ccTheo1{count}=ccTheo(1:n);
ccTheo2_square{count}=ccTheo(n+1:2*n);
ccTheo2_cross{count}=ccTheo(2*n+1:end);
ccSim1{count}=ccSim(1:n);
ccSim2_square{count}=ccSim(n+1:2*n);
ccSim2_cross{count}=ccSim(2*n+1:end);



[ccTheo1BS,ccTheo2_squareBS,ccTheo2_crossBS,ccSim1BS,ccSim2_squareBS,ccSim2_crossBS] = ccSignificance (Stimulus1,Choice1,Orientation1,response1,shuffle,sT,sp,sn);
mu_Theo1{count}=mean(ccTheo1BS);
std_Theo_1{count}=std(ccTheo1BS);
mu_Theo2square{count}=mean(ccTheo2_squareBS);
std_Theo2square{count}=std(ccTheo2_squareBS);
mu_Theo2cross{count}=mean(ccTheo2_crossBS);
std_Theo2cross{count}=std(ccTheo2_crossBS);
mu_Sim1{count}=mean(ccSim1BS);
std_Sim_1{count}=std(ccSim1BS);
mu_Sim2square{count}=mean(ccSim2_squareBS);
std_Sim2square{count}=std(ccSim2_squareBS);
mu_Sim2cross{count}=mean(ccSim2_crossBS);
std_Sim2cross{count}=std(ccSim2_crossBS);

ccSim1BS_All{count}=ccSim1BS;
ccSim2_squareBS_All{count}=ccSim2_squareBS;
ccSim2_crossBS_All{count}=ccSim2_crossBS;


end


%% Filtered CCs and Significance
ccSim2_cross_F=double.empty(1,0);
ccSim2_square_F=double.empty(1,0);
ccSim1_F=double.empty(1,0);
ccTheo2_cross_F=double.empty(1,0);
ccTheo2_square_F=double.empty(1,0);
ccTheo1_F=double.empty(1,0);

ccSim1BS_F=double.empty(1,0);
ccSim2_squareBS_F=double.empty(1,0);
ccSim2_crossBS_F=double.empty(1,0);

for jj=1:count
    SessionData=data{jj};
    response=double(SessionData.Counts_matrix);
    n=size(response,2);

    FilteredSessionNum=jj;
    ccSim2_cross_F=[ccSim2_cross_F,ccSim2_cross{FilteredSessionNum}];
    ccTheo2_cross_F=[ccTheo2_cross_F,ccTheo2_cross{FilteredSessionNum}];
    ccSim2_square_F=[ccSim2_square_F,ccSim2_square{FilteredSessionNum}];
    ccTheo2_square_F=[ccTheo2_square_F,ccTheo2_square{FilteredSessionNum}];
    ccSim1_F=[ccSim1_F,ccSim1{FilteredSessionNum}];
    ccTheo1_F=[ccTheo1_F,ccTheo1{FilteredSessionNum}];
    
    ccSim1BS_F=[ccSim1BS_F,ccSim1BS_All{FilteredSessionNum}];
    ccSim2_squareBS_F=[ccSim2_squareBS_F,ccSim2_squareBS_All{FilteredSessionNum}];
    ccSim2_crossBS_F=[ccSim2_crossBS_F,ccSim2_crossBS_All{FilteredSessionNum}];
end

%% 

pSignificant_Theo1_F=double.empty(1,0);
pSignificant_Theo2sq_F=double.empty(1,0);
pSignificant_Theo2cr_F=double.empty(1,0);
pSignificant_Sim1_F=double.empty(1,0);
pSignificant_Sim2sq_F=double.empty(1,0);
pSignificant_Sim2cr_F=double.empty(1,0);

for jj=1:numel(ccTheo1)
    SessionData=data{jj};
    response=double(SessionData.Counts_matrix);
    n=size(response,2);
    if n<2
        continue;
    end
    FilteredSessionNum=jj;
    pSignificant_Theo1_F=[pSignificant_Theo1_F, 1-normcdf(ccTheo1{FilteredSessionNum},mu_Theo1{FilteredSessionNum},std_Theo_1{FilteredSessionNum})];
    pSignificant_Theo2sq_F =[pSignificant_Theo2sq_F, 1-normcdf(ccTheo2_square{FilteredSessionNum},mu_Theo2square{FilteredSessionNum},std_Theo2square{FilteredSessionNum})];
    pSignificant_Theo2cr_F  = [pSignificant_Theo2cr_F,1-normcdf(ccTheo2_cross{FilteredSessionNum},mu_Theo2cross{FilteredSessionNum},std_Theo2cross{FilteredSessionNum})];
    pSignificant_Sim1_F  = [pSignificant_Sim1_F,1-normcdf(ccSim1{FilteredSessionNum},mu_Sim1{FilteredSessionNum},std_Sim_1{FilteredSessionNum})];
    pSignificant_Sim2sq_F  = [pSignificant_Sim2sq_F,1-normcdf(ccSim2_square{FilteredSessionNum},mu_Sim2square{FilteredSessionNum},std_Sim2square{FilteredSessionNum})];
    pSignificant_Sim2cr_F = [pSignificant_Sim2cr_F,1-normcdf(ccSim2_cross{FilteredSessionNum},mu_Sim2cross{FilteredSessionNum},std_Sim2cross{FilteredSessionNum})];
%     detectNAN(jj)=sum(isnan(pSignificant_Theo1_F));
end

nonNANindexCross=find((isnan(pSignificant_Sim2cr_F)|isnan(pSignificant_Theo2cr_F))==1);
ccSim2_cross_F(nonNANindexCross)=[];
ccTheo2_cross_F(nonNANindexCross)=[];
ccSim2_crossBS_F(:,nonNANindexCross)=[];
pSignificant_Theo2cr_F(nonNANindexCross)=[];
pSignificant_Sim2cr_F(nonNANindexCross)=[];

nonNANindexSquare=find((isnan(pSignificant_Theo2sq_F)|isnan(pSignificant_Sim2sq_F))==1);
ccSim2_square_F(nonNANindexSquare)=[];
ccTheo2_square_F(nonNANindexSquare)=[];
pSignificant_Theo2sq_F(nonNANindexSquare)=[];
pSignificant_Sim2sq_F(nonNANindexSquare)=[];
ccSim2_squareBS_F(:,nonNANindexSquare)=[];

nonNANindexLinear=find((isnan(pSignificant_Sim1_F)|isnan(pSignificant_Theo1_F))==1);
ccSim1_F(nonNANindexLinear)=[];
ccTheo1_F(nonNANindexLinear)=[];
pSignificant_Theo1_F(nonNANindexLinear)=[];
pSignificant_Sim1_F(nonNANindexLinear)=[];

ccSim1BS_F(:,nonNANindexLinear)=[];






end









%% Defined function
function [ccTheo,ccSim,CCpredRatio,zeta,delta,n_half_pred,us_S1,us_S0,OrientationMidMeaningful,Choice_A_Ratio,predRatio] = v1_ccAnalysis (Stimulus1,Choice1,Orientation1,response1,shuffle,sT,sbar,sp,sn)
sum(Stimulus1==Choice1)./numel(Stimulus1)
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
sum(Stimulus1==Choice1)./numel(Stimulus1)
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

% DataAccuracy=sum(Stimulus1==Choice1)./numel(Stimulus1)
% dp_Population = norminv(DataAccuracy,0,1).*2;
% Rb=[response1,R2_square,R2_cross];
% Rbn=Rb( Stimulus1==-1 , :);
% Rbp=Rb(Stimulus1==1 , :);
% dpRb=((mean(Rbp,1)-mean(Rbn,1))./(0.5.*(std(Rbp,1,1).^2+std(Rbn,1,1).^2)).^0.5);
% ccTheo=dpRb./dp_Population;






Corr_s_R=corr(sT',[response1,R2_square,R2_cross]);
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

%%
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
ccdat=[ccTheo;ccSim];
coeff = pca(ccdat,2);
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



%%
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
% ccTheo(find(isnan(ccTheo)==1))=0;

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
% ccSim(find(isnan(ccSim)==1))=0;
end
