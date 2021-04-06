MonkeyNum=2;
if MonkeyNum==1
    SessionTotalNum=59;
    load('monkey1.mat');
else
    SessionTotalNum=71;
    load('monkey2.mat');
end
MonkeySessionList=[1:SessionTotalNum];
contrastAll=1; % flag to denote if we want to pick the all the contrasts or just the highest contrast.
%% Compute Theo and experimental CCs for each session
for SelectedSession=20:20 
    SelectedSession
    if MonkeyNum==1
        SessionData=monkey1{SelectedSession};
    else
        SessionData=monkey2{SelectedSession};
    end
    %% Find same contrast  
    Contrast=SessionData.contrast;
    Union_Contrast=unique(Contrast)
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
    trialnum(SelectedSession)=m;
    Orientation_p=Orientation1(Stimulus1==1);
    Orientation_n=Orientation1(Stimulus1==-1);
    sigma_p=15;sigma_n=3;
    sp=sigma_p^2;
    sn=sigma_n^2;
    thetahat_DB=sqrt((sp+sn)/2);
    DecisionBoundary=sqrt(2*log(sqrt((sp)/(sn)))/(1/(sn)-1/(sp)));
    
     sbar=thetahat_DB^2;
     clearvars sT
     for j=1:m
        if Stimulus1(j)==1
            sT(j)=sp;
        else
            sT(j)=sn;
        end
     end
    %% Seperate according to stimulus and select data around DB
    binsize=3;
    Condition1=(Stimulus1==1)';
    Condition0=(Stimulus1==-1)';
    Condition2=(Orientation1.^2<=(DecisionBoundary+binsize/2).^2);
    Condition3= Orientation1.^2>=(DecisionBoundary-binsize/2).^2 ;
   
    
    
    SelectedIndex_Sp= Condition1 & Condition2 & Condition3;
    SelectedIndex_Sn= Condition0 & Condition2 & Condition3;
    
    Orientation_Sp=Orientation1(SelectedIndex_Sp);
    Orientation_Sn=Orientation1(SelectedIndex_Sn);
    Choice_Sp=Choice1(SelectedIndex_Sp);
    Choice_Sn=Choice1(SelectedIndex_Sn);
    Stim_Sp=Stimulus1(SelectedIndex_Sp);
    Stim_Sn=Stimulus1(SelectedIndex_Sn);
    response1_Sp= response1(SelectedIndex_Sp,:);
    response1_Sn= response1(SelectedIndex_Sn,:);

%     figure
%     histogram(Orientation_n,'Normalization','probability','BinWidth',1); hold on;
%     histogram(Orientation_p,'Normalization','probability','BinWidth',1); hold on;
%     histogram(Orientation_Sn,'Normalization','probability','BinWidth',1); hold on;
%     histogram(Orientation_Sp,'Normalization','probability','BinWidth',1); hold on;
    
    %% Combine
    Stimulus1=[Stim_Sp;Stim_Sn];
    Choice1=[Choice_Sp;Choice_Sn];
    Orientation1=[Orientation_Sp,Orientation_Sn];
    response1=[response1_Sp;response1_Sn];
    n=size(response1,2);
    m=size(response1,1);
    [rmean_Shatn,rmean_Shatp] = R_Centralization_decoded(response1_Sn, response1_Sp, Stimulus1);
    [R1_square_Sp,R1_cross_Sp] = rMomentsGenerate_v1_centra (response1_Sp,response1_Sp, rmean_Shatp,rmean_Shatp);
    [R1_square_Sn,R1_cross_Sn] = rMomentsGenerate_v1_centra (response1_Sn,response1_Sn, rmean_Shatn,rmean_Shatn);
    R1_square1=[R1_square_Sp;R1_square_Sn];
    R1_cross1=[R1_cross_Sp;R1_cross_Sn];
        
    dp_1(SelectedSession,:)= dp_compute(response1_Sp,response1_Sn);
    dp_square(SelectedSession,:)= dp_compute(R1_square_Sp,R1_square_Sn);
    dp_cross(SelectedSession,:)= dp_compute(R1_cross_Sp,R1_cross_Sn);
    
end

Binwidth=0.1;
figure
subplot(131)
histogram(dp_1(SelectedSession,:),'Normalization','probability','BinWidth',Binwidth,'FaceColor','b'); hold on;
axis square
box off
axis([-0.8 0.8 0,0.2])
subplot(132)
histogram(dp_square(SelectedSession,:),'Normalization','probability','BinWidth',Binwidth,'FaceColor','g'); hold on;
axis square
box off
axis([-0.8 0.8 0,0.2])
subplot(133)
histogram(dp_cross(SelectedSession,:),'Normalization','probability','BinWidth',Binwidth,'FaceColor','r'); hold on;
axis square
box off
axis([-0.8 0.8 0,0.2])

range=0.25;
Binwidth=2;
figure
histogram(Orientation_n,'Normalization','probability','BinWidth',Binwidth,'EdgeColor','none'); hold on;
histogram(Orientation_p,'Normalization','probability','BinWidth',Binwidth,'EdgeColor','none'); hold on;
% histogram(Orientation_Sn,'Normalization','probability','BinWidth',1); hold on;
% histogram(Orientation_Sp,'Normalization','probability','BinWidth',1); hold on;
legend('s+','s-')

plot([DecisionBoundary,DecisionBoundary],[0,range],'k-');hold on;
plot([DecisionBoundary+binsize/2,DecisionBoundary+binsize/2],[0,range],'k-');hold on;
plot([DecisionBoundary-binsize/2,DecisionBoundary-binsize/2],[0,range],'k-');hold on;


plot([-DecisionBoundary,-DecisionBoundary],[0,range],'k-');hold on;
plot([-DecisionBoundary+binsize/2,-DecisionBoundary+binsize/2],[0,range],'k-');hold on;
plot([-DecisionBoundary-binsize/2,-DecisionBoundary-binsize/2],[0,range],'k-');hold on;
axis square
box off

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

function [dpRb] = dp_compute(Rbp,Rbn)
    Fb_minus=mean(Rbn,1);
    Fb_plus=mean(Rbp,1);
    dpRb=((mean(Rbp,1)-mean(Rbn,1))./(0.5.*(std(Rbp,1,1).^2+std(Rbn,1,1).^2)).^0.5);
end

