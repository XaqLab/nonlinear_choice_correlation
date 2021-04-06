MonkeyNum=2;
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
%     ResponseAll(SelectedSession,:,:)
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
    
     clearvars sT
     for j=1:m
        if Stimulus1(j)==1
            sT(j)=sp;
        else
            sT(j)=sn;
        end
     end
    
end

% Binwidth=0.1;
% figure
% subplot(131)
% histogram(dp_1(SelectedSession,:),'Normalization','probability','BinWidth',Binwidth,'FaceColor','b'); hold on;
% axis square
% box off
% axis([-0.8 0.8 0,0.2])
% subplot(132)
% histogram(dp_square(SelectedSession,:),'Normalization','probability','BinWidth',Binwidth,'FaceColor','g'); hold on;
% axis square
% box off
% axis([-0.8 0.8 0,0.2])
% subplot(133)
% histogram(dp_cross(SelectedSession,:),'Normalization','probability','BinWidth',Binwidth,'FaceColor','r'); hold on;
% axis square
% box off
% axis([-0.8 0.8 0,0.2])
% 
% range=0.25;
% Binwidth=2;
% figure
% histogram(Orientation_n,'Normalization','probability','BinWidth',Binwidth,'EdgeColor','none'); hold on;
% histogram(Orientation_p,'Normalization','probability','BinWidth',Binwidth,'EdgeColor','none'); hold on;
% % histogram(Orientation_Sn,'Normalization','probability','BinWidth',1); hold on;
% % histogram(Orientation_Sp,'Normalization','probability','BinWidth',1); hold on;
% legend('s+','s-')
% 
% plot([DecisionBoundary,DecisionBoundary],[0,range],'k-');hold on;
% plot([DecisionBoundary+binsize/2,DecisionBoundary+binsize/2],[0,range],'k-');hold on;
% plot([DecisionBoundary-binsize/2,DecisionBoundary-binsize/2],[0,range],'k-');hold on;
% 
% 
% plot([-DecisionBoundary,-DecisionBoundary],[0,range],'k-');hold on;
% plot([-DecisionBoundary+binsize/2,-DecisionBoundary+binsize/2],[0,range],'k-');hold on;
% plot([-DecisionBoundary-binsize/2,-DecisionBoundary-binsize/2],[0,range],'k-');hold on;
% axis square
% box off

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

