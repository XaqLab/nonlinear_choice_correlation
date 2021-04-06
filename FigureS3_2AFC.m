clear all
close all
%% Setting
n=20;     %neuron number
s0=0.5;     %stimulus reference
s_range=0.2
% sampletimes=[1:10];

Times=20;

% fp=1.*randn(n,1);
% fp=3*[-1:2./(n-1):1]';

% fp=1.*randn(n,1);
fp=1*[-1:2./(n-1):1]';
% fp=3*[-1:2./(n-1):1]';

%% Reference responses and statistics




%% Compute theoretical and experimental CC
sampletimes=[1:10];
m_samples=sampletimes.*(n*(n+1)/2+n);   %sample

for kk=1:numel(m_samples)
    A1=randn(n,n);
    for jj=1:Times
        m=m_samples(kk);
        sref=s0.*ones(1,m);

    Stimulus1=[-1.*ones(1,m/2),1.*ones(1,m/2)]';
    s_binary=Stimulus1.*s_range./2+s0;

    [response0] = QuadraticBrain(fp,A1,s0,m);
    [response1_Sp] = QuadraticBrain(fp,A1,s0+s_range./2,m/2);
    [response1_Sn] = QuadraticBrain(fp,A1,s0-s_range./2,m/2);


%% Decode with all statistics
[rmean_Shatn,rmean_Shatp] = R_Centralization_decoded(response1_Sn, response1_Sp, Stimulus1);
[R1_square_Sp,R1_cross_Sp] = rMomentsGenerate_v1_centra (response1_Sp,response1_Sp, rmean_Shatp,rmean_Shatp);
[R1_square_Sn,R1_cross_Sn] = rMomentsGenerate_v1_centra (response1_Sn,response1_Sn, rmean_Shatn,rmean_Shatn);
rmean0=repmat(mean(response0,1),[m,1]);
[R0_square,R0_cross] = rMomentsGenerate_v1_centra (response0,response0, rmean0,rmean0);


R_all_Sp=[response1_Sp,R1_square_Sp,R1_cross_Sp];
R_all_Sn=[response1_Sn,R1_square_Sn,R1_cross_Sn];
R_all_S0=[response0,R0_square,R0_cross];

R_all=[R_all_Sp;R_all_Sn];
[Choice1,Choice0,sbhat,dpRb,wq] = DecodingEngine(R_all,R_all_S0,s0,Stimulus1);
        [ccSim1(jj,:),ccTheo1(jj,:),ccSimSquare(jj,:),ccTheoSquare(jj,:),ccSimCross(jj,:),ccTheoCross(jj,:),slope1(jj),slopeSquare(jj),slopeCross(jj),dp_Population(jj)] = v1_ccAnalysis (Stimulus1,Choice1,Choice0,R_all_Sn,R_all_Sp,R_all_S0)
    
    end
dp_pop_mean(kk)=mean(dp_Population);
ccTheo1_example(kk,:)=ccTheo1(1,:);
ccTheoCross_example(kk,:)=ccTheoCross(1,:);
ccTheoSquare_example(kk,:)=ccTheoSquare(1,:);

ccSim1_example(kk,:)=ccSim1(1,:);
ccSimCross_example(kk,:)=ccSimCross(1,:);
ccSimSquare_example(kk,:)=ccSimSquare(1,:);



slope1std(kk)=std(slope1,1,2);
slopeSquarestd(kk)=std(slopeSquare,1,2);
slopeCrossstd(kk)=std(slopeCross,1,2);

slope1Mean(kk)=mean(slope1,2);
slopeSquareMean(kk)=mean(slopeSquare,2);
slopeCrossMean(kk)=mean(slopeCross,2);

end


dotsize=20;

figure
plot(sampletimes,slope1Mean,'b-');
hold on;
plot(sampletimes,slopeSquareMean,'g-');
hold on;
plot(sampletimes,slopeCrossMean,'r-');
set(gca,'linewidth',1,'fontsize',18,'fontname','CMU Serif');
axis square

figure
plot(sampletimes,slope1std,'b-');
hold on;
plot(sampletimes,slopeSquarestd,'g-');
hold on;
plot(sampletimes,slopeCrossstd,'r-');
set(gca,'linewidth',1,'fontsize',18,'fontname','CMU Serif');
axis square
% 
% axis square
% % legend('2','1');
% % set(gca,'XTick',[-0.5:0.5:0.5])
% % set(gca,'YTick',[-0.5:0.5:0.5])
% % xlabel('Theoretical CC','FontSize',18);
% % ylabel('Simulated CC','FontSize',18);
% set(gca,'linewidth',1,'fontsize',18,'fontname','CMU Serif');
% % title('Informative and decoded')
% % PlotAxisAtOrigin(0,0)
% 
% box off

marksize=10;
plotrange=0.3;
figure
m_plot=[2,5,10];
pp=0;
for jj=1:numel(m_plot)
    kk=m_plot(jj);
    pp=pp+1;
    Pointsize=5;
    subplot(1,3,pp)
    plot(ccTheoCross_example(kk,:),ccSimCross_example(kk,:),'r.','markersize', marksize);hold on;
    plot(ccTheo1_example(kk,:),ccSim1_example(kk,:),'b.','markersize', marksize);hold on;
    plot(ccTheoSquare_example(kk,:),ccSimSquare_example(kk,:),'g.','markersize', marksize);hold on;
%     [slope1_exp] = Fitting_cc (ccTheo1_example(kk,:),ccSim1_example(kk,:));
%     [slopeSquare_exp] = Fitting_cc (ccTheoSquare_example(kk,:),ccSimSquare_example(kk,:));
%     [slopeCross_exp] = Fitting_cc (ccTheoCross_example(kk,:),ccSimCross_example(kk,:));
%     plot([-plotrange,plotrange],slope1_exp.*[-plotrange,plotrange],'b-','markersize', 6);hold on;
%     plot([-plotrange,plotrange],slopeCross_exp.*[-plotrange,plotrange],'r-','markersize', 6);hold on;
%     plot([-plotrange,plotrange],slopeSquare_exp.*[-plotrange,plotrange],'g-','markersize', 6);hold on;
    plot([-plotrange,plotrange],[-plotrange,plotrange] ,'k-','markersize', 6);hold on;
    axis square;
    set(gca,'linewidth',1,'fontsize',18,'fontname','CMU Serif');
end


%% Useful functions

function [slopePCA] = Fitting_cc (ccTheo,ccSim)
ccdat=[ccTheo;ccSim];
coeff = pca(ccdat,2);
slopePCA = abs(coeff(1,2) / coeff(1,1));

end

function [r] = QuadraticBrain(fp,A1,s,m)
    %% Linear Setting
    n=numel(fp);
    u=fp.*s; % Informative linear
    %% Quadratic Setting
    A=A1-A1'; %Skew-symmetric matrix  
    [U,Omg] = eig(A); %U*Omg*U'=A  eigen-decomposition
    % Lambda=diag(1:1/(n-1):2); %positive definite matrix's eigenvalues
    Lambda=diag(1:1/(n-1):2);
    V=expm(A.*s); %rotation matrix, controled by stimulus
    sig=V*Lambda*V'; %V'=inv(V)=expm(-As0) 

    eOmg=expm(Omg.*s); %expm(A)=U*expm(Omg)*U'
    eiOmg=expm(-Omg.*s);

    sigp=U*Omg*eOmg*U'*Lambda*U*eiOmg*U'-...
         U*eOmg*U'*Lambda*U*Omg*eiOmg*U';
    sigp=real(sigp);
    isig=inv(sig);

    sig0=eye(n); 
    sigp0=zeros(n,n);

    r=mvnrnd(u,sig,m);
end

function [Choice,Choice0,sbhat,dpRb,wq] = DecodingEngine(Rb,R0,s0,sb)
m=size(Rb,1);
srange=sb(end)-sb(1);
Rbn=Rb( 1:m/2 , :);
Rbp=Rb( m/2+1:end , :);
Fb_minus=mean(Rbn,1);
Fb_plus=mean(Rbp,1);
fp=(Fb_plus-Fb_minus)./srange;
F0=mean(Rb,1);
Rref=repmat(F0,[m,1]);
wq=pinv(Rb-Rref)*(sb-s0);
wq=wq./(fp*wq);
% shat=(R0-Rref)*wq+s0;
sbhat=(Rb-Rref)*wq+s0;
Choice=((sbhat>s0)-0.5).*2;
dpRb=((mean(Rbp,1)-mean(Rbn,1))./(0.5.*(std(Rbp,1,1).^2+std(Rbn,1,1).^2)).^0.5);

F0=mean(R0,1);
Rref0=repmat(F0,[m,1]);
sbhat0=(R0-Rref0)*wq+s0;
Choice0=((sbhat0>s0)-0.5).*2;

end

%% Estimate Stimulus from linear response, use that to centralize the r to get deltar

%% Optimal decoder

% ij=0;
% for i = 1:n
%     for j = i:n
%         ij=ij+1;
%         kl=0;
%         for k=1:n
%             for l=k:n
%                 kl=kl+1;
%                 gamma(ij,kl) = sig(i,k)*sig(j,l) + sig(j,k)*sig(i,l); % with mean removed, gives cov(sig)
%             end
%         end
%     end
% end
% gamma = (gamma+gamma')/2.; % enforce symmetry
% igamma=inv(gamma);
% 
% J2=Fp'*igamma*Fp %Fisher information
% J1=fp'*isig*fp
% 
%     J=J1+J2;
%     wopt2=igamma*Fp./J;
%     wopt1=isig*fp./J;
% 
% shopt=wopt2'*(R2'-repmat(F,[1,m])) + wopt1'*(r'-repmat(u,[1,m]))  + s;





%%
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

function [ccSim1,ccTheo1,ccSimSquare,ccTheoSquare,ccSimCross,ccTheoCross,slope1,slopeSquare,slopeCross,dp_Population] = v1_ccAnalysis (Stimulus1,Choice1,Choice0,R_all_Sn,R_all_Sp,R_all_S0)
n=20;
DataAccuracy=sum(Stimulus1==Choice1)./numel(Stimulus1);
dp_Population = norminv(DataAccuracy,0,1).*2;
R_all=[R_all_Sp;R_all_Sn];

%% Seperation

% Choice_A=Choice(Stimulus1==-1);
% Choice_B=Choice(Stimulus1==1);

% ccSim_n=corr(Choice_A,R_all_Sn);
% ccSim_p=corr(Choice_B,R_all_Sp);
    



%% Theo
Fb_n=mean(R_all_Sn,1);
Fb_p=mean(R_all_Sp,1);

dpRb=-((Fb_p-Fb_n)./(0.5.*(std(R_all_Sn,1,1).^2+std(R_all_Sp,1,1).^2)).^0.5);
ccTheo=dpRb./dp_Population;

%% Sim

% ccSim=(numel(Choice_A).*ccSim_n+numel(Choice_B).*ccSim_p)./numel(Choice);

ccSim=corr(Choice0,R_all_S0);

% ccSim=ccSim_p;
%% Seperation
ccTheo1=ccTheo(1:n);
ccTheoSquare=ccTheo(n+1:n+n);
ccTheoCross=ccTheo(2*n+1 : end);

ccSim1=ccSim(1:n);
ccSimSquare=ccSim(n+1:n+n);
ccSimCross=ccSim(2*n+1 : end);
[slope1] = Fitting_cc (ccTheo1,ccSim1);
[slopeSquare] = Fitting_cc (ccTheoSquare,ccSimSquare);
[slopeCross] = Fitting_cc (ccTheoCross,ccSimCross);

end
