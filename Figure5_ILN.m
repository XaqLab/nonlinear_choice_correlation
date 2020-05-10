%% This code will generate Figure 5 in the paper: "Revealing nonlinear neural decoding by analyzing choices". 
%% The neural encoding model is a quadratic neural code. Then we build our decoder model as a quadratic decoder.
%% The weights are selected as optimal(Least square) or suboptimal (blind to higher).
%% We compare the choice correlation test for the optimal and suboptimal decoding cases.
%% Then we add information-limiting noise into the sufficient statistics.
%% Then we compare the choice correlation test for the optimal and suboptimal decoding cases under the limited information case.
clear
close all

%% Basic Setting
n=30;     %neuron number
s0=1.5;     %stimulus reference
m=50000;   %sample
A1=randn(n,n);
A=A1-A1'; %Skew-symmetric matrix  
[U,Omg] = eig(A); %U*Omg*U'=A  eigen-decomposition
Lambda=diag(1:3/(n-1):4); %positive definite matrix's eigenvalues
V=expm(A.*s0); %rotation matrix, controled by stimulus
sig=V*Lambda*V'; %V'=inv(V)=expm(-As0) 
eOmg=expm(Omg.*s0); %expm(A)=U*expm(Omg)*U'
eiOmg=expm(-Omg.*s0);
sigp=U*Omg*eOmg*U'*Lambda*U*eiOmg*U'-...
     U*eOmg*U'*Lambda*U*Omg*eiOmg*U'; 
sigp=real(sigp);% covariance derivative.
isig=inv(sig);
%% Generate neural responses and quadratic elements.
u=zeros(1,n);
r=mvnrnd(u,sig,m);
ij=0;
for i=1:n
    for j=1:i
        ij=ij+1;
        fsim=mean(r,1);
        rmean=repmat(fsim,m,1);
        R(:,ij)=(r(:,i)-rmean(:,i)).*(r(:,j)-rmean(:,j)); % compute quadratic elements
        Fp(ij) =sigp(i,j);
        F0(ij) = sig(i,j); 
    end
end

%% Decoding quadratic elements to get estimate of the stimulus.
ij=0;
for i = 1:n
    for j = 1:i
        ij=ij+1;
        kl=0;
        for k=1:n
            for l=1:k
                kl=kl+1;
                gamma0(ij,kl) = sig(i,k)*sig(j,l) + sig(j,k)*sig(i,l); % compute Cov(R1)
            end
        end
    end
end
gamma0 = (gamma0+gamma0')/2.; % enforce symmetry
igamma0=inv(gamma0);
J0=Fp*igamma0*Fp'; %Fisher information
epsilon = 4./J0; % Amplitude of bad noise
J = 1/(1/J0 + epsilon);
ds=sqrt(epsilon).*randn(1,m); % bad noise
s=s0+ds;
Rb=zeros(m,n*(n+1)/2); % R with bad noise
ij=0;
for i=1:n
    for j=1:i
        ij=ij+1;
        Rb(:,ij)=r(:,i).*r(:,j)+ds'.*Fp(ij); % add bad noise into the quadratic elements.
    end
end

wopt=igamma0*Fp'./J0; % The optimal weights
woptb=wopt; % The bad noise doesn't change the weights direction
wsub = sign(Fp');% suboptimal decoding: blind to the correlation
wsub=wsub./(Fp*wsub);
randsign = (binornd(1,0.6,numel(wopt),1)-0.5)*2;
wsubest = randsign.*sign(Fp');% suboptimal decoding: blind to the correlation
wsubest=wsubest./(Fp*wsubest);

shopt=(R-repmat(F0,[m,1]))*wopt + s0;
shsub=(R-repmat(F0,[m,1]))*wsub + s0;
shsubest=(R-repmat(F0,[m,1]))*wsubest + s0;

FI_opt=1./var(shopt);
FI_sub=1./var(shsub);
FI_subest=1./var(shsubest);

shoptb=(Rb-repmat(F0,[m,1]))*woptb + s0;
shsubb=(Rb-repmat(F0,[m,1]))*wsub + s0;
shsubestb=(Rb-repmat(F0,[m,1]))*wsubest + s0;

FI_opt_bad=1./var(shoptb);
FI_sub_bad=1./var(shsubb);
FI_subest_bad=1./var(shsubestb);

FI_sub/FI_opt
FI_subest/FI_opt

FI_sub_bad/FI_opt_bad
FI_subest_bad/FI_opt_bad


figure
subplot(1,2,1);
bar([FI_opt,FI_sub,FI_subest]);
axis square;
subplot(1,2,2);
bar([FI_opt_bad,FI_sub_bad,FI_subest_bad]);
axis square;

%% Compute choice correlation under different conditions.
ij=0;
for i=1:n
    for j=1:i
        ij=ij+1;
        Coptsim(ij) = corr(shopt,R(:,ij));
        Csubsim(ij) = corr(shsub,R(:,ij));
        Csubestsim(ij) = corr(shsubest,R(:,ij)); 
       
        Coptbsim(ij) = corr(shoptb,Rb(:,ij));
        Csubbsim(ij) = corr(shsubb,Rb(:,ij));
        Csubestbsim(ij) = corr(shsubestb,Rb(:,ij));
        
        Coptpred(ij) = Fp(ij)/sqrt(gamma0(ij,ij)) / sqrt(J0);
        Coptbpred(ij) = (Fp(ij)/sqrt(epsilon*Fp(ij)^2+gamma0(ij,ij)))/sqrt(J);
    end
end

%% Plotting
figure;
subplot(1,2,1);
plot(abs(Coptpred),  abs(Csubsim),'r.');hold on;
plot(abs(Coptpred),  abs(Csubestsim),'g.');hold on;
plot(abs(Coptpred),  abs(Coptsim),'b.');hold on;

[slope_sub,RSS_sub] = fit_and_rsqure(Csubsim,Coptpred)
[slope_subest,RSS_subest] = fit_and_rsqure(Csubestsim,Coptpred)
[slope_opt,RSS_opt] = fit_and_rsqure(Coptsim,Coptpred)

plot([0,0.25],[0,0.25],'k-');
axis equal;
title('Extensive information','FontSize',10);
xlabel('Optimal CC');
ylabel('Simulated CC');
% legend('Optimal','Suboptimal','worse Suboptimal');
axis([0,0.25,0,0.25]);

subplot(1,2,2);
plot(abs(Coptbpred),  abs(Csubbsim),'r.');hold on;
plot(abs(Coptbpred),  abs(Csubestbsim),'g.');hold on;
plot(abs(Coptbpred),  abs(Coptbsim),'b.');hold on;

[slope_sub_bad,RSS_sub_bad] = fit_and_rsqure(Csubbsim,Coptbpred)
[slope_subest_bad,RSS_subest_bad] = fit_and_rsqure(Csubestbsim,Coptbpred)
[slope_opt_bad,RSS_opt_bad] = fit_and_rsqure(Coptbsim,Coptbpred)

plot([0,0.6],[0,0.6],'k-');
axis equal
title('Limited information','FontSize',10);
axis([0,0.6,0,0.6]);

function [slope,Rsq] = fit_and_rsqure(CCsim,CCtheo)

[slope] = (abs(CCtheo))' \abs(CCsim)';
CCsim_pred = CCtheo*slope;
Rsq = sum((CCsim - CCsim_pred).^2);


end
