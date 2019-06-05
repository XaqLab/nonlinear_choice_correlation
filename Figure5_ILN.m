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
F0=zeros(n*(n+1)/2,1);
F0b=zeros(n*(n+1)/2,1);
Fp=zeros(n*(n+1)/2,1);
R=zeros(m,n*(n+1)/2);
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
J0=Fp'*igamma0*Fp; %Fisher information
epsilon = 10./J0; % Amplitude of bad noise
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

wopt=igamma0*Fp./J0; % The optimal weights
woptb=wopt; % The bad noise doesn't change the weights direction
wsub = sign(Fp);% suboptimal decoding: blind to the correlation

shopt=wopt'*(R'-repmat(F0,[1,m])) + s0;
shsub=wsub'*(R'-repmat(F0,[1,m])) + s0;
shoptb=woptb'*(Rb'-repmat(F0,[1,m])) + s0;
shsubb=wsub'*(Rb'-repmat(F0,[1,m])) + s0;

%% Compute choice correlation under different conditions.
Coptpred=zeros(n*(n+1)/2,1);
Coptbpred=zeros(n*(n+1)/2,1);
Coptsim=zeros(n*(n+1)/2,1);
Csubsim=zeros(n*(n+1)/2,1);
Coptbsim=zeros(n*(n+1)/2,1);
Csubbsim=zeros(n*(n+1)/2,1);
ij=0;
for i=1:n
    for j=1:i
        ij=ij+1;
        Coptsim(ij) = corr(shopt',R(:,ij));
        Csubsim(ij) = corr(shsub',R(:,ij));
        Coptbsim(ij) = corr(shoptb',Rb(:,ij));
        Csubbsim(ij) = corr(shsubb',Rb(:,ij));
        Coptpred(ij) = Fp(ij)/sqrt(gamma0(ij,ij)) / sqrt(J0);
        Coptbpred(ij) = (Fp(ij)/sqrt(epsilon*Fp(ij)^2+gamma0(ij,ij)))/sqrt(J);
    end
end

%% Plotting
figure;
subplot(1,2,1);
plot(sign(Fp).*Coptpred,  sign(Fp).*Csubsim,'r.');hold on;
plot(sign(Fp).*Coptpred,  sign(Fp).*Coptsim,'b.',[0,0.25],[0,0.25],'k-');hold on;
axis equal;
title('Extensive information','FontSize',10);
xlabel('Optimal CC');
ylabel('Simulated CC');
legend('Suboptimal','Optimal');
axis([0,0.25,0,0.25]);

subplot(1,2,2);
plot(sign(Fp).*Coptbpred,sign(Fp).*Coptbsim,'b.',sign(Fp).*Coptbpred,sign(Fp).*Csubbsim,'r.',[0,0.6],[0,0.6],'k-');
axis equal
title('Limited information','FontSize',10);
axis([0,0.6,0,0.6]);