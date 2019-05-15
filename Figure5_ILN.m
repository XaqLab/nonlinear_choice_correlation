clear
close all

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
sigp=real(sigp);
isig=inv(sig);
%% to Get optimal weights without bad noise
u=zeros(1,n);
r=mvnrnd(u,sig,m);
F0=zeros(n*(n+1)/2,1);
F0b=zeros(n*(n+1)/2,1);
Fp=zeros(n*(n+1)/2,1);
R1=zeros(m,n*(n+1)/2);
ij=0;
for i=1:n
    for j=1:i
        ij=ij+1;
        fsim=mean(r,1);
        rmean=repmat(fsim,m,1);
        R1(:,ij)=(r(:,i)-rmean(:,i)).*(r(:,j)-rmean(:,j));
        Fp(ij) =sigp(i,j);
        F0(ij) = sig(i,j); 
    end
end
ij=0;
for i = 1:n
    for j = 1:i
        ij=ij+1;
        kl=0;
        for k=1:n
            for l=1:k
                kl=kl+1;
                gamma0(ij,kl) = sig(i,k)*sig(j,l) + sig(j,k)*sig(i,l); % with mean removed, gives cov(sig)
            end
        end
    end
end
gamma0 = (gamma0+gamma0')/2.; % enforce symmetry
igamma0=inv(gamma0);
J0=Fp'*igamma0*Fp; %Fisher information
%bad noise
epsilon = 10./J0; % bad noise = 2x variance of optimal estimator for unlimited info case
J = 1/(1/J0 + epsilon);
wopt=igamma0*Fp./J0;
woptb=wopt;
wsub = sign(Fp);% suboptimal decoding

ds=sqrt(epsilon).*randn(1,m); % bad noise
rb = zeros(m,n); % responses corrupted by bad noise
s=s0+ds;
R1b=zeros(m,n*(n+1)/2); % R with bad noise
ij=0;
for i=1:n
    for j=1:i
        ij=ij+1;
        R1b(:,ij)=r(:,i).*r(:,j)+ds'.*Fp(ij);
    end
end

shopt=wopt'*(R1'-repmat(F0,[1,m])) + s0;
shsub=wsub'*(R1'-repmat(F0,[1,m])) + s0;
shoptb=woptb'*(R1b'-repmat(F0,[1,m])) + s0;
shsubb=wsub'*(R1b'-repmat(F0,[1,m])) + s0;
%% Shuffle the choice
shopt_shuffle=shopt;
    for i=1:numel(shopt_shuffle)/2
        temp=shopt_shuffle((2*i));
        shopt_shuffle((2*i))=shopt_shuffle((2*i-1));
        shopt_shuffle((2*i-1))=temp;
    end
    
    for j=1:n
        index_shuffle=randperm(m); % shuffled index
        r_shuffle(:,j)=r([index_shuffle],j);
    end
    [R_shuffle] = rMomentsGenerate_v1_exp (r,r_shuffle);




%%
Coptpredij=zeros(n*(n+1)/2,1);
Coptbpredij=zeros(n*(n+1)/2,1);
Coptsimij=zeros(n*(n+1)/2,1);
Csubsimij=zeros(n*(n+1)/2,1);
Coptbsimij=zeros(n*(n+1)/2,1);
Csubbsimij=zeros(n*(n+1)/2,1);
ij=0;
for i=1:n
    for j=1:i
        ij=ij+1;
        Coptsimij(ij) = corr(shopt',R1(:,ij));
        Coptsim_shuffle(ij) = corr(shopt_shuffle',R1(:,ij));
        Coptsim_shuffleR(ij) = corr(shopt',R_shuffle(:,ij));
        
        Csubsimij(ij) = corr(shsub',R1(:,ij));
        Coptbsimij(ij) = corr(shoptb',R1b(:,ij));
        Csubbsimij(ij) = corr(shsubb',R1b(:,ij));
        Coptpredij(ij) = Fp(ij)/sqrt(gamma0(ij,ij)) / sqrt(J0);
        Coptbpredij(ij) = (Fp(ij)/sqrt(epsilon*Fp(ij)^2+gamma0(ij,ij)))/sqrt(J);
    end
end



%% Plotting
figure;
subplot(1,2,1);
plot(sign(Fp).*Coptpredij,  sign(Fp).*Csubsimij,'r.');hold on;
plot(sign(Fp).*Coptpredij,  sign(Fp).*Coptsimij,'b.',[0,0.25],[0,0.25],'k-');hold on;
axis equal;
title('Extensive information','FontSize',10);
xlabel('Optimal CC');
ylabel('Simulated CC');
legend('Suboptimal','Optimal');
axis([0,0.25,0,0.25]);

subplot(1,2,2);
plot(sign(Fp).*Coptbpredij,sign(Fp).*Coptbsimij,'b.',sign(Fp).*Coptbpredij,sign(Fp).*Csubbsimij,'r.',[0,0.6],[0,0.6],'k-');
axis equal
title('Limited information','FontSize',10);
axis([0,0.6,0,0.6]);


function [R2] = rMomentsGenerate_v1_exp (r1,r2)
m = size(r1,1);
n = size(r1,2);
rmean1=repmat(mean(r1,1),m,1);
rmean2=repmat(mean(r2,1),m,1);

z1=r1-rmean1;
z2=r2-rmean2;

k=0;
for i=1:n
    for j=i:n
        k=k+1;
        R2(:,k)=z1(:,i).*z2(:,j);
    end
end

end

