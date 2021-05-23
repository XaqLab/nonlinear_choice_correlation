%% This code generates Figure S6A: Information decoded from various numbers of linear (blue), quadratic (red), and cubic (green) statistics of downstream neural population.


clear
close all

%% Parameter Setting

n=15;
s0=0;     
srange=0.5;
m=100000;   %sample
gamma=[0.5,1,1];
centralization=0;

% Set Global Parameter for rotationMatrix 
A1=randn(n,n);
A=A1-A1'; %Skew-symmetric matrix  
% Bad Noise
E1=1.*eye(3);
% Expansion
ne=100;
Aexpansion=0.5.*randn(n,ne);

%% downstream neurons after expansion
sref=s0.*randn(1,m);
r0 = MixGauBrain (n,m,sref,gamma,A,E1);
r0e=r0*Aexpansion;
sb=[(s0-srange/2).*ones(1,m/2),(s0+srange/2).*ones(1,m/2)];
rb = MixGauBrain (n,m,sb,gamma,A,E1);
rbe=rb*Aexpansion;

%% Centralize rbe and re
if centralization==1
    rn_e=rbe( 1:m/2 , :);
    rp_e=rbe( m/2+1:end , :);
    Fbn_e=mean(rn_e,1);
    Fbp_e=mean(rp_e,1);
    rmeann_e=repmat(Fbn_e,m/2,1);
    rmeanp_e=repmat(Fbp_e,m/2,1);
    rmeanb_e=[rmeann_e;rmeanp_e];
    zb_e=rbe-rmeanb_e;

    tol=1e-5;
    covn_e=cov(rn_e);
    covp_e=cov(rp_e);
    z3n_e=(rn_e-rmeann_e)*sqrtm(pinv(covn_e,tol));
    z3p_e=(rp_e-rmeanp_e)*sqrtm(pinv(covp_e,tol));
    zb3_e=[z3n_e;z3p_e];



    rmean0_e=repmat(mean(r0e,1),m,1);
    z0_e=r0e-rmean0_e;
    cov0_e=cov(r0e);
    z03_e=(r0e-rmean0_e)*sqrtm(pinv(cov0_e,tol));
else
    z0_e=r0e;
    zb_e=rbe;    
    zb3_e=zb_e;
    z03_e=z0_e;
end

%% Centralize rb and r0
if centralization==1
    rmean0=repmat(mean(r0,1),m,1);
    z0=r0-rmean0;

    rn=rb( 1:m/2 , :);
    rp=rb( m/2+1:end , :);
    Fbn=mean(rn,1);
    Fbp=mean(rp,1);
    rmeann=repmat(Fbn,m/2,1);
    rmeanp=repmat(Fbp,m/2,1);
    rmeanb=[rmeann;rmeanp];
    zb=rb-rmeanb;



    tol=1e-5;
    covn=cov(rn);
    covp=cov(rp);
    z3n=(rn-rmeann)*sqrtm(pinv(covn,tol));
    z3p=(rp-rmeanp)*sqrtm(pinv(covp,tol));
    zb3=[z3n;z3p];

    cov0=cov(r0);
    z03=(r0-rmean0)*sqrtm(pinv(cov0,tol));
else
    z0=r0;
    zb=rb;
    zb3=zb;
    z03=z0;
end
    

%% Decoding total informaion in r0 before expansion. 
combin1=[1:n]';
combin2=combinator(n,2,'c','r') ;
combin3=combinator(n,3,'c','r') ;

moments=0;
[FisherInfo_all,shat_rho_all] = RedundancyDecode_MixGauBrain (r0,z0,z03,rb,zb,zb3,s0,sb,moments,combin1,combin2,combin3, n);

moments=1;
[FisherInfo1_total,shat_rho_1] = RedundancyDecode_MixGauBrain (r0,z0,z03,rb,zb,zb3,s0,sb,moments,combin1,combin2,combin3, n);

moments=2;
[FisherInfo2_total,shat_rho_2]= RedundancyDecode_MixGauBrain (r0,z0,z03,rb,zb,zb3,s0,sb,moments,combin1,combin2,combin3, size(combin2,1));

moments=3;
[FisherInfo3_total,shat_rho_3]= RedundancyDecode_MixGauBrain (r0,z0,z03,rb,zb,zb3,s0,sb,moments,combin1,combin2,combin3, size(combin3,1));

FisherInfo_all


%% Decoding informaion in re after expansion, while controling the decoding number
combin1=[1:ne]';
combin2=combinator(ne,2,'c','r') ;
combin3=combinator(ne,3,'c','r') ;

moments=1
for jj=1:n+5
    [FisherInfo1(jj),ccTheo1{jj},ccSim1{jj}] = RedundancyDecode_CCtest (r0e,z0_e,z03_e,rbe,zb_e,zb3_e,s0,sb,moments,combin1,combin2,combin3, jj,shat_rho_1);
end

figure
plot(FisherInfo1./FisherInfo1_total)
axis square


moments=2
for jj=1:n*(n+1)/2+10
    [FisherInfo2(jj),ccTheo2{jj},ccSim2{jj}] = RedundancyDecode_CCtest (r0e,z0_e,z03_e,rbe,zb_e,zb3_e,s0,sb,moments,combin1,combin2,combin3, jj,shat_rho_2);
end

moments=3
for jj=1:n*(n+1)*(n+2)/6+15
    [FisherInfo3(jj),ccTheo3{jj},ccSim3{jj}] = RedundancyDecode_CCtest (r0e,z0_e,z03_e,rbe,zb_e,zb3_e,s0,sb,moments,combin1,combin2,combin3, jj,shat_rho_3);
end



%% Plot

figure
subplot([131])
plot(FisherInfo1./FisherInfo1_total)
set(gca,'linewidth',1,'fontsize',18,'fontname','CMU Serif');
xlabel('Number of Decoded Units','FontSize',18);
ylabel('Information Ratio','FontSize',18);
axis square


subplot([132])
plot(FisherInfo2./FisherInfo2_total)
set(gca,'linewidth',1,'fontsize',18,'fontname','CMU Serif');
xlabel('Number of Decoded Units','FontSize',18);
ylabel('Information Ratio','FontSize',18);
axis square

subplot([133])
plot(FisherInfo3./FisherInfo3_total)
set(gca,'linewidth',1,'fontsize',18,'fontname','CMU Serif');
xlabel('Number of Decoded Units','FontSize',18);
ylabel('Information Ratio','FontSize',18);
axis square



%% useful functions:

function [FisherInfo,ccTheo,ccSim] = RedundancyDecode_CCtest (r0,z0,z03,rb,zb,zb3,s0,sb,moments,combin1,combin2,combin3, jj, shat_upstream)
% Function used to calculate FI when randomly chose jj units to decode
% optimally

[m,n] = size(r0);
srange=sb(end)-sb(1);

%% Generate Learning and decoding data.
if moments==1
    combinNum=size(combin1,1);
    temp=randperm(combinNum,jj);%produce random jj number as decoded linear units
    indexchose=combin1(temp,:);
    Rb=rb(:,indexchose);
    R0=r0(:,indexchose);
elseif moments==2
    combinNum=size(combin2,1);
    temp=randperm(combinNum,jj);%produce random jj number as decoded quadratic units
    combinchose=combin2(temp,:); % jj by 2
    for i=1:jj
         Rb(:,i)=zb(:,combinchose(i,1)).*zb(:,combinchose(i,2));
         R0(:,i)=z0(:,combinchose(i,1)).*z0(:,combinchose(i,2));
%           Rb(:,i)=zb(:,combin2(i,1)).*zb(:,combin2(i,2));
%           R0(:,i)=z0(:,combin2(i,1)).*z0(:,combin2(i,2));
    end
elseif moments==3
    combinNum=size(combin3,1);
    temp=randperm(combinNum,jj);%produce random jj number as decoded quadratic units
    combinchose=combin3(temp,:); % jj by 3
    for i=1:jj   
         Rb(:,i)=zb3(:,combinchose(i,1)).*zb3(:,combinchose(i,2)).*zb3(:,combinchose(i,3));
         R0(:,i)=z03(:,combinchose(i,1)).*z03(:,combinchose(i,2)).*z03(:,combinchose(i,3));
%           Rb(:,i)=zb(:,combin3(i,1)).*zb(:,combin3(i,2)).*zb(:,combin3(i,3));
%           R0(:,i)=z0(:,combin3(i,1)).*z0(:,combin3(i,2)).*z0(:,combin3(i,3));
    end
elseif moments==0% decode all Rs, only used when we decode rho
    combinNum1=size(combin1,1);
    temp1=randperm(combinNum1,jj);%produce random jj number as decoded linear units
    indexchose1=combin1(temp1,:);
    Rb1=rb(:,indexchose1);
    R01=r0(:,indexchose1);
    N2=jj*(jj+1)/2;
    combinNum2=size(combin2,1);
    temp2=randperm(combinNum2,N2);%produce random jj number as decoded quadratic units
    combinchose2=combin2(temp2,:); % jj by 2
    for i=1:N2
         Rb2(:,i)=zb(:,combinchose2(i,1)).*zb(:,combinchose2(i,2));
         R02(:,i)=z0(:,combinchose2(i,1)).*z0(:,combinchose2(i,2));
    end    
    N3=jj*(jj+1)*(jj+2)/6;
    combinNum3=size(combin3,1);
    temp3=randperm(combinNum3,N3);%produce random jj number as decoded quadratic units
    combinchose3=combin3(temp3,:); % jj by 3
    for i=1:N3   
         Rb3(:,i)=zb3(:,combinchose3(i,1)).*zb3(:,combinchose3(i,2)).*zb3(:,combinchose3(i,3));
         R03(:,i)=z03(:,combinchose3(i,1)).*z03(:,combinchose3(i,2)).*z03(:,combinchose3(i,3));
    end
    Rb=[Rb1,Rb2,Rb3];
    R0=[R01,R02,R03];
end


%% Decoding

F0=mean(Rb,1);
Rref=repmat(F0,[m,1]);
wq=pinv(Rb-Rref)*(sb-s0)';

Rbn=Rb( 1:m/2 , :);
Rbp=Rb( m/2+1:end , :);
Fb_minus=mean(Rbn,1);
Fb_plus=mean(Rbp,1);
fp=(Fb_plus-Fb_minus)./srange;

wq=wq./(fp*wq);
sq=(R0-Rref)*wq+s0;
FisherInfo=1./var(sq);
dpRb=((mean(Rbp,1)-mean(Rbn,1))./(0.5.*(std(Rbp,1,1).^2+std(Rbn,1,1).^2)).^0.5);


ccTheo=dpRb./sqrt(FisherInfo)./srange;
ccSim=corr(shat_upstream,R0);


end

function [FisherInfo,sq] = RedundancyDecode_MixGauBrain (r0,z0,z03,rb,zb,zb3,s0,sb,moments,combin1,combin2,combin3, jj)
% Function used to calculate FI when randomly chose jj units to decode
% optimally

[m,n] = size(r0);
srange=sb(end)-sb(1);

%% Generate Learning and decoding data.
if moments==1
    combinNum=size(combin1,1);
    temp=randperm(combinNum,jj);%produce random jj number as decoded linear units
    indexchose=combin1(temp,:);
    Rb=rb(:,indexchose);
    R0=r0(:,indexchose);
elseif moments==2
    combinNum=size(combin2,1);
    temp=randperm(combinNum,jj);%produce random jj number as decoded quadratic units
    combinchose=combin2(temp,:); % jj by 2
    for i=1:jj
         Rb(:,i)=zb(:,combinchose(i,1)).*zb(:,combinchose(i,2));
         R0(:,i)=z0(:,combinchose(i,1)).*z0(:,combinchose(i,2));
%           Rb(:,i)=zb(:,combin2(i,1)).*zb(:,combin2(i,2));
%           R0(:,i)=z0(:,combin2(i,1)).*z0(:,combin2(i,2));
    end
elseif moments==3
    combinNum=size(combin3,1);
    temp=randperm(combinNum,jj);%produce random jj number as decoded quadratic units
    combinchose=combin3(temp,:); % jj by 3
    for i=1:jj   
         Rb(:,i)=zb3(:,combinchose(i,1)).*zb3(:,combinchose(i,2)).*zb3(:,combinchose(i,3));
         R0(:,i)=z03(:,combinchose(i,1)).*z03(:,combinchose(i,2)).*z03(:,combinchose(i,3));
%           Rb(:,i)=zb(:,combin3(i,1)).*zb(:,combin3(i,2)).*zb(:,combin3(i,3));
%           R0(:,i)=z0(:,combin3(i,1)).*z0(:,combin3(i,2)).*z0(:,combin3(i,3));
    end
elseif moments==0% decode all Rs, only used when we decode rho
    combinNum1=size(combin1,1);
    temp1=randperm(combinNum1,jj);%produce random jj number as decoded linear units
    indexchose1=combin1(temp1,:);
    Rb1=rb(:,indexchose1);
    R01=r0(:,indexchose1);
    N2=jj*(jj+1)/2;
    combinNum2=size(combin2,1);
    temp2=randperm(combinNum2,N2);%produce random jj number as decoded quadratic units
    combinchose2=combin2(temp2,:); % jj by 2
    for i=1:N2
         Rb2(:,i)=zb(:,combinchose2(i,1)).*zb(:,combinchose2(i,2));
         R02(:,i)=z0(:,combinchose2(i,1)).*z0(:,combinchose2(i,2));
    end    
    N3=jj*(jj+1)*(jj+2)/6;
    combinNum3=size(combin3,1);
    temp3=randperm(combinNum3,N3);%produce random jj number as decoded quadratic units
    combinchose3=combin3(temp3,:); % jj by 3
    for i=1:N3   
         Rb3(:,i)=zb3(:,combinchose3(i,1)).*zb3(:,combinchose3(i,2)).*zb3(:,combinchose3(i,3));
         R03(:,i)=z03(:,combinchose3(i,1)).*z03(:,combinchose3(i,2)).*z03(:,combinchose3(i,3));
    end
    Rb=[Rb1,Rb2,Rb3];
    R0=[R01,R02,R03];
end


%% Decoding

F0=mean(Rb,1);
Rref=repmat(F0,[m,1]);
wq=pinv(Rb-Rref)*(sb-s0)';

Rbn=Rb( 1:m/2 , :);
Rbp=Rb( m/2+1:end , :);
Fb_minus=mean(Rbn,1);
Fb_plus=mean(Rbp,1);
fp=(Fb_plus-Fb_minus)./srange;

wq=wq./(fp*wq);
sq=(R0-Rref)*wq+s0;
FisherInfo=1./var(sq);


end

%%
function [FisherInfo,sq,ccTheo,ccSim] = RedundancyDecode_MixGauBrain (r0,z0,z03,rb,zb,zb3,s0,sb,moments,combin1,combin2,combin3, jj)
% Function used to calculate FI when randomly chose jj units to decode
% optimally

[m,n] = size(r0);
srange=sb(end)-sb(1);

%% Generate Learning and decoding data.
if moments==1
    combinNum=size(combin1,1);
    temp=randperm(combinNum,jj);%produce random jj number as decoded linear units
    indexchose=combin1(temp,:);
    Rb=rb(:,indexchose);
    R0=r0(:,indexchose);
elseif moments==2
    combinNum=size(combin2,1);
    temp=randperm(combinNum,jj);%produce random jj number as decoded quadratic units
    combinchose=combin2(temp,:); % jj by 2
    for i=1:jj
         Rb(:,i)=zb(:,combinchose(i,1)).*zb(:,combinchose(i,2));
         R0(:,i)=z0(:,combinchose(i,1)).*z0(:,combinchose(i,2));
%           Rb(:,i)=zb(:,combin2(i,1)).*zb(:,combin2(i,2));
%           R0(:,i)=z0(:,combin2(i,1)).*z0(:,combin2(i,2));
    end
elseif moments==3
    combinNum=size(combin3,1);
    temp=randperm(combinNum,jj);%produce random jj number as decoded quadratic units
    combinchose=combin3(temp,:); % jj by 3
    for i=1:jj   
         Rb(:,i)=zb3(:,combinchose(i,1)).*zb3(:,combinchose(i,2)).*zb3(:,combinchose(i,3));
         R0(:,i)=z03(:,combinchose(i,1)).*z03(:,combinchose(i,2)).*z03(:,combinchose(i,3));
%           Rb(:,i)=zb(:,combin3(i,1)).*zb(:,combin3(i,2)).*zb(:,combin3(i,3));
%           R0(:,i)=z0(:,combin3(i,1)).*z0(:,combin3(i,2)).*z0(:,combin3(i,3));
    end
elseif moments==0% decode all Rs, only used when we decode rho
    combinNum1=size(combin1,1);
    temp1=randperm(combinNum1,jj);%produce random jj number as decoded linear units
    indexchose1=combin1(temp1,:);
    Rb1=rb(:,indexchose1);
    R01=r0(:,indexchose1);
    N2=jj*(jj+1)/2;
    combinNum2=size(combin2,1);
    temp2=randperm(combinNum2,N2);%produce random jj number as decoded quadratic units
    combinchose2=combin2(temp2,:); % jj by 2
    for i=1:N2
         Rb2(:,i)=zb(:,combinchose2(i,1)).*zb(:,combinchose2(i,2));
         R02(:,i)=z0(:,combinchose2(i,1)).*z0(:,combinchose2(i,2));
    end    
    N3=jj*(jj+1)*(jj+2)/6;
    combinNum3=size(combin3,1);
    temp3=randperm(combinNum3,N3);%produce random jj number as decoded quadratic units
    combinchose3=combin3(temp3,:); % jj by 3
    for i=1:N3   
         Rb3(:,i)=zb3(:,combinchose3(i,1)).*zb3(:,combinchose3(i,2)).*zb3(:,combinchose3(i,3));
         R03(:,i)=z03(:,combinchose3(i,1)).*z03(:,combinchose3(i,2)).*z03(:,combinchose3(i,3));
    end
    Rb=[Rb1,Rb2,Rb3];
    R0=[R01,R02,R03];
end


%% Decoding

F0=mean(Rb,1);
Rref=repmat(F0,[m,1]);
wq=pinv(Rb-Rref)*(sb-s0)';

Rbn=Rb( 1:m/2 , :);
Rbp=Rb( m/2+1:end , :);
Fb_minus=mean(Rbn,1);
Fb_plus=mean(Rbp,1);
fp=(Fb_plus-Fb_minus)./srange;

wq=wq./(fp*wq);
sq=(R0-Rref)*wq+s0;
FisherInfo=1./var(sq);
dpRb=((mean(Rbp,1)-mean(Rbn,1))./(0.5.*(std(Rbp,1,1).^2+std(Rbn,1,1).^2)).^0.5);


ccTheo=dpRb./sqrt(FisherInfo)./srange;
ccSim=corr(sq,R0);



end

