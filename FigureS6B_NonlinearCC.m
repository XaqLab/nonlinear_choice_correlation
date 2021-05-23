%% This code generates Figure S6B: Nonlinear choice correlations for recordings from a simulated optimal decoder under the same conditions. 

clear
close all

%% Parameter Setting
n=15;
s0=0;     
srange=0.5;
m=100000;   %sample number
gamma=[0.5,1,1];
centralization=1;

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
    

%% Decoding rho's information from the size more than worthwized number units. 
combin1_all=[1:ne]';

decoding_unit_num=n+2;
temp=randperm(ne,decoding_unit_num);%produce random jj number as decoded linear units
indexchose=combin1_all(temp,:);
ii=0;
for jj=1:decoding_unit_num
    for kk=1:decoding_unit_num
        ii=ii+1;
        combin2(ii,:)=[indexchose(jj),indexchose(kk)];
    end
end

ii=0;
for jj=1:decoding_unit_num
    for kk=1:decoding_unit_num
        for ll=1:decoding_unit_num
            ii=ii+1;
            combin3(ii,:)=[indexchose(jj),indexchose(kk),indexchose(ll)];
        end
    end
end


moments=1;
[FisherInfo1_total,shat_rho_1] = RedundancyDecode_MixGauBrain (r0e,z0_e,z03_e,rbe,zb_e,zb3_e,s0,sb,moments,indexchose);

moments=2;
[FisherInfo2_total,shat_rho_2]= RedundancyDecode_MixGauBrain (r0e,z0_e,z03_e,rbe,zb_e,zb3_e,s0,sb,moments,combin2);

moments=3;
[FisherInfo3_total,shat_rho_3]= RedundancyDecode_MixGauBrain (r0e,z0_e,z03_e,rbe,zb_e,zb3_e,s0,sb,moments,combin3);



%% Decoding informaion in re after expansion, while controling the decoding number

ev_list1=[round(n/3),round(n/2),n];


for kk=1:3
    
    temp=randperm(ne,ev_list1(kk));%produce random jj number as decoded linear units
    indexchose_cc=combin1_all(temp,:);
    ii=0;
    for jj=1:ev_list1(kk)
        for dd=1:ev_list1(kk)
            ii=ii+1;
            combin2_cc(ii,:)=[indexchose_cc(jj),indexchose_cc(dd)];
        end
    end

    ii=0;
    for jj=1:ev_list1(kk)
        for dd=1:ev_list1(kk)
            for ll=1:ev_list1(kk)
                ii=ii+1;
                combin3_cc(ii,:)=[indexchose_cc(jj),indexchose_cc(dd),indexchose_cc(ll)];
            end
        end
    end  
    [FisherInfo1(kk),ccTheo1{kk},ccSim1{kk}] = RedundancyDecode_CCtest (r0e,z0_e,z03_e,rbe,zb_e,zb3_e,s0,sb,1,indexchose_cc,shat_rho_1);
    [FisherInfo2(kk),ccTheo2{kk},ccSim2{kk}] = RedundancyDecode_CCtest (r0e,z0_e,z03_e,rbe,zb_e,zb3_e,s0,sb,2,combin2_cc,shat_rho_2);
    [FisherInfo3(kk),ccTheo3{kk},ccSim3{kk}] = RedundancyDecode_CCtest (r0e,z0_e,z03_e,rbe,zb_e,zb3_e,s0,sb,3,combin3_cc,shat_rho_3);
end



%% Plot
ev_list1=[3,2,1];
ev_list2=ev_list1;
ev_list3=ev_list1;

figure
markersize=20;
subplot(131)
plotscale=1;
plot(ccTheo1{ev_list1(1)},ccSim1{ev_list1(1)},'r.','markersize', markersize);hold on;
plot(ccTheo1{ev_list1(2)},ccSim1{ev_list1(2)},'g.','markersize', markersize);hold on;
plot(ccTheo1{ev_list1(3)},ccSim1{ev_list1(3)},'b.','markersize', markersize);hold on;


hold on;
plot([-plotscale,plotscale],[-plotscale,plotscale],'k-');
axis square
set(gca,'XTick',[-plotscale:0.4:plotscale])
set(gca,'YTick',[-plotscale:0.4:plotscale])
set(gca,'linewidth',1,'fontsize',18,'fontname','CMU Serif');
xlabel('Optimal CC','FontSize',18);
ylabel('Simulated CC','FontSize',18);

for i = 1:numel(ev_list1)
    legend_str{i} = ['recorded unit:' num2str(ev_list1(i))];
end
legend(legend_str)
box off

subplot(132)
plotscale=1;
plot(ccTheo2{ev_list2(1)},ccSim2{ev_list2(1)},'r.','markersize', markersize);hold on;
plot(ccTheo2{ev_list2(2)},ccSim2{ev_list2(2)},'g.','markersize', markersize);hold on;
plot(ccTheo2{ev_list2(3)},ccSim2{ev_list2(3)},'b.','markersize', markersize);hold on;


hold on;
plot([-plotscale,plotscale],[-plotscale,plotscale],'k-');
axis square
set(gca,'XTick',[-plotscale:0.4:plotscale])
set(gca,'YTick',[-plotscale:0.4:plotscale])
set(gca,'linewidth',1,'fontsize',18,'fontname','CMU Serif');
xlabel('Optimal CC','FontSize',18);
ylabel('Simulated CC','FontSize',18);
for i = 1:numel(ev_list1)
    legend_str{i} = ['recorded unit:' num2str(ev_list2(i))];
end
legend(legend_str)

box off

subplot(133)
plotscale=0.6;
plot(ccTheo3{ev_list3(1)},ccSim3{ev_list3(1)},'r.','markersize', markersize);hold on;
plot(ccTheo3{ev_list3(2)},ccSim3{ev_list3(2)},'g.','markersize', markersize);hold on;
plot(ccTheo3{ev_list3(3)},ccSim3{ev_list3(3)},'b.','markersize', markersize);hold on;

hold on;
plot([-plotscale,plotscale],[-plotscale,plotscale],'k-');
axis square
set(gca,'XTick',[-plotscale:0.2:plotscale])
set(gca,'YTick',[-plotscale:0.2:plotscale])
set(gca,'linewidth',1,'fontsize',18,'fontname','CMU Serif');
xlabel('Optimal CC','FontSize',18);
ylabel('Simulated CC','FontSize',18);
for i = 1:numel(ev_list1)
    legend_str{i} = ['recorded unit:' num2str(ev_list3(i))];
end
legend(legend_str)
box off


%% useful functions:
function [slopePCA] = Fitting_cc (ccTheo,ccSim)
ccdat=[ccTheo;ccSim];
coeff = pca(ccdat,2);
slopePCA = abs(coeff(1,2) / coeff(1,1));

end

function [FisherInfo,ccTheo,ccSim] = RedundancyDecode_CCtest (r0,z0,z03,rb,zb,zb3,s0,sb,moments,combin1_chose, shat_upstream)
% Function used to calculate FI when randomly chose jj units to decode
% optimally

[m,n] = size(r0);
srange=sb(end)-sb(1);

%% Generate Learning and decoding data.
if moments==1
    Rb=rb(:,combin1_chose);
    R0=r0(:,combin1_chose);
elseif moments==2
    for i=1:size(combin1_chose,1)
         Rb(:,i)=zb(:,combin1_chose(i,1)).*zb(:,combin1_chose(i,2));
         R0(:,i)=z0(:,combin1_chose(i,1)).*z0(:,combin1_chose(i,2));
    end
elseif moments==3
    for i=1:size(combin1_chose,1)   
         Rb(:,i)=zb3(:,combin1_chose(i,1)).*zb3(:,combin1_chose(i,2)).*zb3(:,combin1_chose(i,3));
         R0(:,i)=z03(:,combin1_chose(i,1)).*z03(:,combin1_chose(i,2)).*z03(:,combin1_chose(i,3));
    end
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
FisherInfo_behavior=1./var(shat_upstream);
FisherInfo=1./var(sq);
dpRb=((mean(Rbp,1)-mean(Rbn,1))./(0.5.*(std(Rbp,1,1).^2+std(Rbn,1,1).^2)).^0.5);


ccTheo=dpRb./sqrt(FisherInfo_behavior)./srange;
ccSim=corr(shat_upstream,R0);


end

function [FisherInfo,sq] = RedundancyDecode_MixGauBrain (r0,z0,z03,rb,zb,zb3,s0,sb,moments,combinchose)
% Function used to calculate FI when randomly chose jj units to decode
% optimally

[m,n] = size(r0);
srange=sb(end)-sb(1);

%% Generate Learning and decoding data.
if moments==1
    Rb=rb(:,combinchose);
    R0=r0(:,combinchose);
elseif moments==2
    for i=1:size(combinchose,1)
         Rb(:,i)=zb(:,combinchose(i,1)).*zb(:,combinchose(i,2));
         R0(:,i)=z0(:,combinchose(i,1)).*z0(:,combinchose(i,2));
    end
elseif moments==3
    for i=1:size(combinchose,1)  
         Rb(:,i)=zb3(:,combinchose(i,1)).*zb3(:,combinchose(i,2)).*zb3(:,combinchose(i,3));
         R0(:,i)=z03(:,combinchose(i,1)).*z03(:,combinchose(i,2)).*z03(:,combinchose(i,3));
    end
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



