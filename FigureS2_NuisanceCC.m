%% This code generates Figure S2: Comparing choice correlations caused by internal and external noise.
clear
close all
clc

%% Cross FI 
J12eqzero=0;
[cc_s1_theo_Jc0,cc_s1s2_sim_Jc0]=NuisanceCCsimulation(J12eqzero);

J12eqzero=1;
[cc_s1_theo_Jc1,cc_s1s2_sim_Jc1]=NuisanceCCsimulation(J12eqzero);

%% Plot
figure
subplot(121)
errorbar(cc_s1_theo_Jc0,mean(cc_s1s2_sim_Jc0,2),std(cc_s1s2_sim_Jc0,1,2),'o','MarkerSize',5); hold on;
plot([-0.4,0.4],[-0.4,0.4],'k-');
axis square
set(gca,'XTick',[-0.4:0.2:0.4])
set(gca,'YTick',[-0.4:0.2:0.4])
title('J_{12}=0');
xlabel('Predicted CC','FontSize',10);
ylabel('Simulated CC; Suboptimal decoding','FontSize',10);
set(gca,'linewidth',1,'fontsize',10,'fontname','CMU Serif');

subplot(122)
errorbar(cc_s1_theo_Jc1,mean(cc_s1s2_sim_Jc1,2),std(cc_s1s2_sim_Jc1,1,2),'o','MarkerSize',5); hold on;
plot([-0.4,0.4],[-0.4,0.4],'k-');
axis square
set(gca,'XTick',[-0.4:0.2:0.4])
set(gca,'YTick',[-0.4:0.2:0.4])
title('J_{12}=16');
xlabel('Predicted CC given s','FontSize',10);
ylabel('Simulated CC; Suboptimal decoding','FontSize',10);
set(gca,'linewidth',1,'fontsize',10,'fontname','CMU Serif');


%% Defined function
function [cc_s1_theo,cc_s1s2_sim]=NuisanceCCsimulation(J12eqzero)
    n=50;
    m=10000;
    RepeatedRounds=10;
    
    if J12eqzero==0
        count=0;
        while (count==0)
            A=rand(n);
            B=A'*A ;
            [P,D] = eig(B) ;
        if ((P'*P - eye(n))>eps) 
                count=0
        else
                f=P(1,1:n) ;
                g=P(2,1:n) ;
                count=1;
        end
        end
    else
        f=1.*rand(1,n);
        g=1.*rand(1,n);
    end
    s1=zeros(1,m);
    nuivar=2;
    s2=sqrt(nuivar).*randn(1,m);
    s2_fix=ones(1,m);
    Sigma=1.*eye(n);
    
    for jj=1:RepeatedRounds
        for j=1:m
            r_s1(j,:) = f.*s1(j)+g.*s2(j)+mvnrnd(zeros(1,n),Sigma,1);
        end

        for j=1:m
            r_s1s2(j,:) = f.*s1(j)+g.*s2_fix(j)+mvnrnd(zeros(1,n),Sigma,1);
        end
        Sigma_tot=Sigma+nuivar.*g'*g;
        J11=f*inv(Sigma)*f';
        J22=g*inv(Sigma)*g';
        J12=f*inv(Sigma)*g';
        J1=J11-J12^2/(J22+1/nuivar);
        wopt_s1=inv(Sigma_tot)*f'/J1;
        wopt_s1s2=inv(Sigma)*f'/J11;
        shat_s1=r_s1*wopt_s1;
        shat_s1s2=r_s1s2*wopt_s1; 
        % Brain doesn't know nuisance is fixed, still use same decoding weights.
        shat_s1s2_opt=r_s1s2*wopt_s1s2; % If brain use the optimal weights?
        sig_shat_sn=sqrt(1/J1-J12^2/(nuivar*J1^2*(1/nuivar+J22)^2));
        cc_s1_sim(:,jj)=corr(r_s1,shat_s1);
        cc_s1s2_sim(:,jj)=corr(r_s1s2,shat_s1s2);

        samedata=0;
        if samedata==1
            cc_s1_theo1=f'./sqrt(diag(Sigma))/sqrt(J11);
        else
            cc_s1_theo=f'./sqrt(diag(Sigma_tot))/sqrt(J1);
        end
        gamma_theo=g'./(sqrt(diag(Sigma)))*J12/(1/nuivar+J22)/J1/sig_shat_sn;
        cc_s1s2_theo=f'./(sqrt(diag(Sigma)))/J1/sig_shat_sn-gamma_theo;
        betta=sqrt(diag(Sigma_tot)./diag(Sigma))./sqrt(J1)./sig_shat_sn;
    end
end