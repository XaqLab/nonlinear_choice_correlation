clear all
close all
%% Setting
n=10;     %neuron number
s=0.5;     %stimulus reference
Informative=12;  %12 means Linear and Quadratic are both informative. 1 means linear informative. 2 means Quadratic informative.
Decoding=12;     %12 means Linear and Quadratic are both decoded. 1 means linear decoded. 2 means Quadratic decoded.
suboptimal=0;    %1 means flip sign of Fp decoding weights. 2 means forget to decode some weights. 
Times=30;

% fp=1.*randn(n,1);
% fp=1.*randn(n,1);
fp=1*[-1:2./(n-1):1]';






%% Linear Setting
u=fp.*s; % Informative linear




%% Decoding

% sampletimes=[2:2:30];
% m_samples=sampletimes.*(n*(n+1)/2+n);   %sample

m_samples=100.*[1:15];

A1=randn(n,n);
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

for kk=1:numel(m_samples)
    %% Quadratic Setting


    for jj=1:Times
        r=mvnrnd(u,sig,m_samples(kk));
        [ccSim1(jj,:),ccTheo1(jj,:),ccSimSquare(jj,:),ccTheoSquare(jj,:),ccSimCross(jj,:),ccTheoCross(jj,:),slope1(jj),slopeSquare(jj),slopeCross(jj),...
         ccSim1_2AFC(jj,:),ccTheo1_2AFC(jj,:),ccSimSquare_2AFC(jj,:),ccTheoSquare_2AFC(jj,:),ccSimCross_2AFC(jj,:),ccTheoCross_2AFC(jj,:),slope1_2AFC(jj),slopeSquare_2AFC(jj),slopeCross_2AFC(jj),...
         J_decoded(jj)] = DecodingEngine_Quadratic(Decoding,suboptimal, u, sig, sigp, fp, r, s);
    end
% ccTheo1Mean=mean((ccTheo1),1);
% ccTheoCrossMean=mean((ccTheoCross),1);
% ccTheoSquareMean=mean((ccTheoSquare),1);
% 
% ccTheo1std=std(ccTheo1,1,1);
% ccTheoCrossstd=std(ccTheoCross,1,1);
% ccTheoSquarestd=std(ccTheoSquare,1,1);
% 
% ccSim1Mean=mean((ccSim1),1);
% ccSimSquareMean=mean((ccSimSquare),1);
% ccSimCrossMean=mean((ccSimCross),1);
% 
% ccSim1std=std(ccSim1,1,1);
% ccSimSquarestd=std(ccSimSquare,1,1);
% ccSimCrossstd=std(ccSimCross,1,1);

ccTheo1_example(kk,:)=ccTheo1(1,:);
ccTheoCross_example(kk,:)=ccTheoCross(1,:);
ccTheoSquare_example(kk,:)=ccTheoSquare(1,:);

ccSim1_example(kk,:)=ccSim1(1,:);
ccSimCross_example(kk,:)=ccSimCross(1,:);
ccSimSquare_example(kk,:)=ccSimSquare(1,:);

slope1_all(:,kk)=slope1';
slopeSquare_all(:,kk)=slopeSquare';
slopeCrosss_all(:,kk)=slopeCross';


ccTheo1_example_2AFC(kk,:)=ccTheo1_2AFC(1,:);
ccTheoCross_example_2AFC(kk,:)=ccTheoCross_2AFC(1,:);
ccTheoSquare_example_2AFC(kk,:)=ccTheoSquare_2AFC(1,:);

ccSim1_example_2AFC(kk,:)=ccSim1_2AFC(1,:);
ccSimCross_example_2AFC(kk,:)=ccSimCross_2AFC(1,:);
ccSimSquare_example_2AFC(kk,:)=ccSimSquare_2AFC(1,:);

slope1_all_2AFC(:,kk)=slope1_2AFC';
slopeSquare_all_2AFC(:,kk)=slopeSquare_2AFC';
slopeCrosss_all_2AFC(:,kk)=slopeCross_2AFC';

J_pop(kk)=mean(J_decoded);

end
CI95 = tinv([0.025 0.975], Times-1);                    

slope1Mean = mean(slope1_all);                                    
slope1_SEM = std(slope1_all)/sqrt(Times);                              
slope1CI = bsxfun(@times, slope1_SEM, CI95(:));              

slopeSquareMean = mean(slopeSquare_all);                                    
slopeSquare_SEM = std(slopeSquare_all)/sqrt(Times);                           
slopeSquareCI = bsxfun(@times, slopeSquare_SEM, CI95(:));             

slopeCrossMean = mean(slopeCrosss_all);                                   
slopeCross_SEM = std(slopeCrosss_all)/sqrt(Times);                            
slopeCrossCI = bsxfun(@times, slopeCross_SEM, CI95(:));             

slope1Mean_2AFC = mean(slope1_all_2AFC); 
slope1_SEM_2AFC = std(slope1_all_2AFC)/sqrt(Times);                              
slope1CI_2AFC = bsxfun(@times, slope1_SEM_2AFC, CI95(:));              

slopeSquareMean_2AFC = mean(slopeSquare_all_2AFC);                                    
slopeSquare_SEM_2AFC = std(slopeSquare_all_2AFC)/sqrt(Times);                           
slopeSquareCI_2AFC = bsxfun(@times, slopeSquare_SEM_2AFC, CI95(:));             

slopeCrossMean_2AFC = mean(slopeCrosss_all_2AFC);                                   
slopeCross_SEM_2AFC = std(slopeCrosss_all_2AFC)/sqrt(Times);                            
slopeCrossCI_2AFC = bsxfun(@times, slopeCross_SEM_2AFC, CI95(:));             




figure
subplot(131)
patch([m_samples fliplr(m_samples)], [slope1CI(1,:)+slope1Mean fliplr(slope1CI(2,:)+slope1Mean)],[0,0.4,0.7]);
hold on;
plot(m_samples,slope1Mean,'b-');
hold on;
patch([m_samples fliplr(m_samples)], [slope1CI_2AFC(1,:)+slope1Mean_2AFC fliplr(slope1CI_2AFC(2,:)+slope1Mean_2AFC)],[0.5,0.4,0.7]);
hold on;
plot(m_samples,slope1Mean_2AFC,'k-');
set(gca,'linewidth',1,'fontsize',18,'fontname','CMU Serif');
axis square
axis([m_samples(1),m_samples(end),0.9,1.9])

subplot(132)
patch([m_samples fliplr(m_samples)], [slopeSquareCI(1,:)+slopeSquareMean fliplr(slopeSquareCI(2,:)+slopeSquareMean)],[0,0.4,0.7]);
hold on;
plot(m_samples,slopeSquareMean,'g-');
hold on;
patch([m_samples fliplr(m_samples)], [slopeSquareCI_2AFC(1,:)+slopeSquareMean_2AFC fliplr(slopeSquareCI_2AFC(2,:)+slopeSquareMean_2AFC)],[0.5,0.4,0.7]);
hold on;
plot(m_samples,slopeSquareMean_2AFC,'k-');
set(gca,'linewidth',1,'fontsize',18,'fontname','CMU Serif');
axis square
axis([m_samples(1),m_samples(end),0.9,1.9])

subplot(133)
patch([m_samples fliplr(m_samples)], [slopeCrossCI(1,:)+slopeCrossMean fliplr(slopeCrossCI(2,:)+slopeCrossMean)],[0,0.4,0.7]);
hold on;
plot(m_samples,slopeCrossMean,'r-');
hold on;
patch([m_samples fliplr(m_samples)], [slopeCrossCI_2AFC(1,:)+slopeCrossMean_2AFC fliplr(slopeCrossCI_2AFC(2,:)+slopeCrossMean_2AFC)],[0.5,0.4,0.7]);
hold on;
plot(m_samples,slopeCrossMean_2AFC,'k-');
set(gca,'linewidth',1,'fontsize',18,'fontname','CMU Serif');
axis square
axis([m_samples(1),m_samples(end),0.9,1.9])

marksize=10;
plotrange=0.3;
figure
m_plot=[2,5,10]-1;
pp=0;
for jj=1:numel(m_plot)
    kk=m_plot(jj);
    pp=pp+1;
    Pointsize=5;
    subplot(1,3,pp)
    plot(ccTheoCross_example_2AFC(kk,:),ccSimCross_example_2AFC(kk,:),'r.','markersize', marksize);hold on;
    plot(ccTheo1_example_2AFC(kk,:),ccSim1_example_2AFC(kk,:),'b.','markersize', marksize);hold on;
    plot(ccTheoSquare_example_2AFC(kk,:),ccSimSquare_example_2AFC(kk,:),'g.','markersize', marksize);hold on;
    plot([-plotrange,plotrange],[-plotrange,plotrange] ,'k-','markersize', 6);hold on;
    axis square;
    set(gca,'linewidth',1,'fontsize',18,'fontname','CMU Serif');
end

%% Useful functions
function [ccSim1,ccTheo1,ccSimSquare,ccTheoSquare,ccSimCross,ccTheoCross,slope1,slopeSquare,slopeCross,...
    ccSim1_2AFC,ccTheo1_2AFC,ccSimSquare_2AFC,ccTheoSquare_2AFC,ccSimCross_2AFC,ccTheoCross_2AFC,slope1_2AFC,slopeSquare_2AFC,slopeCross_2AFC,...
    J_decoded] = DecodingEngine_Quadratic(Decoding,suboptimal, u, sig, sigp, fp, r , s)
% here informative is used to decide if we want to use supoptimal
% weights(non-zero) for non-informative neurons. This is because
% For non-informative neurons, if we decode them with optimal weights w\prop\Sigma^{-1}f?. Because f? is zero, then the weights are zero. Thus the optimal choice correlation should also be zero(assuming non-informative neurons are independent from informative neurons.), which means they are not decoded.
% Maybe here, what we mean by saying decoding is suboptimal decoding. Even though they are not informative, we still assign non-zero weights to them. This will reduce the cc test slope of informative neurons because population information is reduced. Simulated cc for non-informative neurons should be non-zero, whereas theoretical cc is zero.

n=size(r,2);
m=size(r,1);
isig=inv(sig);
z3=r-repmat(u,[1,m])'; %  Get pure quadratic statistics.


F=zeros(n*(n+1)/2,1);
Fp=zeros(n*(n+1)/2,1);
R2=zeros(m,n*(n+1)/2);
ij=0;
for i=1:n
    for j=i:n
        ij=ij+1;
        R2(:,ij)=z3(:,i).*z3(:,j);
        Fp(ij) =sigp(i,j);
        F(ij) = sig(i,j); 
    end
end
  
ij=0;
for i = 1:n
    for j = i:n
        ij=ij+1;
        kl=0;
        for k=1:n
            for l=k:n
                kl=kl+1;
                gamma(ij,kl) = sig(i,k)*sig(j,l) + sig(j,k)*sig(i,l); % with mean removed, gives cov(sig)
            end
        end
    end
end
gamma = (gamma+gamma')/2.; % enforce symmetry
igamma=inv(gamma);

J2=Fp'*igamma*Fp %Fisher information
J1=fp'*isig*fp

nonfliprate=0.8;

forgetrate=0.7;



if Decoding==12
    J=J1+J2;
    wopt2=igamma*Fp./J;
    wopt1=isig*fp./J;
    if suboptimal==1
        randsign1=binornd(1,nonfliprate,n,1).*2-1;
        randsign2=binornd(1,nonfliprate,n*(n+1)/2,1).*2-1;
        wopt1=(fp).*randsign1;
        wopt2=(Fp).*randsign2;
        NormalizeTerm=wopt1'*fp+wopt2'*Fp;
        wopt1=wopt1./NormalizeTerm;
        wopt2=wopt2./NormalizeTerm;
    elseif suboptimal==2
        randsign1=binornd(1,forgetrate,n,1);
        randsign2=binornd(1,forgetrate,n*(n+1)/2,1);
        
        nonzero_entry2=find(randsign2);
        gamma_forget=gamma(nonzero_entry2,nonzero_entry2);
        igamma_f=inv(gamma_forget);
        Fp_f=Fp(nonzero_entry2);
        
        nonzero_entry1=find(randsign1);
        sig_forget=sig(nonzero_entry1,nonzero_entry1);
        isig_f=inv(sig_forget);
        fp_f=fp(nonzero_entry1);
        
        wopt2_f=igamma_f*Fp_f;
        wopt1_f=isig_f*fp_f;
        NormalizeTerm=wopt1_f'*fp_f+wopt2_f'*Fp_f;
        wopt1_f=wopt1_f./NormalizeTerm;
        wopt2_f=wopt2_f./NormalizeTerm;
        
        wopt1=zeros(n,1);
        wopt2=zeros(n*(n+1)/2,1);
        
        wopt1(nonzero_entry1)=wopt1_f;
        wopt2(nonzero_entry2)=wopt2_f;
    end
elseif Decoding==1
    J=J1;
    wopt1=isig*fp./J;
    wopt2=zeros(n*(n+1)/2,1);
elseif Decoding==2
    J=J2;
    wopt1=zeros(n,1);
    wopt2=igamma*Fp./J;
end



shopt=wopt2'*(R2'-repmat(F,[1,m])) + wopt1'*(r'-repmat(u,[1,m]))  + s;
shopt_binary=((shopt>s)-0.5).*2;

J_decoded=1./var(shopt);



gamma_diag=diag(gamma);
sig_diag=diag(sig);
ccSim2 = corr(shopt',R2);
ccTheo2 = Fp./sqrt(gamma_diag)/sqrt(J_decoded);
ccTheo2=ccTheo2';
ccSim1=corr(shopt',r);
ccTheo1=fp./sqrt(sig_diag) /sqrt(J_decoded);
ccTheo1=ccTheo1';


ccSim2_2AFC = corr(shopt_binary',R2);
ccTheo2_2AFC = ccTheo2.*0.76;
ccSim1_2AFC=corr(shopt_binary',r);
ccTheo1_2AFC=ccTheo1.*0.86;

% for j=1:size(R_all_Sn,2)
%     cov_mat_Sp=cov(Choice_Sp,R_all_Sp(:,j));
%     cov_mat_Sn=cov(Choice_Sn,R_all_Sn(:,j));
%     cov_ave=trialNumRatio_Sp.*cov_mat_Sp(1,2)+trialNumRatio_Sn.*cov_mat_Sn(1,2);
%     varR_ave=trialNumRatio_Sp.*var(R_all_Sp(:,j))+trialNumRatio_Sn.*var(R_all_Sn(:,j));    
%     ccSim(j)=(cov_ave)./sqrt(varR_ave.*varShat_ave);
% end




k=0;
p=0;
ij=0;
for i=1:n
    for j=i:n
        ij=ij+1;
        if i==j
            k=k+1;
            ccTheoSquare(k)=ccTheo2(ij);
            ccSimSquare(k)=ccSim2(ij);
            ccTheoSquare_2AFC(k)=ccTheo2_2AFC(ij);
            ccSimSquare_2AFC(k)=ccSim2_2AFC(ij);
        else
            p=p+1;
            ccTheoCross(p)=ccTheo2(ij);
            ccSimCross(p)=ccSim2(ij);
            ccTheoCross_2AFC(p)=ccTheo2_2AFC(ij);
            ccSimCross_2AFC(p)=ccSim2_2AFC(ij);

        end
        
    end
end


[slope1] = Fitting_cc (ccTheo1,ccSim1);
[slopeSquare] = Fitting_cc (ccTheoSquare,ccSimSquare);
[slopeCross] = Fitting_cc (ccTheoCross,ccSimCross);

[slope1_2AFC] = Fitting_cc (ccTheo1_2AFC,ccSim1_2AFC);
[slopeSquare_2AFC] = Fitting_cc (ccTheoSquare_2AFC,ccSimSquare_2AFC);
[slopeCross_2AFC] = Fitting_cc (ccTheoCross_2AFC,ccSimCross_2AFC);


end

function [slopePCA] = Fitting_cc (ccTheo,ccSim)
ccdat=[ccTheo;ccSim];
coeff = pca(ccdat,2);
slopePCA = abs(coeff(1,2) / coeff(1,1));

% dlm = fitlm(ccTheo,ccSim,'Intercept',false);
% slopePCA=dlm.Coefficients.Estimate;

end
