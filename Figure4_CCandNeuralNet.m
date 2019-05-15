clear
close all
clc
%% Basic setting
n=18; %neuron number
s0=1;       
srange=0.5;
m=20000;   %sample
Amplitude=4;
std_train=0.4;
N_Relu=30;
gamma=[0.4,0.3,1];
% Set Global Parameter for rotationMatrix 
A1=randn(6);
A=A1-A1'; %Skew-symmetric matrix  
% Bad Noise
E1=0; % no bad noise
DecodingTpye=1;
% 1: Random expansion, then ReLu, then Linear Regression
% 2: Random expansion, then tanh, then Linear Regression
% E1=2.*eye(3); % no bad noise

%% Reference responses and statistics
sref=s0.*ones(1,m);
r0 = MixGauBrain_SubGroup (n,m,sref,gamma,A,E1);
[R01,R02,R03] = rMomentsGenerate_SubGroup (r0,0);
R0=[R01,R02,R03];
%% Binary input to estimate d prime and ccTheo
sbn=s0-srange/2;
sbp=s0+srange/2;
sb=[sbn.*ones(1,m/2),sbp.*ones(1,m/2)];
rb = MixGauBrain_SubGroup (n,m,sb,gamma,A,E1);
[Rb1,Rb2,Rb3] = rMomentsGenerate_SubGroup (rb,1);
Rb=[Rb1,Rb2,Rb3];
Rbn=Rb( 1:m/2 , :);
Rbp=Rb( m/2+1:end , :);
Fb_minus=mean(Rbn,1);
Fb_plus=mean(Rbp,1);
fp=(Fb_plus-Fb_minus)./srange;

%% Decode with Polynomial until cubic statistics
[sbopt,sopt,~,wopt] = DecodingEngine(R0,Rb,s0,sb);
FisherInfo_opt=1/var(sopt);

%% Decode with random ReLu basis
strain=s0+std_train.*randn(1,m);
rtrain = MixGauBrain_SubGroup (n,m,strain,gamma,A,E1);
net = feedforwardnet([N_Relu,N_Relu]);
net.layers{1}.transferFcn = 'poslin'; % Log-sigmoid transfer function
net.layers{2}.transferFcn = 'poslin'; % Log-sigmoid transfer function

%% Create a neural network estimator.
net_estimator = neuralnet(net);
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
[net, tr] = train(net, rtrain', strain);
strainhat_Relu = net(rtrain');
errors = gsubtract(strain, strainhat_Relu);
% performance = perform(net, targets, outputs)

%%
shat_Relu=net(r0');
shat_Relu=shat_Relu';
FisherInfo=1/var(shat_Relu);
Info_ratio=FisherInfo/FisherInfo_opt;

% Binary input to estimate d prime and ccTheo
Rbn=Rb( 1:m/2 , :);
Rbp=Rb( m/2+1:end , :);
Fb_minus=mean(Rbn,1);
Fb_plus=mean(Rbp,1);
dpRb=((mean(Rbp,1)-mean(Rbn,1))./(0.5.*(std(Rbp,1,1).^2+std(Rbn,1,1).^2)).^0.5);
ccTheo=dpRb./sqrt(FisherInfo)./srange;
ccTheo1=ccTheo(1:size(Rb1,2));
ccTheo2=ccTheo(size(Rb1,2)+1:size(Rb1,2)+size(Rb2,2));
ccTheo3=ccTheo(size(Rb1,2)+size(Rb2,2)+1 : size(Rb1,2)+size(Rb2,2)+size(Rb3,2));
ccTheo_opt=dpRb./sqrt(FisherInfo_opt)./srange;
ccTheo1_opt=ccTheo_opt(1:size(Rb1,2));
ccTheo2_opt=ccTheo_opt(size(Rb1,2)+1:size(Rb1,2)+size(Rb2,2));
ccTheo3_opt=ccTheo_opt(size(Rb1,2)+size(Rb2,2)+1 : size(Rb1,2)+size(Rb2,2)+size(Rb3,2));
ccSim=corr(shat_Relu,R0);
ccSim1=ccSim(1:n);
ccSim2=ccSim(n+1:n+size(R02,2));
ccSim3=ccSim(n+size(R02,2)+1:n+size(R02,2)+size(R03,2));
ccSimopt=corr(sopt,R0);
ccSim1opt=ccSimopt(1:n);
ccSim2opt=ccSimopt(n+1:n+size(R02,2));
ccSim3opt=ccSimopt(n+size(R02,2)+1:n+size(R02,2)+size(R03,2));

%% Plot
plotscale=0.2;
figure
subplot(121)
plot(ccTheo3,ccSim3,'g.','markersize', 10);hold on;
plot(ccTheo2,ccSim2,'r.','markersize', 10);hold on;
plot(ccTheo1,ccSim1,'b.','markersize', 10);
legend('Cubic','Quadratic','Linear');
hold on;
plot([-plotscale,plotscale],[-plotscale,plotscale],'k-');
axis square
set(gca,'XTick',[-plotscale:0.2:plotscale])
set(gca,'YTick',[-plotscale:0.2:plotscale])
set(gca,'linewidth',1,'fontsize',18,'fontname','CMU Serif');
title('ReLu Neural Nets Decoder','FontSize',18);
xlabel('Optimal CC','FontSize',18);
ylabel('Simulated CC','FontSize',18);

subplot(122)
plot(ccTheo3_opt,ccSim3opt,'g.','markersize', 10);hold on;
plot(ccTheo2_opt,ccSim2opt,'r.','markersize', 10);hold on;
plot(ccTheo1_opt,ccSim1opt,'b.','markersize', 10);
legend('Cubic','Quadratic','Linear');
hold on;
plot([-plotscale,plotscale],[-plotscale,plotscale],'k-');
axis square
set(gca,'XTick',[-plotscale:0.2:plotscale])
set(gca,'YTick',[-plotscale:0.2:plotscale])
title('Optimal Polinomial Decoder','FontSize',18);
set(gca,'linewidth',1,'fontsize',18,'fontname','CMU Serif');
xlabel('Optimal CC','FontSize',18);
ylabel('Simulated CC','FontSize',18);



%% Defined function
function r_all = MixGauBrain_SubGroup (n,m,s,gamma,A,E1)
% This function generates m trials of n-d data whose informative moments
% are cubic. 
%With Affine transformation, Quadratic and linear information can also be embedded. 
% Can also be adapted to bad noise case.
% n is neuron number, m is sample number, s is 1 by m vectors. gamma
% controls moments sensitivity.
if E1==0
    s1=s;
    s2=s;
    s3=s;
else
    ds= mvnrnd(zeros(3,1),E1,m)';
    s1=s+ds(1,:);
    s2=s+ds(1,:);
    s3=s+ds(1,:);
end


%% Setting

% Linear Setting
gamma1=gamma(1); % When gamma1 is zero, no linear translation.


% Quadratic Setting
gamma2=gamma(2); % When gamma2 is zero, rotationMatrix is Identity

% Cubic Setting
binary = permn([-1 1],3);
va=binary(prod(binary,2)==1,:);% Generate 4 Lobe directions.
LobeNum=size(va,1);
gamma3=gamma(3);

%% Sampling
% In order to have more neurons, we group 6 neurons in a group and do the transformation only
% within that subgroup. Then we only need (n/6)*6^3 nonlinear units
r_all=double.empty;

for j=1:n/6
    for i=1:m
        nsub=6;
        Lobe = randi([1,LobeNum],1); % Randomly chose one lobe.
        [S,mean] = MixGauMoments (nsub,s3(i),gamma3,va(Lobe,:)); 
        Schose=squeeze(S);
        rcubic(i,:)= mvnrnd(mean,Schose,1);
        Brotate=expm(gamma2.*A.*s2(i)); %rotation matrix, controled by stimulus,gamma2 Control sensitivity
        Bscale=eye(nsub)+gamma2.*s2(i).*diag([1:nsub]./nsub);
        r(i,:)=gamma1.*s1(i).*([1:nsub]./nsub)+rcubic(i,:)*Bscale*Brotate;
    end
    r_all=[r_all,r];
end
end




%% Defined function
function [M, I] = permn(V, N, K)
% PERMN - permutations with repetition
%   Using two input variables V and N, M = PERMN(V,N) returns all
%   permutations of N elements taken from the vector V, with repetitions.
%   V can be any type of array (numbers, cells etc.) and M will be of the
%   same type as V.  If V is empty or N is 0, M will be empty.  M has the
%   size numel(V).^N-by-N. 
%
%   When only a subset of these permutations is needed, you can call PERMN
%   with 3 input variables: M = PERMN(V,N,K) returns only the K-ths
%   permutations.  The output is the same as M = PERMN(V,N) ; M = M(K,:),
%   but it avoids memory issues that may occur when there are too many
%   combinations.  This is particulary useful when you only need a few
%   permutations at a given time. If V or K is empty, or N is zero, M will
%   be empty. M has the size numel(K)-by-N. 
%
%   [M, I] = PERMN(...) also returns an index matrix I so that M = V(I).
% tested in Matlab 2016a
% version 6.1 (may 2016)
% (c) Jos van der Geest
% Matlab File Exchange Author ID: 10584
% email: samelinoa@gmail.com

narginchk(2,3) ;

if fix(N) ~= N || N < 0 || numel(N) ~= 1 ;
    error('permn:negativeN','Second argument should be a positive integer') ;
end
nV = numel(V) ;

if nargin==2, % PERMN(V,N) - return all permutations
    
    if nV==0 || N == 0,
        M = zeros(nV,N) ;
        I = zeros(nV,N) ;
        
    elseif N == 1,
        % return column vectors
        M = V(:) ;
        I = (1:nV).' ;
    else
        % this is faster than the math trick used for the call with three
        % arguments.
        [Y{N:-1:1}] = ndgrid(1:nV) ;
        I = reshape(cat(N+1,Y{:}),[],N) ;
        M = V(I) ;
    end
else % PERMN(V,N,K) - return a subset of all permutations
    nK = numel(K) ;
    if nV == 0 || N == 0 || nK == 0
        M = zeros(numel(K), N) ;
        I = zeros(numel(K), N) ;
    elseif nK < 1 || any(K<1) || any(K ~= fix(K))
        error('permn:InvalidIndex','Third argument should contain positive integers.') ;
    else
        
        V = reshape(V,1,[]) ; % v1.1 make input a row vector
        nV = numel(V) ;
        Npos = nV^N ;
        if any(K > Npos)
            warning('permn:IndexOverflow', ...
                'Values of K exceeding the total number of combinations are saturated.')
            K = min(K, Npos) ;
        end
             
        % The engine is based on version 3.2 with the correction
        % suggested by Roger Stafford. This approach uses a single matrix
        % multiplication.
        B = nV.^(1-N:0) ;
        I = ((K(:)-.5) * B) ; % matrix multiplication
        I = rem(floor(I),nV) + 1 ;
        M = V(I) ;
    end
end
