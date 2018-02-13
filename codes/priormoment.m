

clc
clear
close all


POLICY = 0; % 0 = discretion, 1 = commitment
addpath(genpath('func/'))
%% Prior Specifcations:
%--------------------------------------------------------------------------
% pshape: 0 is point mass, both para and p2 are ignored
%         1 is BETA(mean,stdd)
%         2 is GAMMA(mean,stdd)
%         3 is NORMAL(mean,stdd)
%         4 is INVGAMMA(s^2,nu)
%         5 is UNIFORM [p1,p2]
% para: Theta draws
% p1, p2 are the mean and std. of prior distribution except for uniform
%     pshape    p1      p2
pri=[   1       0.99     0.000001;    % beta ~ BETA
        1       0.45     0.000001;      % alpha ~ BETA
        1       0.6      0.2;      % h ~ BETA      % mean is .86
        2       1     0.5;           % sigma ~ GAMMA
        3       1.5     0.25;   % phi ~ NORMAL
        2       1       0.5;      % eta ~ GAMMA
        1       .5     .2;      % thetaH ~ BETA
        1       .5     .2;      % thetaF ~ BETA
        2       .5     .2;    % a1 ~ BETA
        3       0     0.00000001;    % a2 ~ NORMAL
        3       0    0.000000001;    % a3 ~ NORMAL
        3       0     0.00000001;    % b1 ~ NORMAL
        2       .5     .2;    % b2 ~ BETA
        3       0   0.00000001;    % b3 ~ NORMAL
        3       0    0.00000001;    % c1 ~ NORMAL
        3       0    0.000000001;    % c2 ~ NORMAL
        2       .5     .2;    % c3 ~ BETA
        1       0.0001     0.00002;    % rhoh ~ BETA
        1       0.0001     0.00002;    % rhof ~ BETA
        1       0.5     0.2;    % rhoa ~ BETA
        1       0.5     0.2;    % rhoq ~ BETA
        1       0.5     0.2;    % rhos ~ BETA
        1       0.5     0.2;    % rhor ~ BETA
        2       0.5       0.3;    % mu_q ~ GAMMA   %check genTheta.m, set ne*wTheta(25)=0
        2       0.5       0.3;    % mu_y ~ GAMMA
        2       0.5       0.3;    % mu_r ~ GAMMA
        4       1       .4;    % sigmah ~ GAMMA
        4       1       .4;    % sigmaf ~ GAMMA
        4       1       .4;    % sigmaz ~ GAMMA
        4       1       .4;    % sigmaq ~ GAMMA
        4       1       .4;    % sigmas ~ GAMMA
        4       1       .4;    % sigmapi* ~ GAMMA
        4       1       .4;    % sigmay* ~ GAMMA
        4       1       .4;    % sigmar* ~ GAMMA
        4       1       .4];   % sigmar ~ GAMMA
            
 
                            
% Truncations and bounds
point_b = [-100 100];  % bounds for point mass
bet_b = [0 1];      % bounds for BETA density
gam_b = [0 100];    % bounds for GAMMA density
nor_b = [-10 10];    % bounds for NORMAL density

bounds=zeros(length(pri),2);
for i=1:length(pri)
    if pri(i,1)==0
        bounds(i,:)=point_b;
    elseif pri(i,1)==1
        bounds(i,:)=bet_b;
    elseif pri(i,1)==2
        bounds(i,:)=gam_b;
    elseif pri(i,1)==3
        bounds(i,:)=nor_b;
    else
        bounds(i,:)=[-inf inf];
    end
end
pshape = pri(:,1);
p1 = pri(:,2);
p2 = pri(:,3);
lb = bounds(:,1);
ub = bounds(:,2);
% 
% priordens = priorsim(pshape,p1,p2,N);


% load chain/mh_dis
% load Theta_pr

VARSELECT = [2,3,6,4,9,7,21];
GDPINDEX = 6;

ynames = [  'Consumption     ';
            'Dom. Inflation  ';
            'For. Inflation  ';
            'Real Ex. Rate   ';
            'ToT             ';
            'Output          ';
            'CPI Inflation   ';
            'r_{t} - r_{t-1} ';
            'r_{t}           ';
            '\epsilon_{a}    ';
            '\epsilon_{H}    ';
            '\epsilon_{F}    ';
            '\epsilon_{q}    ';
            '\epsilon_{s}    ';
            '\epsilon_{\pi*} ';
            '\epsilon_{y*}   ';
            '\epsilon_{r*}   ';
            '\epsilon_{r}    ';
			'CPI Inflation -1';      % added Tk, Feb 7, 2006
			'CPI Inflation -2';		 % added Tk, Feb 7, 2006
			'Annual CPI Inf  '];

% Theta_s = priordens;
% NT = size(Theta_s,1);
% Thetamat = Theta_s((NT+rem(NT,2))/2 + 1:end,:);
%     
% % Cartman: "Thin out their numbers!"
% 
% Nthin = 1; % Must be < N
% N = size(Thetamat,1);
% 
% count = N/(Nthin + rem(Nthin,2));
% 
% nselect = count*[1:Nthin];
% 
% Theta_mat = Thetamat(nselect,:);
% NT = size(Theta_mat,1);

h = waitbar(0,'Solving models for thinned parameter MCs...');
NT =200;
n=1;
sigyy0 = zeros(22,22,NT);
sigxx0 = zeros(22,22,NT);
sigxy0 = zeros(22,22,NT);
sigyx0 = zeros(22,22,NT);
sigyy1 = zeros(22,22,NT);
sigxx1 = zeros(22,22,NT);
sigxy1 = zeros(22,22,NT);
sigyx1 = zeros(22,22,NT);
sigyy2 = zeros(22,22,NT);
sigxx2 = zeros(22,22,NT);
sigxy2 = zeros(22,22,NT);
sigyx2 = zeros(22,22,NT);
indet = 0;
ttdraws = 0;
while n <= NT
    Theta=priorsim(pshape,p1,p2,1);
    ttdraws = ttdraws +1;
    % General solution: y(t+1) = AA*y(t) + BB*z(t)
    [ AA,BB,PROBLEM ] = model_solve(Theta,POLICY);
    
    if sum(sum(isnan(AA))) ~=0 | PROBLEM==1
        indet = indet +1;
        
    else
    % Compute covariogram
   
    Lag0 = 0;
    Lag1 = 1;
    Lag2 = 2;
    [sigyy0(:,:,n),sigxx0(:,:,n),sigxy0(:,:,n),sigyx0(:,:,n)] = mom(0,AA,BB*BB',Lag0);
    [sigyy1(:,:,n),sigxx1(:,:,n),sigxy1(:,:,n),sigyx1(:,:,n)] = mom(0,AA,BB*BB',Lag1);
    [sigyy2(:,:,n),sigxx2(:,:,n),sigxy2(:,:,n),sigyx2(:,:,n)] = mom(0,AA,BB*BB',Lag2);
    n=n+1;
    end
    % Compute similar for Covariances
waitbar(n/NT,h,[sprintf('%0.2f',n*100/NT),' percent done solving... '])
end % for n
close(h)


nv = length(VARSELECT);
Nthin=1;

for i = 1:nv
    var_xx(:,i) = sigxx0(VARSELECT(i),VARSELECT(i),:);
    cov_xy(:,i) = sigxx0(VARSELECT(i),GDPINDEX,:);
    autocov_one(:,i) = sigxx1(VARSELECT(i),VARSELECT(i),:);
    autocov_two(:,i) = sigxx2(VARSELECT(i),VARSELECT(i),:);
end

stdv = sqrt(var_xx);
corry = cov_xy./(stdv.*(stdv(:,GDPINDEX)*ones(1,size(stdv,2))));
auto_one = autocov_one./var_xx;
auto_two = autocov_two./var_xx;

% BUSINESS CYCLE STATISTICS:


% S.d. of variables
figure
for i = 1:nv
subplot((nv+rem(nv,2))/2,2,i)
[f,xi] = ksdensity(stdv(:,i),'support',[0, 1000]);
plot(xi,f);
xlim([0 10])
ylabel(ynames(VARSELECT(i),:))
%suptitle2('Standard deviation',18)
end

% 1st order autocorrelation of variables
figure
for i = 1:nv
subplot((nv+rem(nv,2))/2,2,i)
[f,xi] = ksdensity(auto_one(:,i),'support',[-1, 1]);
plot(xi,f);
ylabel(ynames(VARSELECT(i),:))
%suptitle2('First order autocorrelation',18)
end

% 2nd order autocorrelation of variables
figure
for i = 1:nv
subplot((nv+rem(nv,2))/2,2,i)
[f,xi] = ksdensity(auto_two(:,i),'support',[-1, 1]);
plot(xi,f);
ylabel(ynames(VARSELECT(i),:))
%suptitle2('First order autocorrelation',18)
end

% 1st order correlation of variables with GDP
% figure
% for i = 1:nv
% subplot((nv+rem(nv,2))/2,2,i)
% [f,xi] = ksdensity(corry(:,i),'support',[-1, 1]);
% plot(xi,f);
% title(ynames(VARSELECT(i),:))
% end





