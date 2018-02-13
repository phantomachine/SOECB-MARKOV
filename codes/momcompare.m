% momcompare.m
% Compute moments for sample draws from stationary MCMC distribution.
% WARNING: You need to have run COUNTRY_estimate.m first or have the MCMC
% chain saved in the folder "chain".
% 
clc
clear all

CHECKPRIOR = 0;

POLICY = 0; % 0 = discretion, 1 = commitment
addpath(genpath('func/'))
load chain/mh_dis  % load saved posterior parameter chain/distribution

if CHECKPRIOR
    load Theta_pr  % load prior parameter distribution
end

load datamat  % load actual data

VARSELECT = [4,6,7,9];    
DATASELECT = [2,4,5,6];

if length(VARSELECT) ~= length(DATASELECT)
    warning('Make sure you have the same number of variables for VARSELECT and DATASELECT')
    return
end

GDPINDEX = 2; % second position in ynames

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

% Doing moment comparisons for priors (CHECKPRIOR == 1)
if CHECKPRIOR
   Theta_s = Theta_pr;
end

NT = size(Theta_s,1);
Thetamat = Theta_s((NT+rem(NT,2))/2 + 1:end,:);
    
% Cartman: "Thin out their numbers!"

Nthin = 200; % Must be < N
N = size(Thetamat,1);

count = N/(Nthin + rem(Nthin,2));

nselect = count*[1:Nthin];

Theta_mat = Thetamat(nselect,:);
NT = size(Theta_mat,1);

h = waitbar(0,'Solving models for thinned parameter MCs...');
for n = 1:NT
    
    % General solution: y(t+1) = AA*y(t) + BB*z(t)
    [ AA,BB,PROBLEM ] = model_solve(Theta_mat(n,:),POLICY);
    
    % Compute covariogram
    Lag0 = 0;
    Lag1 = 1;
    Lag2 = 2;
    [sigyy0(:,:,n),sigxx0(:,:,n),sigxy0(:,:,n),sigyx0(:,:,n)] = mom(0,AA,BB*BB',Lag0);
    [sigyy1(:,:,n),sigxx1(:,:,n),sigxy1(:,:,n),sigyx1(:,:,n)] = mom(0,AA,BB*BB',Lag1);
    [sigyy2(:,:,n),sigxx2(:,:,n),sigxy2(:,:,n),sigyx2(:,:,n)] = mom(0,AA,BB*BB',Lag2);
    
    % Compute similar for Covariances
waitbar(n/NT,h,[sprintf('%0.2f',n*100/NT),' percent done solving... '])
end % for n
close(h)


nv = length(VARSELECT);

for i = 1:nv
    var(:,i) = reshape(sigxx0(VARSELECT(i),VARSELECT(i),:),Nthin,1);
    cov(:,i) = reshape(sigxx0(VARSELECT(i),GDPINDEX,:),Nthin,1);
    autocov_one(:,i) = reshape(sigxx1(VARSELECT(i),VARSELECT(i),:),Nthin,1);
    autocov_two(:,i) = reshape(sigxx2(VARSELECT(i),VARSELECT(i),:),Nthin,1);
end

stdv = sqrt(var);
%corry = cov./(stdv.*(stdv(:,GDPINDEX)*ones(1,size(stdv,2))));
auto_one = autocov_one./var;
auto_two = autocov_two./var;

% BUSINESS CYCLE STATISTICS:

% S.d. of variables
figure
for i = 1:nv
subplot((nv+rem(nv,2))/2,2,i)
[f,xi] = ksdensity(stdv(:,i),'support',[0, 100]);
plot(xi,f);
ylabel(ynames(VARSELECT(i),:))
hold on
        thetaline = [min(f):.1:max(f)];
        sddata = std(newdata(:,DATASELECT(i)))
        plot(sddata*ones(length(thetaline)),thetaline,'r','LineWidth',0.8) 
        text(sddata,thetaline(end-3),num2str(sddata))
hold off
suptitle2('Standard deviation',6)
end

% 1st order autocorrelation of variables
figure
for i = 1:nv
subplot((nv+rem(nv,2))/2,2,i)
[f,xi] = ksdensity(auto_one(:,i),'support',[-1, 1]);
plot(xi,f);
ylabel(ynames(VARSELECT(i),:))
hold on
        thetaline = [min(f):.1:max(f)];
        autod1 = autocorr(newdata(:,DATASELECT(i)),1);
        plot(autod1(end)*ones(length(thetaline)),thetaline,'r','LineWidth',0.8) 
        text(autod1(end),thetaline(end-3),num2str(sddata))
hold off
suptitle2('First order autocorrelation',6)
end

% 2nd order autocorrelation of variables
figure
for i = 1:nv
subplot((nv+rem(nv,2))/2,2,i)
[f,xi] = ksdensity(auto_two(:,i),'support',[-1, 1]);
plot(xi,f);
ylabel(ynames(VARSELECT(i),:))
hold on
        thetaline = [min(f):.1:max(f)];
        autod2 = autocorr(newdata(:,DATASELECT(i)),2);
        plot(xi,autod2(end)*ones(length(thetaline)),thetaline,'r','LineWidth',0.8) 
        text(autod2(end),thetaline(end-3),num2str(sddata))
hold off
suptitle2('Second order autocorrelation',6)
end

% 1st order correlation of variables with GDP
% figure
% for i = 1:nv
% subplot((nv+rem(nv,2))/2,2,i)
% [f,xi] = ksdensity(corry(:,i),'support',[-1, 1]);
% plot(xi,f);
% title(ynames(VARSELECT(i),:))
% end




