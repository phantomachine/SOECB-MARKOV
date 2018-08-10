%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aus_estimate.m
%
% This program estimates the model in Kam, Lees and Liu for Australian
% data. Data source: IFS,RBA Bulletin,FRED
%
% Estimation using the Metropolis-Hasting Markov Chain Monte Carlo method
% to simulate posterior densities of parameter vector.
%
% Original author: Philip Liu  (gm_mcmc.m code)
% Last tweak: T. Kam, Feb 7, 2006.
% Final tweak: PL and TK April 4, 2006.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script calls:
%     
%     
% Function calls:
%     aus_datadoc.m to manipulate and plot data series
%     gm_likelihood.m for evaluating the likelihood
%     priorln.m for evaluating the prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
tic
warning off
delete output.txt
diary('output.txt')

disp('Estimating for Australia: NK-foreign VAR model & muq not zero')

LOAD_OLD = 0;

%% Policy Choice
POLICY = 0;         % 0 = discretion; 1 = commitment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LC path options - /short
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath(genpath('/home/401/pxl401/Matlab func/'));
% if POLICY==0
%    load_path='/short/x93/dsge/au_dsgemcmc_com/mh_dis';
% else
%     load_path='/short/x93/dsge/au_dsgemcmc_com/mh_com';
% end
% wip_path='/short/x93/dsge/au_dsgemcmc_com/mh_wip';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local path options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('func/'))
addpath(genpath('stats/'))

if LOAD_OLD == 1
    if POLICY==0
       load_path='chain/mh_dis';
    else
       load_path='chain/mh_com';
    end
end
wip_path='chain/mh_wip';

%% Settings options
COUNTRY = 'Australia - ';
PLOT_DATA = 1;      % plot raw and filtered/dummied data
NLAGPOLICY = 1;     % = nlp if nlp #lag of r(t) appears in the state vector. Input: 0, 1, 2 .. etc
mh_scale = .25;      % scale variance to obtain optimal acceptance rate (20%)
N = 500000;          % number MCMC draws
MH_TRIM = 0;
Nburn = MH_TRIM*N;  % number of burn in periods
loadmh = 1;         % 1 load previously saved chain, 0 start new one
ny = 21;            % number of state variables, rows(ynames)
nx = 1;             % number of policy variables
PER_SAVE=1;         % periodic save option
Xsave=5000;         % saving after Xsave draws
% rand('state',34)
disp(['Number of simulations: ', num2str(N), ' draws']);
disp(['mh scale: ', num2str(mh_scale)]);
%% Import COUNTRY data and setting data frequency and estimation startdate.
[data,datadate,raw] = xlsread('data/ausdata.xls');
% save raw_data data datadate raw
%load raw_data
%data = raw;
data(:,3) = -data(:,3); % innvert RER series to log(eP*/P)

interval = [0.1, 0.2, 0.3, 0.4]'; % set to == 0 if annual, else, %[0.1, 0.2, 0.3, 0.4]'
sampldate = 1980;
startdate = 1990;

% Call function to manipulate and plot data
[newdata,dtitle,date_vec] = aus_datadoc(COUNTRY,data,datadate,interval,sampldate,startdate,PLOT_DATA);
newdata(:,1)=[];            % get rid of consumption series

percent_index = [2,3,4,8];
for i = 1:length(percent_index)
    newdata(:,percent_index(i)) = newdata(:,percent_index(i))*100;
end

meandata = ones(size(newdata,1),1)*mean(newdata);
newdata = newdata-meandata;

para_latex;                 % Create matrix of paramater names - text/LateX strings
npara = rows(paraname);

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


pri=[   0       0.99     0.001;    % beta ~ BETA
        0       0.45     0.1;      % alpha ~ BETA
        1       0.6      0.2;      % h ~ BETA      % mean is .86
        2       1     0.5;           % sigma ~ GAMMA
        3       1.5     0.25;   % phi ~ NORMAL
        2       1       0.5;      % eta ~ GAMMA
        1       0.7     0.2        % deltah ~ BETA
        1       0.7     0.2        % deltaf ~ BETA
        1       .5     .2;      % thetaH ~ BETA
        1       .5     .2;      % thetaF ~ BETA
        2       .5     .2;    % a1 ~ BETA
        0       0     0.1;    % a2 ~ NORMAL
        0       0    0.1;    % a3 ~ NORMAL
        0       0     0.1;    % b1 ~ NORMAL
        2       .5     .2;    % b2 ~ BETA
        0       0   0.1;    % b3 ~ NORMAL
        0       0    0.1;    % c1 ~ NORMAL
        0       0    0.1;    % c2 ~ NORMAL
        2       .5     .2;    % c3 ~ BETA
        0       0     0.2;    % rhoh ~ BETA
        0       0     0.2;    % rhof ~ BETA
        1       0.5     0.2;    % rhoa ~ BETA
        1       0.9     0.2;    % rhoq ~ BETA
        1       0.25     0.2;    % rhos ~ BETA
        0       0       0.2;    % rhor ~ BETA
        2       0.5       0.3;    % mu_q ~ GAMMA   %check genTheta.m, set ne*wTheta(25)=0
        2       0.5       0.3;    % mu_y ~ GAMMA
        2       0.5       0.3;    % mu_r ~ GAMMA
        4       .5        .25;    % sigmah ~ INVGAMMA
        4       .5        .25;    % sigmaf ~ INVGAMMA
        4       1       .4;    % sigmaz ~ INVGAMMA
        4       2       .5;    % sigmaq ~ INVGAMMA
        4       1       .4;    % sigmas ~ INVGAMMA
        4       1       .4;    % sigmapi* ~ INVGAMMA
        4       1       .4;    % sigmay* ~ INVGAMMA
        4       1       .4;    % sigmar* ~ INVGAMMA
        4       1       .4];    % sigmar ~ INVGAMMA
%         3       2.5      1;    % c_pif
%         3       2.5     .5;    % c_pi
%         2       6       1;     % c_r
%         3       2.5     .5;    % c_pist
%         2       5       1];    % c_rst

TPOL = [15:17];    % pick index of policy parameters
TPRI = [1:8];      % pick index of private deep parameters
 
                          
% Truncations and bounds
point_b = [-100 100];  % bounds for point mass
bet_b = [0 1];      % bounds for BETA density
gam_b = [0 100];    % bounds for GAMMA density
nor_b = [-10 10];    % bounds for NORMAL density
invgam_b = [0 25];

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
        bounds(i,:)=invgam_b;
    end
end
pshape = pri(:,1);
p1 = pri(:,2);
p2 = pri(:,3);
lb = bounds(:,1);
ub = bounds(:,2);

%% Jumping distribution
%% Load starting values else use prior mean
load mh_init0
if loadmh==0,
    mh_old=0;
    Theta = mh_init ;
    Theta_s =[];
    loglike_s = [];
    logpri_s = [];
else
    load(load_path);
    nn = min(find(loglike_s==0));
    if isempty(nn)
        nn = length(loglike_s)+1;
    end
    Theta_s = Theta_s(1:nn-1,:);
    loglike_s = loglike_s(1:nn-1,:);
    logpri_s = logpri_s(1:nn-1,:);
    mh_old=length(Theta_s);
    Theta=Theta_s(mh_old,:);
    disp(['Number of previous draws loaded: ', num2str(nn-1), ' draws']);
end

% increment random vector, z ~ N(0,variance_normal) for RW-MH algorithm MCMC: 
% Theta(s) = Theta(s-1) + z
if POLICY == 0,
    variance_normal(1)=0.01;
    variance_normal(2)=0.01;
    variance_normal(3)=0.01;
    variance_normal(4)=0.01;
    variance_normal(5)=0.1;
    variance_normal(6)=0.1;
    variance_normal(7)=0.01;
    variance_normal(8)=0.01;
    variance_normal(9)=0.01;
    variance_normal(10)=0.01;
    variance_normal(11)=0.05;
    variance_normal(12)=0.01;
    variance_normal(13)=0.05;
    variance_normal(14)=0.05;
    variance_normal(15)=0.05;
    variance_normal(16)=0.05;
    variance_normal(17)=0.01;
    variance_normal(18)=0.01;
    variance_normal(19)=0.05;
    variance_normal(20)=0.01;
    variance_normal(21)=0.01;
    variance_normal(22)=0.01;
    variance_normal(23)=0.01;
    variance_normal(24)=0.01;
    variance_normal(25)=0.01;
    variance_normal(26)=0.01;
    variance_normal(27)=0.01;
    variance_normal(28)=0.01;
    variance_normal(29)=0.05;
    variance_normal(30)=0.05;
    variance_normal(31)=0.05;
    variance_normal(32)=0.05;
    variance_normal(33)=0.05;
    variance_normal(34)=0.05;
    variance_normal(35)=0.05;
    variance_normal(36)=0.05;
    variance_normal(37)=0.05;  
%     variance_normal(38)=0.05;
%     variance_normal(39)=0.05;
%     variance_normal(40)=0.05;
%     variance_normal(41)=0.05;
%     variance_normal(42)=0.05;
else 
    variance_normal(1)=0.01;
    variance_normal(2)=0.01;
    variance_normal(3)=0.01;
    variance_normal(4)=0.01;
    variance_normal(5)=0.1;
    variance_normal(6)=0.1;
    variance_normal(7)=0.01;
    variance_normal(8)=0.01;
    variance_normal(9)=0.01;
    variance_normal(10)=0.01;
    variance_normal(11)=0.05;
    variance_normal(12)=0.01;
    variance_normal(13)=0.05;
    variance_normal(14)=0.05;
    variance_normal(15)=0.05;
    variance_normal(16)=0.05;
    variance_normal(17)=0.01;
    variance_normal(18)=0.01;
    variance_normal(19)=0.05;
    variance_normal(20)=0.01;
    variance_normal(21)=0.01;
    variance_normal(22)=0.01;
    variance_normal(23)=0.01;
    variance_normal(24)=0.01;
    variance_normal(25)=0.01;
    variance_normal(26)=0.01;
    variance_normal(27)=0.01;
    variance_normal(28)=0.01;
    variance_normal(29)=0.05;
    variance_normal(30)=0.05;
    variance_normal(31)=0.05;
    variance_normal(32)=0.05;
    variance_normal(33)=0.05;
    variance_normal(34)=0.05;
    variance_normal(35)=0.05;
    variance_normal(36)=0.05;
    variance_normal(37)=0.05;
%     variance_normal(38)=0.05;
%     variance_normal(39)=0.05;
%     variance_normal(40)=0.05;
%     variance_normal(41)=0.05;
%     variance_normal(42)=0.05;
end
variance_normal=variance_normal*mh_scale; 

%% Metropolis-Hastings MCMC algorithm
[loglike,PROBLEM1,PROBLEM2,PROBLEM3] = ...
        model_likelihood(Theta,newdata,ny,nx,POLICY,NLAGPOLICY);

disp('Initial model Log likelihood')
disp(loglike)

logpri = priorln(Theta, pshape, p1, p2);

loglikelogpri=loglike+logpri;

save initiallikepriNK loglikelogpri

pau=0;      % initial value for the acceptance rate, pau is acceptance rate 
test1=0;
probs1=0;
probs2=0;
hh = waitbar(0,'Starting RW Metropolis-Hasting');

set(hh,'Name','RWMH-MCMC: Please wait.')
Theta_s = [Theta_s; zeros(N,npara)];
loglike_s = [loglike_s; zeros(N,1)];
logpri_s = [logpri_s; zeros(N,1)];
for j=1:N,
        
        [newTheta] = genTheta(Theta, 0, variance_normal, npara, bounds);      %generate new Theta
	     
        newlogpri = priorln(newTheta, pshape, p1, p2);
        
        if newlogpri>-Inf 
                [newloglike,PROBLEM1,PROBLEM2,PROBLEM3] = ...
                        model_likelihood(newTheta,newdata,ny,nx,POLICY,NLAGPOLICY);
                
                probs1=probs1+PROBLEM1;             % indeterminacy
                probs2=probs2+PROBLEM2+PROBLEM3;    % problem in evaluating the likelihood
                
                if newloglike==-Inf;
                        ratio=0;
                else
                        ratio=exp(newloglike+newlogpri-(loglike+logpri)); %alpha in notes
                        
                        if rand<=ratio              % we accept with prob ratio
                                loglike=newloglike;
                                logpri=newlogpri;
                                Theta=newTheta;
                                pau=pau+1;
                        end
                end
        end
    
	Theta_s(mh_old+j,:)=Theta;
	loglike_s(mh_old+j,1)=loglike;
	logpri_s(mh_old+j,1)=logpri;
    
     %saving per X draws
%    if rem(length(loglike_s),Xsave)==0 & PER_SAVE==1;
%        save(load_path, 'loglike_s', 'logpri_s', 'Theta_s')
%    else
%    end

    prtfrc = j/N;
    if rem(100*j,N)==0;
      if PER_SAVE==1;
        save(load_path, 'loglike_s', 'logpri_s', 'Theta_s')
        disp([num2str(mh_old+j), ' draws saved']);
      else
      end
    toc;
    disp(['percent done: ', num2str(prtfrc*100)]);
    disp(['accept: ', num2str(pau/j)]);
    disp(['percent indeterminancies: ' num2str((probs1+probs2)/j*100)]);
    time=fix(clock)  
    end     
    
    waitbar(prtfrc,hh,sprintf('%f percent done NK-AR-q, accept rate %f',prtfrc,pau/j)); 
end  
toc;  
close(hh)

pau=pau/N*100;
probs1=probs1/N*100;
probs2=probs2/N*100;
ThetaMn=mean(Theta_s);

fid4=fopen('NKsummary.txt','w+');
fprintf(fid4,'%s ','# of draws:');
fprintf(fid4,'%d\n',N);
fprintf(fid4,'%s ','rate of acceptance (in %):');
fprintf(fid4,'%f\n',pau);
fprintf(fid4,'%s ','% of indeterminancies:');
fprintf(fid4,'%f\n',probs1);
fprintf(fid4,'%s ','% of invalid likelihood:');
fprintf(fid4,'%f\n',probs2); 
fclose(fid4);

%% Marginal density
nn = min(find(loglike_s==0));
if isempty(nn)
   nn = length(loglike_s)+1;
end
Theta_s = Theta_s(1:nn-1,:);
loglike_s = loglike_s(1:nn-1,:);
logpri_s = logpri_s(1:nn-1,:);


% getting rid of fix parameters
Theta_marginal = [];
Thetapost = [];
paran = [];

for i=1:npara
     if pshape(i)~=0
        Theta_marginal = [Theta_marginal, Theta_s(:,i)];
        Thetapost  = [Thetapost, Theta_s(Nburn+1:end,i)];
        paran = [paran; paraname(i,:)];
     else
     end
end
marginal = marginal_density(Theta_marginal,logpri_s,loglike_s,MH_TRIM)


save(load_path, 'loglike_s', 'logpri_s', 'Theta_s', 'marginal')
save(wip_path);


%% PLOT PRIOR and POSTERIOR DENSITIES OF PARAMETERS

% You need to modify paran, prior, if you are calibrating some of the
% parameters. Delete the ones from the list of paraname and pri which are
% calibrated.



priordens = priorsim(pshape,p1,p2,N-Nburn);


% % Construct table of posterior estimator statistics
statmat = []; % basic stats
altstatmat = []; % other stats
i = 1;
while i < cols(Thetapost)
    
    prior_up = sort(priordens(:,i));
        thetapri_l = prior_up(round(0.025*(N-Nburn)));
        thetapri_u = prior_up(round(0.975*(N-Nburn)));
        thetapri_mean = mean(prior_up);
        thetapri_std = std(prior_up);
        thetapri_med = median(prior_up);
    
    post_up = sort(Thetapost(:,i));
        thetapost_l = post_up(round(0.025*(N-Nburn)));
        thetapost_u = post_up(round(0.975*(N-Nburn)));
        thetapost_mean = mean(post_up);
        thetapost_std = std(post_up);
        thetapost_med = median(post_up);

statmat(i,:) = [ thetapri_mean, thetapri_l, thetapri_u, ...
                                thetapost_mean, thetapost_l, thetapost_u]; 

altstatmat(i,:) = [ thetapri_mean, thetapri_med, thetapri_std, ...
                                thetapost_mean, thetapost_med, thetapost_std];

i = i+1;   
end

% Now write statmat and altstatmat into a TeX Table
collabel = {'Prior mean','5\%','95\%','Posterior mean','5\%','95\%'};        
for i=1:rows(paran) 
    rowlabel{:,i} = paran(i,:); 
end
matrix2latex2(statmat, 'AUSVAR1stat1.tex', 'rowLabels', rowlabel, ...
            'columnLabels', collabel, 'alignment', 'c', 'format', '%-6.2f','size', 'tiny');        


collabel = {'Prior mean','median','Std','Posterior mean','median','Std'};        
matrix2latex2(altstatmat, 'AUSVAR1stat2.tex', 'rowLabels', rowlabel, ...
            'columnLabels', collabel, 'alignment', 'c', 'format', '%-6.2f','size', 'tiny');        

        
% PLOT PRIOR AND POSTERIOR DENSITIES
optim = 0;  % kernel estimation bandwidth parameter

plot_density(Thetapost,priordens,optim,paran,TPOL,TPRI);

% PLOT POSTERIOR CHAINS INDIVIDUALLY
figure
TPLOT = length(Thetapost);
for i = 1:cols(Thetapost)
    subplot(5,6,i)
    plot(Thetapost(end-TPLOT+1:end,i))
    title(paran(i,:),'FontSize',18)
end
print -depsc MCMCq
hgsave('MCMCq')
toc;
diary off;

disp('DONE estimating for Australia: NK-foreign VAR model & muq not zero')
