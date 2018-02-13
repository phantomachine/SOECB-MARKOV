% redpolicy.m
%
%-------------------------------------------------------------------------%
% Calculate reduced form optimal linear policy rule from posterior estimates
% Author: T.Kam, (c) 2008
% last update: 15/05/08
%-------------------------------------------------------------------------%

clear all
clc
%warning off

%-------------------------------------------------------------------------%
% STEP 0: Settings
%-------------------------------------------------------------------------%

loadOLD = 1;  % 0 = no; 1 = yes
  
muq = 'muq';  % 'muq0'

POLICY = 0; %0 = discretion; 1 = commitment
COUNTRY = 1; % 1 = Australia; 2 = Canada; 3 = New Zealand



    if COUNTRY == 1
        country = 'Australia';
    elseif COUNTRY == 2
        country = 'Canada';
    elseif COUNTRY == 3
        country = 'NZ';
    end
    
fraction = 0.01;

variables;

% Counting number of endogenous, control and exogenous variables
ny = size(ynames,1);
nx = size(xnames,1);
nz = size(znames,1);
        
% List of variables names as arguments in loss function: (q,ygap,rt-rt-1,pibar)
%targets_y = [ynames(4,:); ynames(6,:); ynames(8,:); ynames(21,:)];

% Do not touch below:
addpath(genpath('../func/'))
addpath(genpath('../chain/'))

% Thetapost = xlsread('postmeantheta.xls', 1, 'B2:D38');
% Theta = Thetapost(:,COUNTRY);

% Load save MCMC chain of parameters

if loadOLD == 0
    
    load mh_dis;
    NOBS = size(Theta_s,1);

    Theta_s = Theta_s(ceil(NOBS*(1-fraction)):end,:);

    [NEXP,NVAR] = size(Theta_s);



    P1 = zeros(NEXP,ny);
    P2 = zeros(NEXP,nz);

    h = waitbar(0,'Starting ...');

    for n = 1:NEXP
    %---------------------------------------------------------------------%
    % STEP 1: Find Markov equilibrium solution: 
    % y(t) endog., x(t) policy, z(t) exog.
    % 
    % y(t) = H1*y(t-1) + H2*z(t)
    % x(t) = F1*y(t-1) + F2*z(t)
    %---------------------------------------------------------------------%
    
        [F1,F2,H1,H2,AA,BB,xnames,ynames,znames,PROBLEM] = ...
                                        model_solve2(Theta_s(n,:),POLICY);

        if PROBLEM == 1;
        disp(['Maximum number of iterations reached! Problem in model_solve.'; 
          'You may not get proper convergence to fixed point solutions. ' ])
        end

    % STEP 2: Extract policy rule
    %---------------------------------------------------------------------%
    % x(t) = (F1*inv(H1))*y(t) + (F2 - inv(H1)*H2)*z(t)
    %---------------------------------------------------------------------%
    
        iH1 = pinv(H1);
        P1(n,:) = F1*iH1;
        P2(n,:) = F2(:,1:nz) - F1*iH1*H2(:,1:nz);

        waitbar(n/NEXP,h,[sprintf('%0.1f',n*100/NEXP),' percent done.'])
    end % ENDFOR
    close(h)

    if COUNTRY == 1 && POLICY == 0
        save Oz-muq-disc P1 P2
    elseif COUNTRY == 1 && POLICY == 1
        save Oz-muq-comm P1 P2
    elseif COUNTRY == 2 && POLICY == 0
        save Can-muq-disc P1 P2
    elseif COUNTRY == 2 && POLICY == 1
        save Can-muq-comm P1 P2
    elseif COUNTRY == 3 && POLICY == 0
        save NZ-muq-disc P1 P2
    elseif COUNTRY == 3 && POLICY == 1
        save NZ-muq-comm P1 P2
    end

else
    if COUNTRY == 1 && POLICY == 0
        load Oz-muq-disc
    elseif COUNTRY == 1 && POLICY == 1
        load Oz-muq-comm
    elseif COUNTRY == 2 && POLICY == 0
        load Can-muq-disc
    elseif COUNTRY == 2 && POLICY == 1
        load Can-muq-comm
    elseif COUNTRY == 3 && POLICY == 0
        load NZ-muq-disc
    elseif COUNTRY == 3 && POLICY == 1
        load NZ-muq-comm
    end
    
    NEXP = size(P1,1);
end

% Calculate statistics

P1 = sort(P1);
P2 = sort(P2);

P1mean = mean(P1,1);
P2mean = mean(P2,1);
Pmean = [P1mean, P2mean];

P1std = std(P1,1);
P2std = std(P2,1);
Pstd = [P1std, P2std];

P1med = median(P1,1);
P2med = median(P2,1);
Pmed = [P1med, P2med];

P1_l = P1(ceil(0.025*NEXP),:);
P1_u = P1(ceil(0.975*NEXP),:);

P2_l = P2(ceil(0.025*NEXP),:);
P2_u = P2(ceil(0.975*NEXP),:);

Plow = [P1_l, P2_l];
Phi = [P1_u, P2_u];


% % STEP 3: Report policy rule in Table form

yznames = [ynames; znames];
for i = 1:size(yznames,1)
    rowlabel{:,i} = yznames(i,:);
end

collabel = {'Mean','Median','2.5\%','97.5\%','Std'};

texfile = sprintf('%s',[country,num2str(POLICY),muq,'-redpolicy.tex']);
matrix2latex2([Pmean', Pmed', Plow', Phi', Pstd'], texfile, 'rowLabels', rowlabel, ...
            'columnLabels', collabel, 'alignment', 'c', 'format', '%-6.3f','size', 'small');





      