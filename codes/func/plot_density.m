function density_plot = plot_density( Theta, Prior,optimal,paraname,TPOL,TPRI)

% Original author: Philip Liu
% This function produces plot of the posterior and prior density
% Needs: mh_optimal_bandwidth.m, kernel_density_estimate.m

nn=length(Theta);
mm=length(Prior);
kk='gaussian';
[draws npara] = size(Theta);

%%  Calculates the optimal bandwidth parameter of a kernel estimator 
%%  used to estimate a posterior univariate density from realisations of a 
%%  Metropolis-Hastings algorithm.  
%% 
%%  * M. Skold and G.O. Roberts [2003], "Density estimation for the Metropolis-Hastings algorithm". 
%%  * Silverman [1986], "Density estimation for statistics and data analysis". 
%%
%%  data            :: a vector with n elements.
%%  bandwidth       :: a scalar equal to 0,-1 or -2. For a value different from 0,-1 or -2 the
%%                     function will return optimal_bandwidth = bandwidth.
%%  kernel_function :: 'gaussian','uniform','triangle','epanechnikov',
%%                     'quartic','triweight','cosinus'.
%%  hh              :: returns optimal bandwidth            
%--------------------------------------------------------------------------
if optimal==0;
    bb=0;
    display('Rule of thumb bandwidth parameter');
                    %  Rule of thumb bandwidth parameter (Silverman [1986] corrected by 
                    %  Skold and Roberts [2003] for Metropolis-Hastings). 
elseif optimal==-1;
    bb=-1;
    display('Plug-in estimation of the optimal bandwidth');
                    % Adaptation of the Sheather and Jones [1991] plug-in estimation of the optimal bandwidth 
                    % parameter for metropolis hastings algorithm.
elseif optimal==-2;
    bb=-2;
    display('Bump killing to smooth long rejecting periods');
                    % Bump killing... We construct local bandwith parameters in order to remove 
                    % spurious bumps introduced by long rejecting periods.   
elseif  optimal>0;
    bb=optimal;   
    display('User specified');
                    % User specified.
end

%Posterior bandwidth:
obandp=[];
i=1;
while i<=npara;
     hh=mh_optimal_bandwidth(Theta(:,i),nn,bb,kk);
     obandp=[obandp; hh];
     i=i+1;
end

%Prior bandwidth:
obandpr=[];
i=1;
while i<=npara;
     hh=mh_optimal_bandwidth(Prior(:,i),mm,bb,kk);
     obandpr=[obandpr; hh];
     i=i+1;
end

%%  Estimating a continuous density. A kernel density 
%%  estimator is used (see Silverman [1986]). 
%% 
%%  * Silverman [1986], "Density estimation for statistics and data analysis". 
%%
%%  The code is adapted from DYNARE TOOLBOX. 

grid=2^9;
% Posterior Density:
%==========================================================================
pdens=[];
ffp=[];
    i=1;
    while i<=npara;
        [pden f]=kernel_density_estimate(Theta(:,i),grid,obandp(i),kk);
        pdens=[pdens; pden'];
        ffp=[ffp; f'];
        i=i+1;
    end
% Prior Density:
%==========================================================================
prdens=[];
ffpr=[];
    i=1;
    while i<=npara;
        [prden f]=kernel_density_estimate(Prior(:,i),grid,obandpr(i),kk);
        prdens=[prdens; prden'];
        ffpr=[ffpr; f'];
        i=i+1;
    end
    

% Plots:
%--------------------------------------------------------------------------

%s1 = strvcat('h', '\sigma', '\eta', '\phi', '\theta_H', '\theta_F',  '\phi_1','\phi_2','\delta', '\rho_{r^*}', '\rho_a', '\lambda_1'); 

% if rem(npara,2) ~= 0
%     prows = 6;
%     pcols = (npara+1)/prows;
% else
%     prows = 6;
%     pcols = npara/prows;
% end
    
TDEEP = [TPOL,TPRI];
TEXO = 1:npara; 
TEXO(TDEEP) = [];  % pick index of exogenous processes parameters

figure 
    i=1;
    while i <= length(TDEEP);
        subplot(4,3,i,'align');
        ymin = min(min(ffp(TDEEP(i),:)),min(ffpr(TDEEP(i),:)));
        ymax = max(max(ffp(TDEEP(i),:)),max(ffpr(TDEEP(i),:)));
        xmin = min(min(pdens(TDEEP(i),:)),min(prdens(TDEEP(i),:)));
        xmax = max(max(pdens(TDEEP(i),:)),max(prdens(TDEEP(i),:)));
        plot(pdens(TDEEP(i),:), ffp(TDEEP(i),:) , 'k', prdens(TDEEP(i),:), ffpr(TDEEP(i),:), '--', 'LineWidth', 0.8);
        hold on
        thetaline = ymin:0.1:ymax;
        meanTh = mean(Theta(:,TDEEP(i)));
        plot(meanTh*ones(length(thetaline)),thetaline,'r','LineWidth',0.8)    
        text(meanTh,thetaline(end-3),num2str(meanTh))
        hold off
        xlim([xmin,xmax])
        ylim([ymin,ymax])
        title(paraname(TDEEP(i),:), 'fontsize', 18);
        i=i+1;
    end
    print -depsc 01PriorPostDeep
    hgsave('01PriorPostDeep')
    
figure 
    i=1;
    while i <= length(TEXO);
        subplot(5,3,i,'align');
        plot(pdens(TEXO(i),:), ffp(TEXO(i),:) , 'k', prdens(TEXO(i),:), ffpr(TEXO(i),:), '--', 'LineWidth', 0.8);
        ymin = min(min(ffp(TEXO(i),:)),min(ffpr(TEXO(i),:)));
        ymax = max(max(ffp(TEXO(i),:)),max(ffpr(TEXO(i),:)));
        xmin = min(min(pdens(TEXO(i),:)),min(prdens(TEXO(i),:)));
        xmax = max(max(pdens(TEXO(i),:)),max(prdens(TEXO(i),:)));
        hold on
        thetaline = ymin:0.1:ymax;
        meanTh = mean(Theta(:,TEXO(i)));
        plot(meanTh*ones(length(thetaline)),thetaline,'r','LineWidth',0.8) 
        text(meanTh,thetaline(end-3),num2str(meanTh))
        hold off
        xlim([xmin,xmax])
        ylim([ymin,ymax])
        title(paraname(TEXO(i),:), 'fontsize', 18);
        i=i+1;
    end
    print -depsc 02PriorPostExog
    hgsave('02PriorPostExog')

