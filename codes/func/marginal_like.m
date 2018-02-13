function marginal = marginal_like(Theta_s, loglike_s,logpri_s,MH_TRIM)
% marginal_like.m
% 
% Marginal Likelihood estimator
% This function computes the marginal likelihood of a model given some data
% following the Harmonic Mean procedure proposed by Gelfand and
% Dey (1994). The weighting function follows Geweke (1998). 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orginal author: Juan Rubio-Ramirez [New York and Minneapolis, 6-21-2001]
%                 and DYNARE
% Modified by: Philip Liu, 20 Feb 2006
%
% References:
%    Gelfand, A.E. and D.K. Dey (1994), "Bayesian Model Choice:
%    Asymptotics and Exact Calculations"
%    Journal of the Royal Statistical Society B, 56, pp. 501-514.
%    Geweke, J. (1998), "Using Simulation Methods for Bayesian 
%    Econometric Models: Inference, Development, and Communication" 
%    Staff Report 249, Federal Reserve Bank of Minneapolis.
%=========================================================================
% input:    Theta_s = contains the posterior theta (after eliminating fixed
%                     parameters)
%           loglike_s = log of the likelihood
%           logpri_s = log of prior
%           MH_TRIM = % of initial draws to be dropped
%
% output:   marginal = [values of p,
%                       log of the marginal likelihood p(Data|model)]
%                       
% Function calls:     
%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nruns = size(Theta_s,1);
npara = size(Theta_s,2);
logpost=logpri_s+loglike_s;
T = nruns-round(MH_TRIM*nruns)+1;
N = T;
simulations = zeros(N,npara);
lposterior   = zeros(N,1);
% Trimming the chain
simulations = Theta_s(round(MH_TRIM*nruns):nruns,:);
lposterior = logpost(round(MH_TRIM*nruns):nruns);

lpost_mode = max(lposterior);

MU = mean(simulations)';
SIGMA = zeros(npara,npara);
for i=1:N;
    SIGMA = SIGMA + (simulations(i,:)'-MU)*(simulations(i,:)'-MU)';
end;
SIGMA = SIGMA/N;

DetSIGMA = det(SIGMA);
InvSIGMA = inv(SIGMA);
marginal = [];

% calculate the marginal density for serveral values of p; smaller value of
% p will result in better behavior of P(Y|A) over the domain of Theta but
% great simulation error due to smaller of theta in Theta
for p = 0.1:0.1:0.9;
    critval = chi2inv(p,npara);
    tmp = 0;
    % Eliminates thetas according to geweke (1998) eq. 4.32
    j=0;
    for i = 1:N;
        deviation  = (simulations(i,:)-MU')*InvSIGMA*((simulations(i,:)-MU'))';
        if deviation <= critval;
           lftheta = -log(p)-(npara*log(2*pi)+log(DetSIGMA)+deviation)/2;
           % Calculate f(theta)/posterior 
           tmp = tmp + exp(lftheta - lposterior(i));
           j=j+1;
        end;    
    end;
    marginal = cat(1,marginal,[p,real(-log(tmp/N))]); 
end;    
    



