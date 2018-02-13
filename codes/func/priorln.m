function lnprior = priorln(para, pshape, p1, p2)
% Orginal author: Frank Schorfheide 
% Modified by: Philip Liu, july 2005
% Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

% This procedure computes a prior density for
% the structural parameters of the DSGE models
% pshape: 0 is point mass, both para and p2 are ignored
%         1 is BETA(mean,stdd)
%         2 is GAMMA(mean,stdd)
%         3 is NORMAL(mean,stdd)
%         4 is INVGAMMA(s^2,nu)
%         5 is UNIFORM [p3,p4]
% p1=mean   p2=std

lnprior = 0;
nprio = length(pshape);

i = 1;
while i <=  nprio;
a = 0;
b = 0;

   if pshape(i) == 1;     %  BETA Prior 
     a = (1-p1(i))*p1(i)^2/p2(i)^2 - p1(i);
     b = a*(1/p1(i) - 1);
     lnprior = lnprior + lpdfbeta(para(i),a,b);   
     
   elseif pshape(i) == 2; % GAMMA PRIOR 
     b = p2(i)^2/p1(i);
     a = p1(i)/b;
     lnprior = lnprior + lpdfgam(para(i),a,b);
     
   elseif pshape(i) == 3; % GAUSSIAN PRIOR 
     lnprior = lnprior + lpdfnorm(para(i),p1(i),p2(i));
     
   elseif pshape(i) == 4; % INVGAMMA1 PRIOR 
     lnprior = lnprior + lpdfig1(para(i),p1(i),p2(i));
     
   elseif pshape(i) == 5; % UNIFORM PRIOR 
     lnprior = lnprior + log(1/(p2(i)-p1(i)));
     
   end;

  i = i+1;
end;