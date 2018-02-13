function prior = priorsim( pshape, p1, p2, ndraws)

% priorsim.m
% This function simulates ndraws from the prior distribution used for
% generating the prior density plot
% The code comes with no performance guarantees.
% However, if you find this code useful then please cite the paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% July 2005: orginal author Philip Liu
%
%==========================================================================
% Input:
%       pshape = parameter distributions of size n
%       p1 = mean
%       p2 = std deviations
%       ndraws = number of draws generated
%
% Return: 
%       prior = n x ndraws matrix containing prior draws
%       
%==========================================================================
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

npara = size(pshape,1);
prior=[];

ii=1;
while ii <= npara
    a=0;
    b=0;
        
    if pshape(ii) == 1;     %  BETA Prior 
     a = (1-p1(ii))*p1(ii)^2/p2(ii)^2 - p1(ii);
     b = a*(1/p1(ii) - 1);
     draws=betarnd(a,b,ndraws,1);
     prior=[prior draws];
     
   elseif pshape(ii) == 2; % GAMMA PRIOR 
     b = p2(ii)^2/p1(ii);
     a = p1(ii)/b;
     draws=gamrnd(a,b,ndraws,1);
     prior=[prior draws];
     
   elseif pshape(ii) == 3; % GAUSSIAN PRIOR 
     draws=normrnd(p1(ii),p2(ii),ndraws,1);
     prior=[prior draws];
     
   elseif pshape(ii) == 4; % INVGAMMA1 PRIOR 
     b = p2(ii)^2/p1(ii);
     a = p1(ii)/b;
     draws=1./gamrnd(a,b,ndraws,1);
     prior=[prior draws];
     
   elseif pshape(ii) == 5; % UNIFORM PRIOR 
     draws= a + (b-a) * rand(ndraws);
     prior=[prior draws];
   end;
   ii = ii+1;
end;
