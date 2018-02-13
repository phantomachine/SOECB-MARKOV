function [Sigma, PROBLEM] = kalman_initial(A, G, V)

% This function initializes Sigma for the Kalman algortim using the
% doubling algorithm. Please refer to kfilter for more general setup for
% Sigma0 and K0.
% Original author: Philip Liu 25 October 2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		x[t+1] = Ax[t] + Cw[t+1]
%       y[t] = Gx[t] + Dv[t] 
%       where:  E C[w(t+1)] [w(t+1)']C' = V;
%       and     E [v(t+1)] [v(t+1)'] = 0;
%               E [w(t+1)] [v(t+1)] = 0;
% Riccati equation:
% S = A*S*A' + V - A*S*G'*inv(G*S*G')*G*S*A';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PROBLEM=0;
m=max(size(A));
s0=.01*eye(m);
dd=1;
it=1;
maxit=2500;
  while (dd>1e-4 && it<=maxit);
    s1 = A*s0*A' + V - A*s0*G'*inv(G*s0*G')*G*s0*A';
    dd=max(max(abs(s1-s0)));
    it=it+1;
    s0=s1;
  end;
Sigma=s0;

if it>=maxit; 
    disp(['WARNING: Iteration limit of ', num2str(maxit),' reached in kalman_initial.m']);
    PROBLEM=1;
end;