function  [F1,F2,H1,H2,PROBLEM] = discretion(A0, A1, A2, A3, A4, A5, n, beta, Q, W)
%   disc	Solves model for discretion rule and companion representation.
%	Sourced directly from Richard Dennis' Gauss algorithm


H1 = A1;
H2 = A5;
F1 = zeros(1,n);
F2 = zeros(1,n);
M = zeros(n,n);
PHIY = eye(n);

% solving for feedback matrices H1, H2, F1, F2 
maxiter = 1000;

lenj = 1;
j = 0;
while lenj >= 1e-6 
  leni = 1;
  i = 0;
while leni >= 1e-6  
    MN = W+beta*F1'*Q*F1+beta*H1'*M*H1;
    leni = max(abs(MN-M));    % d_max metric, max| mn - m |
    M = MN;
    i = i + 1;
    if i == maxiter,
        PROBLEM=1;
        return;
    else
        prob=0;
    end
end
  D = A0-A2*H1-A4*F1;
  Dinv = inv(D);
  CC = -(Q+A3'*Dinv'*M*Dinv*A3)\A3'*Dinv'*M*Dinv;
  F1n = CC*A1;
  F2  = CC*A5;
  H1 = Dinv*(A1+A3*F1n);
  H2 = Dinv*(A5+A3*F2);
  lenj = max(abs(F1n'-F1'));
  F1 = F1n;
    
  if j == maxiter
      PROBLEM = 1;
      return
  else
      PROBLEM = 0;
  end
  j = j + 1;
end;
