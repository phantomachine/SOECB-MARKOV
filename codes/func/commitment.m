function [H1,H2,PROBLEM]=RD_COMM(a0,a1,a2,a3,a4,a5,n,beta,q,w)

%This code was written by Richard Dennis to implement the
%optimal pre-commitment policy rule section in 'Solving for
%Optimal Policy in Rational-Expectations Models: New Solution
%Algorithms', Federal Reserve Bank of San Francisco Working
%Paper #01-09.  The code comes with no performance guarantees.
%However, if you find this code useful then please cite the paper.

%The stochastic constraints take the form

%A0yt = A1yt-1 + A2Etyt+1 + A3xt + A4Etxt+1 + A5vt

%The solution takes the form:
%[lt' yt' xt']' = H1*[lt-1' yt-1' xt-1']' + H2*[vt' 0 0]'
%where lt-1 is a vector of pre-determined Lagrange multipliers.

TOL = 1e-6;
maxiter = 5000;

onn = zeros(n,n);
onp = zeros(n,cols(a3));
opp = zeros(cols(a3),cols(a3));
enn = eye(n);

a = [ onn a0 -a3;
      a0' w  onp;
   (-a3') onp' q];

b = [ onn a1 onp;
   (a2'/beta) onn onp;
   (a4'/beta) onp' opp];
   
c = [onn a2 a4;
   (beta*a1') onn onp;
   onp' onp' opp];

d = [a5 onn onp;
    onn onn onp;
    onp' onp' opp];
   
% solve the rational expectations model using Binder-Pesaran (1995) 
       
h = b;
g = d;
leni = 1;
i = 0;
while leni >= TOL;
  hn = inv(a-c*h)*b;
  leni = max(max(abs(hn-h)));  % d_infinity metric, sup| hn - h |
  h = hn;
  if i == maxiter
      PROBLEM = 1;
      return
  else
      PROBLEM = 0;
  end
  i = i + 1;
end;
g = inv(a-c*h)*d;      


% print results 

H1 = makezero(h);
H2 = makezero(g);
%uncvar=diag(RD_makezero(sigma));

