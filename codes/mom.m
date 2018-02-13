function [sigyyJ,sigxxJ,sigxyJ,sigyxJ]=mom(gx,hx,varshock,J)
%[sigyyJ,sigxxJ,sigxyJ,sigyxJ]=mom(gx,hx,varshock,J)
%Let w(t)=[y(t) x(t)] Then the program computes:
% sigxxJ=E{x(t)*x(t+J)'}
% sigyyJ=E{y(t)*y(t+J)'}
% sigxyJ=E{x(t)*y(t+J)'}
% sigyxJ=E{y(t)*x(t+J)'}
% where x(t) and y(t) evolve as
% x(t+1) = hx x(t) + e(t+1)
% y(t) = gx x(t)
%and Ee(t)e(t)'=varshock
% (c) Stephanie Schmitt-Grohe and Martin Uribe, January 24, 1994
if nargin<4
J=0;
end

[X,D]=eig(hx);
if rank(X)<size(X,1)

%Doubling algorithm
hx_old=hx;
sig_old=varshock;
sigxx_old=eye(size(hx));
diferenz=.1;
while diferenz>1e-25;
sigxx=hx_old*sigxx_old*hx_old'+sig_old;
diferenz = max(max(abs(sigxx-sigxx_old)));
sig_old=hx_old*sig_old*hx_old'+sig_old;
hx_old=hx_old*hx_old;
sigxx_old=sigxx;
end    %while diferenz

else

%Algebraic method
%Get the variance of inv(X)*x
vare=X\varshock *(X\eye(size(X)))';

for i=1:size(hx,1)
for j=1:size(hx,1)
sig(i,j) = 1/(1-D(i,i)*D(j,j)')*vare(i,j);
end   %for j
end   %for i

%Get the variance of x
sigxx=X*sig*X';

end   %if rank


%Get E{x(t)*x(t+J)'}
sigxxJ=hx^(-min(0,J))*sigxx*(hx')^(max(0,J));

%Get E{y(t)*y(t+J)'}
sigyyJ=gx*sigxxJ*gx';


%Get E{y(t)*x(t+J)'}
sigyxJ=gx*sigxxJ;


%Get E{x(t)*y(t+J)'}
sigxyJ=sigxxJ*gx';

sigxxJ = real(sigxxJ);
sigyyJ = real(sigyyJ);
sigyxJ = real(sigyxJ);
sigxyJ = real(sigxyJ);