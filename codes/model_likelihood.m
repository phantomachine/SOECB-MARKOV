function [ L, PROBLEM1,PROBLEM2,PROBLEM3 ] = model_likelihood(Theta,dy,ny,nx,POLICY,nlp)

% model_likelihood.m
%
% This function computes the likelihood function of a model described by
% (Theta,policy) with respect to the dataset dy of dimension T x Ny.
% nlp is no. of lagged policy in state vector yt
%--------------------------------------------------------------------------
% Orginal author: Philip Liu, October 2005 gm_likelihood.m
% Nipped and tucked by: T.Kam, Feb 2006.
% Last: Feb 6, P.Liu. Feb 7, T.Kam
%--------------------------------------------------------------------------
% input: dy = data
%        Theta = vector of variables
%        n = number of endogenous variables
% output: L is the likelihood(dy, Theta)
%
% Function calls:     
% model-solve.m     solves RE model with optimal commitment or discretion policy 
%                   and write solution in state space form
%                       
% doublek.m         to find stationary kalman gain and MSE
%                   matrix of a matrix Riccati equation
%--------------------------------------------------------------------------                   

%% Setting up the preliminaries for the S.S
[T,nz]=size(dy);        % T=no of obs, Ny=ns (no of variables must be equal to the number shocks)
PROBLEM1=0;             % problem in the solution of the model
PROBLEM2=0;             % problem in initializing the kalmand filter
PROBLEM3=0;             % reduced rank var-cov matrix in likelihood loop

% it has to be a way of calculating the steady-state of the state-variables
% given some Theta, x0=initial state for kalman filter

   
if POLICY == 1;pl
    x0=zeros(2*ny+nx-nlp,1);  % order for x(t)=[lambda(t) yt pit qt pict gt ut rst it]
else                    
    x0=zeros(ny+nx-nlp,1);    % order for x(t)=[yt pit qt pict gt ut rst it]
end

Ny = length(x0);


D = zeros(nz);          % D is the serial correlation matrix for n(t)

% H = [Theta(38); 0; 0; 0; Theta(39:41)'; 0; Theta(42)];  % constants

% cc = length(H) - length(find(~H));

% variance covariance
Omega = diag(Theta(end-nz+1:end));

% Check for parameter constraints upon the model

% npara=length(Theta);
% for i=1:npara,
%     if Theta(i) < 0
%         L=-Inf;
%         PROBLEM1 = 1;
%         return;
%     end
% end

%% solving the model

[PP,QQ,PROBLEM1] = model_solve(Theta,POLICY);

% Delete last nlp # state variable in the solution y' = PPy + QQz which is the
% repeated equation for the policy instrument/state to avoid covariance
% matrix singularity in the Kalman gain calculation
    PP = PP(1:end-nlp,1:end-nlp);
    QQ = QQ(1:end-nlp,:);
    
if PROBLEM1 == 1
        L=-Inf;
        return;
end

%% rewritting the state-space system into innovation representation 
% Reference: hansen and sargent manuscript, chapter 9
%==========================================================================
% state space: X(t) = AoX(t-1) + Cw(t); Y(t) = G X(t) + v(t)
%              where w(t) and v(t) are Martingale difference sequences
%              with E(w(t)w(t)')=I, E(v(t)v(t)')=R
%              and  E(X(t)X(t)')=Sigma
%--------------------------------------------------------------------------
% innovation: Xhat(t) = AoXhat(t-1) + K(t)a(t); Y(t) = G Xhat(t) + a(t)
%             where a(t) = Y(t) - GXhat(t), with
%             E(a(t)a(t)')=Omega=G*Sigma*G'+R and K(t) is the Kalman gain
%==========================================================================
C=QQ*Omega;
V=C*C';
Ao=PP;
R= diag(Theta(end-nz+1:end));        %zeros(nz,nz);         % E(vtvt')

% POLICY=0, state vector = x(t)=[lambda(t) yt pit qt pict gt ut rst it]
% POLICY=1, state vector = x(t)=[yt pit qt pict gt ut rst it]

G = zeros(nz,Ny); % observation equation - picks out observables from state vector

if POLICY == 1            %commitment
    nlm = ny;
else                      %discretionary
    nlm = 0;
end

% G(1,1+nlm) = 1;           % Consumption == 1
G(1,3+nlm) = 4;           % pi_F == 3
G(2,4+nlm) = 1;           % q
G(3,5+nlm) = 1;           % s
G(4,6+nlm) = 1;           % y
G(5,7+nlm) = 4;           % pi
G(6,9+nlm) = 4;           % r
G(7,15+nlm) = 4;          % pi*
G(8,16+nlm) = 1;          % y*
G(9,17+nlm) = 4;         % r*


%% Initializing the Kalman filter
V1=V;               % V1=CC'
V2=R;               % V2 = E(vtvt') zero in this case

% use the limiting matrix of K(t) and Sigma(t) as initial values
% more general verison 

%[K, Sigma, PROBLEM2]=kfilter(Ao, G, V1, V2);

[Sigma, PROBLEM2] = kalman_initial(Ao, G, V);

%[K, Sigma,PROBLEM2] = doublek(Ao,G,V1,V2);
if PROBLEM2==1;
       L=-Inf;
       return;
end

% vec(Sigma1|0)=(I-FOF)^-1 vec(CC')
%  exog=size(Ao,1);
%  Sigma     =inv(eye((exog^2))-kron(Ao,Ao))*reshape(V,(exog^2),1);    
%  Sigma     =reshape(Sigma,exog,exog);

%save Sigma

%% likelihood function using the Kalman filter
% L=-0.5(sum_0^T{ Ny*ln(2pi) + ln|Omega(t)| + a(t)inv(Omega(t))a(t)

xt        = x0';                    % initial x(0)
innov     = dy(1,:)-xt*G';          % initial a(0)  % constant here: May 10
L         = 0;
for i=2:T;
    
    Omega   = R + G*Sigma*G';                                 %Omega(t) = G*Sigma(t)*G'+R
%           save Omega Omega G Sigma
    if det(Omega) <= 0; %| rank(Omega) < length(Omega);       % OR conditional statement, TK Feb 6, 2006
        L=-Inf;
        PROBLEM3 = 1;
        return;
    end
    

    L       = L + log(det(Omega))+innov*inv(Omega)*innov';
    K       = Ao*Sigma*G'*inv(G*Sigma*G'+R);                 %K(t) = Ao*Sigma(t)*G'*inv(G*Sigma(t)*G'+R); (eqn:9.1.9)
    Sigma   = Ao*Sigma*Ao'+V-K*G*Sigma*Ao';                    %Sigma(t+1)= Ao*Sigma(t)*Ao'+C*C'-K*G*Sigma(t)*Ao'; (eqn:9.1.11)
    xt      = xt*Ao'+innov*K';                                 %current state
    innov   = dy(i,:)-xt*G' ;                                   %current innovation a_t % constant here: May 10
end
    
L              = -0.5*(L + log(det(Omega)) + innov*inv(Omega)*innov' + Ny*T*log(2*pi)); %add last observation T

% Reference: Hansen, McGrattan and Sargent, Mechanics of forming and estimating 
% dynamic linear economies, Minneapolis Fed Staff report 182
% Given the model of the form:
% x(t+1) = Ao x(t) + C w(t+1) [State equation]
%         under commitment (where policy=0)
%               x(t) = AAx(t-1) + BB w(t);
%               AA=(n+1)x(n+1), BB=(n+1)xns
%               x(t)=[lambdat yt xt]'; w(t)=[eq, eg, eu, ers]';                
%         under discretion (where policy=1)
%               x(t) = AAx(t-1) + BB w(t);    
%               AA=(2n+1)x(2n+1), BB=(2n+1)xns
%               x(t)=[yt xt]'; w(t)=[eq, eg, eu, ers]';
% z(t) = G x(t) + v(t) [Measurement equation]
% where G picks out the observed variables
% v(t) = D v(t-1) + n(t), E[n(t)n(t)']=R, 
% E[w(t+1)n(s)]=0 [Serially correlated shocks]
