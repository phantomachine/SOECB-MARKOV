function r = randnv(mu,sigma,m,n);
% randnv Random arrays from the normal distribution.
%   R = NORMRND(MU,SIGMA) returns an array of random numbers chosen from a
%   normal distribution with mean MU and standard deviation SIGMA.  The size
%   of R is the common size of MU and SIGMA if both are arrays.  If either
%   parameter is a scalar, the size of R is the size of the other
%   parameter.
%
%   R = NORMRND(MU,SIGMA,M,N,...) or R = NORMRND(MU,SIGMA,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   See also NORMCDF, NORMFIT, NORMINV, NORMLIKE, NORMPDF, NORMSTAT,
%   RANDOM, RANDN.

%   NORMRND uses Marsaglia's ziggurat method.

sigma(sigma < 0) = NaN;
sizeOut = [m,n];
r = randn(sizeOut) .* sigma + mu;
