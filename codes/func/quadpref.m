function [Q,W] = quadpref(mu_vec,targets_yvec,ynames,ny,nx,varargin)

% This function specifies the elements of the matrix quadratic preferences
% of the policy controller. Note dimensions and variable ordering follows
% the setup in the main file KLL06model_run.m.
%
% Loss function: L = 0.5*(x'Qx + y'Wy)
%


% Input:
% mu_vec                vector of (nty + ntx) x 1 parameters of loss function parameters
% targets_yvec          nty x 1 string vector of loss function argument names
% ynames                list of y variables
% targets_xvec          ntx x 1 string vector of loss function argument
%                       names: optional
% xnames                list of x's: optional

%%%%%%%%%% Use i = strmatch('Dom. Inflation  ',ynames)

nty = size(targets_yvec,1);
W = zeros(ny,ny);
Q = zeros(nx,nx);

if nargin > 5
    ntx = size(targets_xvec,1);    
else
    ntx = 0;
end


for i = 1:(nty + ntx)
    INDEX = strmatch(targets_yvec(i,:),ynames);
    W(INDEX,INDEX) = mu_vec(i);
    if i > nty
        INDEX = strmatch(targets_xvec(i-nty,:),xnames);
        Q(INDEX,INDEX) = mu_vec(i);
    end    
end