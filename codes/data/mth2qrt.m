function [yq] = mth2qrt(ym)

% mth2qrt.m
% 
% function [yq] = mth2qrt(ym)
%
% Function takes in monthly time series, ym, and produces quarterly
% average, yq
%
% T.Kam, 2006. Use and abuse freely subject to GNU GPL spirit.

yq = []; 
i=1; 
while i < rows(ym) 
    yq = [yq; sum(ym(i:i+2))/length(ym(i:i+2))]; 
i = i+3; 
end