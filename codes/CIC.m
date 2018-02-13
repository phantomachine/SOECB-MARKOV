function [Lower,Upper,CIC2] = CIC(DIST1, DIST2, p);
% -------------------------------------------------------------------------
% Calculate Confidence Interal Criterion (CIC) statistic for DIST1 and
% DIST2.  This measures the overlap of two distributions.
% CIC = 1/(1-w) x integral of DIST2 over the 0.5p to (1-0.5)p confidence 
% interval of DIST1, where p is the confidence level.
% DIST1 and DIST2 are assumed to be nx1 column vectors.
% -------------------------------------------------------------------------
[ndraws1,aa] = size(DIST1);
[ndraws2,bb] = size(DIST2);

nbins   = 100;
[n1,x1] = hist(DIST1,nbins);
[n2,x2] = hist(DIST2,x1);

% calculate lower and upper bounds of DIST1:

Lower = 0;
prob_L = 0;
bin=0;
while prob_L < (p/2)
    bin = bin+1;
    prob_L = prob_L + n1(bin)/ndraws1;
end
Lower = bin;
   
Upper = nbins;
prob_U = 1;
bin = nbins+1;
while prob_U > (1-p/2)
    bin = bin-1;
    prob_U = prob_U - n1(bin)/ndraws1;
end
Upper = bin;

% calculate overlap
CIC2 = 0;
CIC1 = 0;
for i = (Lower+1) : (Upper-1)
    CIC2 = CIC2 + n2(i);
    CIC1 = CIC1 + n1(i);
end;

CIC2 = CIC2/(1-p)/ndraws2;
CIC1 = CIC1/(1-p)/ndraws1;  %check ~ 1


