function marginal = marginal_density(x2,logpri2,loglike2,MH_TRIM)

nruns = size(x2,1);
npara = size(x2,2);
logpost2=logpri2+loglike2;
B = 1;
T = nruns-round(MH_TRIM*nruns)+1;
N = T;
simulations = zeros(N,npara);
lposterior   = zeros(N,1);
nblck=1;
while B <= nblck;
    simulations((B-1)*T+1:B*T,:) = x2(round(MH_TRIM*nruns):nruns,:,B);
    lposterior((B-1)*T+1:B*T) = logpost2(round(MH_TRIM*nruns):nruns,B);
    B = B + 1;
end;    

lpost_mode = max(lposterior);

MU = mean(simulations)';
SIGMA = zeros(npara,npara);
for i=1:N;
    SIGMA = SIGMA + (simulations(i,:)'-MU)*(simulations(i,:)'-MU)';
end;
SIGMA = SIGMA/N;

DetSIGMA = det(SIGMA);
InvSIGMA = inv(SIGMA);
marginal = [];

for p = 0.1:0.1:0.9;
    critval = qchisq(p,npara);
    tmp = 0;
    for i = 1:N;
        deviation  = (simulations(i,:)-MU')*InvSIGMA*((simulations(i,:)-MU'))';
        if deviation <= critval;
	  lftheta = -log(p)-(npara*log(2*pi)+log(DetSIGMA)+deviation)/2;
	  tmp = tmp + exp(lftheta - lposterior(i)+lpost_mode);
        end;    
    end;
    marginal = cat(1,marginal,[p,lpost_mode-log(tmp/N)]); 
end;    
    



