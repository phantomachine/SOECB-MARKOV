% dataprocess.m
%
% Import and manipulate data for file 001data.xls
%


% import and convert monthly AUS POP age 15+ data (RBA G07hist.xls) into
% quarterly avg
aupop15 = xlsread('data/auspop15monthly.xls');
[aupop] = mth2qrt(aupop15(:,2));

% Import australian consumption data
aucons0 = xlsread('data/auscons.xls');
aucons = log(aucons0(:,2)./aupop);  % log per capita consumption


% Import AUS GDP
augdp0 = xlsread('data/ausgdp.xls');
augdp = log(augdp0(:,2)./aupop);

% import and convert monthly AUS 90day irate to quarterly average
auirm = xlsread('data/ausirm.xls');
[auir] = mth2qrt(auirm);
auir = auir/100;
auir = auir(35:end);

% import and convert quarterly AUS CPI into quarterly inflation rate
aup = xlsread('data/auscpi.xls');
aupi = log(aup(2:end))-log(aup(1:end-1));

% AUS terms of trade data
[autot,titles] = xlsread('data/autot.xls');
Pf = autot(:,3);
Ph = autot(:,2);
autot = log(Pf./Ph); autot = autot(2:end); % log(Pf/Ph)
aupif = log(Pf(2:end)) - log(Pf(1:end-1)); % pi_F

% AUS log RER monthly to quarterly avg.
aurerm = xlsread('data/ausrerm.xls');
[aurer] = mth2qrt(aurerm);

% US cpi inflation
uscpi0 = xlsread('data/uscpi.xls');
[uscpi] = mth2qrt(uscpi0);
uspi = uscpi(2:end,:)-uscpi(1:end-1,:);

% US population
uspop0 = xlsread('data/uspop.xls');
[uspop] = mth2qrt(uspop0);

% US GDP
usgdp = xlsread('data/usgdp.xls');
usgdp = log(usgdp*1000./uspop);
usgdp = usgdp(2:end,:);

% US FFR
usffr = xlsread('data/usffr.xls');
usffr = mth2qrt(usffr)/100;
usffr = usffr(2:end,:);

% I copied this data matrix by hand into data/001ausdata.xls

datamat = [aucons, aupif, aurer, autot, augdp, aupi, auir, uspi, usgdp, usffr];

