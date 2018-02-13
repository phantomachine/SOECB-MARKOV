function [newdata,dtitle,date_vec] = aus_datadoc(COUNTRY,data,datadate,interval,sampldate,startdate,PLOT_DATA)

% aus_datadoc.m
%
% This script manipulates data imported into aus_estimate.m. It also plots
% raw and doctored data.
% 
% T.Kam and K.Lees, Feb 2, 2006.
%
% Log; Feb 2, 2006: Cruel and unnatural things done to the data - HP
% filtered GDP/person, ToT, RER. Also removed GST spike in inflation data
% using least squares regression.
%
% Functions called:
%   multiplot2D
%   suptitle2
%   hpfilter
%   ols
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Option
fontsizeplus = 8; % increase default font size in function suptitle2

date = sampldate;
frequency = rows(interval);
date_vec = [];
i = 1;
while i < rows(datadate)-1
    date_vec = [date_vec; date*ones(rows(interval),1) + interval];
    date = date+1;
i = i + frequency;
end

date_vec = date_vec(1:end-1);
dtitle = datadate(1,2:end);

data(:,2) = data(:,2)*4; % annualize inflation rate of QonQ imports inflation
data(:,6) = data(:,6)*4; % annualize inflation rate of QonQ CPI inflation
data(:,8) = data(:,8)*4; % annualize inflation rate of QonQ US CPI inflation

if PLOT_DATA
    multiplot2D(data,dtitle,date_vec);
    suptitle2([COUNTRY,'Raw Data ', num2str(date_vec(1)),' to ', num2str(date_vec(end))],fontsizeplus)
end

% Pass data through HP filter to get s, smoothed series for output, TOT and
% RER
lambda = 1600; % smoothing parameter

[ctrend,desvabs] = hpfilter(data(:,1),lambda);
[ytrend,desvabs] = hpfilter(data(:,5),lambda);
[strend,desvabs] = hpfilter(data(:,4),lambda);
[qtrend,desvabs] = hpfilter(data(:,3),lambda);
[usytrend,desvab] = hpfilter(data(:,9),lambda);

cgap = data(:,1)-ctrend;
ygap = data(:,5)-ytrend;
sgap = data(:,4)-strend;
qgap = data(:,3)-qtrend;
usygap = data(:,9)-usytrend;

newdata = data;
newdata(:,1) = cgap;
newdata(:,3) = qgap;
newdata(:,4) = sgap;
newdata(:,5) = ygap;
newdata(:,9) = usygap;

nd = (startdate-sampldate)*rows(interval)+1;
newdata = newdata(nd:end,:);
date_vec = date_vec(nd:end);

% remove GST effect on CPI inflation in 2001
inf = newdata(:,6);
GSTindex = find(inf==3.6566);
inf(GSTindex) = (inf(GSTindex-1)+inf(GSTindex+1))/2;
newdata(:,6) = inf;
% infc=ones(rows(inf),1);
% gstdum=zeros(rows(inf),1);
% gstdum(83)=1; % we know that inf is highest on this GST date
% y=inf;
% x=[infc gstdum];
% results = ols(y,x);  % Need to install LeSage's jplv7 Econometrics Toolbox
% inf_gst = results.resid;
% newdata(:,6) = inf_gst;

if PLOT_DATA
    multiplot2D(newdata,dtitle,date_vec);
    suptitle2([COUNTRY,'Doctored Data ', num2str(date_vec(1)),' to ', num2str(date_vec(end))],fontsizeplus)
end