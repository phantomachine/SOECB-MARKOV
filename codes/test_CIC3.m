% Calculates CIC for policy parameters
clear 
clc

nbins=25;

load mh_dis_diag_muq_zero_aus;

thetabig=d2(:,:);
theta=thetabig(:,[3:11 15 19 22:24 27:37]);

aus_muy=theta(:,15);
aus_mur=theta(:,16);

clear d1;
clear d2;
clear theta;
clear thetabig;

load mh_dis_diag_muq_zero_can;

thetabig=d2(:,:);
theta=thetabig(:,[3:11 15 19 22:24 27:37]);

can_muy=theta(:,15);
can_mur=theta(:,16);

clear d1;
clear d2;
clear theta;
clear thetabig;

load mh_dis_diag_muq_zero_nz;

% New Zealand Case
%d2=Theta_s(700001:1200000,:);
thetabig=d2;

%Remove point parameters
theta=thetabig(:,[3:11 15 19 22:24 27:37]);

nz_muy=theta(:,15);
nz_mur=theta(:,16);

clear Theta_s;
clear loglike_s;
clear logpri_s;
clear theta;
clear d2;


[fay,xay]=ksdensity(aus_muy,'kernel','triangle','width',0.05);
[far,xar]=ksdensity(aus_mur,'kernel','triangle','width',0.05);
[fcy,xcy]=ksdensity(can_muy,'kernel','triangle','width',0.05);
[fcr,xcr]=ksdensity(can_mur,'kernel','triangle','width',0.05);
[fny,xny]=ksdensity(nz_muy,'kernel','triangle','width',0.05);
[fnr,xnr]=ksdensity(nz_mur,'kernel','triangle','width',0.05);

set(gcf,'Color','w');

subplot(2,1,1);
plot(xay,fay,'-',xcy,fcy,'--',xny,fny,':','LineWidth',1.5);
legend('Australia','Canada','New Zealand');
title('Output preference, \mu_{y}', 'fontsize', 12);
xlim([0 2]);

subplot(2,1,2);
plot(xar,far,'-',xcr,fcr,'--',xnr,fnr,':','LineWidth',1.5);
legend('Australia','Canada','New Zealand');
title('Smoothing preference, \mu_{r}', 'fontsize', 12);
xlim([0 2]);

[lyanz,uyanz,y_aus_nz] = CIC(aus_muy, nz_muy, 0.1);
[lycnz,uycnz,y_can_nz] = CIC(can_muy, nz_muy, 0.1);
[lyac,uyac,y_aus_can] = CIC(aus_muy, can_muy, 0.1);

[lranz,uranz,r_aus_nz] = CIC(aus_mur, nz_mur, 0.1);
[lrcnz,urcnz,r_can_nz] = CIC(can_mur, nz_mur, 0.1);
[lrac,urac,r_aus_can] = CIC(aus_mur, can_mur, 0.1);

AUS=[aus_muy aus_mur];
NZ=[nz_muy nz_mur];
CAN=[can_muy can_mur];

% Joint test Australia and New Zealand
[ANX ACX] = hist3(AUS,[nbins nbins]);
AC1=ACX{1};
AC2=ACX{2};

% Need to move from bin centres to edges
AC1=[AC1(1:end) AC1(end)+(AC1(1,2)-AC1(1,1))];
AC2=[AC2(1:end) AC2(end)+(AC2(1,2)-AC2(1,1))];

ACXedge1=AC1-(AC1(1,2)-AC1(1,1))/2;
ACXedge2=AC2-(AC2(1,2)-AC2(1,1))/2;
ACXedge={ ACXedge1 ACXedge2};

AUSNZ=hist3(NZ, 'Edges', ACXedge);
SUMAUSNZ=sumc(sumc(AUSNZ));

% Joint test Canada and New Zealand
[CNX CCX] = hist3(CAN,[nbins nbins]);
CC1=CCX{1};
CC2=CCX{2};

% Need to move from bin centres to edges

CC1=[CC1(1:end) CC1(end)+(CC1(1,2)-CC1(1,1))];
CC1=[CC1(1:end) CC1(end)+(CC1(1,2)-CC1(1,1))];

CCXedge1=CC1-(CC1(1,2)-CC1(1,1))/2;
CCXedge2=CC2-(CC2(1,2)-CC2(1,1))/2;
CCXedge={ CCXedge1 CCXedge2};


CANNZ=hist3(NZ, 'Edges', CCXedge);
SUMCANNZ=sumc(sumc(CANNZ));

% Joint test Australia and Canada
AUSCAN=hist3(CAN, 'Edges', ACXedge);
SUMAUSCAN=sumc(sumc(AUSCAN));

[NNX NCX] = hist3(NZ,[nbins nbins]);

figure;
subplot(2,2,1);
title('Australia');
hist3(AUS,[nbins nbins])

subplot(2,2,2);
title('New Zealand');
hist3(NZ,[nbins nbins])

subplot(2,2,3);
title('Overlap');
hist3(NZ, 'Edges', ACXedge);


