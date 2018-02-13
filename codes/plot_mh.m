clear;
load mh_dis;
Theta_s(:,12:16)=[];
Theta_s(:,6)=[];
Theta_s(:,3)=[];


chain=Theta_s;
likelihood=loglike_s;
logprior=logpri_s;
   
trim=1/3;
ntrim=trim*length(chain);
chaintr=chain(ntrim:end,:);
chain_ci=prctile(chaintr, [2.5 50 97.5]);


s1 = strvcat('\sigma', '\eta', '\theta','\phi', '\rho_g', '\rho_{u}', '\rho_{r^*}', 'W_y', 'q'); 
figure;
chaintrd=chaintr;
i=1;
while i<=9;
   subplot(3,3,i)
     plot(chaintrd(:,i)) ;  
     axis([1 length(chaintrd) min(chaintrd(:,i)) max(chaintrd(:,i))]);
     title(s1(i,:), 'fontsize', 12);
   i=i+1;
end

