function [newTheta] = genTheta(Theta, mean_Theta, var_Theta, npara, bounds)

%this function generates new draws and check again bounds

lenb=1;

while lenb==1,

    newTheta = Theta + normrnd(mean_Theta,var_Theta,1,npara);  

    check=sum(newTheta'<=bounds(:,1)| newTheta'>=bounds(:,2));

    if check ~= 0;

        lenb=1;

    else

        lenb=0;

    end

end



% If you are calibrating certain parameters fix them here. 
newTheta(1)=0.99; % pshape = 0
newTheta(2)=0.45; %pshape = 0

% newTheta(26)= 0; % Set mu_q =0 with Prob(1), pshape = 0

% newTheta(7) = newTheta(8); % delta_h and delta_f

newTheta(12) = 0;      % pi-output a2
newTheta(13) = 0;     % pi-interest a3

newTheta(14) = 0;      % output-pi b1
newTheta(16) = 0;    % output-interest b3

newTheta(17) = 0;     % interest-pi c1
newTheta(18) = 0;     % interest-output c2
newTheta(25) = 0;       % rhor

newTheta(20:21) = 0;   % turn off pho_h and rho_f
