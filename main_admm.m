function [beta, zeta, err]= main_admm(Gamma, Z, gamma, zeta, beta, n, p, J, h, tau, w, Eta, lam_max, YY, D, lDA, DA, eta, A, sigma)

err = [];

for ijk = 1:100  
    
    Gammai = Gamma;
    Zi = Z;
    gammai = gamma ;
    zetai = zeta;
    betai =  beta;

    
%Step 1: update beta

beta = []; 

for k = 1:p    
  T = -D((2*J*(k-1) + 1) : (2*J*k), :) * Gamma + A * gamma(:,k) - ...
    eta * (D((2*J*(k-1) + 1) : (2*J*k),:) * (Z - YY) - A * zeta(:,k) - lDA((2*J*(k-1) + 1) : (2*J*k), :) * betai(:));
  
  no_val = norm(T);
  if no_val <=  w(k) / sqrt(J)
      beta = [beta, repmat(0, 2*J, 1)];
  end
  if no_val >  w(k) / sqrt(J)
        uk = (no_val - w(k) / sqrt(J)) / (eta * lam_max);
       beta = [beta,  T / (w(k)/ (uk * sqrt(J)) + eta * lam_max)];
  end
end
  

err_beta = norm(betai-beta, 'fro')^2;  
  

 
% Step 2: update zeta

zeta = []; 
for k = 1:p
  T = - gamma(:,k)  + eta * A * beta(:,k);
  no_val = norm(T);
  if no_val <=  Eta(k) / sqrt(J)
      zeta = [zeta, repmat(0, 2*J, 1)];
  end
  if  no_val >  Eta(k) / sqrt(J)
        uk = (no_val - Eta(k) / sqrt(J)) / (eta);
       zeta = [zeta, T / (Eta(k)/ (uk * sqrt(J)) + eta) ];
  end
end
  
err_zeta = norm(zetai-zeta, 'fro')^2;  



% Step 3: update Z (nJ-dim)

DB = D' *  beta(:);

Z = zeros([n*J, 1]);
for i = 1:(n*J)   
  uk = -Gamma(i) - eta * (-YY(i) + DB(i));
  Z(i) = (uk > tau/J) *  (uk - tau/J) / eta  +  (uk < (tau-1)/J) *  (uk - (tau-1)/J) / eta;
end


err_Z = norm(Zi - Z)^2;


% Step 4: update Gamma (nJ by 1) and gamma (2J by p)

Gamma = Gamma + sigma * eta * (Z - YY + DB);
gamma = gamma + sigma * eta * (zeta - A * beta);

err_Gamma = norm(Gammai - Gamma)^2;
err_gamma = norm(gammai - gamma)^2;


err0 = err_beta + err_zeta  + err_Z + err_Gamma + err_gamma;
err_dual = norm(zeta - A * beta, 'fro')^2 ;
err_Z = norm(Z - YY + D' * beta(:))^2 ;
err1 = sqrt(err0 + err_dual + err_Z) ;

err = [err, err1];

if err1 < 10^(-5)
    break
end

end

