
addpath(genpath(pwd))


rng(30)


ISE_total = [];
index_total=[];
toc_time = [];

parfor reprep = 1:100

%% model setting

% sim1
n = 400;  %400
p = 7;
J = 100;
h = n^(-1/5);
tau = 0.25;

% eta is penalty parameter 
eta = 0.5;
sigma =1;


% U vector
U = rand(n,1);
% u vector
u = 1/J : 1/J : 1;

% X design n by p matrix
X = normrnd(0,1, [n,p]);

% sim1
cov = zeros(6,6);
for i = 1:6
    for j = 1:6
        cov(i,j) = 0.5^(abs(i-j));
    end
end
X = mvnrnd(repmat(0,6,1),cov,n);
X = [repmat(1,n,1), X];



%Y sim1         V,V,Z,IV,IV,Z,Z
Y = 2*sin(2*pi*U) + 8* U.* (1-U).* X(:,2) + 2*X(:,4) + 3*X(:,5) + normrnd(0,1, [n,1]);
%Y = 2*sin(2*pi*U) + 8* U.* (1-U).* X(:,2) + 2*X(:,4) + 3*X(:,5) + trnd(3,n,1);

beta_true = zeros(2*J, p);
beta_true(1:2:(2*J), 1) = 2*sin(2*pi*u) + norminv(tau,0,1);   %tinv(tau,3);     %  norminv(tau,0,1);
beta_true(2:2:(2*J), 1) = 4*pi*cos(2*pi*u);
beta_true(1:2:(2*J), 2) = 8*u .* (1-u);
beta_true(2:2:(2*J), 2) = 8- 16*u;
beta_true(1:2:(2*J), 4) = 2;
beta_true(1:2:(2*J), 5) = 3;


% for scaling
X = X/(n);
Y = Y/(n);

% YY is length-nJ vecctor
YY = [];

for i = 1:n
    YY = [YY;  (Y(i) * keru(u - U(i),h))'];
end


% c = (1,0,... 1,0)' is 2J-dim vector
c = repmat([1,0]', [J,1]);

% A is 2J by 2J matrix
A = eye(2*J) -  c * c' / J;


% construct D'
D = [];
for i = 1:n
    di = zeros([J, 2*J]);
    for l = 1:J
        di(l, [2*l-1, 2*l]) = [1, (U(i) - u(l))/h] * keru(u(l) - U(i),h);       
    end
    D = [D ; kron(X(i,:), di)];   
end

%D = D' ; 


DA = D' * D + kron(eye(p), A' * A);
[~, lamm] = eig(DA);
lam_max = max(max(lamm));
lDA = lam_max * eye(2*J*p) - DA;


% We update     beta, zeta, Z, Gamma, gamma
% Z, Gamma : nJ by 1 
% beta, zeta, gamma : 2J by p

% w, Eta: length-p, weight vectors
w = repmat(1,[p,1]);
Eta = w;

Gamma = zeros([n*J,1]);
Z = Gamma ;

gamma = zeros([2*J, p]);
zeta = ones([2*J, p]);


%% initial estimate

est_beta_ini = ini_real(J,p,D,YY,tau);
est_zeta_ini = A * est_beta_ini;



%% Regular quantile estimator (initial one)
eta = 0.5;
sigma =1;
w = repmat(0.0000001, [p,1]);
Eta = repmat(0.0000001, [p,1]);


lam_set = 0.2:0.1:0.6

%% Lasso estimator (our estimator)

beta_ini = est_beta_ini;
zeta_ini = est_zeta_ini;

BIC_set =[];   BIC_set_QKLASSO=[];
BIC_coef = cell(1,8);  BIC_coef_QKLASSO = cell(1,8);
BIC_index = [];  BIC_index_QKLASSO = [];

iii=0;
for lam = lam_set
   iii = iii+1; 
w = [];
for ii = 1:p
  w = [w, scad(norm(beta_ini(:,ii)) / sqrt(J), lam)];
end

Eta = [];
for ii = 1:p
  Eta = [Eta, scad(norm(A * beta_ini(:,ii)) / sqrt(J), lam)];
end

[est_beta, est_zeta, err]= main_admm(Gamma, Z, gamma, A*beta_ini, beta_ini, n, p, J, h, tau, w, Eta, lam_max, YY, D', lDA, DA, eta, A, sigma);

Eta = repmat(0.000000001, [p,1]);

[QKLASSO_beta, QKLASSO_zeta, err]= main_admm(Gamma, Z, gamma, A*beta_ini, beta_ini, n, p, J, h, tau, w, Eta, lam_max, YY, D', lDA, DA, eta, A, sigma);


% V(2),V(2),Z(0),IV(1),IV(1),Z(0),Z(0)
ind_VIZ = [];  ind_VIZ_QKLASSO = []; 
for jj = 1:p
    if norm(est_beta(1:2:end,jj)) < 10^(-3)
        est_beta(1:2:end,jj) = 0;
        ind_VIZ = [ind_VIZ, 0];
    elseif std(est_beta(1:2:end,jj)) < 10^(-3)  &  norm(est_beta(1:2:end,jj)) >= 10^(-3)
        est_beta(1:2:end,jj) = mean(est_beta(1:2:end,jj));
        ind_VIZ = [ind_VIZ, 1];
    else
        ind_VIZ = [ind_VIZ, 2];
    end
    
    if norm(QKLASSO_beta(1:2:end,jj)) < 10^(-3)
        QKLASSO_beta(1:2:end,jj) = 0;
        ind_VIZ_QKLASSO = [ind_VIZ_QKLASSO, 0];
    elseif std(QKLASSO_beta(1:2:end,jj)) < 10^(-3)  &  norm(QKLASSO_beta(1:2:end,jj)) >= 10^(-3)
        QKLASSO_beta(1:2:end,jj) = mean(QKLASSO_beta(1:2:end,jj));
        ind_VIZ_QKLASSO = [ind_VIZ_QKLASSO, 1];
    else
        ind_VIZ_QKLASSO = [ind_VIZ_QKLASSO, 2];
    end    
end
  

tt= YY(:) - D * est_beta(:);
BIC_set = [BIC_set, log(sum((tau * tt.* (tt>0) + (tau-1) * tt.* (tt<0))/J)) + sum(ind_VIZ == 2)*log(n*h)/(2*n*h) + sum(ind_VIZ == 1)*log(n)/(2*n)];       

tt= YY(:) - D * QKLASSO_beta(:);
BIC_set_QKLASSO = [BIC_set_QKLASSO, log(sum((tau * tt.* (tt>0) + (tau-1) * tt.* (tt<0))/J)) + sum(ind_VIZ_QKLASSO == 2)*log(n*h)/(2*n*h) + sum(ind_VIZ_QKLASSO == 1)*log(n)/(2*n)];       

BIC_coef_QKLASSO{iii} = QKLASSO_beta(1:2:end,:);
BIC_coef{iii} = est_beta(1:2:end,:);

BIC_index_QKLASSO = [BIC_index_QKLASSO; ind_VIZ_QKLASSO];
BIC_index = [BIC_index; ind_VIZ];
end

ind_our = find(BIC_set == min(BIC_set));
ind_QKLASSO = find(BIC_set_QKLASSO == min(BIC_set_QKLASSO));


est_beta = BIC_coef{ind_our};
QKLASSO_beta = BIC_coef_QKLASSO{ind_QKLASSO};

%  V(2),V(2),Z(0),IV(1),IV(1),Z(0),Z(0)  (2201100)

est_ind_our = BIC_index(ind_our,:);
est_ind_QKLASSO = BIC_index_QKLASSO(ind_QKLASSO,:);


%% Compute ISE measures
ISE = zeros(3, p);
for ii = 1:p
    ISE(1, ii) = sum((beta_ini(1:2:end,ii) - beta_true(1:2:end,ii)).^2)/J;
    ISE(2, ii) = sum((QKLASSO_beta(:,ii) - beta_true(1:2:end,ii)).^2)/J;
    ISE(3, ii) = sum((est_beta(:,ii) - beta_true(1:2:end,ii)).^2)/J;
end
ISE = [ISE, sum(ISE,2)];
ISE_total = [ISE_total;ISE];


ind_ini = [];
for jj = 1:p
    if norm(beta_ini(1:2:end,jj)) < 10^(-3)
        ind_ini = [ind_ini, 0];
    elseif std(beta_ini(1:2:end,jj)) < 10^(-3)  &  norm(beta_ini(1:2:end,jj)) >= 10^(-3)
        ind_ini = [ind_ini, 1];
    else
        ind_ini = [ind_ini, 2];
    end
end

index_total=[index_total; [ind_ini;  est_ind_QKLASSO; est_ind_our]];
end

save('ttsim1_n400_normal_tau025.mat', 'index_total', 'ISE_total')













