function [S] = sspace(params,G)

%% Shocks

sigma_r = params(23); % shock, regular
sigma_n = params(24); % shock, non-regular
sigma_i = params(25); % shock, unemployed

[e, wt] = GaussHermite(G.Ne);
eps_r = sqrt(2)*e*sigma_r; % error vector
eps_n = sqrt(2)*e*sigma_n; % error vector
eps_i = sqrt(2)*e*sigma_i; % error vector

% vcv = [sigma_eps^2,0;0,sigma_v^2]; % DON'T NEED THIS NOW, THEY'RE INDEPENDENT
% detV = det(vcv);
% detV = det(vcv);
% detR = det(R);

shocks_i = kron(eps_i,ones(length(eps_n)*length(eps_r),1));
shocks_r = repmat(kron(eps_r,ones(length(eps_n),1)),[length(eps_i) 1]);
shocks_n = repmat(eps_n,[length(eps_i)*length(eps_r) 1]);
weight = kron(wt, kron(wt,wt)); % 27x1 % kron(wt,wt);

% Basis for Income Shocks
zeps_r= 2*(eps_r-eps_r(1))/(eps_r(G.Ne,1)-eps_r(1))-1; 
zeps_n= 2*(eps_n-eps_n(1))/(eps_n(G.Ne,1)-eps_n(1))-1; 
zeps_i= 2*(eps_i-eps_i(1))/(eps_i(G.Ne,1)-eps_i(1))-1; 
Teps_r=chebpoly_base(G.Ne-1,zeps_r);
Teps_n=chebpoly_base(G.Ne-1,zeps_n);
Teps_i=chebpoly_base(G.Ne-1,zeps_i);
T2eps_r = diag(Teps_r'*Teps_r);
T2eps_n = diag(Teps_n'*Teps_n);
T2eps_i = diag(Teps_i'*Teps_i);   

%% Discrete Variables

matstat = [1 1 0]; %[1 1 0]
children = [1 2 0]; %[1 2 0]
workexp = [0:9];

%% Continuous Variables

assets_lb  = 0;
assets_int = 8000;
assets_ub  = 28000;
assets1 = linspace(assets_lb,assets_int,G.n_assets-3);
assets2 = linspace(assets_int,assets_ub,4);
%assets = [assets1,assets2([2:4])];

%% Chevyshev Approximation

[assets,nA,extmin_A,extmax_A,d_A,T_A,T2_A] = cheby_values(G.n_assets,assets_ub,assets_lb);

%% SS for chevyshev

SS_A = assets;
SS_X = repmat(workexp, [1 3]);
SS_N = kron(children, ones([1, length(workexp)]));
SS_M = kron(matstat, ones([1, length(workexp)]));

%% output

S = struct(...
    'Teps_i',Teps_i,'Teps_r',Teps_r,'Teps_n',Teps_n,'T2eps_i',T2eps_i,'T2eps_r',T2eps_r,'T2eps_n',T2eps_n,...
    'shocks_i',shocks_i,'shocks_r',shocks_r,'shocks_n',shocks_n,'weight',weight,...
    'nA',nA,'extmin_A',extmin_A,'extmax_A',extmax_A,'d_A',d_A,'T_A',T_A,'T2_A',T2_A,...
    'SS_A',SS_A,'SS_X',SS_X,'SS_M',SS_M,'SS_N',SS_N,'eps_i',eps_i,'eps_r',eps_r,'eps_n',eps_n);

end
