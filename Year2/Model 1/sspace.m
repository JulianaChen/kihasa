function [S] = sspace(params,G)

%% Shocks

[e, wt] = GaussHermite(G.Ne);
eps_r = sqrt(2)*e*params(19); % error vector
eps_n = sqrt(2)*e*params(20); % error vector
eps_i = sqrt(2)*e*params(21); % error vector

% DON'T NEED THIS NOW, THEY'RE INDEPENDENT
% vcv = [sigma_eps^2,0;0,sigma_v^2];
% detV = det(vcv);
% detV = det(vcv);
% detR = det(R);

shocks_i = kron(eps_i,ones(length(eps_n)*length(eps_r),1));
shocks_r = repmat(kron(eps_r,ones(length(eps_n),1)),[length(eps_i) 1]);
shocks_n = repmat(eps_n,[length(eps_i)*length(eps_r) 1]);
weight = kron(wt,wt); %kron(wt, kron(wt,wt)); % 27x1

% Basis for Income Shocks
zeps_r= 2*(eps_r-eps_r(1))/(eps_r(G.Ne,1)-eps_r(1))-1; 
zeps_n= 2*(eps_n-eps_n(1))/(eps_n(G.Ne,1)-eps_n(1))-1; 
%zeps_i= 2*(eps_i-eps_i(1))/(eps_i(Ne,1)-eps_i(1))-1; 
Teps_r=chebpoly_base(G.Ne-1,zeps_r);
Teps_n=chebpoly_base(G.Ne-1,zeps_n);
%Teps_i=chebpoly_base(Ne-1,zeps_u);
T2eps_r = diag(Teps_r'*Teps_r);
T2eps_n = diag(Teps_n'*Teps_n);
%T2eps_i = diag(Teps_i'*Teps_i);   

%% State Space

matstat = [1 1 0]; %[1 1 0]
children = [1 2 0]; %[1 2 0]
workexp = [0:9];

assets_lb  = -5;
assets_int = 100;
assets_ub  = 500;
assets1 = linspace(assets_lb,assets_int,G.n_assets-3);
assets2 = linspace(assets_int,assets_ub,4);
assets=[assets1,assets2([2:4])];
    
%% SS for linear (only 1 continuous)

SS_A = assets;
SS_X = repmat(workexp, [1 3]);
SS_N = kron(children, ones([1, length(workexp)]));
SS_M = kron(matstat, ones([1, length(workexp)]));

%% output

S = struct(...
    'Teps_r',Teps_r,'Teps_n',Teps_n,'T2eps_r',T2eps_r,'T2eps_n',T2eps_n,...
    'shocks_i',shocks_i,'shocks_r',shocks_r,'shocks_n',shocks_n,'weight',weight,...
    'SS_A',SS_A,'SS_X',SS_X,'SS_M',SS_M,'SS_N',SS_N,'eps_r',eps_r,'eps_n',eps_n);

end
