function [S] = sspace(G,P)

% shocks
[e, wt] = GaussHermite(G.Ne);
eps = sqrt(2)*e*P.sigma_w; 
shocks = eps;
weight = wt;

% basis for shocks
zeps = 2*(eps-eps(1))/(eps(G.Ne,1)-eps(1))-1; 
Teps=chebpoly_base(G.Ne-1,zeps);
T2eps = diag(Teps'*Teps);

% assets
assets_lb  = 0;
assets_int = 8000;
assets_ub  = 40000;

% values for polynomial approximation
[assets,nA,extmin_A,extmax_A,d_A,T_A,T2_A] = cheby_values(G.n_assets,assets_ub,assets_lb);

% asset vector
assets1 = linspace(assets_lb,assets_int,G.n_assets-3);
assets2 = linspace(assets_int,assets_ub,4);
SS_A = [assets1,assets2([2:4])];

% experience
SS_X = [0:9];

S = struct('Teps',Teps,'T2eps',T2eps,'shocks',shocks,'weight',weight,'nA',nA,...
    'extmin_A',extmin_A,'extmax_A',extmax_A,'d_A',d_A,'T_A',T_A,'T2_A',T2_A,...
    'SS_A',SS_A,'SS_X',SS_X);

end
