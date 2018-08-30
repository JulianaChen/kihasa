function [S] = sspace(G,P)

% shocks
[e, wt] = GaussHermite(G.Ne);
eps_w = sqrt(2)*e*P.sigma_w; 
shocks = eps_w;
weight = wt;

% basis for shocks
zeps_w = 2*(eps_w-eps_w(1))/(eps_w(G.Ne,1)-eps_w(1))-1; 
Teps_w = chebpoly_base(G.Ne-1,zeps_w);
T2eps_w = diag(Teps_w'*Teps_w);

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

S = struct('Teps_w',Teps_w,'T2eps_w',T2eps_w,'shocks',shocks,'weight',weight,'eps_w',eps_w,...
    'nA',nA,'extmin_A',extmin_A,'extmax_A',extmax_A,'d_A',d_A,'T_A',T_A,'T2_A',T2_A,...
    'SS_A',SS_A,'SS_X',SS_X);

end
