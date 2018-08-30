clear all; clc;

%% Parameters

% General

Ne = 3; % Gauss-Hermite Points
beta = 0.8; % Discount Rate
r = 0.07; % Interest Rate

n_shocks = 3;
n_period = 20;
n_assets = 15;
n_cons = 30;
n_exp = 10;
n_pop = 100;

Eps=randn(3,n_pop,n_period);

G = struct('Ne',Ne,'beta',beta,'r',r,'Eps',Eps,'n_exp',n_exp,'n_pop',n_pop,...
'n_shocks',n_shocks,'n_period',n_period,'n_assets',n_assets,'n_cons',n_cons);

% Estimated
[eparams, enames] = xlsread('param_vector.xlsx', 'eParams');

% Calbrated
[cparams, cnames] = xlsread('param_vector.xlsx', 'cParams');

params = [eparams;cparams];
names = [enames;cnames];

for k = 1:length(params)
     CurrVarname=cell2mat(names(k));
     CurrValue=params(k);
     P.(CurrVarname)=CurrValue;
end

%% State Space

S = sspace(G,P);

%% Solution

% Model 1: only consumption and assets
[c_star1, v_star1] = solution1(G,P,S);

% Model 2: work, consumption and assets
[c_star2, v_star2, w_star2] = solution2(G,P,S);

% reshape for simulation
for t=1:1:G.n_period-1
	for x = 1:1:G.n_exp
        C(:,x,t) = reshape(c_star2(:,:,x,t),[],1);
        W(:,x,t) = reshape(w_star2(:,:,x,t),[],1);
    end
end


% Model 3: labor, consumption and assets

% Model 4:

% Model 5:
