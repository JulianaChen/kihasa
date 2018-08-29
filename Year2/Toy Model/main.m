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

G = struct('Ne',Ne,'beta',beta,'r',r,...
'n_shocks',n_shocks,'n_period',n_period,'n_assets',n_assets,'n_cons',n_cons,...
'n_exp',n_exp);

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

% Model 3: labor, consumption and assets
% Model 4:
% Model 5: