%% Sensitivity Analysis
clear all; clc;

% set version
version = 'July26';
namefile = strcat('results',version,'.xls');

%% Set up Parameters
run Setup_Parameters.m

%% Set up State Space
S = sspace(params0,G);

%% Select 1 Type and N < 3,000
z = 1;
G.n_pop = 10; % only 10 people

%% Solution (for a single type)
[C,M,R,N,U,Ar_out,An_out,Au_out,wh_aux,w_j_r_aux,w_j_n_aux]= solution(G,types(z,1),types(z,2),S,params0);

%% Simulation (for a single type)

% Draw Types (for 1 type) - only 10
type = ones(G.n_pop,1)*z;
abi = ones(G.n_pop,1)*types(z,1);
edu = ones(G.n_pop,1)*types(z,2);

[c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,exp_s] = simulation(params0,G,S,abi,edu,type,C,M,R,N,U);

%% Plots
plots(G,c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,abi,edu,type)

% print parameters used for plots
xlswrite(namefile,params);
