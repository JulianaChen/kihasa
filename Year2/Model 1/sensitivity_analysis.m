%% Sensitivity Analysis
clear all; clc;

% set version
version = 'young_02';
paramfile = strcat('params_',version,'.xlsx');

%% Set up Parameters
run Setup_Parameters.m

%% Set up State Space
params=params0;
S = sspace(params0,G);

%% Select 1 Type and N < 3,000
z = 1;
G.n_pop = 10; % only 10 people

%% Solution (for a single type)
%[C,M,R,N,U,Ar_out,An_out,Au_out,wh_aux,w_j_r_aux,w_j_n_aux]= solution(G,types(z,1),types(z,2),S,params0);
[C,M,R,N,U,Ar_out,An_out,Au_out,wh_aux,w_j_r_aux,w_j_n_aux]= solution_cheb(G,types(z,1),types(z,2),S,params0);

%% Simulation (for a single type)

% Draw Types (for 1 type) - only 10
type = ones(G.n_pop,1)*z;
abi = ones(G.n_pop,1)*types(z,1);
edu = ones(G.n_pop,1)*types(z,2);

[c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,exp_s] = simulation(params0,G,S,abi,edu,type,C,M,R,N,U);

%% Plots
plots(G,c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,abi,edu,type)

%% simulation output
for n = 1:1:G.n_pop
    test = [c_s(n,1:19)',r_s(n,1:19)',n_s(n,1:19)',u_s(n,1:19)',m_s(n,1:19)',ch_s(n,1:19)',inv_s(n,1:19)',a_s(n,1:19)',exp(wh_s(n,1:19)'),wr_s(n,1:19)',wn_s(n,1:19)',exp_s(n,1:19)'];
    sheetname = strcat('Sheet',num2str(n));
    xlswrite(strcat('results_',version,'.xlsx'),test,sheetname)
end
xlswrite(strcat('results_',version,'.xlsx'),params0,'params')
