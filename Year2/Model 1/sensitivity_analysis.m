%% Sensitivity Analysis
clear all; clc;

% set version
% version = 'new_marr_04_2';
% paramfile = strcat('params_',version,'.xlsx');

%% Set up Parameters
run Setup_Parameters.m

%% Set up State Space
%params0(34)=4*params0(34); 
S = sspace3(params0,G);
params=params0;

%% Select 1 Type and N < 3,000
z = 1;
abi=types(z,1);
edu=types(z,2);

G.n_pop = 100; % 1,000 people

%% Solution (for a single type)
%[C,M,R,N,U,Ar_out,An_out,Au_out,wh_aux,w_j_r_aux,w_j_n_aux]= solution(G,types(z,1),types(z,2),S,params0);
%[C,M,R,N,U,Ar_out,An_out,Au_out,wh_aux,w_j_r_aux,w_j_n_aux]= solution_cheb(G,types(z,1),types(z,2),S,params0); %% marriage prob - 3 sectors
[C,M,R,N,U,Ar_out,An_out,Au_out,wh_aux]= solution(G,types(z,1),types(z,2),S,params0); %% marriage prob added (conditional on simulated asset)

M2=squeeze(m_func(:,1,1,21:30,:));
C1=squeeze(c_func(:,1,1,1:10,:));
C2=squeeze(c_func(:,1,1,11:20,:));
C3=squeeze(c_func(:,1,1,21:30,:));
Lr1=squeeze(lr_func(:,1,1,21:30,:));
Ln1=squeeze(ln_func(:,1,1,21:30,:));
Lu1=squeeze(lu_func(:,1,1,21:30,:));

M2=squeeze(M(:,1,1,21:30,:));
R2=squeeze(R(:,1,1,21:30,:));
N2=squeeze(N(:,1,1,21:30,:));
U2=squeeze(U(:,1,1,21:30,:));

C1=squeeze(C(:,1,1,1:10,:));
C2=squeeze(C(:,1,1,11:20,:));
C3=squeeze(C(:,1,1,21:30,:));

%save solution9.mat

%% Simulation (for a single type)

% Draw Types (for 1 type) - only 10
type = ones(G.n_pop,1)*z;
abi = ones(G.n_pop,1)*types(z,1);
edu = ones(G.n_pop,1)*types(z,2);

%[c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,exp_s] = simulation(params0,G,S,abi,edu,type,C,M,R,N,U);
[c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,exp_s] = simulation(params0,G,S,abi,edu,type,C,M,R,N,U); % fixed assets,inv
%[c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,exp_s] = simulation_test(params0,G,S,abi,edu,type,C,M,R,N,U); % fixed assets,inv

plot(sum(c_s)/G.n_pop)
hold on
plot(sum(a_s)/G.n_pop)
hold off
saveas(gcf,'cons.png')

plot(sum(m_s)/G.n_pop)
saveas(gcf,'married.png')

% %% Plots
plots(G,c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,abi,edu,type)
 
% %% simulation output
% for n = 1:1:G.n_pop
%     test = [c_s(n,1:19)',r_s(n,1:19)',n_s(n,1:19)',u_s(n,1:19)',m_s(n,1:19)',ch_s(n,1:19)',inv_s(n,1:19)',a_s(n,1:19)',exp(wh_s(n,1:19)'),wr_s(n,1:19)',wn_s(n,1:19)',exp_s(n,1:19)'];
%     sheetname = strcat('Sheet',num2str(n));
%     xlswrite(strcat('results_',version,'.xlsx'),test,sheetname)
% end
% xlswrite(strcat('results_',version,'.xlsx'),params0,'params')
