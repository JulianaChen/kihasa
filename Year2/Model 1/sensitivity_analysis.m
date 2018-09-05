%% Sensitivity Analysis
clear all; clc;

% set version
%version = 'new_marr_04_2';
% paramfile = strcat('params_',version,'.xlsx');

%% Set up Parameters
run Setup_Parameters.m

%% Set up State Space
%params0(34)=4*params0(34); 
S = sspace2(params0,G);
params=params0;
%G.beta=0.85;

%% Select 1 Type and N < 3,000
z = 1;
G.n_pop = 1000; % 1,000 people
abi=types(z,1);
edu=types(z,2);

%% Solution (for a single type)
%[C,M,R,N,U,Ar_out,An_out,Au_out,wh_aux,w_j_r_aux,w_j_n_aux]= solution(G,types(z,1),types(z,2),S,params0);
%[C,M,R,N,U,Ar_out,An_out,Au_out,wh_aux,w_j_r_aux,w_j_n_aux]= solution_cheb(G,types(z,1),types(z,2),S,params0); %% marriage prob - 3 sectors
%[C,M,R,N,U,Ar_out,An_out,Au_out,wh_aux]= solution_cheb2_output(G,types(z,1),types(z,2),S,params0); %% marriage prob added (conditional on simulated asset)
[C,M,R,N,U,Ar_out,An_out,Au_out,wh] = solution_cheb2_output(G,abi,edu,S,params);

M2=squeeze(M(:,1,1,21:30,:));
C1=squeeze(C(:,1,1,1:10,:));
C2=squeeze(C(:,1,1,21:30,:));

%% Simulation (for a single type)

% Draw Types (for 1 type) - only 10
type = ones(G.n_pop,1)*z;
abi = ones(G.n_pop,1)*types(z,1);
edu = ones(G.n_pop,1)*types(z,2);

[c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,exp_s] = simulation(params0,G,S,abi,edu,type,C,M,R,N,U);
%[c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,exp_s] = simulation2(params0,G,S,abi,edu,type,C,M,R,N,U); % fixed assets,inv
%[c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,exp_s] = simulation3(params0,G,S,abi,edu,type,C,M,R,N,U); % fixed assets,inv

%% Plots
%plots(G,c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,abi,edu,type)
t=[1:19]
plot(t,mean(c_s,1),t,mean(a_s(:,(1:19)),1))

plot(sum(c_s)/G.n_pop)
hold on
plot(sum(a_s)/G.n_pop)
hold off
saveas(gcf,'cons.png')

plot(sum(m_s)/G.n_pop)
% %% simulation output
% for n = 1:1:G.n_pop
%     test = [c_s(n,1:19)',r_s(n,1:19)',n_s(n,1:19)',u_s(n,1:19)',m_s(n,1:19)',ch_s(n,1:19)',inv_s(n,1:19)',a_s(n,1:19)',exp(wh_s(n,1:19)'),wr_s(n,1:19)',wn_s(n,1:19)',exp_s(n,1:19)'];
%     sheetname = strcat('Sheet',num2str(n));
%     xlswrite(strcat('results_',version,'.xlsx'),test,sheetname)
% end
% xlswrite(strcat('results_',version,'.xlsx'),params0,'params')
