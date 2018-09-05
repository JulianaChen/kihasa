%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% MATLAB CODE MARRIAGE AND FERTILITY MODEL %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

% set version
version = '_vector';
paramfile = strcat('params',version,'.xlsx');

%% Set up Parameters
run Setup_Parameters.m

%% Select 1 Type and N < 3,000 (for testing):
z=1;
G.n_pop = 1000; % 1,000 people
abi=types(z,1);
edu=types(z,2);
params=params0;

%% Set up State Space
S = sspace(params0,G); % optimal assets for poly
%S = sspace2(params0,G); % kinked asset vector, 3 spots
%S = sspace3(params0,G); % kinked asset vector, 15 spots 
%S = sspace4(params0,G); % linspace low assets

%% Test Solution (only 1 type)
[C,M,R,N,U,Ar_out,An_out,Au_out,wh,W,W2] = solution(G,abi,edu,S,params);
[C,M,R,N,U,Ar_out,An_out,Au_out,wh] = solution_both(G,abi,edu,S,params);

%% Simulation (only 1 type)

% Draw Types (for 1 type)
type = ones(G.n_pop,1)*z;
abi = ones(G.n_pop,1)*types(z,1);
edu = ones(G.n_pop,1)*types(z,2);

[c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,exp_s] = simulation(params0,G,S,abi,edu,type,C,M,R,N,U);

%% Plots
%plots(G,c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,abi,edu,type)
t=[1:19]
plot(t,mean(c_s,1),t,mean(a_s(:,(1:19)),1))

plot(sum(c_s)/G.n_pop)
hold on
plot(sum(a_s)/G.n_pop)
hold off

plot(sum(m_s)/G.n_pop)

%% Solution (loop all 6 types)
tic;
parfor z=1:1:G.n_incond
    z
    [C(:,:,:,:,:,z),M(:,:,:,:,:,z),R(:,:,:,:,:,z),N(:,:,:,:,:,z),U(:,:,:,:,:,z)]= solution(G,types(z,1),types(z,2),S,params0);
    toc
end

%save solution.mat

tic;
parfor z=1:1:G.n_incond
    z
    [C(:,:,:,:,:,z),M(:,:,:,:,:,z),R(:,:,:,:,:,z),N(:,:,:,:,:,z),U(:,:,:,:,:,z)]= solution_cheb2(G,types(z,1),types(z,2),S,params0);
    toc
end

save solution_cheb_newasset.mat

%% Drawing Types

%temporary 
G.n_pop=1000

for n=1:1:G.n_pop
        seed(n)=rand;
        if seed(n)<0.16
            type(n,1)=1;
            abi(n,1)=1;
            edu(n,1)=1;
        elseif seed(n)<0.189 && seed(n)>=0.16
            type(n,1)=2;
            abi(n,1)=1;
            edu(n,1)=2;
        elseif seed(n)<0.215 && seed(n)>=0.189
            type(n,1)=3;
            abi(n,1)=1;
            edu(n,1)=3;
        elseif seed(n)<0.579 && seed(n)>=0.215
            type(n,1)=4;
            abi(n,1)=2;
            edu(n,1)=1;
        elseif seed(n)<0.73 && seed(n)>=0.579
            type(n,1)=5;
            abi(n,1)=2;
            edu(n,1)=2;
        elseif seed(n)<=1 && seed(n)>=0.73
            type(n,1)=6;
            abi(n,1)=2;
            edu(n,1)=3;
        end
end

%% Simulation
tic;
[c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,exp_s] = simulation(params0,G,S,abi,edu,type,C,M,R,N,U);
toc;

%save simulation.mat

tic;
[c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,exp_s] = simulation2(params0,G,S,abi,edu,type,C,M,R,N,U);
toc;

save simulation_cheb_newasset.mat

%% Load Previous Functions & Simulated Data
%load('simulation_July10.mat')

%% Plots
plots(G,c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,abi,edu,type)

xlswrite(strcat('results_',version,'.xlsx'),params0,'params')

%% simulation output
for n = 1:1:10 %only first 10 individuals
    test = [c_s(n,1:19)',r_s(n,1:19)',n_s(n,1:19)',u_s(n,1:19)',m_s(n,1:19)',ch_s(n,1:19)',inv_s(n,1:19)',a_s(n,1:19)',exp(wh_s(n,1:19)'),wr_s(n,1:19)',wn_s(n,1:19)',exp_s(n,1:19)'];
    sheetname = strcat('Sheet',num2str(n));
    xlswrite(strcat('results_',version,'.xlsx'),test,sheetname)
end

%% Estimation (WE NEED TO DISCUSS THE ESTIMATION ALGORITHM)

options = optimoptions('fmincon','Algorithm','interior-point','Display','iter');
[theta,fval,exit,output,lambda]=fminsearch(@(theta) mdf(theta,G,S),params0,[],[],[],[],bounds(:,1),bounds(:,2),[],options);
theta

exit;
