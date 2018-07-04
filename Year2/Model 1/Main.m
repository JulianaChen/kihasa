%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% MATLAB CODE MARRIAGE AND FERTILITY MODEL %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

%% Set up Parameters
run Setup_Parameters.m

%% Temporary:
z=1;
abi=types(z,1);
edu=types(z,2);
params=params0;

%% Set up State Space
S = sspace(params0,G);

%% Testing Solution
tic;
for z=1:1:G.n_incond
    z
    [C(:,:,:,:,:,z),M(:,:,:,:,:,z),R(:,:,:,:,:,z),N(:,:,:,:,:,z),U(:,:,:,:,:,z)]= solution(G,types(z,1),types(z,2),S,params0);
    toc
end

save solution.mat

%% Drawing Types
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
[c_s,r_s,n_s,u_s,m_s,a_s,wh_s,inv_s,wr_s,wn_s] = simulation(params0,G,S,abi,edu,type,C,M,R,N,U);
toc;

save simulation.mat
