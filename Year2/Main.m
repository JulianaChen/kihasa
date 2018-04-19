%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% MATLAB CODE MARRIAGE AND FERTILITY MODEL %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

%% Structural Parameters

% EDUCATION STAGE

% %Family Background
% famb1=0.3; %delta1
% famb2=0.3; %delta2
% 
% %Parental Income Shocks by Age 14
% pincshocks=0.3; %delta3
% 
% %Private Expenditures in Schooling
% pexpschool=4; %delta4
% 
% %SD of Shocks
% sd_shock=sqrt(0.5);

% ADULT STAGE

%Disutility of work by Sector (1-2)
psi_r=-4;
psi_n=-4;

%Value of Marriage in HH Production (3-5)
theta1_r=-3;
theta1_n=-3;
theta1_u=0;

%Value of Child HC in HH Production (6-8)
theta3_r=0.5;
theta3_n=0.5;
theta3_u=1;
theta=[theta1_r;theta1_n;theta1_u;theta3_r;theta3_n;theta3_u];

%Child HC Production Function [table 3, CGST_Oct82015] (9-10)
gamma1=0.5595;
phi=0.4282;

% %Female Share of Consumption
% delta1=0.3;
% delta2=0.3;

%Family Background Types (11-14)
alpha01_r=4.246;
alpha01_n=3.921;
alpha02_r=0.0995;
alpha02_n=0.0115;

%Return to 2yr College (15-16)
alpha11_n=0.116;
alpha11_r=0.016;

%Return to 4yrs college(17-18)
alpha12_n=0.474;
alpha12_r=0.174;

%Return to Recent Sector Experience (19-20)
alpha2_r=0.437;
alpha2_n=0.302;
%unemployed?
alpha=[alpha01_r;alpha01_n;alpha02_r;alpha02_n;alpha11_r;alpha11_n;alpha12_r;alpha12_n;alpha2_r;alpha2_n];

% Shocks (21-23)
sigma_r = 0.43; %shock to regular
sigma_n = 0.43; %shock to non-regular
sigma_i = 0; %245; %shock to hh income

% Probability of marriage (24-28)
omega0_w  =  0.3349; 
omega0_u  =  0.401; 
omega11 = - 0.0581;
omega12 = - 0.0904;
omega2 = 0.005;
omega=[omega0_w;omega0_u;omega11;omega12;omega2];

% Terminal Value function (29-31)
lambda1=1;
lambda2=1;
lambda3=1;
lambda=[lambda1;lambda2;lambda3];
%Vector of Initial Parameters
params0 = [psi_r;psi_n;gamma1;phi;theta;alpha;sigma_r;sigma_n;sigma_i;omega;lambda];

% %Bounds
% alpha_b=[0.1 0.8];
% delta_b=[0.1 0.8];
% gamma_b=[0.2 0.9];
% theta_b=[0.1 0.9];
% bounds = [alpha_b,delta_b,gamma_b,theta_b];

%% Initial Conditions

edu_levels = [1:3];
abi_levels = [1 2];
types = [kron(abi_levels',ones(length(edu_levels),1)) repmat(edu_levels',[length(abi_levels) 1])];

%% General Parameters
Ne = 3; % Gauss-Hermite Points
%Utility of Consumption
sigma=0.5;
%Discount rate
beta=0.95;
%Interest rate
r=0.07;
%Investment in Children
Inv=3;
% state parameters
n_incond = 1; length(types);
n_shocks = 9; %27;
n_period = 20;
n_pop = 3000;
n_cons = 10; %20;
n_wrkexp = 10;
n_matstat = 2;
n_assets = 10;
n_hwages = 3;
n_childK = 3; 
n_SS = n_assets*n_hwages*n_childK ; %500;
% simulation parameters
Eps=randn(3,n_pop,n_period);

G = struct('Ne',Ne,'sigma',sigma,'beta',beta,'r',r,'Inv',Inv,'Eps',Eps,...
    'n_incond',n_incond,'n_period',n_period,'n_shocks',n_shocks,'n_SS',n_SS,...
    'n_pop',n_pop,'n_cons',n_cons,'n_wrkexp',n_wrkexp,'n_matstat',n_matstat,...
    'n_assets',n_assets,'n_hwages',n_hwages,'n_childK',n_childK);


% %% Drawing Types
% 
% for n=1:1:G.n_pop
%         seed(n)=rand;
%         if seed(n)<0.16
%             type(n,1)=1;
%             abi(n,1)=1;
%             edu(n,1)=1;
%         elseif seed(n)<0.189 && seed(n)>=0.16
%             type(n,1)=2;
%             abi(n,1)=1;
%             edu(n,1)=2;
%         elseif seed(n)<0.215 && seed(n)>=0.189
%             type(n,1)=3;
%             abi(n,1)=1;
%             edu(n,1)=3;
%         elseif seed(n)<0.579 && seed(n)>=0.215
%             type(n,1)=4;
%             abi(n,1)=2;
%             edu(n,1)=1;
%         elseif seed(n)<0.73 && seed(n)>=0.579
%             type(n,1)=5;
%             abi(n,1)=2;
%             edu(n,1)=2;
%         elseif seed(n)<=1 && seed(n)>=0.73
%             type(n,1)=6;
%             abi(n,1)=2;
%             edu(n,1)=3;
%         end
% end
% 
% % husband wages
% wh_s = 1 + (20-1).*rand(G.n_pop,1);

%% Test Functions
S = sspace_small_newcgrid(params0,G);
%S = sspace(params0,G);
tic;
for z=1:1:n_incond
    z
    [V_lin(:,:,:,:,:,:,:,z),C(:,:,:,z),C_lin(:,:,:,:,:,:,:,z),M(:,:,:,z),M_lin(:,:,:,:,:,:,:,z),R(:,:,:,z),R_lin(:,:,:,:,:,:,:,z),N(:,:,:,z),N_lin(:,:,:,:,:,:,:,z),U(:,:,:,z),U_lin(:,:,:,:,:,:,:,z)] = solution_newcgrid_rsp_linear(G,types(z,1),types(z,2),S,params0);
    toc
end
toc;
save solution_linear2
tic;
for z=1:1:n_incond
    z;
    [V_lin(:,:,:,:,:,:,:,z),C(:,:,:,z),C_lin(:,:,:,:,:,:,:,z),M(:,:,:,z),M_lin(:,:,:,:,:,:,:,z),R(:,:,:,z),R_lin(:,:,:,:,:,:,:,z),N(:,:,:,z),N_lin(:,:,:,:,:,:,:,z),U(:,:,:,z),U_lin(:,:,:,:,:,:,:,z)] = solution_newcgrid_rsp(G,types(z,1),types(z,2),S,params0);
    toc
end
toc;
save solution_polynomial2
 
tic;
for z=1:n_incond
    z
    %[alpC(:,:,:,z),alpR(:,:,:,z),alpN(:,:,:,z),alpU(:,:,:,z),alpM(:,:,:,z)]=polfunc_approx(C(:,:,:,z),R(:,:,:,z),N(:,:,:,z),U(:,:,:,z),M(:,:,:,z),S,G);
    [alpC(:,:,:,z),alpR(:,:,:,z),alpN(:,:,:,z),alpU(:,:,:,z),alpM(:,:,:,z)]=polfunc_approx_small(C(:,:,:,z),R(:,:,:,z),N(:,:,:,z),U(:,:,:,z),M(:,:,:,z),S,G);
    toc
end
toc;

tic;
[c_slin,r_slin,n_slin,u_slin,m_slin,a_slin,k_slin,wr_slin,wn_slin] = simulation(params0,alpC,alpR,alpN,alpU,alpM,G,S,abi,edu,type,wh_s,C_lin,M_lin,R_lin,N_lin,U_lin);
toc;
save simulation_smallNov27_linapprox2
%% save output
%save solutiontest.mat;