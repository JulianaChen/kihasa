%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% MATLAB CODE MARRIAGE AND FERTILITY MODEL %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

%% Structural Parameters

%Family Background
famb1=0.3; %delta1
famb2=0.3; %delta2

%Parental Income Shocks by Age 14
pincshocks=0.3; %delta3

%Private Expenditures in Schooling
pexpschool=4; %delta4

%SD of Shocks
sd_shock=sqrt(0.5);

%Utility of Consumption
sigma=0.4;

%Disutility of work by Sector
psi_r=0.4;
psi_n=0.4;
%unemployed?

%MRS of HH production wrt consumption/leisure
kappa=0.5;

%Value of Marriage in HH Production
theta1=0.7;

%Value of Number of Children in HH Production
theta2=0.3;

%Value of Child HC in HH Production
theta3=0.4;

%Child HC Production Function
gamma1=0.3;
gamma2=0.2;
gamma3=0.2;
phi=0.4;
rho=0.3;

%Female Share of Consumption
delta1=0.3;
delta2=0.3;

%Variance of HH Income Shocks
var_hhincshock=0.2;

%Family Background Types
alpha_f1=4;
alpha_f2=4;

%Return to College
alpha1=0.3;

%Return to General Experiencce
alpha2=0.3;

%Return to Recent Sector Experience
alpha3_r=0.4;
alpha3_n=0.4;
%unemployed?

%Interest rate
r=0.07;

%Discount rate
beta=0.05;

%Investment in Children
inv=1;

%Vector of Initial Parameters
Params0 = [famb1,famb2,pincshocks,pexpschool,sd_shock,...
           sigma,psi_r,psi_n,kappa,theta1,theta2,theta3,...
           gamma1,gamma2,gamma3,phi,rho,delta1,delta2,...
           var_hhincshock,alpha_f1,alpha_f2,alpha1,alpha2,...
           alpha3_r,alpha3_n];

%Bounds
alpha_b=[0.1 0.8];
delta_b=[0.1 0.8];
gamma_b=[0.2 0.9];
theta_b=[0.1 0.9];
bounds = [alpha_b,delta_b,gamma_b,theta_b];

P = struct();

%% Initial Conditions

edu_levels = [1:3];
abi_levels = [1 2];
types = [kron(abi_levels',ones(length(edu_levels),1)) repmat(edu_levels',[length(abi_levels) 1])];

%% Shocks

Ne = 3; % Gauss-Hermite Points
[e, wt] = GaussHermite(Ne);
sigma_r = sqrt(0.05); %shock to regular
sigma_n = sqrt(0.05); %shock to non-regular
sigma_i = sqrt(0.05); %shock to hh income

eps_r = sqrt(2)*e*sigma_r; % error vector
eps_n = sqrt(2)*e*sigma_n; % error vector
eps_i = sqrt(2)*e*sigma_i; % error vector

% DON'T NEED THIS NOW, THEY'RE INDEPENDENT
% vcv = [sigma_eps^2,0;0,sigma_v^2];
% detV = det(vcv);
% detV = det(vcv);
% detR = det(R);

shocks_i = kron(eps_i,ones(length(eps_n)*length(eps_r),1));
shocks_r = repmat(kron(eps_r,ones(length(eps_n),1)),[length(eps_i) 1]);
shocks_n = repmat(eps_n,[length(eps_i)*length(eps_r) 1]);

shocks = [shocks_i shocks_r shocks_n];
weight = kron(wt, kron(wt,wt)); % 27x1

%% General Parameters

n_incond = length(types);
n_shocks = length(shocks);
n_period = 20;
n_pop = 1000;

G = struct('n_incond',n_incond,'n_period',n_period,'n_shocks',n_shocks,'n_pop',n_pop);

%% QUESTIONS:

%PARAMETERS BY AGE AND EDUCATION
%PROBABILITY FUNCTIONS

%% State Space

%Endogenous

assets_lb = 1;
assets_ub = 5;
n_assets = 20;
assets = linspace(assets_lb,assets_ub,n_assets);

matstat = [1 0];
n_matstat = length(matstat);

workexp = [1:10];
workexp_r = [1:3];
workexp_n = [1:3];
n_wrkexp = length(workexp);

sector = [1:3];

% c_min = 0;
% c_max = 5;
c_n = 20;
% c_vector = linspace(c_min,c_max,c_n);

%Exogenous

children = [1 0];

hwages_lb = 1;
hwages_ub = 5;
n_hwages = 5;
hwages = linspace(hwages_lb,hwages_ub,n_hwages);

childK_lb = 1.1;
childK_ub = 4.5;
n_childK = 5;
childK = linspace(childK_lb,childK_ub,n_childK);

%% SS for linear

% SS_K = repmat(childHC',[length(assets)*length(hearnings)*length(workexp)*length(matstat) 1]);
% SS_A = repmat(kron(assets',ones(length(childHC),1)),[length(hearnings)*length(workexp)*length(matstat) 1]);
% SS_H = repmat(kron(hearnings',ones(length(assets)*length(childHC),1)),[length(workexp)*length(matstat) 1]);
% SS_X = repmat(kron(workexp',ones(length(hearnings)*length(assets)*length(childHC),1)),[length(matstat) 1]);
% SS_N = kron(children',ones([length(childHC)*length(assets)*length(hearnings)*length(workexp),1]));
% SS_M = kron(matstat',ones([length(childHC)*length(assets)*length(hearnings)*length(workexp),1]));

%SS = [SS_M SS_N SS_X SS_H SS_A SS_K];

%% SS for chevyshev

SS_K = repmat(childK',[length(assets)*length(hwages) 1]);
SS_A = repmat(kron(assets',ones(length(childK),1)),[length(hwages) 1]);
SS_H = repmat(kron(hwages',ones(length(assets)*length(childK),1)), 1 );
SS_rows = [SS_H SS_A SS_K]; % rows

SS_X = repmat(workexp, [1 length(matstat)]);
SS_M = kron(matstat, ones([1, length(workexp)]));
SS_N = kron(children, ones([1, length(workexp)]));
SS_cols = [SS_M' SS_N' SS_X']; % columns

%% Chevyshev Approximation

[z_A,ext_A,extmin_A,extmax_A,d_A,vector_A,T_A,T2_A] = cheby_values(n_assets,assets_ub,assets_lb);
[z_H,ext_H,extmin_H,extmax_H,d_H,vector_H,T_H,T2_H] = cheby_values(n_hwages,hwages_ub,hwages_lb);
[z_K,ext_K,extmin_K,extmax_K,d_K,vector_K,T_K,T2_K] = cheby_values(n_childK,assets_ub,assets_lb);
    
%% save output
%save solutiontest.mat;
