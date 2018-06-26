%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% MATLAB CODE MARRIAGE AND FERTILITY MODEL %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Structural Parameters

% ADULT STAGE

%Disutility of work by Sector (1-2) (EST)
psi_r=-4;
psi_n=-5;

%Value of Marriage in HH Production (3-5) (EST)
theta1_r=-3;
theta1_n=-3;
theta1_u=0;

%Value of Child HC in HH Production (6-8) (EST)
theta3_r=0.5;
theta3_n=0.5;
theta3_u=1;
theta=[theta1_r;theta1_n;theta1_u;theta3_r;theta3_n;theta3_u];

% Female Share of Consumption (CAL)
% delta=0.5;

%Family Background Types (9-12) (EST)
alpha01_r=4.246;
alpha01_n=3.921;
alpha02_r=0.0995;
alpha02_n=0.0115;

%Return to 2yr College (13-14) (EST)
alpha11_r=0.116;
alpha11_n=0.016;

%Return to 4yrs college(15-16) (EST)
alpha12_r=0.474;
alpha12_n=0.174;

%Return to General Work Experience (17-18) (EST)
alpha2_r=0.437;
alpha2_n=0.302;

alpha=[alpha01_r;alpha01_n;alpha02_r;alpha02_n;alpha11_r;alpha11_n;alpha12_r;alpha12_n;alpha2_r;alpha2_n];

% Shocks (19-21) (EST)
sigma_r = 0.43; %shock to regular
sigma_n = 0.43; %shock to non-regular
sigma_i = 0; %245; %shock to hh income

% Probability of marriage (22-26) (CAL)
omega0_w  =  -0.066; 
omega0_u  =  0.401; 
omega11 = - 0.0581;
omega12 = - 0.0904;
omega2 = 0.015;

omega=[omega0_w;omega0_u;omega11;omega12;omega2];

% Terminal Value function (27-30) (EST)
lambda1=1;
lambda2=1;
lambda3=1;
lambda4=1;

lambda=[lambda1;lambda2;lambda3;lambda4];

% husband wages (regression coefficients on means/sd) (31-40) (CAL)
eta01=6.685;
eta02=0.00000000825;
eta11=0.218;
eta12=0.297;
eta2=0.036;

eta03=0.0113;
eta04=0.0000214;
eta21=0.0097;
eta22=0.0058;
eta3=0.0000349;

eta=[eta01;eta02;eta11;eta12;eta2;eta03;eta04;eta21;eta22;eta3];

% child investment (regression coefficients on means/sd) (41-50)(CAL)
iota01=2.426;
iota02=0.00000000185; 
iota11=0.155;
iota12=0.177;
iota2=0.0971;

iota03=0.0202;
iota04=0.0000751; 
iota21=0.0139;
iota22=0.0084;
iota3=-0.0000297;

iota=[iota01;iota02;iota11;iota12;iota2;iota03;iota04;iota21;iota22;iota3];

% child human capital (51-55) (CAL) 
kappa01=-1.279;
kappa02=-0.0288;
kappa03=0.0586;
kappa04=0.0711;
kappa05=0.00655;

kappa=[kappa01;kappa02;kappa03;kappa04;kappa05];

% probability of losing a regular job (56-60) (CAL)
tau10=2.282;
tau11=-0.155;
tau12=-0.297;
tau13=-0.405;
tau14=-0.815;

% probability of getting a regular job (60-65) (CAL)
tau20=-2.105;
tau21=0.209;
tau22=0.0424;
tau23=-0.324;
tau24=1.271;
  
tau=[tau10;tau11;tau12;tau13;tau14;tau20;tau21;tau22;tau23;tau24];

% probabilities of having a second child (66-77)
phi10 = 3.679;
phi11 = -2.89;
phi12 = -3.197;
phi13 = 1.121;
phi20 = 8.569;
phi21 = -2.528;
phi22 = -4.114;
phi23 = 0.52;
phi30 = 5.692;
phi31 = -0.898;
phi32 = -1.69;
phi33 = -0.379;

phi=[phi10;phi11;phi12;phi13;phi20;phi21;phi22;phi23;phi30;phi31;phi32;phi33];

% Vector of Initial Parameters
params0 = [psi_r;psi_n;theta;alpha;sigma_r;sigma_n;sigma_i;omega;lambda;eta;iota;kappa;tau;phi];

%% Initial Conditions

edu_levels = [1:3];
abi_levels = [1 2];
types = [kron(abi_levels',ones(length(edu_levels),1)) repmat(edu_levels',[length(abi_levels) 1])];

%% General Parameters

%Gauss-Hermite Points
Ne = 3;

%Utility of Consumption
sigma=0.5;

%Discount rate
beta=0.95;

%Interest rate
r=0.07;

% state parameters
n_incond = length(types);
n_shocks = 9; %27;
n_period = 20;
n_pop = 3000;
n_cons = 10; %20;
n_wrkexp = 10;
n_matstat = 3;
n_assets = 10;

% simulation parameters
Eps=randn(3,n_pop,n_period);

%% Save General Parameters

G = struct('Ne',Ne,'sigma',sigma,'beta',beta,'r',r,'Eps',Eps,'n_assets',n_assets,...
    'n_incond',n_incond,'n_period',n_period,'n_shocks',n_shocks,'n_pop',n_pop,...
    'n_cons',n_cons,'n_wrkexp',n_wrkexp,'n_matstat',n_matstat);