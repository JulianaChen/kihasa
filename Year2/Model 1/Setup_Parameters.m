%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% MATLAB CODE MARRIAGE AND FERTILITY MODEL %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Structural Parameters

% ADULT STAGE

% %Disutility of work by Sector (1-2) (EST)
% psi_r=-1.4;
% psi_n=-1.5;
% 
% %Value of Marriage in HH Production (3-5) (EST)
% theta1_r=5; %error: gamma1=0.5595, old: -3;
% theta1_n=5; %error: phi=0.4282, old: -3;
% theta1_u=0; %error: theta1_r, old: 0;
% 
% %Value of Child HC in HH Production (6-8) (EST)
% theta3_r=3; %error: theta1_n=-3, old: 0.5;
% theta3_n=3; %error: theta1_u=0, old: 0.5;
% theta3_u=5; %error: theta3_r=0.5, old: 1;
% 
% theta=[theta1_r;theta1_n;theta1_u;theta3_r;theta3_n;theta3_u];
% 
% % Female Share of Consumption (CAL)
% delta=0.5;
% 
% %Family Background Types (9-12) (EST)
% alpha01_r=4.246;
% alpha01_n=3.921;
% alpha02_r=0.0995;
% alpha02_n=0.0115;
% 
% %Return to 2yr College (13-14) (EST)
% alpha11_r=0.116;
% alpha11_n=0.016;
% 
% %Return to 4yrs college(15-16) (EST)
% alpha12_r=0.474;
% alpha12_n=0.174;
% 
% %Return to General Work Experience (17-18) (EST)
% alpha2_r=0.437;
% alpha2_n=0.302;
% 
% %Family Background Types (9-12) (EST)
% alpha01_r=6.887;
% alpha01_n=6.54;
% alpha02_r=0.059;
% alpha02_n=0.103;
% 
% %Return to 2yr College (13-14) (EST)
% alpha11_r=0.115;
% alpha11_n=0.0062;
% 
% %Return to 4yrs college(15-16) (EST)
% alpha12_r=0.391;
% alpha12_n=0.194;
% 
% %Return to General Work Experience (17-18) (EST)
% alpha2_r=0.181;
% alpha2_n=0.246;
%  
% alpha=[alpha01_r;alpha01_n;alpha02_r;alpha02_n;alpha11_r;alpha11_n;alpha12_r;alpha12_n;alpha2_r;alpha2_n];
% 
% % Shocks (19-21) (EST)
% sigma_r=0.43; %shock to regular
% sigma_n=0.43; %shock to non-regular
% sigma_i=0; %245; %shock to hh income
%  
% % Probability of marriage (22-26) (CAL) - younger cohort only
% omega0_w=-1.427; 
% omega0_u=-15.30; 
% omega11=-0.411;
% omega12=-1.254;
% omega2=0.548;

% % Probability of marriage (22-26) (CAL)
% omega0_w=-1.198; old
% omega0_u=-4.657; old
% omega11=0.865; old
% omega12=-1.179; old
% omega2=0.211; old
% 
% omega0_r=-68.53
% omega11=-0.566
% omega12=-0.795
% omega13=20.22
% omega14=-0.152
% 
% omega0_n=-73.30
% omega21=-1.437
% omega22=-2.299
% omega23=21.20
% omega24=0.0625
% 
% omega0_u=-62.76
% omega31=-0.11
% omega32=-2.324
% omega33=19.08
% omega34=-0.124
%     
% omega=[omega0_r;omega11;omega12;omega13;;omega14;omega0_n;omega21;omega22;omega23;omega24;omega0_u;omega31;omega32;omega33;omega34];
% 
% % Terminal Value function (27-30) (EST)
% lambda1=1;
% lambda2=1;
% lambda3=1;
% lambda4=1;
% 
% lambda=[lambda1;lambda2;lambda3;lambda4];

% husband wages (regression coefficients on means/sd) (31-40) (CAL)
% eta01=6.685;
% eta02=0.00000000825;
% eta11=0.218;
% eta12=0.297;
% eta2=0.036;
% 
% eta03=0.0113;
% eta04=0.0000214;
% eta21=0.0097;
% eta22=0.0058;
% eta3=0.0000349;
% 
% eta01=3.747;
% eta02=0.0296;
% eta11=0.209;
% eta12=0.231;
% eta2=3.747;
% 
% eta03=-0.00391;
% eta04=-0.00679;
% eta21=0.0171;
% eta22=0.00909;
% eta3=0.0413;
% 
% eta=[eta01;eta02;eta11;eta12;eta2;eta03;eta04;eta21;eta22;eta3];

% % child investment (regression coefficients on means/sd) (41-50)(CAL)
% iota01=2.426;
% iota02=0.00000000185; 
% iota11=0.155;
% iota12=0.177;
% iota2=0.0971;
% 
% iota03=0.0202;
% iota04=0.0000751; 
% iota21=0.0139;
% iota22=0.0084;
% iota3=-0.0000297;
% 
% iota01=3.018;
% iota02=0.0936; 
% iota11=0.309;
% iota12=0.423;
% iota2=1.386;
% 
% iota03=0.0172;
% iota04=-0.00567; 
% iota21=0.12;
% iota22=0.00625;
% iota3=0.000375;
% 
% iota=[iota01;iota02;iota11;iota12;iota2;iota03;iota04;iota21;iota22;iota3];

% % child human capital (51-55) (CAL) 
% kappa01=-1.275;
% kappa02=-0.019;
% kappa03=0.0555;
% kappa04=0.115;
% kappa05=0.00769;
% 
% kappa=[kappa01;kappa02;kappa03;kappa04;kappa05];
% 
% % probability of losing a regular job (56-60) (CAL)
% tau10=1.535;
% tau11=-0.15;
% tau12=-0.277;
% tau13=-0.0241;
% tau14=-0.823;

% % probability of getting a regular job (60-65) (CAL)
% tau20=-0.979;
% tau21=0.415;
% tau22=-0.0738;
% tau23=-0.0436;
% tau24=0.602;
%   
% tau=[tau10;tau11;tau12;tau13;tau14;tau20;tau21;tau22;tau23;tau24];
% 
% % probabilities of having a second child (66-77) - younger cohort only
% phi10=0.305;
% phi11=-1.022;
% phi12=-1.256;
% phi13=-0.32;
%  
% phi20=-0.287;
% phi21=0.0994;
% phi22=-0.69;
% phi23=-0.587;
%  
% phi30=-0.194;
% phi31=-0.0853;
% phi32=-0.245;
% phi33=-0.173;

% % probabilities of having a second child (66-77)
% phi10=0.789;
% phi11=-3.452;
% phi12=-3.535;
% phi13=0.544;
% 
% phi20=2.916;
% phi21=-1.544;
% phi22=-2.047;
% phi23=-0.386;
% 
% phi30=0.569;
% phi31=-0.603;
% phi32=-0.849;
% phi33=-0.0719;
% 
% phi=[phi10;phi11;phi12;phi13;phi20;phi21;phi22;phi23;phi30;phi31;phi32;phi33];
% 
% % Vector of Initial Parameters
% params0 = [psi_r;psi_n;theta;alpha;sigma_r;sigma_n;sigma_i;omega;lambda;eta;iota;kappa;tau;phi];

%% Estimated and Calibrated Parameters (import from excel)

% [params0, names] = xlsread(paramfile);
% for k = 1:length(params0)
%      CurrVarname=cell2mat(names(k));
%      CurrValue=params0(k);
%      P.(CurrVarname)=CurrValue;
% end

[eparams, enames] = xlsread('param_vector.xlsx', 'eParams');
[cparams, cnames] = xlsread('param_vector.xlsx', 'cParams');

params0 = [eparams;cparams];
names0 = [enames;cnames];

for k = 1:length(params0)
     CurrVarname=cell2mat(names0(k));
     CurrValue=params0(k);
     P.(CurrVarname)=CurrValue;
end

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
beta=0.8;

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
