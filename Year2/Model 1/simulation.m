function [c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,exp_s] = simulation(params,G,S,abi,edu,type,C,M,R,N,U)

%% Parameters

% alpha01_r=params(9); % wage return for the ability type, regular
% alpha01_n=params(10); % wage return for the ability type, non-regular
% alpha02_r=params(11); % additional return for high ability type, regular
% alpha02_n=params(12); % additional return for high ability type, non-regular
% alpha11_r=params(13); % wage return to 2yr college, regular
% alpha11_n=params(14); % wage return to 2yr college, non-regular
% alpha12_r=params(15); % wage return to 4yr college, regular
% alpha12_n=params(16); % wage return to 4yr college, non-regular
% alpha2_r=params(17); % wage return to general work experience, regular
% alpha2_n=params(18); % wage return to general work experience, non-regular
% 
% sigma_r = params(19); % shock, regular
% sigma_n = params(20); % shock, non-regular
% sigma_i = params(21); % shock, unemployed

alpha01_r=params(13);
alpha01_n=params(14);
alpha02_r=params(15);
alpha02_n=params(16);
alpha11_r=params(17);
alpha11_n=params(18);
alpha12_r=params(19);
alpha12_n=params(20);
alpha2_r=params(21);
alpha2_n=params(22);
sigma_r=params(23);
sigma_n=params(24);
sigma_i=params(25);

%     phi10 = params(76); % probability of second child
%     phi11 = params(77); % probability of second child
%     phi12 = params(78); % probability of second child
%     phi13 = params(79); % probability of second child
%     phi20 = params(80); % probability of second child
%     phi21 = params(81); % probability of second child
%     phi22 = params(82); % probability of second child
%     phi23 = params(83); % probability of second child
%     phi30 = params(84); % probability of second child
%     phi31 = params(85); % probability of second child
%     phi32 = params(86); % probability of second child
%     phi33 = params(87); % probability of second child
    
% phi10 = params(66); %3.679;
% phi11 = params(67); %-2.89;
% phi12 = params(68); %-3.197;
% phi13 = params(69); %1.121;
% phi20 = params(70); %8.569;
% phi21 = params(71); %-2.528;
% phi22 = params(72); %-4.114;
% phi23 = params(73); %0.52;
% phi30 = params(74); %5.692;
% phi31 = params(75); %-0.898;
% phi32 = params(76); %-1.69;
% phi33 = params(77); %-0.379;

phi10=params(80);
phi11=params(81);
phi12=params(82);
phi13=params(83);
phi14=params(84);
phi20=params(85);
phi21=params(86);
phi22=params(87);
phi23=params(88);
phi24=params(89);
phi30=params(90);
phi31=params(91);
phi32=params(92);
phi33=params(93);
phi34=params(94);

%     eta01 = params(41); % husgand's wage return to low ability type (mean)
%     eta02 = params(42); % husgand's wage return to high ability type (mean)
%     eta11 = params(43); % husgand's wage return to 2yr college (mean)
%     eta12 = params(44); % husgand's wage return to 4yr college (mean)
%     eta2 = params(45); % husgand's wage return to women's age (mean)
%     eta03 = params(46); % husgand's wage return to low ability type (sd)
%     eta04 = params(47); % husgand's wage return to high ability type (sd)
%     eta21 = params(48); % husgand's wage return to 2yr college (sd)
%     eta22 = params(49); % husgand's wage return to 4yr college (sd)
%     eta3 = params(50); % husgand's wage return to women's age (sd)
%      
%     iota01 = params(51); % child investment of low ability type (mean)
%     iota02 = params(52); % child investment of high ability type (mean)
%     iota11 = params(53); % child investment of 2yr college (mean)
%     iota12 = params(54); % child investment of 4yr college (mean)
%     iota2 = params(55); % child investment by women's age (mean)
%     iota03 = params(56); % child investment of low ability type (sd)
%     iota04 = params(57); % child investment of high ability type (sd)
%     iota21 = params(58); % child investment of 2yr college (sd)
%     iota22 = params(59); % child investment of 4yr college (sd)
%     iota3 = params(60); % child investment by women's age (sd)

% eta01 = params(31); % husgand's wage return to low ability type (mean)
% eta02 = params(32); % husgand's wage return to high ability type (mean)
% eta11 = params(33); % husgand's wage return to 2yr college (mean)
% eta12 = params(34); % husgand's wage return to 4yr college (mean)
% eta2 = params(35); % husgand's wage return to women's age (mean)
% eta03 = params(36); % husgand's wage return to low ability type (sd)
% eta04 = params(37); % husgand's wage return to high ability type (sd)
% eta21 = params(38); % husgand's wage return to 2yr college (sd)
% eta22 = params(39); % husgand's wage return to 4yr college (sd)
% eta3 = params(40); % husgand's wage return to women's age (sd)
%      
% iota01 = params(41); % child investment of low ability type (mean)
% iota02 = params(42); % child investment of high ability type (mean)
% iota11 = params(43); % child investment of 2yr college (mean)
% iota12 = params(44); % child investment of 4yr college (mean)
% iota2 = params(45); % child investment by women's age (mean)
% iota03 = params(46); % child investment of low ability type (sd)
% iota04 = params(47); % child investment of high ability type (sd)
% iota21 = params(48); % child investment of 2yr college (sd)
% iota22 = params(49); % child investment of 4yr college (sd)
% iota3 = params(50); % child investment by women's age (sd)

eta01=params(45);
eta02=params(46);
eta11=params(47);
eta12=params(48);
eta2=params(49);
eta03=params(50);
eta04=params(51);
eta21=params(52);
eta22=params(53);
eta3=params(54);
iota01=params(55);
iota02=params(56);
iota11=params(57);
iota12=params(58);
iota2=params(59);
iota03=params(60);
iota04=params(61);
iota21=params(62);
iota22=params(63);
iota3=params(64);

% rho01	=	8.275;
% rho02	=	0.353;
% rho11	=	0.389;
% rho12	=	0.734;
% rho03	=	0.254;
% rho04	=	-0.135;
% rho21	=	0.0044;
% rho22	=	-0.00975;

% % option 1: rescaled asset (1/100) - DOESN'T WORK (ASSETS IS OUT OF BOUNDS)
% rho01 = 35.15; % assets at age 18-20 (mean)
% rho02 = 48.29; % assets at age 18-20 (mean)
% rho11 = 27.81; % assets at age 18-20 (mean)
% rho12 = 96.12; % assets at age 18-20 (mean)
% rho03 = 15.94; % assets at age 18-20 (sd)
% rho04 = -8.921; % assets at age 18-20 (sd)
% rho21 = 0.577; % assets at age 18-20 (sd)
% rho22 = -0.549; % assets at age 18-20 (sd)
 
% % option 2: log asset - LOOKS A BIT BETTER
% rho01 = 8.316; % assets at age 18-20 (mean)
% rho02 = 0.327; % assets at age 18-20 (mean)
% rho11 = 0.437; % assets at age 18-20 (mean)
% rho12 = 0.702; % assets at age 18-20 (mean)
% rho03 = 0.141; % assets at age 18-20 (sd)
% rho04 = -0.0771; % assets at age 18-20 (sd)
% rho21 = 0.00122; % assets at age 18-20 (sd)
% rho22 = -0.00944; % assets at age 18-20 (sd)

    % option 3: real asset 
    rho01 = 5112; % assets at age 18-20 (mean)
    rho02 = 2464; % assets at age 18-20 (mean)
    rho11 = 2133; % assets at age 18-20 (mean)
    rho12 = 8671; % assets at age 18-20 (mean)
    rho03 = 2628; % assets at age 18-20 (sd)
    rho04 = -1479; % assets at age 18-20 (sd)
    rho21 = 158.9; % assets at age 18-20 (sd)
    rho22 = -16.16; % assets at age 18-20 (sd)

%% Initial Conditions

exp_s = zeros(G.n_pop,G.n_period); 
m_s = zeros(G.n_pop,G.n_period); 
ch_s = zeros(G.n_pop,G.n_period); 
a_s = zeros(G.n_pop,G.n_period); 

%Run regressions of mean and sd of assets at age 18-20 on education/ability
a_mean = rho01 + rho02*(abi==2) + rho11*(edu==2) + rho12*(edu==3);
a_sd = rho03 + rho04*(abi==2) + rho21*(edu==2) + rho22*(edu==3);
a_s(:,1) = 0; %normrnd(a_mean, a_sd); % draw from distribution

%Check initial assets are within the bounds
sum(a_s(:,1)<min(S.SS_A))
sum(a_s(:,1)>max(S.SS_A))
for n=1:G.n_pop
    if a_s(n,1)<min(S.SS_A)
        a_s(n,1)= min(S.SS_A);
    end
    if a_s(:,1)>max(S.SS_A)
        a_s(:,1)= max(S.SS_A);
    end
end
sum(a_s(:,1)<min(S.SS_A))
sum(a_s(:,1)>max(S.SS_A))

%Wider Assets Vector
A_wide = S.SS_A;
A_wide(1) = min(S.SS_A)-1000;
A_wide(G.n_assets)= max(S.SS_A)+5000;

%% Loop over all periods/individuals
for n=1:G.n_pop
for t=1:1:G.n_period-1
    %for n=1:n%1:G.n_pop
    
    % Obtain sector shocks
    epssim_r(n,t)=sqrt(2)*G.Eps(1,n,t)'*sigma_r;
    epssim_n(n,t)=sqrt(2)*G.Eps(2,n,t)'*sigma_n;
    
    if epssim_r(n,t)<S.eps_r(1) || epssim_r(n,t)>S.eps_r(3)
        eps_rg=epssim_r(n,t)*[-1;0;1];
    else 
        eps_rg=S.eps_r;
    end
    if epssim_n(n,t)<S.eps_n(1) || epssim_n(n,t)>S.eps_n(3)
        eps_ng=epssim_n(n,t)*[-1;0;1];
    else 
        eps_ng=S.eps_n;
    end
    
    % Calculate wages
    wr_s(n,t) = exp(alpha01_r + alpha02_r*(abi(n)==2) + alpha11_r*(edu(n)==2) + alpha12_r*(edu(n)==3) + alpha2_r*log(1+exp_s(n,t)) + epssim_r(n,t)); 
    wn_s(n,t) = exp(alpha01_n + alpha02_n*(abi(n)==2) + alpha11_n*(edu(n)==2) + alpha12_n*(edu(n)==3) + alpha2_n*log(1+exp_s(n,t)) + epssim_n(n,t));

    % Locate in the experience/marriage/children vector   
    if m_s(n,t)==0 
        x=min(30, 20 + exp_s(n,t)+1);  
    elseif m_s(n,t)==1
        if ch_s(n,t)==1
           x=min(10, exp_s(n,t)+1); 
        elseif ch_s(n,t)==2
           x=min(20, 10 + exp_s(n,t)+1);  
        end
    end
                            
    % Optimal Choices (using linear interpolation)
    cc_s(n,t)= interpn(A_wide, eps_rg, eps_ng, C(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t));
    rr_s(n,t)= interpn(A_wide, eps_rg, eps_ng, R(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t)); 
    nn_s(n,t)= interpn(A_wide, eps_rg, eps_ng, N(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t)); 
    uu_s(n,t)= interpn(A_wide, eps_rg, eps_ng, U(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t)); 
       
    [v, Ind] = max([rr_s(n,t), nn_s(n,t), uu_s(n,t)]);
    if Ind==1  
     r_s(n,t)=1;
     n_s(n,t)=0;
     u_s(n,t)=0;
     exp_s(n,t+1)=exp_s(n,t)+1;
    elseif Ind==2
     r_s(n,t)=0;
     n_s(n,t)=1;
     u_s(n,t)=0;
     exp_s(n,t+1)=exp_s(n,t)+1;
    else
     r_s(n,t)=0;
     n_s(n,t)=0;
     u_s(n,t)=1;
     exp_s(n,t+1)=exp_s(n,t);
    end
    
    % Probability of Children
    prob_2kids_r = normcdf(phi10 + phi11*(edu(n)==2) + phi12*(edu(n)==3) + phi13*exp_s(n,t));
    prob_2kids_n = normcdf(phi20 + phi21*(edu(n)==2) + phi22*(edu(n)==3) + phi23*exp_s(n,t));
    prob_2kids_u = normcdf(phi30 + phi31*(edu(n)==2) + phi32*(edu(n)==3) + phi33*exp_s(n,t));
       
    if m_s(n,t)==0 
       marr(n,t)=interpn(A_wide, eps_rg, eps_ng,M(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t));
       if marr(n,t)<0.5 
          m_s(n,t+1)=0;
       else
          m_s(n,t+1)=1;
          ch_s(n,t+1)=1;
       end
    else
        m_s(n,t+1)=m_s(n,t);
        if ch_s(n,t)==1 && r_s(n,t)==1
           if prob_2kids_r<0.5
               ch_s(n,t+1)=1;
           else
               ch_s(n,t+1)=2;
           end
        elseif ch_s(n,t)==1 && n_s(n,t)==1
           if prob_2kids_n<0.5
               ch_s(n,t+1)=1;
           else
               ch_s(n,t+1)=2;
           end
        elseif ch_s(n,t)==1 && u_s(n,t)==1
           if prob_2kids_u<0.5
               ch_s(n,t+1)=1;
           else
               ch_s(n,t+1)=2;
           end
        elseif ch_s(n,t)==2
           ch_s(n,t+1)=2;       
        end
    end
         
    % Find age of woman
    age = 18*(edu(n)==1) + 20*(edu(n)==2) + 22*(edu(n)==3) + t;

    % Draw husband wages
	wh_mean = eta01 + eta02*(abi(n)==2) + eta11*(edu(n)==2) + eta12*(edu(n)==3) + eta2*age;
	wh_sd = 0.7022257; %eta03 + eta04*(abi(n)==2) + eta21*(edu(n)==2) + eta22*(edu(n)==3) + eta3*age;
	wh_s(n,t) = normrnd(wh_mean,wh_sd);

    % Draw child investments
	Inv_mean = iota01 + iota02*(abi(n)==2) + iota11*(edu(n)==2) + iota12*(edu(n)==3) + iota2*age;
	Inv_sd = 0.9270494; %iota03 + iota04*(abi(n)==2) + iota21*(edu(n)==2) + iota22*(edu(n)==3) + iota3*age;
	inv_s(n,t) = normrnd(Inv_mean,Inv_sd);

    % Optimal simulated consumption with validations
    a_tmw(n,t)= (1+G.r)*(a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*exp(wh_s(n,t)) - ch_s(n,t)*exp(inv_s(n,t)) - cc_s(n,t));
    
    if a_tmw(n,t)<S.SS_A(1) % or make into -10
       c_s(n,t)=max(0,a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*exp(wh_s(n,t)) - ch_s(n,t)*exp(inv_s(n,t))- S.SS_A(1)/(1+G.r));  
      
    elseif a_tmw(n,t)>S.SS_A(G.n_assets)
       c_s(n,t)=max(0,a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*exp(wh_s(n,t)) - ch_s(n,t)*exp(inv_s(n,t)) - S.SS_A(G.n_assets)/(1+G.r));    
        
    else       
       c_s(n,t)=cc_s(n,t);
    end

    % Transition for assets (Budget Constraint)    
    a_s(n,t+1)= (1+G.r)*(a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*exp(wh_s(n,t)) - ch_s(n,t)*exp(inv_s(n,t)) - c_s(n,t));
   
    n    
    %end
    t
end
% test = [c_s(n,1:19)',r_s(n,1:19)',n_s(n,1:19)',u_s(n,1:19)',m_s(n,1:19)',ch_s(n,1:19)',inv_s(n,1:19)',a_s(n,1:19)',exp(wh_s(n,1:19)'),wr_s(n,1:19)',wn_s(n,1:19)',exp_s(n,1:19)'];
% sheetname = strcat('Sheet',num2str(n));
% xlswrite('simulation_July24_newparams.xls',test,sheetname)
end
end