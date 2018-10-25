function [c_func,m_func,ch_func,lr_func,ln_func,lu_func,Ar_out,An_out,Au_out,wh,W,W2] = solution(G,abi,edu,S,params)

%% index for parameters:

% estimated
psi_r=params(1);
psi_n=params(2);
psi_u=params(3);
theta1_r=params(4);
theta1_n=params(5);
theta1_u=params(6);
theta2_r=params(7);
theta2_n=params(8);
theta2_u=params(9);
theta3_r=params(10);
theta3_n=params(11);
theta3_u=params(12);
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
lambda1=params(26);
lambda2=params(27);
lambda3=params(28);
lambda4=params(29);

% calibrated    
omega0_r=params(30);
omega11=params(31);
omega12=params(32);
omega13=params(33);
omega14=params(34);
omega0_n=params(35);
omega21=params(36);
omega22=params(37);
omega23=params(38);
omega24=params(39);
omega0_u=params(40);
omega31=params(41);
omega32=params(42);
omega33=params(43);
omega34=params(44);
eta01=params(45);
eta02=params(46);
eta11=params(47);
eta12=params(48);
eta2=params(49);
% eta03=params(50);
% eta04=params(51);
% eta21=params(52);
% eta22=params(53);
% eta3=params(54);
% iota01=params(55);
% iota02=params(56);
% iota11=params(57);
% iota12=params(58);
% iota2=params(59);
% iota03=params(60);
% iota04=params(61);
% iota21=params(62);
% iota22=params(63);
% iota3=params(64);
kappa01=params(65);
kappa02=params(66);
kappa03=params(67);
kappa04=params(68);
kappa05=params(69);
tau10=params(70);
tau11=params(71);
tau12=params(72);
tau13=params(73);
tau14=params(74);
tau20=params(75);
tau21=params(76);
tau22=params(77);
tau23=params(78);
tau24=params(79);
% phi10=params(80);
% phi11=params(81);
% phi12=params(82);
% phi13=params(83);
% phi14=params(84);
% phi20=params(85);
% phi21=params(86);
% phi22=params(87);
% phi23=params(88);
% phi24=params(89);
% phi30=params(90);
% phi31=params(91);
% phi32=params(92);
% phi33=params(93);
% phi34=params(94);
% chi01=params(95);
% chi02=params(96);
% chi11=params(97);
% chi12=params(98);
% chi03=params(99);
% chi04=params(100);
% chi21=params(101);
% chi22=params(102);

%% new child costs & probabilities parameters:

iota01=-279.3;
iota02=14.50;
iota03=-1.223;
iota04=22.55;
iota05=10.76;

iota11=-601.1;
iota12=22.12;
iota13=61.02;
iota14=115.4;
iota15=21.24;

iota21=-706.3;
iota22=26.82;
iota23=226.2;
iota24=273.1;
iota25=24.74;

iota31=-1388;
iota32=-70.53;
iota33=155.4;
iota34=156.3;
iota35=49.59;

phi01=6.197;
phi02=0.201;
phi03=0.306;
phi04=0.462;
phi05=-0.210;

phi11=4.531;
phi12=-0.0317;
phi13=0.0987;
phi14=0.0931;
phi15=-0.132;

phi21=-2.053;
phi22=-0.0807;
phi23=-0.0415;
phi24=-0.0439;
phi25=0.0413;

phi31=-15.22;
phi32=-0.149;
phi33=-0.611;
phi34=-0.748;
phi35=0.350;

omi01=3.121;
omi02=0.187;
omi03=0.261;
omi04=0.335;
omi05=-0.131;

omi11=-0.00928;
omi12=0.123;
omi13=0.312;
omi14=0.131;
omi15=-0.0484;

omi21=-3.851;
omi22=0.173;
omi23=0.165;
omi24=0.0189;
omi25=0.0478;

omi31=-8.215;
omi32=-0.128;
omi33=0.188;
omi34=-0.239;
omi35=0.145;

omi41=-1.357;
omi42=-0.188;
omi43=-0.251;
omi44=-0.102;
omi45=-0.00468;

omi51=4.303;
omi52=-0.252;
omi53=-0.134;
omi54=-0.252;
omi55=-0.120;

omi61=-2.768;
omi62=-0.181;
omi63=-0.216;
omi64=-0.287;
omi65=0.0628;

omi71=-12.88;
omi72=-0.0946;
omi73=-0.602;
omi74=-0.655;
omi75=0.291;

%% other parameters:

% consumption
chh_min = 50; % minimun consumption
delta = 0.5; % Female Share of Consumption (CAL)

% expanded assets vector for linear interpolation
A_wide = S.SS_A;
A_wide(1) = S.extmin_A;
A_wide(G.n_assets) = S.extmax_A;

% taxes
basetax = 0.1;
marrtax = 0.05;
kidtax = 0.05;
w_min = 100; % minimum wage for child tax credits

%% Terminal Value Function:

% age at TVF
age_TVF = 18*(edu==1) + 20*(edu==2) + 22*(edu==3) + G.n_period;

% assets
assets = S.SS_A; 

% husband wages
% wh_mean = eta01 + eta02*(abi==2) + eta11*(edu==2) + eta12*(edu==3) + eta2*age_TVF;
% wh_sd = 0; %0.7022257; %eta03 + eta04*(abi==2) + eta21*(edu==2) + eta22*(edu==3) + eta3*age_TVF;
% wh_TVF = normrnd(wh_mean,wh_sd);

% child human capital 
% Inv_mean = iota01 + iota02*(abi==2) + iota11*(edu==2) + iota12*(edu==3) + iota2*age_TVF;
% Inv_sd = 0; %0.9270494; %iota03 + iota04*(abi==2) + iota21*(edu==2) + iota22*(edu==3) + iota3*age_TVF;
% Inv_TVF = normrnd(Inv_mean,Inv_sd);
% K_TVF = exp(kappa01 + kappa02*(abi==2) + kappa03*(edu==2) + kappa04*(edu==3) + kappa05*(Inv_TVF));

% child investments (by type and age)
Inv_mean_05 = iota01 + iota02*(abi==2) + iota03*(edu==2) + iota04*(edu==3) + iota05*age_TVF;
Inv_sd_05 = 233.7977;
Inv_mean_611 = iota11 + iota12*(abi==2) + iota13*(edu==2) + iota14*(edu==3) + iota15*age_TVF;
Inv_sd_611 = 320.3355;
Inv_mean_1217 = iota21 + iota22*(abi==2) + iota23*(edu==2) + iota24*(edu==3) + iota25*age_TVF;
Inv_sd_1217 = 375.8709;
Inv_mean_18 = iota31 + iota32*(abi==2) + iota33*(edu==2) + iota34*(edu==3) + iota35*age_TVF;
Inv_sd_18 = 508.7068;

% draw investments
Inv05 = normrnd(Inv_mean_05,Inv_sd_05); %investment years 0 to 5
Inv611 = normrnd(Inv_mean_611,Inv_sd_611); %investment years 6 to 11
Inv1217 = normrnd(Inv_mean_1217,Inv_sd_1217); %investment years 12 to 17
Inv18 = normrnd(Inv_mean_18,Inv_sd_18); %investment 18+

% 1 child prob:
prob_1k_05 = normcdf(omi01 + omi02*(abi==2) + omi03*(edu==2) + omi04*(edu==3) + omi05*age_TVF);
prob_1k_611 = normcdf(omi11 + omi12*(abi==2) + omi13*(edu==2) + omi14*(edu==3) + omi15*age_TVF);
prob_1k_1217 = normcdf(omi21 + omi22*(abi==2) + omi23*(edu==2) + omi24*(edu==3) + omi25*age_TVF);
prob_1k_18 = normcdf(omi31 + omi32*(abi==2) + omi33*(edu==2) + omi34*(edu==3) + omi35*age_TVF);

% 2 child prob:
prob_2k_05 = normcdf(omi41 + omi42*(abi==2) + omi43*(edu==2) + omi44*(edu==3) + omi45*age_TVF);
prob_2k_611 = normcdf(omi51 + omi52*(abi==2) + omi53*(edu==2) + omi54*(edu==3) + omi55*age_TVF);
prob_2k_1217 = normcdf(omi61 + omi62*(abi==2) + omi63*(edu==2) + omi64*(edu==3) + omi65*age_TVF);
prob_2k_18 = normcdf(omi71 + omi72*(abi==2) + omi73*(edu==2) + omi74*(edu==3) + omi75*age_TVF);

% Expected Investment => turn into log(Inv) vector for the TVF!
Inv_1k_TVF = prob_1k_05*(0.5*Inv05+Inv611+Inv1217+Inv1217) + ...
    + prob_1k_611*(0.5*Inv611+Inv1217+Inv18)...
    + prob_1k_1217*(0.5*Inv1217+Inv18) + prob_1k_18*(0.5*Inv18); % 1 child
Inv_2k_TVF = prob_2k_05*(0.5*Inv05+Inv611+Inv1217+Inv1217) + ...
    + prob_2k_611*(0.5*Inv611+Inv1217+Inv18)...
    + prob_2k_1217*(0.5*Inv1217+Inv18) + prob_2k_18*(0.5*Inv18); % 2 child
Inv_TVF = [repmat(Inv_1k_TVF,1,10),repmat(2*Inv_2k_TVF,1,10),zeros(1,10)];

% TVF 
TVF = repmat(real(lambda1*(assets).^(1-G.sigma)/(1-G.sigma))',1,30) ...
    + repmat(lambda2*(S.SS_X.^(1-G.sigma))/(1-G.sigma),G.n_assets,1) ...
    + lambda3*(S.SS_M.^(1-G.sigma))/(1-G.sigma) ...
    + lambda4*(Inv_TVF.^(1-G.sigma))/(1-G.sigma);

tic
% loop for time (20):
for t = G.n_period-1:-1:1
    t
    toc
    
    % Coefficients for Chebyshev Approximation
    if t==G.n_period-1
        Emax = TVF;
        for x = 1:1:(G.n_matstat*G.n_wrkexp)
            Num(x,:) = Emax(:,x)'*S.T_A;
            Den = S.T2_A;
            coeff(x,:) = Num(x,:)./Den';
            coeff2(x,:) = coeff(x,:);
        end
    else
        Emax = W(:,:,t+1);
        Emax2 = W2(:,:,t+1);
        for x = 1:1:(G.n_matstat*G.n_wrkexp)
            Num(x,:) = Emax(:,x)'*S.T_A;
            Num2(x,:) = Emax2(:,x)'*S.T_A;
            Den = S.T2_A;
            coeff(x,:) = Num(x,:)./Den';
            coeff2(x,:) = Num2(x,:)./Den';        
        end
    end
    
    % age
    age = 18*(edu==1) + 20*(edu==2) + 22*(edu==3) + t;
    
    % draw husband wage
    wh_mean = eta01 + eta02*(abi==2) + eta11*(edu==2) + eta12*(edu==3) + eta2*age;
    wh_sd = 0.7022257; %eta03 + eta04*(abi==2) + eta21*(edu==2) + eta22*(edu==3) + eta3*age;
    wh(t) = normrnd(wh_mean,wh_sd);
    
    % draw investments
%     Inv_mean = iota01 + iota02*(abi==2) + iota11*(edu==2) + iota12*(edu==3) + iota2*age;
%     Inv_sd = 0.9270494; %iota03 + iota04*(abi==2) + iota21*(edu==2) + iota22*(edu==3) + iota3*age;
%     Inv(t) = normrnd(Inv_mean,Inv_sd);
    
    % child investments (by type and age)
    Inv_mean_05 = iota01 + iota02*(abi==2) + iota03*(edu==2) + iota04*(edu==3) + iota05*age_TVF;
    Inv_sd_05 = 233.7977;
    Inv_mean_611 = iota11 + iota12*(abi==2) + iota13*(edu==2) + iota14*(edu==3) + iota15*age_TVF;
    Inv_sd_611 = 320.3355;
    Inv_mean_1217 = iota21 + iota22*(abi==2) + iota23*(edu==2) + iota24*(edu==3) + iota25*age_TVF;
    Inv_sd_1217 = 375.8709;
    Inv_mean_18 = iota31 + iota32*(abi==2) + iota33*(edu==2) + iota34*(edu==3) + iota35*age_TVF;
    Inv_sd_18 = 508.7068;

    % draw investments
    Inv05 = normrnd(Inv_mean_05,Inv_sd_05); %investment years 0 to 5
    Inv611 = normrnd(Inv_mean_611,Inv_sd_611); %investment years 6 to 11
    Inv1217 = normrnd(Inv_mean_1217,Inv_sd_1217); %investment years 12 to 17
    Inv18 = normrnd(Inv_mean_18,Inv_sd_18); %investment 18+

    % Probabilty of child ages - just like other probabiliy, by type & age
    prob05 = normcdf(phi01 + phi02*(abi==2) + phi03*(edu==2) + phi04*(edu==3) + phi05*age);
    prob611 = normcdf(phi11 + phi12*(abi==2) + phi13*(edu==2) + phi14*(edu==3) + phi15*age);
    prob1217 = normcdf(phi21 + phi22*(abi==2) + phi23*(edu==2) + phi24*(edu==3) + phi25*age);
    prob18 = normcdf(phi31 + phi32*(abi==2) + phi33*(edu==2) + phi34*(edu==3) + phi35*age);

    % expected investment
    Inv(t) = prob05*Inv05 + prob611*Inv611 + prob1217*Inv1217 + prob18*Inv18;
    
    % child human capital 
    K(t) = exp(kappa01 + kappa02*(abi==2) + kappa03*(edu==2) + kappa04*(edu==3) + kappa05*(Inv(t)));
    
    % loop for work experience and marital status (30):
    for x = 1:1:(G.n_matstat*G.n_wrkexp)
        x;
        
        % current state discrete variables:
        m_j = S.SS_M(x);  % marital status
        n_j = S.SS_N(x);  % children
        X_j = S.SS_X(x);  % experience
        
        % job probabilities
        prob_lamba = normcdf(tau10 + tau11*(edu==2) + tau12*(edu==3) + tau13*age + tau14*X_j); % losing a regular job
        prob_pi = normcdf(tau20 + tau21*(edu==2) + tau22*(edu==3) + tau23*age + tau24*X_j); % getting a regular job

        % 2nd child probabilities
%         prob_2kids_r = normcdf(phi10 + phi11*(edu==2) + phi12*(edu==3) + phi13*X_j);
%         prob_2kids_n = normcdf(phi20 + phi21*(edu==2) + phi22*(edu==3) + phi23*X_j);
%         prob_2kids_u = normcdf(phi30 + phi31*(edu==2) + phi32*(edu==3) + phi33*X_j);

        % loop for shocks (27):
        for i = 1:1:G.n_shocks
            i;
            
            % shocks
            shock_i = S.shocks_i(i);
            shock_r = S.shocks_r(i);
            shock_n = S.shocks_n(i);
            
            % sector-specific wages
            PT_w_j_r = exp(alpha01_r + alpha02_r*(abi==2) + alpha11_r*(edu==2) + alpha12_r*(edu==3) + alpha2_r*log(1+X_j) + shock_r);
            PT_w_j_n = exp(alpha01_n + alpha02_n*(abi==2) + alpha11_n*(edu==2) + alpha12_n*(edu==3) + alpha2_n*log(1+X_j) + shock_n);
            w_j_u = 0; % unemployed women don't have earnings
            
            % post-tax wages
            w_j_r = (1 - basetax + m_j*marrtax)*PT_w_j_r + n_j*kidtax*PT_w_j_r*(PT_w_j_r<= w_min);
            w_j_n = (1 - basetax + m_j*marrtax)*PT_w_j_n + n_j*kidtax*PT_w_j_n*(PT_w_j_n<= w_min);
            
            % loop over assets (10):
            for j = 1:1:G.n_assets
                j;
                
                % HH's assets
                A_j = S.SS_A(j); 
                
                % marriage probabilities
                prob_marr_r = normcdf(omega0_r + omega11*(edu==2) + omega12*(edu==3) + omega13*log(1+age) + omega14*log(10+A_j));
                prob_marr_n = normcdf(omega0_n + omega21*(edu==2) + omega22*(edu==3) + omega23*log(1+age) + omega24*log(10+A_j));
                prob_marr_u = normcdf(omega0_u + omega31*(edu==2) + omega32*(edu==3) + omega33*log(1+age) + omega34*log(10+A_j));

                % consumption vector
                chh_max = A_j + max(w_j_r,w_j_n) + exp(wh(t));
                c_vector = linspace(chh_min,chh_max,G.n_cons);
                cr_vector = c_vector;
                cn_vector = c_vector;
                cu_vector = c_vector;

                % loop over consumption (30):
                for k = 1:1:G.n_cons
                    k;
                    
                    % HH consumption
                    chh_r = cr_vector(k);
                    chh_n = cn_vector(k);
                    chh_u = cu_vector(k);
                    
                    % validate child investments
                     if n_j*exp(Inv(t)) > A_j + w_j_r + exp(wh(t))*m_j + shock_i - chh_r 
                        inv_r = A_j + w_j_r + exp(wh(t))*m_j + shock_i - chh_r;
                     else
                        inv_r = n_j*exp(Inv(t));
                     end
                     % validate child investments
                     if n_j*exp(Inv(t)) > A_j + w_j_n + exp(wh(t))*m_j + shock_i - chh_n 
                        inv_n = A_j + w_j_n + exp(wh(t))*m_j + shock_i - chh_n;
                     else
                        inv_n = n_j*exp(Inv(t));
                     end
                     % validate child investments
                     if n_j*exp(Inv(t)) > A_j + w_j_u + exp(wh(t))*m_j + shock_i - chh_u 
                        inv_u = A_j + w_j_u + exp(wh(t))*m_j + shock_i - chh_u;
                     else
                        inv_u = n_j*exp(Inv(t));
                     end
                    
                    % woman's consumption
                    cw_r = delta*chh_r;
                    cw_n = delta*chh_n;
                    cw_u = delta*chh_u;
                    
                    % Sector-Specific Utility
                    u_r(k) = (cw_r^(1-G.sigma))/(1-G.sigma) + psi_r + theta1_r*log(1+m_j) + theta2_r*log(1+n_j) + theta3_r*log(K(t))*n_j;
                    u_n(k) = (cw_n^(1-G.sigma))/(1-G.sigma) + psi_n + theta1_n*log(1+m_j) + theta2_n*log(1+n_j) + theta3_n*log(K(t))*n_j;
                    u_u(k) = (cw_u^(1-G.sigma))/(1-G.sigma) + psi_u + theta1_u*log(1+m_j) + theta2_u*log(1+n_j) + theta3_u*log(K(t))*n_j;

                    % Maried with 1 kid:
                    if x <= 10
                        
                        % regular job
                        A_next = (1+G.r) * (A_j + (w_j_r + exp(wh(t))*m_j + shock_i) - chh_r - inv_r); %n_j*exp(Inv(t))); % eq. 8
                        A_next = S.extmax_A*(A_next>S.extmax_A) + A_next*(A_next<=S.extmax_A);
                        x_next = x + 1;
                        if x_next == 11
                            x_next = 10;
                        end

                        % polynomial approximation of VF
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vm_r_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vm_r_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        Vm_r_next_2 = sum(coeff(x_next+10,:).*Base,2); %cheby_approx
                        Vm_r_next2_2 = sum(coeff2(x_next+10,:).*Base,2); %cheby_approx
                        % save output
                        Amr_next(k)=A_next;
                        if A_next < S.extmin_A
                            Vmr_next(k)=NaN;
                            Vmr_next2(k)=NaN;
                            Vmr_next_2(k)=NaN;
                            Vmr_next2_2(k)=NaN;
                        else
                            Vmr_next(k)=Vm_r_next;
                            Vmr_next2(k)=Vm_r_next2;
                            Vmr_next_2(k)=Vm_r_next_2;
                            Vmr_next2_2(k)=Vm_r_next2_2;
                        end
                                            
                        % non-regular job
                        A_next = (1+G.r) * (A_j + (w_j_n + exp(wh(t))*m_j + shock_i) - chh_n - inv_n); %n_j*exp(Inv(t)));
                        A_next = S.extmax_A*(A_next>S.extmax_A) + A_next*(A_next<=S.extmax_A);
                        x_next = x + 1;
                        if x_next == 11
                            x_next = 10;
                        end

                        % polynomial approximation of VF
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vm_n_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vm_n_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        Vm_n_next_2 = sum(coeff(x_next+10,:).*Base,2); %cheby_approx
                        Vm_n_next2_2 = sum(coeff2(x_next+10,:).*Base,2); %cheby_approx
                        % save outpu
                        Amn_next(k)=A_next;
                        if A_next < S.extmin_A
                            Vmn_next(k)=NaN;
                            Vmn_next2(k)=NaN;
                            Vmn_next_2(k)=NaN;
                            Vmn_next2_2(k)=NaN;
                        else
                            Vmn_next(k)=Vm_n_next;
                            Vmn_next2(k)=Vm_n_next2;
                            Vmn_next_2(k)=Vm_n_next_2;
                            Vmn_next2_2(k)=Vm_n_next2_2;
                        end
                        
                        % unemployed
                        A_next = (1+G.r) * (A_j + (w_j_u + exp(wh(t))*m_j + shock_i) - chh_u - inv_u); %n_j*exp(Inv(t)));
                        A_next = S.extmax_A*(A_next>S.extmax_A) + A_next*(A_next<=S.extmax_A);
                        x_next = x;

                        % polynomial approximation of VF                      
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vm_u_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vm_u_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        Vm_u_next_2 = sum(coeff(x_next+10,:).*Base,2); %cheby_approx
                        Vm_u_next2_2 = sum(coeff2(x_next+10,:).*Base,2); %cheby_approx
                        % save output
                        Amu_next(k)=A_next;
                        if A_next < S.extmin_A
                            Vmu_next(k)=NaN;
                            Vmu_next2(k)=NaN;
                            Vmu_next_2(k)=NaN;
                            Vmu_next2_2(k)=NaN;
                        else
                            Vmu_next(k)=Vm_u_next;
                            Vmu_next2(k)=Vm_u_next2;
                        	Vmu_next_2(k)=Vm_u_next_2;
                            Vmu_next2_2(k)=Vm_u_next2_2;
                        end
                        
                        % Sector-Specific Value Functions (1 child) 
                          Vm1_r(k) = u_r(k) + G.beta * ((prob_lamba*Vmr_next(k))+(1-prob_lamba)*Vmr_next2(k));
                          Vm1_n(k) = u_n(k) + G.beta * ((prob_pi*Vmn_next(k))+(1-prob_pi)*Vmn_next2(k));
                          Vm1_u(k) = u_u(k) + G.beta * ((prob_pi*Vmu_next(k))+(1-prob_pi)*Vmu_next2(k));
                        
                        % Sector-Specific Value Functions (2 child, maybe)
                          Vm1_r_2(k) = u_r(k) + G.beta * ((prob_lamba*Vmr_next_2(k))+(1-prob_lamba)*Vmr_next2_2(k));
                          Vm1_n_2(k) = u_n(k) + G.beta * ((prob_pi*Vmn_next_2(k))+(1-prob_pi)*Vmn_next2_2(k));
                          Vm1_u_2(k) = u_u(k) + G.beta * ((prob_pi*Vmu_next_2(k))+(1-prob_pi)*Vmu_next2_2(k));

                        % Sector-Specific Value Fucntions (married)
%                           Vm_r(k) = Vm1_r(k)*(1-prob_2kids_r) + Vm1_r_2(k)*(prob_2kids_r);
%                           Vm_n(k) = Vm1_n(k)*(1-prob_2kids_n) + Vm1_n_2(k)*(prob_2kids_n);
%                           Vm_u(k) = Vm1_u(k)*(1-prob_2kids_u) + Vm1_u_2(k)*(prob_2kids_u);
                        
                        % save marriage 1 child (for marriage decision)
                        Vm_r_aux(k,x,j,i) = Vm1_r(k); %ERROR: before no j,i
                        Vm_n_aux(k,x,j,i) = Vm1_n(k); %ERROR: before no j,i
                        Vm_u_aux(k,x,j,i) = Vm1_u(k); %ERROR: before no j,i
                        
                    % Married with 2 kids: 
                    elseif x <= 20
                        
                        % regular job
                        A_next = (1+G.r) * (A_j + (w_j_r + exp(wh(t))*m_j + shock_i) - chh_r - inv_r); %n_j*exp(Inv(t)));% eq. 8
                        A_next = S.extmax_A*(A_next>S.extmax_A) + A_next*(A_next<=S.extmax_A);
                        x_next = x + 1;
                        if x_next == 21
                            x_next = 20;
                        end

                        % Chebyshev Approximation
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vm_r_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vm_r_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        % save output
                        Am2r_next(k)=A_next;
                        if A_next < S.extmin_A
                            Vm2r_next(k)=NaN;
                            Vm2r_next2(k)=NaN;
                        else
                            Vm2r_next(k)=Vm_r_next;
                            Vm2r_next2(k)=Vm_r_next2;
                        end
                          
                        % non-regular job
                        A_next = (1+G.r) * (A_j + (w_j_n + exp(wh(t))*m_j + shock_i) - chh_n - inv_n); %n_j*exp(Inv(t)));
                        A_next = S.extmax_A*(A_next>S.extmax_A) + A_next*(A_next<=S.extmax_A);
                        x_next = x + 1;
                        if x_next == 21
                            x_next = 20;
                        end

                        % Chebyshev Approximation
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vm_n_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vm_n_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        % save output
                        Am2n_next(k)=A_next;
                        if A_next < S.extmin_A
                            Vm2n_next(k)=NaN;
                            Vm2n_next2(k)=NaN;
                        else
                            Vm2n_next(k)=Vm_n_next;
                            Vm2n_next2(k)=Vm_n_next2;
                        end

                        % unemployed
                        A_next = (1+G.r) * (A_j + (w_j_u + exp(wh(t))*m_j + shock_i) - chh_u - inv_u); %n_j*exp(Inv(t)));
                        A_next = S.extmax_A*(A_next>S.extmax_A) + A_next*(A_next<=S.extmax_A);
                        x_next = x;

                        % Chebyshev Approximation
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vm_u_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vm_u_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        % save output
                        Am2u_next(k)=A_next;
                        if A_next < S.extmin_A
                            Vm2u_next(k)=NaN;
                            Vm2u_next2(k)=NaN;
                        else
                            Vm2u_next(k)=Vm_u_next;
                            Vm2u_next2(k)=Vm_u_next2;
                        end

                        % Sector-Specific Value Functions (2 child)
                        Vm2_r(k) = u_r(k) + G.beta * ((prob_lamba*Vm2r_next(k))+(1-prob_lamba)*Vm2r_next2(k));
                        Vm2_n(k) = u_n(k) + G.beta * ((prob_pi*Vm2n_next(k))+(1-prob_pi)*Vm2n_next2(k));
                        Vm2_u(k) = u_u(k) + G.beta * ((prob_pi*Vm2u_next(k))+(1-prob_pi)*Vm2u_next2(k));
                        
                    % Single
                    else
                        
                        % Regular
                        A_next = (1+G.r) * (A_j + (w_j_r + exp(wh(t))*m_j + shock_i) - chh_r - inv_r); %n_j*exp(Inv(t)));
                        A_next = S.extmax_A*(A_next>S.extmax_A) + A_next*(A_next<=S.extmax_A);
                        x_next = x + 1;
                        if x_next == 31
                            x_next = 30;
                        end

                        % polynomial approximation of VF
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vs_r_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vs_r_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        % save output
                        Asr_next(k)=A_next;
                        if A_next < S.extmin_A
                            Vsr_next(k)=NaN;
                            Vsr_next2(k)=NaN;
                        else
                            Vsr_next(k)=Vs_r_next;
                            Vsr_next2(k)=Vs_r_next2;
                        end

                        % Non-regular
                        A_next = (1+G.r) * (A_j + (w_j_n + exp(wh(t))*m_j + shock_i) - chh_n - inv_n); %n_j*exp(Inv(t)));
                        A_next = S.extmax_A*(A_next>S.extmax_A) + A_next*(A_next<=S.extmax_A);
                        x_next = x + 1;
                        if x_next == 31
                            x_next = 30;
                        end

                        % polynomial approximation of VF
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vs_n_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vs_n_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        % save output
                        Asn_next(k)=A_next;
                        if A_next < S.extmin_A
                            Vsn_next(k)=NaN;
                            Vsn_next2(k)=NaN;
                        else
                            Vsn_next(k)=Vs_n_next;
                            Vsn_next2(k)=Vs_n_next2;
                        end

                        % Unemployed
                        A_next = (1+G.r) * (A_j + (w_j_u + exp(wh(t))*m_j + shock_i) - chh_u - inv_u); %n_j*exp(Inv(t)));
                        A_next = S.extmax_A*(A_next>S.extmax_A) + A_next*(A_next<=S.extmax_A);
                        x_next = x;

                        % linear approximation of VF
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vs_u_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vs_u_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        % save output
                        Asu_next(k)=A_next;
                        if A_next < S.extmin_A
                            Vsu_next(k)=NaN;
                            Vsu_next2(k)=NaN;
                        else
                            Vsu_next(k)=Vs_u_next;
                            Vsu_next2(k)=Vs_u_next2;
                        end

                        % Sector-Specific Value Functions (single)
                        Vs_r(k) = u_r(k) + G.beta * ((prob_lamba*Vsr_next(k))+(1-prob_lamba)*Vsr_next2(k));
                        Vs_n(k) = u_n(k) + G.beta * ((prob_pi*Vsn_next(k))+(1-prob_pi)*Vsn_next2(k));
                        Vs_u(k) = u_u(k) + G.beta * ((prob_pi*Vsu_next(k))+(1-prob_pi)*Vsu_next2(k));
                        Vsm_r(k) = prob_marr_r*Vm_r_aux(k,x-20,j,i) + (1-prob_marr_r)*Vs_r(k);
                        Vsm_n(k) = prob_marr_n*Vm_n_aux(k,x-20,j,i) + (1-prob_marr_n)*Vs_n(k);
                        Vsm_u(k) = prob_marr_u*Vm_u_aux(k,x-20,j,i) + (1-prob_marr_u)*Vs_u(k);
                    end
                end

                % optimal consumption and max VF
                if x <= 10
                    % check
                    Vm1_r(Amr_next < min(A_wide)) = NaN;
                    Vm1_r(Amr_next > max(A_wide)) = NaN;
                    Vm1_n(Amn_next < min(A_wide)) = NaN;
                    Vm1_n(Amn_next > max(A_wide)) = NaN;
                    Vm1_u(Amu_next < min(A_wide)) = NaN;
                    Vm1_u(Amu_next > max(A_wide)) = NaN;
                    Vm1_r_2(Amr_next < min(A_wide)) = NaN;
                    Vm1_r_2(Amr_next > max(A_wide)) = NaN;
                    Vm1_n_2(Amn_next < min(A_wide)) = NaN;
                    Vm1_n_2(Amn_next > max(A_wide)) = NaN;
                    Vm1_u_2(Amu_next < min(A_wide)) = NaN;
                    Vm1_u_2(Amu_next > max(A_wide)) = NaN;
                    % save optimal
                    [Vm1_r_star, index_m1r_k] = max(Vm1_r);
                    [Vm1_n_star, index_m1n_k] = max(Vm1_n);
                    [Vm1_u_star, index_m1u_k] = max(Vm1_u);
                    [Vm1_r_2_star, index_m1r_2_k] = max(Vm1_r_2);
                    [Vm1_n_2_star, index_m1n_2_k] = max(Vm1_n_2);
                    [Vm1_u_2_star, index_m1u_2_k] = max(Vm1_u_2);
                    cm1_r_star = cr_vector(index_m1r_k);
                    cm1_n_star = cn_vector(index_m1n_k);
                    cm1_u_star = cu_vector(index_m1u_k);
                    cm1_r_2_star = cr_vector(index_m1r_2_k);
                    cm1_n_2_star = cn_vector(index_m1n_2_k);
                    cm1_u_2_star = cu_vector(index_m1u_2_k);
                    cm_star_aux = [cm1_r_star, cm1_n_star, cm1_u_star, cm1_r_2_star, cm1_n_2_star, cm1_u_2_star];
                    [Vm_star, lm_index] = max([Vm1_r_star,Vm1_n_star,Vm1_u_star, Vm1_r_2_star, Vm1_n_2_star, Vm1_u_2_star]); % 3 job options
                    [Vm_star2] = max([Vm1_n_star,Vm1_u_star, Vm1_n_2_star, Vm1_u_2_star]); % emax2 (2 job options)
                elseif x <= 20
                    % check
                    Vm2_r(Am2r_next < min(A_wide)) = NaN;
                    Vm2_r(Am2r_next > max(A_wide)) = NaN;
                    Vm2_n(Am2n_next < min(A_wide)) = NaN;
                    Vm2_n(Am2n_next > max(A_wide)) = NaN;
                    Vm2_u(Am2u_next < min(A_wide)) = NaN;
                    Vm2_u(Am2u_next > max(A_wide)) = NaN;
                    % save optimal
                    [Vm2_r_star, index_m2r_k] = max(Vm2_r); %ERROR was Vm_r
                    [Vm2_n_star, index_m2n_k] = max(Vm2_n); %ERROR was Vm_n
                    [Vm2_u_star, index_m2u_k] = max(Vm2_u); %ERROR was Vm_u
                    cm2_r_star = cr_vector(index_m2r_k);
                    cm2_n_star = cn_vector(index_m2n_k);
                    cm2_u_star = cu_vector(index_m2u_k);
                    cm2_star_aux = [cm2_r_star, cm2_n_star, cm2_u_star];
                    [Vm2_star, lm2_index] = max([Vm2_r_star,Vm2_n_star,Vm2_u_star]); % 3 job options
                    [Vm2_star2] = max([Vm2_n_star,Vm2_u_star]); % emax2 (2 job options)
                else
                    % check
                    Vs_r(Asr_next < min(A_wide)) = NaN;
                    Vs_r(Asr_next > max(A_wide)) = NaN;
                    Vs_n(Asn_next < min(A_wide)) = NaN;
                    Vs_n(Asn_next > max(A_wide)) = NaN;
                    Vs_u(Asu_next < min(A_wide)) = NaN;
                    Vs_u(Asu_next > max(A_wide)) = NaN;
                    % save optimal
                    [Vs_r_star, index_sr_k] = max(Vs_r);
                    [Vs_n_star, index_sn_k] = max(Vs_n);
                    [Vs_u_star, index_su_k] = max(Vs_u);
                    [Vsm_r_star, index_smr_k] = max(Vsm_r);
                    [Vsm_n_star, index_smn_k] = max(Vsm_n);
                    [Vsm_u_star, index_smu_k] = max(Vsm_u);
                    cs_r_star = cr_vector(index_sr_k);
                    cs_n_star = cn_vector(index_sn_k);
                    cs_u_star = cu_vector(index_su_k);
                    csm_r_star = cr_vector(index_smr_k);
                    csm_n_star = cn_vector(index_smn_k);
                    csm_u_star = cu_vector(index_smu_k);
                    cs_star_aux = [csm_r_star, csm_n_star, csm_u_star, cs_r_star, cs_n_star, cs_u_star];
                    % check if V_star are the same for single and married
                    max_m = max([Vsm_r_star, Vsm_n_star, Vsm_u_star]);
                    max_s = max([Vs_r_star, Vs_n_star, Vs_u_star]);
                    % if they are the same, choose from the single options
                    if max_m == max_s
                        [Vs_star, ls_index] =  max([Vs_r_star, Vs_n_star, Vs_u_star]);
                        [Vs_star2] = max([Vs_n_star, Vs_u_star]);
                        ls_index = ls_index + 3;
                    else
                        [Vs_star, ls_index] = max([Vsm_r_star, Vsm_n_star, Vsm_u_star, Vs_r_star, Vs_n_star, Vs_u_star]);
                        [Vs_star2] = max([Vsm_n_star, Vsm_u_star, Vs_n_star, Vs_u_star]);
                    end
                end

            % save choice:
            if x <= 10
                c_star(j, i, x, t) = cm_star_aux(lm_index);
                l_star(j, i, x, t) = lm_index;
                m_star(j, i, x, t) = 1; % married
                if lm_index<=3
                    ch_star(j,i,x,t) = 0; % 1 child
                else
                    ch_star(j,i,x,t) = 1; % 2 children
                end
                V_star(j, i, x, t) = Vm_star;
                V2_star(j,i, x, t) = Vm_star2;
            elseif x <= 20
                c_star(j, i, x, t) = cm2_star_aux(lm2_index);
                l_star(j, i, x, t) = lm2_index;
                m_star(j, i, x, t) = 1; % married
                ch_star(j,i, x, t) = 2; % 2 children
                V_star(j, i, x, t) = Vm2_star;
                V2_star(j,i, x, t) = Vm2_star2;
            else
                c_star(j, i, x, t) = cs_star_aux(ls_index);
                l_star(j, i, x, t) = ls_index;
                m_star(j, i, x, t) = (ls_index<=3); % married OR single (0)
                ch_star(j,i, x, t) = (ls_index<=3); % 1 child OR 0 children
                V_star(j, i, x, t) = Vs_star;
                V2_star(j,i, x, t) = Vs_star2;
            end
            
            % save the number assets outside grid
            if x <= 10
                Ar_out(j,i,x,t) = sum(Amr_next < min(A_wide)) + sum(Amr_next > max(A_wide));
                An_out(j,i,x,t) = sum(Amn_next < min(A_wide)) + sum(Amn_next > max(A_wide));
                Au_out(j,i,x,t) = sum(Amu_next < min(A_wide)) + sum(Amu_next > max(A_wide));
            elseif x <= 20
                Ar_out(j,i,x,t) = sum(Am2r_next < min(A_wide)) + sum(Am2r_next > max(A_wide));
                An_out(j,i,x,t) = sum(Am2n_next < min(A_wide)) + sum(Am2n_next > max(A_wide));
                Au_out(j,i,x,t) = sum(Am2u_next < min(A_wide)) + sum(Am2u_next > max(A_wide));
            else
                Ar_out(j,i,x,t) = sum(Asr_next < min(A_wide)) + sum(Asr_next > max(A_wide));
                An_out(j,i,x,t) = sum(Asn_next < min(A_wide)) + sum(Asn_next > max(A_wide));
                Au_out(j,i,x,t) = sum(Asu_next < min(A_wide)) + sum(Asu_next > max(A_wide));
            end
            end
        end
        
        % Integrate over shocks
        W(:,x,t) = pi^(-1/2)*V_star(:,:,x,t)*S.weight;
        W2(:,x,t)= pi^(-1/2)*V2_star(:,:,x,t)*S.weight;
        
        % reshape policy func
        c_func(:,:,:,:,x,t) = reshape(c_star(:,:,x,t), [G.n_assets,3,3,3]);
        l_func(:,:,:,:,x,t) = reshape(l_star(:,:,x,t), [G.n_assets,3,3,3]);
        m_func(:,:,:,:,x,t) = reshape(m_star(:,:,x,t), [G.n_assets,3,3,3]);
        ch_func(:,:,:,:,x,t)= reshape(ch_star(:,:,x,t), [G.n_assets,3,3,3]);
    end
end

% three labor functions (as 0 or 1)
lr_func = l_func == 1 | l_func == 4;
ln_func = l_func == 2 | l_func == 5;
lu_func = l_func == 3 | l_func == 6;

end
