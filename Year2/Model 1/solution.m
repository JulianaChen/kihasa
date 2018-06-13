function [c_func,m_func,lr_func,ln_func,lu_func] = solution(G,abi,edu,S,params)

% Index for parameters
% params0 = [psi_r;psi_n;theta;alpha;sigma_r;sigma_n;sigma_i;omega;lambda;eta;iota;kappa;tau];
    
    psi_r=params(1); % disutility of work by sector (regular)
    psi_n=params(2); % disutility of work by sector (non-regular)
    
    theta1_r=params(3); % value of marriage in HH production, regular
    theta1_n=params(4); % value of marriage in HH production, non-regular
    theta1_u=params(5); % value of marriage in HH production, unemployed
    theta3_r=params(6); % value of child HC in HH production, regular
    theta3_n=params(7); % value of child HC in HH production, non-regular
    theta3_u=params(8); % value of child HC in HH production, unemployed
    
    alpha01_r=params(9); % wage return for the ability type, regular
    alpha01_n=params(10); % wage return for the ability type, non-regular
    alpha02_r=params(11); % additional return for high ability type, regular
    alpha02_n=params(12); % additional return for high ability type, non-regular
    alpha11_r=params(13); % wage return to 2yr college, regular
    alpha11_n=params(14); % wage return to 2yr college, non-regular
    alpha12_r=params(15); % wage return to 4yr college, regular
    alpha12_n=params(16); % wage return to 4yr college, non-regular
    alpha2_r=params(17); % wage return to general work experience, regular
    alpha2_n=params(18); % wage return to general work experience, non-regular
    
    sigma_r = params(19); % shock, regular
    sigma_n = params(20); % shock, non-regular
    sigma_i = params(21); % shock, unemployed
    
    omega0_w = params(22); % probability of marriage for workers
    omega0_u = params(23); % probability of marriage for non-workers
    omega11  = params(24); % probability of marriage for 2yr college
    omega12  = params(25); % probability of marriage for 4yr college
    omega2 = params(26); % probability of marriage for age
    
    lambda1 = params(27); % terminal value function (for assets)
    lambda2 = params(28); % terminal value function (for HH income)
    lambda3 = params(29); % terminal value function (for HH production)
    lambda4 = params(30); % terminal value function (for work exp)    
    
    eta01 = params(31); % husgand's wage return to low ability type (mean)
    eta02 = params(32); % husgand's wage return to high ability type (mean)
    eta11 = params(33); % husgand's wage return to 2yr college (mean)
    eta12 = params(34); % husgand's wage return to 4yr college (mean)
    eta2 = params(35); % husgand's wage return to women's age (mean)
    eta03 = params(36); % husgand's wage return to low ability type (sd)
    eta04 = params(37); % husgand's wage return to high ability type (sd)
    eta21 = params(38); % husgand's wage return to 2yr college (sd)
    eta22 = params(39); % husgand's wage return to 4yr college (sd)
    eta3 = params(40); % husgand's wage return to women's age (sd)
     
    iota01 = params(41); % child investment of low ability type (mean)
    iota02 = params(42); % child investment of high ability type (mean)
    iota11 = params(43); % child investment of 2yr college (mean)
    iota12 = params(44); % child investment of 4yr college (mean)
    iota2 = params(45); % child investment by women's age (mean)
    iota03 = params(46); % child investment of low ability type (sd)
    iota04 = params(47); % child investment of high ability type (sd)
    iota21 = params(48); % child investment of 2yr college (sd)
    iota22 = params(49); % child investment of 4yr college (sd)
    iota3 = params(50); % child investment by women's age (sd)
    
    kappa01 = params(51); % child human capital for low ability type 
    kappa02 = params(52); % child human capital for high ability type
    kappa03 = params(53); % child human capital for 2yr college
    kappa04 = params(54); % child human capital for 4yr college
    kappa05 = params(55); % child human capital for cihld investment level
    
    tau10 = params(56); % probability of losing a regular job
    tau11 = params(57); % probability of losing a regular job (2yr college)
    tau12 = params(58); % probability of losing a regular job (4yr college)
    tau13 = params(59); % probability of losing a regular job (age)
    tau14 = params(60); % probability of losing a regular job (work exp)
    tau20 = params(61); % probability of finding a regular job
    tau21 = params(62); % probability of losing a regular job (2yr college)
    tau22 = params(63); % probability of losing a regular job (4yr college)
    tau23 = params(64); % probability of losing a regular job (age)
    tau24 = params(65); % probability of losing a regular job (work exp)
    % move to params
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
    
%% Terminal Value Function:

% age at TVF
age_TVF = 18*(edu==1) + 20*(edu==2) + 22*(edu==3) + G.n_period;

% assets
assets = S.SS_A; 
    
% husband wages
wh_mean = eta01 + eta02*(abi==2) + eta11*(edu==2) + eta12*(edu==3) + eta2*age_TVF;
wh_sd = eta03 + eta04*(abi==2) + eta21*(edu==2) + eta22*(edu==3) + eta3*age_TVF;
wh = normrnd(wh_mean,wh_sd);
    
% child human capital 
Inv_mean = iota01 + iota02*(abi==2) + iota11*(edu==2) + iota12*(edu==3) + iota2*age_TVF;
Inv_sd = iota03 + iota04*(abi==2) + iota21*(edu==2) + iota22*(edu==3) + iota3*age_TVF;
Inv = normrnd(Inv_mean,Inv_sd);

K = kappa01 + kappa02*(abi==2) + kappa03*(edu==2) + kappa04*(edu==3) + kappa05*(Inv); 
% issues: K is negative, can it be? Change TVF?

% TVF 
TVF = lambda1*(1-exp(-assets)) + lambda2*(unique(S.SS_X).^(1-G.sigma))/(1-G.sigma) + lambda3*(wh^(1-G.sigma))/(1-G.sigma) + lambda4*(1-exp(-K));
% issues: (1) check assets function is concave (2) need lambda4 for K (3)
% check all lambdas 
% note: changed lambda4*(K^(1-G.sigma))/(1-G.sigma) to lambda1*(1-exp(-K))
% bc of negative K?

tic
% loop for time (20):
for t = G.n_period-1:-1:1
    t
    toc
    
    % age
    age = 18*(edu==1) + 20*(edu==2) + 22*(edu==3) + t;
    
    % draw husband wage
    wh_mean = eta01 + eta02*(abi==2) + eta11*(edu==2) + eta12*(edu==3) + eta2*age; %function of women's age, education and ability types
    wh_sd = eta03 + eta04*(abi==2) + eta21*(edu==2) + eta22*(edu==3) + eta3*age; %same as above
    wh = normrnd(wh_mean,wh_sd); 
    
    % draw investments
    Inv_mean = iota01 + iota02*(abi==2) + iota11*(edu==2) + iota12*(edu==3) + iota2*age;
    Inv_sd = iota03 + iota04*(abi==2) + iota21*(edu==2) + iota22*(edu==3) + iota3*age;
    Inv = normrnd(Inv_mean,Inv_sd); 
    
    % child human capital 
    K = kappa01 + kappa02*(abi==2) + kappa03*(edu==2) + kappa04*(edu==3) + kappa05*(Inv); %function of women's ability, education
  
    % marriage probabilities 
    prob_marr_w = omega0_w + omega11*(edu==2) + omega12*(edu==3) + omega2*age;
    prob_marr_u = omega0_u + omega11*(edu==2) + omega12*(edu==3) + omega2*age;
    
    % loop for work experience and marital status (20): %%%%% change to 30
    for x = 1:1:(G.n_matstat*G.n_wrkexp)
        x;
    
        % current state discrete variables:
        m_j = S.SS_M(x);  % marital status
        n_j = S.SS_N(x);  % children
        X_j = S.SS_X(x);  % experience
    
        % job probabilities %% how to make sure the prob is in [0 1] - use
        % probit 
        prob_lamba = tau10 + tau11*(edu==2) + tau12*(edu==3) + tau13*age + tau14*X_j; % probability of losing a regular job
        prob_pi = tau20 + tau21*(edu==2) + tau22*(edu==3) + tau23*age + tau24*X_j; % probability of getting a regular job
    
        % child probabilities  %% how to make sure the prob is in [0 1]
        prob_2kids_r = phi10 + phi11*(edu==2) + phi12*(edu==3) + phi13*X_j;
        prob_2kids_n = phi20 + phi21*(edu==2) + phi22*(edu==3) + phi23*X_j;
        prob_2kids_u = phi30 + phi31*(edu==2) + phi32*(edu==3) + phi33*X_j;
        
        % loop for shocks (9):
        for i = 1:1:G.n_shocks
            i;
        
            shock_i = S.shocks_i(i);
            shock_r = S.shocks_r(i);
            shock_n = S.shocks_n(i);
            
            % sector-specific state variables: - take out of assets
            w_j_r = exp(alpha01_r + alpha02_r*(abi==2) + alpha11_r*(edu==2) + alpha12_r*(edu==3) + alpha2_r*log(1+X_j) + shock_r);
            w_j_n = exp(alpha01_n + alpha02_n*(abi==2) + alpha11_n*(edu==2) + alpha12_n*(edu==3) + alpha2_n*log(1+X_j) + shock_n);  
            w_j_u = 0; % unemployed women don't have earnings
              
            % loop over assets
            for j = 1:1:G.n_SS
                j;
            
                % HH's assets
                A_j = S.SS_A(j); 
            
                % consumption vector
                chh_min = 0.1; %  move outside
                chh_r = w_j_r + wh*m_j + A_j;
                chh_n = w_j_n + wh*m_j + A_j;
                chh_u = w_j_u + wh*m_j + A_j; % + transfer
                chh_r_max = max(chh_min,chh_r);
                chh_n_max = max(chh_min,chh_n);
                chh_u_max = max(chh_min,chh_u);
                cr_vector = linspace(chh_min,chh_r_max,G.n_cons);
                cn_vector = linspace(chh_min,chh_n_max,G.n_cons);
                cu_vector = linspace(chh_min,chh_u_max,G.n_cons);

                % wide assets vector for linear interpolation % move -10 to params
                A_wide = linspace(-10,max(S.SS_A),length(unique(S.SS_A))); % need?
                
                % loop over consumption
                for k = 1:1:G.n_cons
                    k;
                    
                    % HH consumption
                    chh_r = cr_vector(k);
                    chh_n = cn_vector(k);
                    chh_u = cu_vector(k);
                    % woman's consumption
                    cw_r = 0.5*chh_r;
                    cw_n = 0.5*chh_n;
                    cw_u = 0.5*chh_u;
                    
                    % Sector-Specific Utility: %% complex number bc -K!!!!!
                    u_r(k) = (cw_r^(1-G.sigma))/(1-G.sigma) + psi_r + theta1_r*log(1+m_j) + theta3_r*log(K)*n_j;
                    u_n(k) = (cw_n^(1-G.sigma))/(1-G.sigma) + psi_n + theta1_n*log(1+m_j) + theta3_n*log(K)*n_j;
                    u_u(k) = (cw_u^(1-G.sigma))/(1-G.sigma) + theta1_u*log(1+m_j) + theta3_u*log(K)*n_j;
                    
                    % Maried with 1 kid:
                    if x <= 10
                        
                        % regular job:
                        A_next = (1+G.r) * (A_j + (w_j_r + wh*m_j + shock_i) - chh_r - n_j*Inv); % eq. 8
                        x_next = x + 1;
                        if x_next == 11
                            x_next = 10;
                        end
                        % value function:
                        if t==G.n_period-1
                            Emax = TVF;
                            Emax2 = TVF;
                        else
                            Emax = W(:,x_next,t+1);
                            Emax2 = W2(:,x_next,t+1);
                        end
                        % linear approximation of VF
                        Vm_r_next_linear = interpn(A_wide,Emax,A_next);
                        Vm_r_next_linear2 = interpn(A_wide,Emax2,A_next);
                        Amr_next(k)=A_next;
                        Vmr_next_linear(k,x)=Vm_r_next_linear;
                        Vmr_next_linear2(k,x)=Vm_r_next_linear2;
                        
                        % non-regular job:
                        A_next = (1+G.r) * (A_j + (w_j_n + wh*m_j + shock_i) - chh_n - n_j*Inv);
                        x_next = x + 1;
                        if x_next == 11
                            x_next = 10;
                        end
                        % value function:
                        if t==G.n_period-1
                            Emax = TVF;
                            Emax2 = TVF;
                        else
                            Emax = W(:,x_next,t+1);
                            Emax2 = W2(:,x_next,t+1);
                        end
                        % linear approximation of VF
                        Vm_n_next_linear = interpn(A_wide,Emax,A_next);
                        Vm_n_next_linear2 = interpn(A_wide,Emax2,A_next);                                            
                        Amn_next(k)=A_next;
                        Vmn_next_linear(k,x)=Vm_n_next_linear;
                        Vmn_next_linear2(k,x)=Vm_n_next_linear2;
                        
                        % unemployed:
                        A_next = (1+G.r) * (A_j + (w_j_u + wh*m_j + shock_i) - chh_u - n_j*Inv);
                        x_next = x;
                        % value function:
                        if t==G.n_period-1
                            Emax = TVF;
                            Emax2 = TVF;
                        else
                            Emax = W(:,x_next,t+1);
                            Emax2 = W2(:,x_next,t+1);
                        end
                        % linear approximation of VF
                        Vm_u_next_linear = interpn(A_wide,Emax,A_next); 
                        Vm_u_next_linear2 = interpn(A_wide,Emax2,A_next);                                                
                        Amu_next(k)=A_next;
                        Vmu_next_linear(k,x)=Vm_u_next_linear;
                        Vmu_next_linear2(k,x)=Vm_u_next_linear2;
                        
                        % Sector-Specific Value Functions
                        Vm_r(k) = u_r(k) + G.beta * ((prob_lamba*Vm_r_next_linear)+(1-prob_lamba)*Vm_r_next_linear2);
                        Vm_n(k) = u_n(k) + G.beta * ((prob_pi*Vm_n_next_linear)+(1-prob_pi)*Vm_n_next_linear2);
                        Vm_u(k) = u_u(k) + G.beta * ((prob_pi*Vm_u_next_linear)+(1-prob_pi)*Vm_u_next_linear2);

                    % Married with 2 kids: %% Julie will add this part
                    elseif x <= 20
                        
                        % regular job:
                        A_next = (1+G.r) * (A_j + (w_j_r + wh*m_j + shock_i) - chh_r - n_j*Inv); % eq. 8
                        x_next = x + 1;
                        if x_next == 21
                            x_next = 20;
                        end
                        % value function:
                        if t==G.n_period-1
                            Emax = TVF;
                            Emax2 = TVF;
                        else
                            Emax = W(:,x_next,t+1);
                            Emax2 = W2(:,x_next,t+1);
                        end
                        % linear approximation of VF
                        Vm2_r_next_linear = interpn(A_wide,Emax,A_next);
                        Vm2_r_next_linear2 = interpn(A_wide,Emax2,A_next);
                        Am2r_next(k)=A_next;
                        Vm2r_next_linear(k,x)=Vm2_r_next_linear;
                        Vm2r_next_linear2(k,x)=Vm2_r_next_linear2;
                        
                        % non-regular job:
                        A_next = (1+G.r) * (A_j + (w_j_n + wh*m_j + shock_i) - chh_n - n_j*Inv);
                        x_next = x + 1;
                        if x_next == 21
                            x_next = 20;
                        end
                        % value function:
                        if t==G.n_period-1
                            Emax = TVF;
                            Emax2 = TVF;
                        else
                            Emax = W(:,x_next,t+1);
                            Emax2 = W2(:,x_next,t+1);
                        end
                        % linear approximation of VF
                        Vm2_n_next_linear = interpn(A_wide,Emax,A_next); 
                        Vm2_n_next_linear2 = interpn(A_wide,Emax2,A_next);
                        Am2n_next(k)=A_next;
                        Vm2n_next_linear(k,x)=Vm2_n_next_linear;
                        Vm2n_next_linear2(k,x)=Vm2_n_next_linear2;
                        
                        % unemployed:
                        A_next = (1+G.r) * (A_j + (w_j_u + wh*m_j + shock_i) - chh_u - n_j*Inv);
                        x_next = x;
                        % value function:
                        if t==G.n_period-1
                            Emax = TVF;
                            Emax2 = TVF;
                        else
                            Emax = W(:,x_next,t+1);
                            Emax2 = W2(:,x_next,t+1);
                        end
                        % linear approximation of VF
                        Vm2_u_next_linear = interpn(A_wide,Emax,A_next); 
                        Vm2_u_next_linear2 = interpn(A_wide,Emax2,A_next);
                        Am2u_next(k)=A_next;
                        Vm2u_next_linear(k,x)=Vm2_u_next_linear;
                        Vm2u_next_linear2(k,x)=Vm2_u_next_linear2;
                        
                        % Sector-Specific Value Functions
                        Vm2_r(k) = u_r(k) + G.beta * ((prob_lamba*Vm2_r_next_linear)+(1-prob_lamba)*Vm2_r_next_linear2);
                        Vm2_n(k) = u_n(k) + G.beta * ((prob_pi*Vm2_n_next_linear)+(1-prob_pi)*Vm2_n_next_linear2);
                        Vm2_u(k) = u_u(k) + G.beta * ((prob_pi*Vm2_u_next_linear)+(1-prob_pi)*Vm2_u_next_linear2);
                        % save marriage values (for marriage decision)
                        Vm_r_aux(k,x) = Vm2_r(k)*(prob_2kids_r) + Vm_r(k)*(1-prob_2kids_r);
                        Vm_n_aux(k,x) = Vm2_n(k)*(prob_2kids_n) + Vm_n(k)*(1-prob_2kids_n);
                        Vm_u_aux(k,x) = Vm2_u(k)*(prob_2kids_u) + Vm_u(k)*(1-prob_2kids_u);
                    
                    % Single:
                    else
                        
                        % Regular:
                        A_next = (1+G.r) * (A_j + (w_j_r + wh*m_j + shock_i) - chh_r - n_j*Inv);
                        x_next = x + 1;
                        if x_next == 21
                            x_next = 20;
                        end
                        % value function:
                        if t==G.n_period-1
                            Emax = TVF;
                            Emax2 = TVF;
                        else
                            Emax = W(:,x_next,t+1);
                            Emax2 = W2(:,x_next,t+1);
                        end
                        % linear approximation of VF
                        Vs_r_next_linear = interpn(A_wide,Emax,A_next);
                        Vs_r_next_linear2 =interpn(A_wide,Emax2,A_next);
                        Asr_next(k)=A_next;
                        Vsr_next_linear(k,x)=Vs_r_next_linear;
                        Vsr_next_linear2(k,x)=Vs_r_next_linear2;
                        
                        % Non-regular:
                        A_next = (1+G.r) * (A_j + (w_j_n + wh*m_j + shock_i) - chh_n - n_j*Inv);
                        x_next = x + 1;
                        if x_next == 21
                            x_next = 20;
                        end
                        % value function:
                        if t==G.n_period-1
                            Emax = TVF;
                            Emax2 = TVF;
                        else
                            Emax = W(:,x_next,t+1);
                            Emax2 = W2(:,x_next,t+1);
                        end
                        % linear approximation of VF
                        Vs_n_next_linear = interpn(A_wide,Emax,A_next);
                        Vs_n_next_linear2 = interpn(A_wide,Emax2,A_next);
                        Asn_next(k)=A_next;
                        Vsn_next_linear(k,x)=Vs_n_next_linear;
                        Vsn_next_linear2(k,x)=Vs_n_next_linear2;
                        
                        % Unemployed:
                        A_next = (1+G.r) * (A_j + (w_j_u + wh*m_j + shock_i) - chh_u - n_j*Inv);
                        x_next = x;
                        % value function:
                        if t==G.n_period-1
                            Emax = TVF;
                            Emax2 = TVF;
                        else
                            Emax = W(:,x_next,t+1);
                            Emax2 = W2(:,x_next,t+1);
                        end
                        % linear approximation of VF
                        Vs_u_next_linear = interpn(A_wide,Emax,A_next);
                        Vs_u_next_linear2 = interpn(A_wide,Emax2,A_next);
                        Asu_next(k)=A_next;
                        Vsu_next_linear(k,x)=Vs_u_next_linear;
                        Vsu_next_linear2(k,x)=Vs_u_next_linear2;
                        
                        % Sector-Specific Value Functions
                        Vs_r(k) = u_r(k) + G.beta * ((prob_lamba*Vs_r_next_linear)+(1-prob_lamba)*Vs_r_next_linear2);
                        Vs_n(k) = u_n(k) + G.beta * ((prob_pi*Vs_n_next_linear)+(1-prob_pi)*Vs_n_next_linear2);
                        Vs_u(k) = u_u(k) + G.beta * ((prob_pi*Vs_u_next_linear)+(1-prob_pi)*Vs_u_next_linear2);
                        Vsm_r(k) = prob_marr_w*Vm_r_aux(k,x-10) + (1-prob_marr_w)*Vs_r(k);
                        Vsm_n(k) = prob_marr_w*Vm_n_aux(k,x-10) + (1-prob_marr_w)*Vs_n(k);
                        Vsm_u(k) = prob_marr_u*Vm_u_aux(k,x-10) + (1-prob_marr_u)*Vs_u(k);
                    end
                end

                % optimal consumption and max Emax
                if x <= 10
                    % check
                    Vm_r(Amr_next < -10) = NaN;
                    Vm_r(Amr_next > max(S.SS_A)) = NaN;
                    Vm_n(Amn_next < -10) = NaN;
                    Vm_n(Amn_next > max(S.SS_A)) = NaN;
                    Vm_u(Amu_next < -10) = NaN;
                    Vm_u(Amu_next > max(S.SS_A)) = NaN;
                    % save optimal
                    [Vm_r_star, Index_mr_k] = max(Vm_r);
                    [Vm_n_star, Index_mn_k] = max(Vm_n);
                    [Vm_u_star, Index_mu_k] = max(Vm_u);
                    cm_r_star = cr_vector(Index_mr_k);
                    cm_n_star = cn_vector(Index_mn_k);
                    cm_u_star = cu_vector(Index_mu_k);
                    cm_star_aux = [cm_r_star, cm_n_star, cm_u_star];
                    [Vm_star, Index_lm] = max([Vm_r_star,Vm_n_star,Vm_u_star]); % with 3 job options
                    [Vm_star2, Index_lm2] = max([Vm_n_star,Vm_u_star]); % is this emax2 (2 job options?)
                else %%%% add 2 kids?
                    % check
                    Vs_r(Asr_next < -10) = NaN;
                    Vs_r(Asr_next > max(S.SS_A)) = NaN;
                    Vs_n(Asn_next < -10) = NaN;
                    Vs_n(Asn_next > max(S.SS_A)) = NaN;
                    Vs_u(Asu_next < -10) = NaN;
                    Vs_u(Asu_next > max(S.SS_A)) = NaN;
                    % save optimal
                    [Vs_r_star, Index_sr_k] = max(Vs_r);
                    [Vs_n_star, Index_sn_k] = max(Vs_n);
                    [Vs_u_star, Index_su_k] = max(Vs_u);
                    [Vsm_r_star, Index_smr_k] = max(Vsm_r);
                    [Vsm_n_star, Index_smn_k] = max(Vsm_n);
                    [Vsm_u_star, Index_smu_k] = max(Vsm_u);
                    cs_r_star = cr_vector(Index_sr_k);
                    cs_n_star = cn_vector(Index_sn_k);
                    cs_u_star = cu_vector(Index_su_k);
                    csm_r_star = cr_vector(Index_smr_k);
                    csm_n_star = cn_vector(Index_smn_k);
                    csm_u_star = cu_vector(Index_smu_k);
                    cs_star_aux = [csm_r_star, csm_n_star, csm_u_star, cs_r_star, cs_n_star, cs_u_star];
                    [Vs_star, Index_ls] = max([Vsm_r_star, Vsm_n_star, Vsm_u_star, Vs_r_star, Vs_n_star, Vs_u_star]);
                    [Vs_star2, Index_ls2] = max([Vsm_n_star, Vsm_u_star, Vs_n_star, Vs_u_star]); % is this emax2?
                end
            % save choice:
            if x <= 10 %20
                c_star(j, i, x, t) = cm_star_aux(Index_lm);
                l_star(j, i, x, t) = Index_lm; % which L choice?
                V_star(j, i, x, t) = Vm_star;
                V2_star(j,i, x, t) = Vm_star2;
            else %%%% add 2 kids?
                c_star(j, i, x, t) = cs_star_aux(Index_ls);
                l_star(j, i, x, t) = Index_ls; % which L choice?
                V_star(j, i, x, t) = Vs_star;
                V2_star(j,i, x, t) = Vs_star2;
            end
            
            % save the number assets outside grid
            if x <= 10
                Ar_out(j,i,x,t) = sum(Amr_next < -10) + sum(Amr_next > max(S.SS_A));
                An_out(j,i,x,t) = sum(Amn_next < -10) + sum(Amn_next > max(S.SS_A));
                Au_out(j,i,x,t) = sum(Amu_next < -10) + sum(Amu_next > max(S.SS_A));
            else %%%% add 2 kids?
                Ar_out(j,i,x,t) = sum(Asr_next < -10) + sum(Asr_next > max(S.SS_A));
                An_out(j,i,x,t) = sum(Asn_next < -10) + sum(Asn_next > max(S.SS_A));
                Au_out(j,i,x,t) = sum(Asu_next < -10) + sum(Asu_next > max(S.SS_A));
            end
            end
        end
        
        % Integrate over shocks
        W(:,x,t) = pi^(-1/2)*V_star(:,:,x,t)*S.weight;
        W2(:,x,t)= pi^(-1/2)*V2_star(:,:,x,t)*S.weight;
        
        % reshape policy func
        v_func(:,x,t) = reshape(V_star(:,:,x,t),[],1);
        v2_func(:,x,t)= reshape(V2_star(:,:,x,t),[],1);
        c_func(:,x,t) = reshape(c_star(:,:,x,t),[],1);
        l_func(:,x,t) = reshape(l_star(:,:,x,t),[],1);
    end
end

% three labor functions (as 0 or 1)
lr_func = l_func == 1 | l_func == 4;
ln_func = l_func == 2 | l_func == 5;
lu_func = l_func == 3 | l_func == 6;

% marriage function for single women: equals 1 if labor choice is 1, 2, or 3
m_func = l_func == 1 | l_func == 2 | l_func == 3;
