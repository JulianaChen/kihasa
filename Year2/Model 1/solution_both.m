%function [c_func,m_func,lr_func,ln_func,lu_func,Ar_out,An_out,Au_out,wh_aux,w_j_r_aux,w_j_n_aux] = solution_both(G,abi,edu,S,params)

%% index for parameters:

% params0 = [psi_r;psi_n;theta;alpha;sigma_r;sigma_n;sigma_i;omega;lambda;eta;iota;kappa;tau];
    
    psi_r=params(1); % disutility of work by sector (regular)
    psi_n=params(2); % disutility of work by sector (non-regular)
    psi_u= -1;       % disutility of work by sector (unemployment)
    
    theta1_r=params(3); % value of marriage in HH production, regular
    theta1_n=params(4); % value of marriage in HH production, non-regular
    theta1_u=params(5); % value of marriage in HH production, unemployed
    theta2_r=1.5;       % value of children in HH production, regular
    theta2_n=1.5;       % value of children in HH production, non-regular
    theta2_u=2;         % value of children in HH production, unemployed
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
    
    phi10 = params(66); % probability of second child
    phi11 = params(67); % probability of second child
    phi12 = params(68); % probability of second child
    phi13 = params(69); % probability of second child
    phi20 = params(70); % probability of second child
    phi21 = params(71); % probability of second child
    phi22 = params(72); % probability of second child
    phi23 = params(73); % probability of second child
    phi30 = params(74); % probability of second child
    phi31 = params(75); % probability of second child
    phi32 = params(76); % probability of second child
    phi33 = params(77); % probability of second child

%% other parameters:

chh_min = 0.1; % minimun consumption
delta = 0.5; % Female Share of Consumption (CAL)

% expanded assets vector for linear interpolation
% A_min = min(S.SS_A)-20;
% A_max = max(S.SS_A)+1000;
A_min=S.extmin_A;
A_max=S.extmax_A;
A_wide = S.SS_A;
A_wide(1) = A_min;
A_wide(10) = A_max;

%% Terminal Value Function:

% age at TVF
age_TVF = 18*(edu==1) + 20*(edu==2) + 22*(edu==3) + G.n_period;

% assets
assets = S.SS_A; 

% husband wages
wh_mean = eta01 + eta02*(abi==2) + eta11*(edu==2) + eta12*(edu==3) + eta2*age_TVF;
wh_sd = 0.7022257; %eta03 + eta04*(abi==2) + eta21*(edu==2) + eta22*(edu==3) + eta3*age_TVF;
wh = normrnd(wh_mean,wh_sd);

% child human capital 
Inv_mean = iota01 + iota02*(abi==2) + iota11*(edu==2) + iota12*(edu==3) + iota2*age_TVF;
Inv_sd = 0.9270494; %iota03 + iota04*(abi==2) + iota21*(edu==2) + iota22*(edu==3) + iota3*age_TVF;
Inv = normrnd(Inv_mean,Inv_sd);

K = exp(kappa01 + kappa02*(abi==2) + kappa03*(edu==2) + kappa04*(edu==3) + kappa05*(Inv));

% TVF 
TVF = repmat(real(lambda1*(assets).^(1-G.sigma)/(1-G.sigma))',1,30) + repmat(lambda2*(S.SS_X.^(1-G.sigma))/(1-G.sigma),10,1) ...
+ lambda3*(wh.^(1-G.sigma))/(1-G.sigma) + lambda4*(K.^(1-G.sigma))/(1-G.sigma);

tic
% loop for time (20):
for t = G.n_period-1:-1:1
    t
    %toc
    
    % Coefficients for Chebyshev Approximation
    
    if t==G.n_period-1
        Emax = TVF;
        Emax2 = TVF;
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
    wh = normrnd(wh_mean,wh_sd); 
    wh_aux(t)=wh; %%% SAVE FOR LATER
    % draw investments
    Inv_mean = iota01 + iota02*(abi==2) + iota11*(edu==2) + iota12*(edu==3) + iota2*age;
    Inv_sd = 0.9270494; %iota03 + iota04*(abi==2) + iota21*(edu==2) + iota22*(edu==3) + iota3*age;
    Inv = normrnd(Inv_mean,Inv_sd); 
    
    % child human capital 
    K = exp(kappa01 + kappa02*(abi==2) + kappa03*(edu==2) + kappa04*(edu==3) + kappa05*(Inv));
    
    % marriage probabilities 
    prob_marr_w = normcdf(omega0_w + omega11*(edu==2) + omega12*(edu==3) + omega2*age);
    prob_marr_u = normcdf(omega0_u + omega11*(edu==2) + omega12*(edu==3) + omega2*age);
    
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
        
        % child probabilities
        prob_2kids_r = normcdf(phi10 + phi11*(edu==2) + phi12*(edu==3) + phi13*X_j);
        prob_2kids_n = normcdf(phi20 + phi21*(edu==2) + phi22*(edu==3) + phi23*X_j);
        prob_2kids_u = normcdf(phi30 + phi31*(edu==2) + phi32*(edu==3) + phi33*X_j);
        
        % loop for shocks (9):
        for i = 1:1:G.n_shocks
            i;
            
            % shocks
            shock_i = S.shocks_i(i);
            shock_r = S.shocks_r(i);
            shock_n = S.shocks_n(i);
            
            % sector-specific state variables
            w_j_r = exp(alpha01_r + alpha02_r*(abi==2) + alpha11_r*(edu==2) + alpha12_r*(edu==3) + alpha2_r*log(1+X_j) + shock_r);
            w_j_n = exp(alpha01_n + alpha02_n*(abi==2) + alpha11_n*(edu==2) + alpha12_n*(edu==3) + alpha2_n*log(1+X_j) + shock_n);  
            w_j_u = 0; % unemployed women don't have earnings
            w_j_r_aux(i,x,t)=w_j_r; %%% save for later
            w_j_n_aux(i,x,t)=w_j_n; %%% save for later
            
            % loop over assets (10):
            for j = 1:1:G.n_assets
                j;
                
                % HH's assets
                A_j = S.SS_A(j); 
                
                % consumption vector
                chh_r = w_j_r + exp(wh)*m_j + A_j;
                chh_n = w_j_n + exp(wh)*m_j + A_j;
                chh_u = w_j_u + exp(wh)*m_j + A_j;
                chh_r_max = max(chh_min,chh_r);
                chh_n_max = max(chh_min,chh_n);
                chh_u_max = max(chh_min,chh_u);
                cr_vector = linspace(chh_min,chh_r_max,G.n_cons);
                cn_vector = linspace(chh_min,chh_n_max,G.n_cons);
                cu_vector = linspace(chh_min,chh_u_max,G.n_cons);
                
                % loop over consumption (10):
                for k = 1:1:G.n_cons
                    k;
                    
                    % HH consumption
                    chh_r = cr_vector(k);
                    chh_n = cn_vector(k);
                    chh_u = cu_vector(k);
                    
                    % woman's consumption
                    cw_r = delta*chh_r;
                    cw_n = delta*chh_n;
                    cw_u = delta*chh_u;
                    
                    % Sector-Specific Utility
                    u_r(k) = (cw_r^(1-G.sigma))/(1-G.sigma) + psi_r + theta1_r*log(1+m_j) + theta2_r*log(1+n_j) + theta3_r*log(K)*n_j;
                    u_n(k) = (cw_n^(1-G.sigma))/(1-G.sigma) + psi_n + theta1_n*log(1+m_j) + theta2_n*log(1+n_j) + theta3_n*log(K)*n_j;
                    u_u(k) = (cw_u^(1-G.sigma))/(1-G.sigma) + psi_u + theta1_u*log(1+m_j) + theta2_u*log(1+n_j) + theta3_u*log(K)*n_j;
                    
                    % Maried with 1 kid:
                    if x <= 10
                        
                        % regular job
                        A_next = (1+G.r) * (A_j + (w_j_r + exp(wh)*m_j + shock_i) - chh_r - n_j*Inv); % eq. 8
                        x_next = x + 1;
                        if x_next == 11
                            x_next = 10;
                        end

                        % linear approximation of VF
                        Vm_r_next_linear = interpn(A_wide,Emax(:,x_next),A_next);
                        Vm_r_next_linear2 = interpn(A_wide,Emax2(:,x_next),A_next);
                        Vm_r_next_linear_2 = interpn(A_wide,Emax(:,x_next+10),A_next);
                        Vm_r_next_linear2_2 = interpn(A_wide,Emax2(:,x_next+10),A_next);
                                    
                        % polynomial approx
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vm_r_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vm_r_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        Vm_r_next_2 = sum(coeff(x_next+10,:).*Base,2); %cheby_approx
                        Vm_r_next2_2 = sum(coeff2(x_next+10,:).*Base,2); %cheby_approx
                      
                        % save approximations
                        Amr_next(k)=A_next;
                        Vmr_next(k,x)=Vm_r_next;
                        Vmr_next2(k,x)=Vm_r_next2;
                        Vmr_next_2(k,x)=Vm_r_next_2;
                        Vmr_next2_2(k,x)=Vm_r_next2_2;
                        Vmr_next_linear(k,x)=Vm_r_next_linear;
                        Vmr_next_linear2(k,x)=Vm_r_next_linear2;
                        Vmr_next_linear_2(k,x)=Vm_r_next_linear_2;
                        Vmr_next_linear2_2(k,x)=Vm_r_next_linear2_2;
                                            
                        % non-regular job
                        A_next = (1+G.r) * (A_j + (w_j_n + exp(wh)*m_j + shock_i) - chh_n - n_j*Inv);
                        x_next = x + 1;
                        if x_next == 11
                            x_next = 10;
                        end

                        % linear approximation of VF
                        Vm_n_next_linear = interpn(A_wide,Emax(:,x_next),A_next);
                        Vm_n_next_linear2 = interpn(A_wide,Emax2(:,x_next),A_next);
                        Vm_n_next_linear_2 = interpn(A_wide,Emax(:,x_next+10),A_next);
                        Vm_n_next_linear2_2 = interpn(A_wide,Emax2(:,x_next+10),A_next); 
                        
                        % polynomial approx
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vm_n_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vm_n_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        Vm_n_next_2 = sum(coeff(x_next+10,:).*Base,2); %cheby_approx
                        Vm_n_next2_2 = sum(coeff2(x_next+10,:).*Base,2); %cheby_approx
                        
                        % save approx
                        Amn_next(k)=A_next;
                        Vmn_next(k,x)=Vm_n_next;
                        Vmn_next2(k,x)=Vm_n_next2;
                        Vmn_next_2(k,x)=Vm_n_next_2;
                        Vmn_next2_2(k,x)=Vm_n_next2_2;
                        Vmn_next_linear(k,x)=Vm_n_next_linear;
                        Vmn_next_linear2(k,x)=Vm_n_next_linear2;
                        Vmn_next_linear_2(k,x)=Vm_n_next_linear_2;
                        Vmn_next_linear2_2(k,x)=Vm_n_next_linear2_2;
                        
                        % unemployed
                        A_next = (1+G.r) * (A_j + (w_j_u + exp(wh)*m_j + shock_i) - chh_u - n_j*Inv);
                        x_next = x;

                        % linear approximation of VF
                        Vm_u_next_linear = interpn(A_wide,Emax(:,x_next),A_next); 
                        Vm_u_next_linear2 = interpn(A_wide,Emax2(:,x_next),A_next);
                        Vm_u_next_linear_2 = interpn(A_wide,Emax(:,x_next+10),A_next); 
                        Vm_u_next_linear2_2 = interpn(A_wide,Emax2(:,x_next+10),A_next);   
                        
                        % polynomial approx
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vm_u_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vm_u_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        Vm_u_next_2 = sum(coeff(x_next+10,:).*Base,2); %cheby_approx
                        Vm_u_next2_2 = sum(coeff2(x_next+10,:).*Base,2); %cheby_approx

                        % save output
                        Amu_next(k)=A_next;
                        Vmu_next(k,x)=Vm_u_next;
                        Vmu_next2(k,x)=Vm_u_next2;
                        Vmu_next_2(k,x)=Vm_u_next_2;
                        Vmu_next2_2(k,x)=Vm_u_next2_2;
                        Vmu_next_linear(k,x)=Vm_u_next_linear;
                        Vmu_next_linear2(k,x)=Vm_u_next_linear2;
                        Vmu_next_linear_2(k,x)=Vm_u_next_linear_2;
                        Vmu_next_linear2_2(k,x)=Vm_u_next_linear2_2;
                        
                        % Sector-Specific Value Functions (1 child) 
%                         Vm1_r(k) = u_r(k) + G.beta * ((prob_lamba*Vm_r_next_linear)+(1-prob_lamba)*Vm_r_next_linear2);
%                         Vm1_n(k) = u_n(k) + G.beta * ((prob_pi*Vm_n_next_linear)+(1-prob_pi)*Vm_n_next_linear2);
%                         Vm1_u(k) = u_u(k) + G.beta * ((prob_pi*Vm_u_next_linear)+(1-prob_pi)*Vm_u_next_linear2);    
                          Vm1_r(k) = u_r(k) + G.beta * ((prob_lamba*Vm_r_next)+(1-prob_lamba)*Vm_r_next2);
                          Vm1_n(k) = u_n(k) + G.beta * ((prob_pi*Vm_n_next)+(1-prob_pi)*Vm_n_next2);
                          Vm1_u(k) = u_u(k) + G.beta * ((prob_pi*Vm_u_next)+(1-prob_pi)*Vm_u_next2);
                        
                        % Sector-Specific Value Functions (2 child, maybe)
%                         Vm1_r_2(k) = u_r(k) + G.beta * ((prob_lamba*Vm_r_next_linear_2)+(1-prob_lamba)*Vm_r_next_linear2_2);
%                         Vm1_n_2(k) = u_n(k) + G.beta * ((prob_pi*Vm_n_next_linear_2)+(1-prob_pi)*Vm_n_next_linear2_2);
%                         Vm1_u_2(k) = u_u(k) + G.beta * ((prob_pi*Vm_u_next_linear_2)+(1-prob_pi)*Vm_u_next_linear2_2);
                          Vm1_r_2(k) = u_r(k) + G.beta * ((prob_lamba*Vm_r_next_2)+(1-prob_lamba)*Vm_r_next2_2);
                          Vm1_n_2(k) = u_n(k) + G.beta * ((prob_pi*Vm_n_next_2)+(1-prob_pi)*Vm_n_next2_2);
                          Vm1_u_2(k) = u_u(k) + G.beta * ((prob_pi*Vm_u_next_2)+(1-prob_pi)*Vm_u_next2_2);

                        % Sector-Specific Value Fucntions (married)
                          Vm_r(k) = Vm1_r(k)*(1-prob_2kids_r) + Vm1_r_2(k)*(prob_2kids_r);
                          Vm_n(k) = Vm1_n(k)*(1-prob_2kids_n) + Vm1_n_2(k)*(prob_2kids_n);
                          Vm_u(k) = Vm1_u(k)*(1-prob_2kids_u) + Vm1_u_2(k)*(prob_2kids_u);
                        
                        % save marriage 1 child (for marriage decision)
                        Vm_r_aux(k,x) = Vm1_r(k);
                        Vm_n_aux(k,x) = Vm1_n(k);
                        Vm_u_aux(k,x) = Vm1_u(k);
                        
                    % Married with 2 kids: 
                    elseif x <= 20
                        
                        % regular job
                        A_next = (1+G.r) * (A_j + (w_j_r + exp(wh)*m_j + shock_i) - chh_r - n_j*Inv); % eq. 8
                        x_next = x + 1;
                        if x_next == 21
                            x_next = 20;
                        end

                        % linear approximation of VF
                        Vm2_r_next_linear = interpn(A_wide,Emax(:,x_next),A_next);
                        Vm2_r_next_linear2 = interpn(A_wide,Emax2(:,x_next),A_next);
                        
                        % Chebyshev Approximation
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vm2_r_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vm2_r_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx

                        % save output
                        Am2r_next(k)=A_next;
                        Vmr_next(k,x)=Vm2_r_next;
                        Vmr_next2(k,x)=Vm2_r_next2;
                        Vmr_next_linear(k,x)=Vm2_r_next_linear;
                        Vmr_next_linear2(k,x)=Vm2_r_next_linear2;
                          
                        % non-regular job
                        A_next = (1+G.r) * (A_j + (w_j_n + exp(wh)*m_j + shock_i) - chh_n - n_j*Inv);
                        x_next = x + 1;
                        if x_next == 21
                            x_next = 20;
                        end

                        % linear approximation of VF
                        Vm2_n_next_linear = interpn(A_wide,Emax(:,x_next),A_next); 
                        Vm2_n_next_linear2 = interpn(A_wide,Emax2(:,x_next),A_next);
                        
                        % Chebyshev Approximation
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vm2_n_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vm2_n_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx

                        % save approx
                        Am2n_next(k)=A_next;
                        Vmn_next(k,x)=Vm2_n_next;
                        Vmn_next2(k,x)=Vm2_n_next2;
                        Vmn_next_linear(k,x)=Vm2_n_next_linear;
                        Vmn_next_linear2(k,x)=Vm2_n_next_linear2;

                        % unemployed
                        A_next = (1+G.r) * (A_j + (w_j_u + exp(wh)*m_j + shock_i) - chh_u - n_j*Inv);
                        x_next = x;

                        % linear approximation of VF
                        Vm2_u_next_linear = interpn(A_wide,Emax(:,x_next),A_next); 
                        Vm2_u_next_linear2 = interpn(A_wide,Emax2(:,x_next),A_next);

                        % Chebyshev Approximation
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vm2_u_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vm2_u_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx

                        % save output
                        Am2u_next(k)=A_next;
                        Vmu_next(k,x)=Vm2_u_next;
                        Vmu_next2(k,x)=Vm2_u_next2;
                        Vmu_next_linear(k,x)=Vm2_u_next_linear;
                        Vmu_next_linear2(k,x)=Vm2_u_next_linear2;

                        % Sector-Specific Value Functions (2 child)
                        Vm2_r(k) = u_r(k) + G.beta * ((prob_lamba*Vm2_r_next)+(1-prob_lamba)*Vm2_r_next2);
                        Vm2_n(k) = u_n(k) + G.beta * ((prob_pi*Vm2_n_next)+(1-prob_pi)*Vm2_n_next2);
                        Vm2_u(k) = u_u(k) + G.beta * ((prob_pi*Vm2_u_next)+(1-prob_pi)*Vm2_u_next2);
                        
                    % Single
                    else
                        
                        % Regular
                        A_next = (1+G.r) * (A_j + (w_j_r + exp(wh)*m_j + shock_i) - chh_r - n_j*Inv);
                        x_next = x + 1;
                        if x_next == 31
                            x_next = 30;
                        end

                        % linear approximation of VF
                        Vs_r_next_linear = interpn(A_wide,Emax(:,x_next),A_next);
                        Vs_r_next_linear2 =interpn(A_wide,Emax2(:,x_next),A_next);
                      
                        % polynomial approx
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vs_r_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vs_r_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        
                        % save output
                        Asr_next(k)=A_next;
                        Vsr_next(k,x)=Vs_r_next;
                        Vsr_next2(k,x)=Vs_r_next2;
                        Vsr_next_linear(k,x)=Vs_r_next_linear;
                        Vsr_next_linear2(k,x)=Vs_r_next_linear2;

                        % Non-regular
                        A_next = (1+G.r) * (A_j + (w_j_n + exp(wh)*m_j + shock_i) - chh_n - n_j*Inv);
                        x_next = x + 1;
                        if x_next == 31
                            x_next = 30;
                        end

                        % linear approximation of VF
                        Vs_n_next_linear = interpn(A_wide,Emax(:,x_next),A_next);
                        Vs_n_next_linear2 = interpn(A_wide,Emax2(:,x_next),A_next);

                        % polynomial approx
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vs_n_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vs_n_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        
                        % save output
                        Asn_next(k)=A_next;
                        Vsn_next(k,x)=Vs_n_next;
                        Vsn_next2(k,x)=Vs_n_next2;
                        Vsn_next_linear(k,x)=Vs_n_next_linear;
                        Vsn_next_linear2(k,x)=Vs_n_next_linear2;

                        % Unemployed
                        A_next = (1+G.r) * (A_j + (w_j_u + exp(wh)*m_j + shock_i) - chh_u - n_j*Inv);
                        x_next = x;

                        % linear approximation of VF
                        Vs_u_next_linear = interpn(A_wide,Emax(:,x_next),A_next);
                        Vs_u_next_linear2 = interpn(A_wide,Emax2(:,x_next),A_next);
                        
                        % polynomial approx
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vs_u_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vs_u_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx

                        % save output
                        Asu_next(k)=A_next;
                        Vsu_next(k,x)=Vs_u_next;
                        Vsu_next2(k,x)=Vs_u_next2;
                        Vsu_next_linear(k,x)=Vs_u_next_linear;
                        Vsu_next_linear2(k,x)=Vs_u_next_linear2;

                        % Sector-Specific Value Functions (single)
%                         Vs_r(k) = u_r(k) + G.beta * ((prob_lamba*Vs_r_next_linear)+(1-prob_lamba)*Vs_r_next_linear2);
%                         Vs_n(k) = u_n(k) + G.beta * ((prob_pi*Vs_n_next_linear)+(1-prob_pi)*Vs_n_next_linear2);
%                         Vs_u(k) = u_u(k) + G.beta * ((prob_pi*Vs_u_next_linear)+(1-prob_pi)*Vs_u_next_linear2);
                        Vs_r(k) = u_r(k) + G.beta * ((prob_lamba*Vs_r_next)+(1-prob_lamba)*Vs_r_next2);
                        Vs_n(k) = u_n(k) + G.beta * ((prob_pi*Vs_n_next)+(1-prob_pi)*Vs_n_next2);
                        Vs_u(k) = u_u(k) + G.beta * ((prob_pi*Vs_u_next)+(1-prob_pi)*Vs_u_next2);
                        
                        % Sector-Specific Value Functions (marriage, maybe)
                        Vsm_r(k) = prob_marr_w*Vm_r_aux(k,x-20) + (1-prob_marr_w)*Vs_r(k);
                        Vsm_n(k) = prob_marr_w*Vm_n_aux(k,x-20) + (1-prob_marr_w)*Vs_n(k);
                        Vsm_u(k) = prob_marr_u*Vm_u_aux(k,x-20) + (1-prob_marr_u)*Vs_u(k);
                    end
                end
                        
                % optimal consumption and max VF
                if x <= 10
                    % check
                    Vm_r(Amr_next < min(A_wide)) = NaN;
                    Vm_r(Amr_next > max(A_wide)) = NaN;
                    Vm_n(Amn_next < min(A_wide)) = NaN;
                    Vm_n(Amn_next > max(A_wide)) = NaN;
                    Vm_u(Amu_next < min(A_wide)) = NaN;
                    Vm_u(Amu_next > max(A_wide)) = NaN;
                    % save optimal
                    [Vm_r_star, index_mr_k] = max(Vm_r);
                    [Vm_n_star, index_mn_k] = max(Vm_n);
                    [Vm_u_star, index_mu_k] = max(Vm_u);
                    cm_r_star = cr_vector(index_mr_k);
                    cm_n_star = cn_vector(index_mn_k);
                    cm_u_star = cu_vector(index_mu_k);
                    cm_star_aux = [cm_r_star, cm_n_star, cm_u_star];
                    [Vm_star, lm_index] = max([Vm_r_star,Vm_n_star,Vm_u_star]); % 3 job options
                    [Vm_star2] = max([Vm_n_star,Vm_u_star]); % emax2 (2 job options)
                elseif x <= 20
                    % check
                    Vm2_r(Am2r_next < min(A_wide)) = NaN;
                    Vm2_r(Am2r_next > max(A_wide)) = NaN;
                    Vm2_n(Am2n_next < min(A_wide)) = NaN;
                    Vm2_n(Am2n_next > max(A_wide)) = NaN;
                    Vm2_u(Am2u_next < min(A_wide)) = NaN;
                    Vm2_u(Am2u_next > max(A_wide)) = NaN;
                    % save optimal
                    [Vm2_r_star, index_m2r_k] = max(Vm_r);
                    [Vm2_n_star, index_m2n_k] = max(Vm_n);
                    [Vm2_u_star, index_m2u_k] = max(Vm_u);
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
                    [Vs_star, ls_index] = max([Vsm_r_star, Vsm_n_star, Vsm_u_star, Vs_r_star, Vs_n_star, Vs_u_star]);
                    [Vs_star2] = max([Vsm_n_star, Vsm_u_star, Vs_n_star, Vs_u_star]);
                end
            
            % save choice:
            if x <= 10
                c_star(j, i, x, t) = cm_star_aux(lm_index);
                l_star(j, i, x, t) = lm_index;
                V_star(j, i, x, t) = Vm_star;
                V2_star(j,i, x, t) = Vm_star2;
            elseif x <= 20
                c_star(j, i, x, t) = cm2_star_aux(lm2_index);
                l_star(j, i, x, t) = lm2_index;
                V_star(j, i, x, t) = Vm2_star;
                V2_star(j,i, x, t) = Vm2_star2;
            else
                c_star(j, i, x, t) = cs_star_aux(ls_index);
                l_star(j, i, x, t) = ls_index;
                V_star(j, i, x, t) = Vs_star;
                V2_star(j,i, x, t) = Vs_star2;
            end
            
            % save the number assets outside grid
            if x <= 10
                Ar_out(j,i,x,t) = sum(Amr_next < A_min) + sum(Amr_next > max(A_wide));
                An_out(j,i,x,t) = sum(Amn_next < A_min) + sum(Amn_next > max(A_wide));
                Au_out(j,i,x,t) = sum(Amu_next < A_min) + sum(Amu_next > max(A_wide));
            elseif x <= 20
                Ar_out(j,i,x,t) = sum(Am2r_next < A_min) + sum(Am2r_next > max(A_wide));
                An_out(j,i,x,t) = sum(Am2n_next < A_min) + sum(Am2n_next > max(A_wide));
                Au_out(j,i,x,t) = sum(Am2u_next < A_min) + sum(Am2u_next > max(A_wide));
            else
                Ar_out(j,i,x,t) = sum(Asr_next < A_min) + sum(Asr_next > max(A_wide));
                An_out(j,i,x,t) = sum(Asn_next < A_min) + sum(Asn_next > max(A_wide));
                Au_out(j,i,x,t) = sum(Asu_next < A_min) + sum(Asu_next > max(A_wide));
            end
            
            %% output to compare
            Amr_next_aux(:,x,j,t)=Amr_next';
            Amn_next_aux(:,x,j,t)=Amn_next';
            Amu_next_aux(:,x,j,t)=Amu_next';
            Vu_next_linear2(:,:,j,t) = [Vmu_next_linear2,Vsu_next_linear2(:,21:30)];
            Vn_next_linear2(:,:,j,t) = [Vmn_next_linear2,Vsn_next_linear2(:,21:30)];
            Vr_next_linear2(:,:,j,t) = [Vmr_next_linear2,Vsr_next_linear2(:,21:30)];
            Vu_next_linear(:,:,j,t) = [Vmu_next_linear,Vsu_next_linear(:,21:30)];
            Vn_next_linear(:,:,j,t) = [Vmn_next_linear,Vsn_next_linear(:,21:30)];
            Vr_next_linear(:,:,j,t) = [Vmr_next_linear,Vsr_next_linear(:,21:30)];
            Vu_next2(:,:,j,t) = [Vmu_next2,Vsu_next2(:,21:30)];
            Vn_next2(:,:,j,t) = [Vmn_next2,Vsn_next2(:,21:30)];
            Vr_next2(:,:,j,t) = [Vmr_next2,Vsr_next2(:,21:30)];
            Vu_next(:,:,j,t) = [Vmu_next,Vsu_next(:,21:30)];
            Vn_next(:,:,j,t) = [Vmn_next,Vsn_next(:,21:30)];
            Vr_next(:,:,j,t) = [Vmr_next,Vsr_next(:,21:30)];
            end
        end
        
        % Integrate over shocks
        W(:,x,t) = pi^(-1/2)*V_star(:,:,x,t)*S.weight;
        W2(:,x,t)= pi^(-1/2)*V2_star(:,:,x,t)*S.weight;
        
        % reshape policy func
        c_func(:,:,:,x,t) = reshape(c_star(:,:,x,t), [G.n_assets,3,3]);
        l_func(:,:,:,x,t) = reshape(l_star(:,:,x,t), [G.n_assets,3,3]);
    end
end

% three labor functions (as 0 or 1)
lr_func = l_func == 1 | l_func == 4;
ln_func = l_func == 2 | l_func == 5;
lu_func = l_func == 3 | l_func == 6;

% marriage function for single women: equals 1 if labor choice is 1, 2, or 3
m_func = l_func == 1 | l_func == 2 | l_func == 3;
