function [v_func_rsp,c_func,c_func_rsp,m_func,m_func_rsp,lr_func,lr_func_rsp,ln_func,ln_func_rsp,lu_func,lu_func_rsp] = solution_linear(G,abi,edu,S,params)

% Index for parameters

    psi_r=params(1);
    psi_n=params(2);
    theta1_r=params(3);
    theta1_n=params(4);
    theta1_u=params(5);
    theta3_r=params(6);
    theta3_n=params(7);
    theta3_u=params(8);
    gamma1=params(9);
    phi=params(10);
    alpha01_r=params(11);
    alpha01_n=params(12);
    alpha02_r=params(13);
    alpha02_n=params(14);
    alpha11_n=params(15);
    alpha11_r=params(16);
    alpha12_n=params(17);
    alpha12_r=params(18);
    alpha2_r=params(19);
    alpha2_n=params(20);
    sigma_r = params(21);
    sigma_n = params(22);
    sigma_i = params(23); 
    omega0_w = params(24); 
    omega0_u = params(25); 
    omega11  = params(26);
    omega12  = params(27);
    omega2 = params(28);
    lambda1 = params(29);
    lambda2 = params(30);
    lambda3 = params(31);
    eta01 = params(32);
    eta02 = params(33);
    eta11 = params(34);
    eta12 = params(35);
    eta2 = params(36);
    eta03 = params(37);
    eta04 = params(38);
    eta21 = params(39);
    eta22 = params(40);
    eta3 = params(41);
    kappa01 = params(42);
    kappa02 = params(43);
    kappa11 = params(44);
    kappa12 = params(45);
    kappa2 = params(46);
    kappa03 = params(47);
    kappa04 = params(48);
    kappa21 = params(49);
    kappa22 = params(50);
    kappa3 = params(51);
 
% Terminal Value Function: TVF = A_T + W_T + Q_T = assets + wages + HH_prod
age_TVF = 18*(edu==1) + 20*(edu==2) + 22*(edu==3) + G.n_period;

    % assets
    assets = S.SS_A; 
    
    % husband wages
    wh_mean = eta01 + eta02*(abi==2) + eta11*(edu==2) + eta12*(edu==3) + eta2*age_TVF;
    wh_sd = eta03 + eta04*(abi==2) + eta21*(edu==2) + eta22*(edu==3) + eta3*age_TVF;
    wh = normrnd(wh_mean,wh_sd);
    
    % child human capital
    Inv_mean = kappa01 + kappa02*(abi==2) + kappa11*(edu==2) + kappa12*(edu==3) + kappa2*age_TVF;
    Inv_sd = kappa03 + kappa04*(abi==2) + kappa21*(edu==2) + kappa22*(edu==3) + kappa3*age_TVF;
    Inv = normrnd(Inv_mean,Inv_sd);
    K_j = 5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% what is this?
    K = (gamma1*K_j^phi + (1-gamma1)*Inv^phi)^(1/phi);

TVF = lambda1*assets + lambda2*unique(S.SS_X) + lambda3*wh + lambda3*K;

tic    
% loop for time (20):
for t = G.n_period-1:-1:1
    t
    toc
   
    age = 18*(edu==1) + 20*(edu==2) + 22*(edu==3) + t;
    
    % draw husband wage
    wh_mean = eta01 + eta02*(abi==2) + eta11*(edu==2) + eta12*(edu==3) + eta2*age; %function of women's age, education and ability types
    wh_sd = eta03 + eta04*(abi==2) + eta21*(edu==2) + eta22*(edu==3) + eta3*age; %same as above
    wh_next = normrnd(wh_mean,wh_sd);
    
    % draw investments
    Inv_mean = kappa01 + kappa02*(abi==2) + kappa11*(edu==2) + kappa12*(edu==3) + kappa2*age;
    Inv_sd = kappa03 + kappa04*(abi==2) + kappa21*(edu==2) + kappa22*(edu==3) + kappa3*age;
    Inv = normrnd(Inv_mean,Inv_sd);
    K_next = (gamma1*K_j^phi + (1-gamma1)*G.Inv^phi)^(1/phi); 
    
    % marriage probabilities: should be function of age? - take out of assets
    prob_marr_w = omega0_w + omega11*(edu==2) + omega12*(edu==3) + omega2*t; %or switch it to age?
    prob_marr_u = omega0_u + omega11*(edu==2) + omega12*(edu==3) + omega2*t; %or switch it to age?
    
    % terminal value/Emax
    if t==G.n_period-1
        Emax = TVF;
        %Emax2 = TVF;
    else
        Emax = W(:,:,t+1);
        %Emax2 = W2(:,;,t+1);
    end
    
    % reshape value function?
%     for x = 1:1:(G.n_matstat*G.n_wrkexp)
%         Emax_rsp(:,:,:,x) = reshape(Emax(:,x),[G.n_childK,G.n_assets,G.n_hwages]);
%         %Emax2_rsp
%     end
    
    % loop for work experience and marital status (20): %%%%% change to 30?
    for x = 1:1:(G.n_matstat*G.n_wrkexp)
        x;
        
        % current state discrete variables:
        m_j = S.SS_M(x);  % marital status
        n_j = S.SS_N(x);  % children
        X_j = S.SS_X(x);  % experience
    
        % loop for shocks (9):
        for i = 1:1:G.n_shocks
            i;
        
            shock_i = S.shocks_i(i);
            shock_r = S.shocks_r(i);
            shock_n = S.shocks_n(i);
            
            % sector-specific state variables: - take out of assets
            w_j_r = exp(alpha01_r + alpha02_r*(abi==2) + alpha11_r*(edu==2) + alpha12_r*(edu==3) + alpha2_r*log(1+X_j) + shock_r); % same
            w_j_n = exp(alpha01_n + alpha02_n*(abi==2) + alpha11_n*(edu==2) + alpha12_n*(edu==3) + alpha2_n*log(1+X_j) + shock_n); % same  
            w_j_u = 0; % unemployed women don't have earnings
              
            % loop over assets
            for j = 1:1:G.n_SS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change to 10?
                j;
            
                % HH's assets
                A_j = S.SS_A(j); 
            
                % consumption vector
                chh_min = 0.1;
                chh_max = max(0,A_j);
                cr_vector = linspace(chh_min,chh_max,G.n_cons);                
                cn_vector=cr_vector;
                cu_vector=cr_vector;

                % wide assets vector for linear interpolation (move up)
                A_wide = linspace(-10,max(S.SS_A),length(unique(S.SS_A)));
                
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
                    
                    % Sector-Specific Utility:
                    u_r(k) = (cw_r^(1-G.sigma))/(1-G.sigma) + psi_r + theta1_r*log(1+m_j) + theta3_r*log(K)*n_j;
                    u_n(k) = (cw_n^(1-G.sigma))/(1-G.sigma) + psi_n + theta1_n*log(1+m_j) + theta3_n*log(K)*n_j;
                    u_u(k) = (cw_u^(1-G.sigma))/(1-G.sigma) + theta1_u*log(1+m_j) + theta3_u*log(K)*n_j;
                    
                    % married AND 1 KID:
                    if x <= 10
                    
                        % regular job:
                        A_next = (1+G.r) * (A_j + (w_j_r + wh*m_j + shock_i) - chh_r - n_j*Inv); % eq. 8
                        x_next = x + 1;
                        if x_next == 11
                            x_next = 10;
                        end
                        % value function:
                        Vm_r_next_linear = interpn(A_wide,Emax,A_next); %interpn(unique(S.SS_K),A_wide,unique(S.SS_H),Emax_rsp(:,:,:,x),K_next,A_next,wh_next);
%%%%%%%%%%%%%%%%%%%%%%%%Vm_r_next_linear2 = interpn(unique(S.SS_K),A_wide,unique(S.SS_H),Emax2_rsp(:,:,:,x),K_next,A_next,wh_next);                        
                        Amr_next(k)=A_next;
                        Vmr_next_linear(k,x)=Vm_r_next_linear;
                        % non-regular job:
                        A_next = (1+G.r) * (A_j + (w_j_n + wh*m_j + shock_i) - chh_n - n_j*Inv);
                        x_next = x + 1;
                        if x_next == 11
                            x_next = 10;
                        end
                        % value function:
                        Vm_n_next_linear = interpn(A_wide,Emax,A_next); %interpn(unique(S.SS_K),A_wide,unique(S.SS_H),Emax_rsp(:,:,:,x),K_next,A_next,wh_next);
%%%%%%%%%%%%%%%%%%%%%%%%Vm_n_next_linear2 = interpn(unique(S.SS_K),A_wide,unique(S.SS_H),Emax2_rsp(:,:,:,x),K_next,A_next,wh_next);                                            
                        Amn_next(k)=A_next;
                        Vmn_next_linear(k,x)=Vm_n_next_linear;
                        % unemployed:
                        A_next = (1+G.r) * (A_j + (w_j_u + wh*m_j + shock_i) - chh_u - n_j*Inv);
                        x_next = x;
                        % value function:
                        Vm_u_next_linear = interpn(A_wide,Emax,A_next); %interpn(unique(S.SS_K),A_wide,unique(S.SS_H),Emax_rsp(:,:,:,x),K_next,A_next,wh_next);
%%%%%%%%%%%%%%%%%%%%%%%%Vm_u_next_linear2 = interpn(unique(S.SS_K),A_wide,unique(S.SS_H),Emax2_rsp(:,:,:,x),K_next,A_next,wh_next);                                                
                        Amu_next(k)=A_next;
                        Vmu_next_linear(k,x)=Vm_u_next_linear;
                        % Sector-Specific Value Functions
                        Vm_r(k) = u_r(k) + G.beta * Vm_r_next_linear; %lambda*V1 + (1-lamba*V2)
                        Vm_n(k) = u_n(k) + G.beta * Vm_n_next_linear; %lambda*V1 + (1-lamba*V2)
                        Vm_u(k) = u_u(k) + G.beta * Vm_u_next_linear; %lambda*V1 + (1-lamba*V2)
                        % save marriage values (for marriage decision)
                        Vm_r_aux(k,x) = Vm_r(k);
                        Vm_n_aux(k,x) = Vm_n(k);
                        Vm_u_aux(k,x) = Vm_u(k);
                    % Married with 2 kids    
                    % Single:
                    else
                        
                        % Regular:
                        A_next = (1+G.r) * (A_j + (w_j_r + wh*m_j + shock_i) - chh_r - n_j*Inv);
                        x_next = x + 1;
                        if x_next == 21
                            x_next = 20;
                        end
                        Vs_r_next_linear = interpn(A_wide,Emax,A_next); %interpn(unique(S.SS_K),A_wide,unique(S.SS_H),Emax_rsp(:,:,:,x),K_next,A_next,wh_next);
                        %Vs_r_next_linear2 (Emax2)
                        Asr_next(k)=A_next;
                        Vsr_next_linear(k,x)=Vs_r_next_linear;
                        % Non-regular:
                        A_next = (1+G.r) * (A_j + (w_j_n + wh*m_j + shock_i) - chh_n - n_j*Inv);
                        x_next = x + 1;
                        if x_next == 21
                            x_next = 20;
                        end
                        Vs_n_next_linear = interpn(A_wide,Emax,A_next); %interpn(unique(S.SS_K),A_wide,unique(S.SS_H),Emax_rsp(:,:,:,x),K_next,A_next,wh_next);
                        %Vs_n_next_linear2 (Emax2)
                        Asn_next(k)=A_next;
                        Vsn_next_linear(k,x)=Vs_n_next_linear;
                        % Unemployed:
                        A_next = (1+G.r) * (A_j + (w_j_u + wh*m_j + shock_i) - chh_u - n_j*Inv);
                        x_next = x;
                        Vs_u_next_linear = interpn(A_wide,Emax,A_next); %interpn(unique(S.SS_K),A_wide,unique(S.SS_H),Emax_rsp(:,:,:,x),K_next,A_next,wh_next);
                        %Vs_r_next_linear2 (Emax2)
                        Asu_next(k)=A_next;
                        Vsu_next_linear(k,x)=Vs_u_next_linear;
                        % Sector-Specific Value Functions
                        Vs_r(k) = u_r(k) + G.beta * Vs_r_next_linear; %pi*Emax1 + (1-pi)*Emax2
                        Vs_n(k) = u_n(k) + G.beta * Vs_n_next_linear; %pi*Emax1 + (1-pi)*Emax2
                        Vs_u(k) = u_u(k) + G.beta * Vs_u_next_linear; %pi*Emax1 + (1-pi)*Emax2
                        Vsm_r(k) = prob_marr_w*Vm_r_aux(k,x-10) + (1-prob_marr_w)*Vs_r(k); %Vm_r_aux CHANGE To Vm_r_expected: prob*E[m,1)+(1-prob)*E(m,2)
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
                    [Vm_star, Index_lm] = max([Vm_r_star,Vm_n_star,Vm_u_star]);
                else
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
                end
            % save choice:
            if x <= 10
                c_star(j, i, x, t) = cm_star_aux(Index_lm);
                l_star(j, i, x, t) = Index_lm;
                V_star(j, i, x, t) = Vm_star;
            else
                c_star(j, i, x, t) = cs_star_aux(Index_ls);
                l_star(j, i, x, t) = Index_ls;
                V_star(j, i, x, t) = Vs_star;
            end
            
            % save the number assets outside grid
            if x <= 10
                Ar_out(j,i,x,t) = sum(Amr_next < S.extmin_A) + sum(Amr_next > S.extmax_A);
                An_out(j,i,x,t) = sum(Amn_next < S.extmin_A) + sum(Amn_next > S.extmax_A);
                Au_out(j,i,x,t) = sum(Amu_next < S.extmin_A) + sum(Amu_next > S.extmax_A);
            else
                Ar_out(j,i,x,t) = sum(Asr_next < S.extmin_A) + sum(Asr_next > S.extmax_A);
                An_out(j,i,x,t) = sum(Asn_next < S.extmin_A) + sum(Asn_next > S.extmax_A);
                Au_out(j,i,x,t) = sum(Asu_next < S.extmin_A) + sum(Asu_next > S.extmax_A);
            end
            end
        end
        
        % Integrate over shocks
        W(:,x,t) = pi^(-1/2)*V_star(:,:,x,t)*S.weight;
        %W2(:,x,t) = pi^(-1/2)*V2_start(:,:,x,t)*S.weight; 
        % reshape policy func
        v_func(:,x,t) = reshape(V_star(:,:,x,t),[],1);
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

%% reshaped policy functions

for t=1:1:G.n_period-1
    for x=1:1:G.n_matstat*G.n_wrkexp
        v_func_rsp(:,:,:,:,:,x,t) = reshape(v_func(:,x,t), [G.n_childK,G.n_assets,G.n_hwages,3,3]);
        c_func_rsp(:,:,:,:,:,x,t) = reshape(c_func(:,x,t), [G.n_childK,G.n_assets,G.n_hwages,3,3]);
        c_func_rsp9(:,:,:,:,x,t) = reshape(c_func(:,x,t),[G.n_childK,G.n_assets,G.n_hwages,G.n_shocks]);
        m_func_rsp(:,:,:,:,:,x,t) = reshape(m_func(:,x,t), [G.n_childK,G.n_assets,G.n_hwages,3,3]);
        m_func_rsp9(:,:,:,:,x,t) = reshape(m_func(:,x,t),[G.n_childK,G.n_assets,G.n_hwages,G.n_shocks]);
        lr_func_rsp(:,:,:,:,:,x,t) = reshape(lr_func(:,x,t), [G.n_childK,G.n_assets,G.n_hwages,3,3]);
        lr_func_rsp9(:,:,:,:,x,t) = reshape(lr_func(:,x,t),[G.n_childK,G.n_assets,G.n_hwages,G.n_shocks]);
        ln_func_rsp(:,:,:,:,:,x,t) = reshape(ln_func(:,x,t), [G.n_childK,G.n_assets,G.n_hwages,3,3]);
        ln_func_rsp9(:,:,:,:,x,t) = reshape(ln_func(:,x,t),[G.n_childK,G.n_assets,G.n_hwages,G.n_shocks]);
        lu_func_rsp(:,:,:,:,:,x,t) = reshape(lu_func(:,x,t), [G.n_childK,G.n_assets,G.n_hwages,3,3]);
        lu_func_rsp9(:,:,:,:,x,t) = reshape(lu_func(:,x,t),[G.n_childK,G.n_assets,G.n_hwages,G.n_shocks]);
    end
end

% % Save in a structure
% 
% Sol = struct(...
%     'c_func',c_func,'c_func_rsp',c_func_rsp,'c_func_rsp9',c_func_rsp9,...
%     'm_func',m_func,'m_func_rsp',m_func_rsp,'m_func_rsp9',m_func_rsp9,...
%     'lr_func',lr_func,'lr_func_rsp',lr_func_rsp,'lr_func_rsp9',lr_func_rsp9,...
%     'ln_func',ln_func,'ln_func_rsp',ln_func_rsp,'ln_func_rsp9',ln_func_rsp9,...
%     'lu_func',lu_func,'lu_func_rsp',lu_func_rsp,'lu_func_rsp9',lu_func_rsp9);
    