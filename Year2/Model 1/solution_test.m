function [c_func,m_func,lr_func,ln_func,lu_func,Ar_out,An_out,Au_out] = solution(G,S,P)

%% index for parameters:

% estimated
psi_r=params(1);
psi_n=params(2);
psi_u=params(3);
theta1_r=0; %params(4);
theta1_n=0; %params(5);
theta1_u=0; %params(6);
theta2_r=0; %params(7);
theta2_n=0; %params(8);
theta2_u=0; %params(9);
theta3_r=0; %params(10);
theta3_n=0; %params(11);
theta3_u=0; %params(12);
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
lambda1=10; %params(26);
lambda2=0; %params(27);
lambda3=0; %params(28);
lambda4=0; %params(29);

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
chi01=params(95);
chi02=params(96);
chi11=params(97);
chi12=params(98);
chi03=params(99);
chi04=params(100);
chi21=params(101);
chi22=params(102);

%% other parameters:

chh_min = 0.1; % minimun consumption
delta = 0.5; % Female Share of Consumption (CAL)

% expanded assets vector for linear interpolation
A_min=S.extmin_A;
A_max=S.extmax_A;
A_wide = S.SS_A;
A_wide(1) = A_min;
A_wide(10) = A_max;

%% Terminal Value Function:

% % age at TVF
% age_TVF = 18*(edu==1) + 20*(edu==2) + 22*(edu==3) + G.n_period;

% assets
assets = S.SS_A; 

% % husband wages
% wh_mean = eta01 + eta02*(abi==2) + eta11*(edu==2) + eta12*(edu==3) + eta2*age_TVF;
% wh_sd = 0.7022257; %eta03 + eta04*(abi==2) + eta21*(edu==2) + eta22*(edu==3) + eta3*age_TVF;
% wh_TVF = normrnd(wh_mean,wh_sd);
% 
% % child human capital 
% Inv_mean = iota01 + iota02*(abi==2) + iota11*(edu==2) + iota12*(edu==3) + iota2*age_TVF;
% Inv_sd = 0.9270494; %iota03 + iota04*(abi==2) + iota21*(edu==2) + iota22*(edu==3) + iota3*age_TVF;
% Inv_TVF = normrnd(Inv_mean,Inv_sd);
% 
% K_TVF = exp(kappa01 + kappa02*(abi==2) + kappa03*(edu==2) + kappa04*(edu==3) + kappa05*(Inv_TVF));

% TVF 
TVF = repmat(real(lambda1*(assets).^(1-G.sigma)/(1-G.sigma))',1,30); %+ repmat(lambda2*(S.SS_X.^(1-G.sigma))/(1-G.sigma),10,1) ...
%+ lambda3*(wh_TVF.^(1-G.sigma))/(1-G.sigma) + lambda4*(K_TVF.^(1-G.sigma))/(1-G.sigma);

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
            coeff2(x,:)=coeff(x,:);
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
     
%     % draw husband wage
%     wh_mean = eta01 + eta02*(abi==2) + eta11*(edu==2) + eta12*(edu==3) + eta2*age;
%     wh_sd = 0.7022257; %eta03 + eta04*(abi==2) + eta21*(edu==2) + eta22*(edu==3) + eta3*age;
%     wh(t) = 0; %normrnd(wh_mean,wh_sd);
%     
%     % draw investments
%     Inv_mean = iota01 + iota02*(abi==2) + iota11*(edu==2) + iota12*(edu==3) + iota2*age;
%     Inv_sd = 0.9270494; %iota03 + iota04*(abi==2) + iota21*(edu==2) + iota22*(edu==3) + iota3*age;
%     Inv(t) = 0; %normrnd(Inv_mean,Inv_sd);
%     
%     % child human capital 
%     K(t) = 0; %exp(kappa01 + kappa02*(abi==2) + kappa03*(edu==2) + kappa04*(edu==3) + kappa05*(Inv(t)));
    
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
                chh_r = w_j_r + A_j; %+ exp(wh(t))*m_j + A_j;
                chh_n = w_j_n + A_j; %+ exp(wh(t))*m_j + A_j;
                chh_u = w_j_u + A_j; %+ exp(wh(t))*m_j + A_j;
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
                    u_r(k) = (cw_r^(1-G.sigma))/(1-G.sigma); %+ psi_r; %+ theta1_r*log(1+m_j) + theta2_r*log(1+n_j) + theta3_r*log(K(t))*n_j;
                    u_n(k) = (cw_n^(1-G.sigma))/(1-G.sigma); %+ psi_n; %+ theta1_n*log(1+m_j) + theta2_n*log(1+n_j) + theta3_n*log(K(t))*n_j;
                    u_u(k) = (cw_u^(1-G.sigma))/(1-G.sigma); %+ psi_u; %+ theta1_u*log(1+m_j) + theta2_u*log(1+n_j) + theta3_u*log(K(t))*n_j;
%                     uC_r(k)=(cw_r^(1-G.sigma))/(1-G.sigma)
%                     uC_n(k)=(cw_n^(1-G.sigma))/(1-G.sigma)
%                     uC_u(k)=(cw_u^(1-G.sigma))/(1-G.sigma)
%                     uM_r(k)=theta1_r*log(1+m_j)
%                     uM_n(k)=theta1_n*log(1+m_j)
%                     uM_u(k)=theta1_u*log(1+m_j)
%                     uN_r(k)=theta2_r*log(1+n_j)
%                     uN_n(k)=theta2_n*log(1+n_j)
%                     uN_u(k)=theta2_u*log(1+n_j)
%                     uK_r(k)=theta3_r*log(K(t))*n_j
%                     uK_n(k)=theta3_n*log(K(t))*n_j
%                     uK_u(k)=theta3_u*log(K(t))*n_j
                    % Maried with 1 kid:
                    if x <= 10
                        
                        % regular job
                        A_next = (1+G.r) * (A_j + (w_j_r + shock_i) - chh_r); % eq. 8
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
                        Vmr_next(k)=Vm_r_next;
                        Vmr_next2(k)=Vm_r_next2;
                        Vmr_next_2(k)=Vm_r_next_2;
                        Vmr_next2_2(k)=Vm_r_next2_2;
                                            
                        % non-regular job
                        A_next = (1+G.r) * (A_j + (w_j_n + shock_i) - chh_n);
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
                        Vmn_next(k)=Vm_n_next;
                        Vmn_next2(k)=Vm_n_next2;
                        Vmn_next_2(k)=Vm_n_next_2;
                        Vmn_next2_2(k)=Vm_n_next2_2;
                        
                        % unemployed
                        A_next = (1+G.r) * (A_j + (w_j_u + shock_i) - chh_u);
                        x_next = x;

                        % polynomial approximation of VF                      
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vm_u_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vm_u_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        Vm_u_next_2 = sum(coeff(x_next+10,:).*Base,2); %cheby_approx
                        Vm_u_next2_2 = sum(coeff2(x_next+10,:).*Base,2); %cheby_approx
                        % save output
                        Amu_next(k)=A_next;
                        Vmu_next(k)=Vm_u_next;
                        Vmu_next2(k)=Vm_u_next2;
                        Vmu_next_2(k)=Vm_u_next_2;
                        Vmu_next2_2(k)=Vm_u_next2_2;
                        
                        % Sector-Specific Value Functions (1 child) 
                          Vm1_r(k) = u_r(k) + G.beta * ((prob_lamba*Vm_r_next)+(1-prob_lamba)*Vm_r_next2);
                          Vm1_n(k) = u_n(k) + G.beta * ((prob_pi*Vm_n_next)+(1-prob_pi)*Vm_n_next2);
                          Vm1_u(k) = u_u(k) + G.beta * ((prob_pi*Vm_u_next)+(1-prob_pi)*Vm_u_next2);
                        
                        % Sector-Specific Value Functions (2 child, maybe)
                          Vm1_r_2(k) = u_r(k) + G.beta * ((prob_lamba*Vm_r_next_2)+(1-prob_lamba)*Vm_r_next2_2);
                          Vm1_n_2(k) = u_n(k) + G.beta * ((prob_pi*Vm_n_next_2)+(1-prob_pi)*Vm_n_next2_2);
                          Vm1_u_2(k) = u_u(k) + G.beta * ((prob_pi*Vm_u_next_2)+(1-prob_pi)*Vm_u_next2_2);

                        % Sector-Specific Value Fucntions (married)
                          Vm_r(k) = Vm1_r(k)*(1-prob_2kids_r) + Vm1_r_2(k)*(prob_2kids_r);
                          Vm_n(k) = Vm1_n(k)*(1-prob_2kids_n) + Vm1_n_2(k)*(prob_2kids_n);
                          Vm_u(k) = Vm1_u(k)*(1-prob_2kids_u) + Vm1_u_2(k)*(prob_2kids_u);
                        
                        % save marriage 1 child (for marriage decision)
                        Vm_r_aux(k,x,j,i) = Vm1_r(k);
                        Vm_n_aux(k,x,j,i) = Vm1_n(k);
                        Vm_u_aux(k,x,j,i) = Vm1_u(k);
                        
                    % Married with 2 kids: 
                    elseif x <= 20
                        
                        % regular job
                        A_next = (1+G.r) * (A_j + (w_j_r + shock_i) - chh_r); % eq. 8
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
                        Vm2r_next(k)=Vm_r_next;
                        Vm2r_next2(k)=Vm_r_next2;
                          
                        % non-regular job
                        A_next = (1+G.r) * (A_j + (w_j_n + shock_i) - chh_n);
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
                        Vm2n_next(k)=Vm_n_next;
                        Vm2n_next2(k)=Vm_n_next2;

                        % unemployed
                        A_next = (1+G.r) * (A_j + (w_j_u + shock_i) - chh_u);
                        x_next = x;

                        % Chebyshev Approximation
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vm_u_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vm_u_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        % save output
                        Am2u_next(k)=A_next;
                        Vm2u_next(k)=Vm_u_next;
                        Vm2u_next2(k)=Vm_u_next2;

                        % Sector-Specific Value Functions (2 child)
                        Vm2_r(k) = u_r(k) + G.beta * ((prob_lamba*Vm_r_next)+(1-prob_lamba)*Vm_r_next2);
                        Vm2_n(k) = u_n(k) + G.beta * ((prob_pi*Vm_n_next)+(1-prob_pi)*Vm_n_next2);
                        Vm2_u(k) = u_u(k) + G.beta * ((prob_pi*Vm_u_next)+(1-prob_pi)*Vm_u_next2);
                        
                    % Single
                    else
                        
                        % Regular
                        A_next = (1+G.r) * (A_j + (w_j_r + shock_i) - chh_r);
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
                        Vsr_next(k)=Vs_r_next;
                        Vsr_next2(k)=Vs_r_next2;

                        % Non-regular
                        A_next = (1+G.r) * (A_j + (w_j_n + shock_i) - chh_n);
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
                        Vsn_next(k)=Vs_n_next;
                        Vsn_next2(k)=Vs_n_next2;

                        % Unemployed
                        A_next = (1+G.r) * (A_j + (w_j_u + shock_i) - chh_u);
                        x_next = x;

                        % linear approximation of VF
                        Base=chebpoly_base(S.nA+1, S.d_A*(A_next - S.extmin_A) - 1);
                        Vs_u_next = sum(coeff(x_next,:).*Base,2); %cheby_approx
                        Vs_u_next2 = sum(coeff2(x_next,:).*Base,2); %cheby_approx
                        % save output
                        Asu_next(k)=A_next;
                        Vsu_next(k)=Vs_u_next;
                        Vsu_next2(k)=Vs_u_next2;

                        % Sector-Specific Value Functions (single)
                        Vs_r(k) = u_r(k) + G.beta * ((prob_lamba*Vs_r_next)+(1-prob_lamba)*Vs_r_next2);
                        Vs_n(k) = u_n(k) + G.beta * ((prob_pi*Vs_n_next)+(1-prob_pi)*Vs_n_next2);
                        Vs_u(k) = u_u(k) + G.beta * ((prob_pi*Vs_u_next)+(1-prob_pi)*Vs_u_next2);
                        Vsm_r(k) = prob_marr_r*Vm_r_aux(k,x-20,j,i) + (1-prob_marr_r)*Vs_r(k);
                        Vsm_n(k) = prob_marr_n*Vm_n_aux(k,x-20,j,i) + (1-prob_marr_n)*Vs_n(k);
                        Vsm_u(k) = prob_marr_u*Vm_u_aux(k,x-20,j,i) + (1-prob_marr_u)*Vs_u(k);
                    end
                end
%                 output1 = [cr_vector;cn_vector;cu_vector;u_r;u_n;u_u;Amr_next;Amn_next;Amu_next]
%                 output2 = [Vmr_next;Vmr_next2;Vmr_next_2;Vmr_next2_2;Vmn_next;Vmn_next2;Vmn_next_2;Vmn_next2_2;Vmu_next;Vmu_next2;Vmu_next_2;Vmu_next2_2]      
%                 output4 = [cr_vector;cn_vector;cu_vector;u_r;uK_r;u_n;uK_n;u_u;uK_u;Am2r_next;Am2n_next;Am2u_next]
%                 output5 = [Vm2r_next;Vm2r_next2;Vmn_next;Vmn_next2;Vmu_next;Vmu_next2]
%                 output7 = [cr_vector;cn_vector;cu_vector;u_r;uK_r;u_n;uK_n;u_u;uK_u;Asr_next;Asn_next;Asu_next]
%                 output8 = [Vsr_next;Vsr_next2;Vsn_next;Vsn_next2;Vsu_next;Vsu_next2]
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
                    [Vs_star, ls_index] = max([Vsm_r_star, Vsm_n_star, Vsm_u_star, Vs_r_star, Vs_n_star, Vs_u_star]);
                    [Vs_star2] = max([Vsm_n_star, Vsm_u_star, Vs_n_star, Vs_u_star]);
                end
%             output3 = [Vm1_r;Vm1_n;Vm1_u;Vm1_r_2;Vm1_n_2;Vm1_u_2;Vm_r;Vm_n;Vm_u]
%             lm_index
%             Vm_star
%             Vm_star2
%             cm_star_aux(lm_index)
%             output6 = [Vm2_r;Vm2_n;Vm2_u]
%             lm2_index
%             Vm2_star
%             Vm2_star2
%             cm2_star_aux(lm2_index)
%             output9 = [Vs_r;Vs_n;Vs_u;Vsm_r;Vsm_n;Vsm_u]
%             ls_index
%             Vs_star
%             Vs_star2
%             cs_star_aux(ls_index)
            % save choice:
            if x <= 10
                c_star(j, i, x, t) = cm_star_aux(lm_index);
                l_star(j, i, x, t) = lm_index;
                m_star(j, i, x, t) = (lm_index<=3);
                V_star(j, i, x, t) = Vm_star;
                V2_star(j,i, x, t) = Vm_star2;
            elseif x <= 20
                c_star(j, i, x, t) = cm2_star_aux(lm2_index);
                l_star(j, i, x, t) = lm2_index;
                m_star(j, i, x, t) = (lm_index<=3);
                V_star(j, i, x, t) = Vm2_star;
                V2_star(j,i, x, t) = Vm2_star2;
            else
                c_star(j, i, x, t) = cs_star_aux(ls_index);
                l_star(j, i, x, t) = ls_index;
                m_star(j, i, x, t) = (ls_index<=3);
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
            end
        end
        
        % Integrate over shocks
        W(:,x,t) = pi^(-1/2)*V_star(:,:,x,t)*S.weight;
        W2(:,x,t)= pi^(-1/2)*V2_star(:,:,x,t)*S.weight;
        
        % reshape policy func
        c_func(:,:,:,x,t) = reshape(c_star(:,:,x,t), [G.n_assets,3,3]);
        l_func(:,:,:,x,t) = reshape(l_star(:,:,x,t), [G.n_assets,3,3]);
        m_func(:,:,:,x,t) = reshape(m_star(:,:,x,t), [G.n_assets,3,3]);
    end
    %check = squeeze(c_star(:, 1, :, t))
end
%squeeze(V_star(:,1,[1,11,21],t))
% three labor functions (as 0 or 1)
lr_func = l_func == 1 | l_func == 4;
ln_func = l_func == 2 | l_func == 5;
lu_func = l_func == 3 | l_func == 6;

% marriage function for single women: equals 1 if labor choice is 1, 2, or 3
m_func2 = l_func == 1 | l_func == 2 | l_func == 3;
end