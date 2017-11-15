% loop for initial conditions:
% for n = 1:1:G.n_incond
n=1;
     abi = types(n,1);
     edu = types(n,2);

% Terminal Value Function:
% TVF = A_T + W_T + Q_T = assets + wages + HH_production

assets = SS_A; 
wages = exp(alpha1*SS_X + alpha2*edu + abi);
hhprod = (SS_M+1).^theta1 .* (SS_N+1).^theta2 .* SS_K;
% matrix of J=10x5x5=500 rows and 10x2=20 cols
TVF = assets + wages + hhprod;

tic
% loop for time (25):
for t = G.n_period-1:-1:1
    t
    toc
    if t==G.n_period-1
        Emax = TVF;
    else
        Emax = W;
    end
    
    % Chevyshev Approximation
    Num = Emax'*kron(T_A, kron(T_H,T_K)); % numerator (bases*function) 
    Den = kron(T2_A, kron(T2_H,T2_K)); % square of T multiplied
    for x = 1:1:(n_matstat*n_wrkexp)
        alpha(x,:) = Num(x,:)./Den';
    end
    % alpha contains 20 rows of 19x4x4 = 304 coefficients  
    
    % loop for shocks (27):
    for i = 1:1:G.n_shocks % 27 x 3
        i
        
        shock_i = shocks_i(i);
        shock_r = shocks_r(i);
        shock_n = shocks_n(i);
        
        % loop for work experience and marital status (20):
        for x = 1:1:(n_matstat*n_wrkexp)
            x;
            
            % current state variables:
            m_j = SS_M(x);  % marital status
            n_j = SS_N(x);  % children
            X_j = SS_X(x);  % experience
         
        % loop over states (20 assets x 5 child HC x 5 husband wages = 500):
        for j = 1:1:length(SS_rows)
            j;
            
            % current state variables:
            wh_j = SS_H(j); % husband's wage
            A_j = SS_A(j);  % HH's assets
            K_j = SS_K(j);  % HC of child
            
            % sector-specific state variables:
            w_j_r = exp(alpha1*X_j + alpha2*edu + abi + shock_r); % same
            w_j_n = exp(alpha1*X_j + alpha2*edu + abi + shock_n); % same
            w_j_u = 0; % unemployed women don't have earnings?
              
            % sector-specific probabilities: % CHECK THIS WHEN AGE VARIES
            prob_marr_r = (edu + abi + 35 + t)/100; % only change with time & type
            prob_marr_n = (edu + abi + 30 + t)/100; % only change with time & type
            prob_marr_u = (edu + abi + 25 + t)/100; % only change with time & type
            %probability = a0 + a1*col4 + a2*col2 + a3*work + a4*t;
            
            % transitions for exogenous variables:
            K_next = K_j + inv; % make a function (CES)
            wh_next = wh_j; % no transition
            
            % consumption vector
            c_vector = unique(assets); % consumption vector (unique assets, 20)
            %%%need to change this, so that if A_j=1, cons go from 0 to 1
            %%%you also have income and lowest shock, 0 to max chh based on
            %%%A_next function
            
            % loop over consumption
            for k = 1:1:length(c_vector)
                k;
                
                chh = c_vector(k); % HH consumption
                cw = 0.5*chh; % woman's consumption
                
                % single & married, get all six :
                if x <= 10

                    % regular job:
                    A_next = (1+r) * (A_j + (w_j_r + wh_j*m_j + shock_i) - chh - n_j*inv); % eq. 8
                    x_next = x + 1;
                    if x_next == 11
                        x_next = 10;
                    end
                    %%%% validate A_next compared to extmin_A and extmax_A
                    % value function:
                    Base=kron(chebpoly_base(18+1, d_A*(A_next - extmin_A) - 1),kron(chebpoly_base(3+1, d_H*(wh_next - extmin_H) - 1),chebpoly_base(3+1, d_K*(K_next - extmin_K) - 1)));
                    Vs_r_next = sum(alpha(x_next,:).*Base,2); %cheby_approx(coeffs,ncheby,extminA,extminB,extminC,dA,dB,dC,A_next,B_next,C_next)

                    % non-regular job:
                    A_next = (1+r) * (A_j + (w_j_n + wh_j*m_j + shock_i) - chh - n_j*inv);
                    x_next = x + 1;
                    if x_next == 11
                        x_next = 10;
                    end
                    % value function:
                    Base=kron(chebpoly_base(18+1, d_A*(A_next - extmin_A) - 1),kron(chebpoly_base(3+1, d_H*(wh_next - extmin_H) - 1),chebpoly_base(3+1, d_K*(K_next - extmin_K) - 1)));
                    Vs_n_next = sum(alpha(x_next,:).*Base,2); 
                    
                    % unemployed:
                    A_next = (1+r) * (A_j + (w_j_u + wh_j*m_j + shock_i) - chh - n_j*inv);
                    x_next = x;
                    % value function:
                    Base=kron(chebpoly_base(18+1, d_A*(A_next - extmin_A) - 1),kron(chebpoly_base(3+1, d_H*(wh_next - extmin_H) - 1),chebpoly_base(3+1, d_K*(K_next - extmin_K) - 1)));
                    Vs_u_next = sum(alpha(x_next,:).*Base,2); % extmin_A was NOT low enough (0.99, A_next = -0.4144), can we change each approx?
                       
                % Married:
                else
                    % Regular:
                    A_next = (1+r) * (A_j + (w_j_r + wh_j*m_j + shock_i) - chh - n_j*inv);
                    x_next = x + 1;
                    if x_next == 21
                        x_next = 20;
                    end
                    Base=kron(chebpoly_base(18+1, d_A*(A_next - extmin_A) - 1),kron(chebpoly_base(3+1, d_H*(wh_next - extmin_H) - 1),chebpoly_base(3+1, d_K*(K_next - extmin_K) - 1)));
                    Vm_r_next = sum(alpha(x_next,:).*Base,2);
                    
                    % Non-regular:
                    A_next = (1+r) * (A_j + (w_j_n + wh_j*m_j + shock_i) - chh - n_j*inv);
                    x_next = x + 1;
                    if x_next == 21
                        x_next = 20;
                    end
                    Base=kron(chebpoly_base(18+1, d_A*(A_next - extmin_A) - 1),kron(chebpoly_base(3+1, d_H*(wh_next - extmin_H) - 1),chebpoly_base(3+1, d_K*(K_next - extmin_K) - 1)));
                    Vm_n_next = sum(alpha(x_next,:).*Base,2);
                    
                    % Unemployed:
                    A_next = (1+r) * (A_j + (w_j_u + wh_j*m_j + shock_i) - chh - n_j*inv);
                    x_next = x;
                    Base=kron(chebpoly_base(18+1, d_A*(A_next - extmin_A) - 1),kron(chebpoly_base(3+1, d_H*(wh_next - extmin_H) - 1),chebpoly_base(3+1, d_K*(K_next - extmin_K) - 1)));
                    Vm_u_next = sum(alpha(x_next,:).*Base,2);
                    %%%%save A_next and consumption
                end
                
                % Sector-Specific Utility:
                u_r(k) = (cw^(1-sigma))/(1-sigma) + psi_r + kappa*(1+theta1*m_j+theta2*n_j+theta3*K_j);
                u_n(k) = (cw^(1-sigma))/(1-sigma) + psi_n + kappa*(1+theta1*m_j+theta2*n_j+theta3*K_j);
                u_u(k) = (cw^(1-sigma))/(1-sigma) + kappa*(1+theta1*m_j+theta2*n_j+theta3*K_j);
                
                if x <= 10
                    % Sector-Specific Value Functions
                    Vs_r(k) = u_r(k) + beta * ((1-prob_marr_r)*Vs_r_next + prob_marr_r*Vm_r_next);
                    Vs_n(k) = u_n(k) + beta * ((1-prob_marr_n)*Vs_n_next + prob_marr_n*Vm_n_next);
                    Vs_u(k) = u_u(k) + beta * ((1-prob_marr_u)*Vs_u_next + prob_marr_u*Vm_u_next); 
                else
                    % Sector-Specific Value Functions
                    Vm_r(k) = u_r(k) + beta * Vm_r_next;
                    Vm_n(k) = u_n(k) + beta * Vm_n_next;
                    Vm_u(k) = u_u(k) + beta * Vm_u_next;
                end
            end
            
            %%% validate that consumption is feasible
            %%% check that A_next is within extmin_A and extmin_B
                % save optimal V* & c*
                [Vs_r_star, Index_sr_k] = max(Vs_r);
                [Vs_n_star, Index_sn_k] = max(Vs_n);
                [Vs_u_star, Index_su_k] = max(Vs_u);
                [Vm_r_star, Index_mr_k] = max(Vm_r);
                [Vm_n_star, Index_mn_k] = max(Vm_n);
                [Vm_u_star, Index_mu_k] = max(Vm_u);
                cs_r_star = c_vector(Index_sr_k);
                cs_n_star = c_vector(Index_sn_k);
                cs_u_star = c_vector(Index_su_k);
                cm_r_star = c_vector(Index_mr_k);
                cm_n_star = c_vector(Index_mn_k);
                cm_u_star = c_vector(Index_mu_k);
                c_star_vector = [c_r_star, c_n_star, c_u_star];
                
                % save labor choice:
                [V_star, Index_l] = max([V_r_star,V_n_star,V_u_star]);
                c_star(j, i, t) = c_star_vector(Index_l); %c_star(j,i,s,t)
                l_star(j, i, t) = Index_l;
                V_star_aux(j, i, t) = V_star; 
                
        end
        end
    end
    
    % Integrate over shocks
%     W(:,t)=pi^(-1/2)*V(:,:,t)*G.wt; 
%     Em(:,t,k) = detR*detV^(-1/2)*pi^(-(size(R,1))/2)*G.w(i)*G.w(j)*V + Em(:,t,k);
    W = pi^(-1/2)*V_star_aux(:,:,t)*weight;   
end
%%% time this shit
% revers order of matstat, to do 10 rows of married first and then 10 rows
% of single, and then when single i can compare 6 value functions