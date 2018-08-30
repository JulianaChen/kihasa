function [c_star, V_star, l_star, A_out, A_out_u] = solution2(G,P,S)

% expanded assets
A_min=S.extmin_A;
A_max=S.extmax_A;
A_wide = S.SS_A;
A_wide(1) = A_min;
A_wide(G.n_assets)= A_max;

% Terminal Value Function:
%TVF = P.lambda*(S.SS_A).^(1-P.sigma)/(1-P.sigma);
TVF = repmat((P.lambda1*(S.SS_A).^(1-P.sigma)/(1-P.sigma))',1,G.n_exp) + ...
      repmat((P.lambda2*(S.SS_X).^(1-P.sigma)/(1-P.sigma)),G.n_assets,1);

tic
% loop for time (20):
for t = G.n_period-1:-1:1
    t
    toc
    
    % Coefficients for Polynomial Approximation
    if t==G.n_period-1
        Emax = TVF;
        Num = Emax'*S.T_A;
        Den = S.T2_A;
        coeff = Num./Den';
    else
        Emax = W(:,:,t+1);
        Num = Emax'*S.T_A;
        Den = S.T2_A;
        coeff = Num./Den';
    end
    
    % age
    age(t) = 22 + t;
    
    % loop for shocks (3):
    for i = 1:1:G.n_shocks
        i;
        
        % shocks
        shock = S.shocks(i);
                
        % loop over assets (15):
        for j = 1:1:G.n_assets
            j;
            
            % household assets
            A_j = S.SS_A(j);
            
            % loop over experience (10):
            for x = 1:1:G.n_exp
                x;
                
                % work experience
                X_j = S.SS_X(x);
                
                % probability of work %% ASK: probability of keeping a job?
                prob_work = normcdf(P.tau10 + P.tau13*age(t) + P.tau14*X_j);
                
                % wage
                wage = P.alpha01 + P.alpha02*log(1+X_j) + shock;
                wage_u = 0;
                
                % consumption vector (non-linear)
                c_max = max(P.c_min, exp(wage) + A_j);
                cu_max = max(P.c_min, exp(wage_u) + A_j);
                c1 = linspace(P.c_min,0.75*c_max,G.n_cons-5);
                c2 = linspace(0.75*c_max,c_max,6);
                c1u = linspace(P.c_min,0.75*cu_max,G.n_cons-5);
                c2u = linspace(0.75*cu_max,cu_max,6);
                c_vector = [c1, c2(2:6)];
                cu_vector = [c1u, c2u(2:6)];
                if A_j == 0
                    cu_vector = repmat(50,1,G.n_cons);
                end
            
                % loop for consumtion (30):
                for k = 1:1:G.n_cons
                    k;

                    % consumption
                    c = c_vector(k);
                    c_u = cu_vector(k);
                    
                    % utility
                    u(k) = (c^(1-P.sigma))/(1-P.sigma) - P.psi;
                    u_u(k) = (c_u^(1-P.sigma))/(1-P.sigma) - P.psi_u;
                    
                    % work
                    % assets
                    A_next(k) = (1+G.r) * (A_j + exp(wage) - c);
                    % experience
                    x_next = x+1;
                    if x_next == 11
                        x_next = 10;
                    end
                    % approx %%% ASK ITALO: DO I NEED TO APPROX FOR NO WORK
                    Base = chebpoly_base(S.nA+1, S.d_A*(A_next(k) - S.extmin_A) - 1);
                    V_next(k) = sum(coeff(x_next,:).*Base,2);
                    
                    % no work
                    % assets
                    A_next_u(k) = (1+G.r) * (A_j + exp(wage_u) - c_u);
                    % experience
                    x_next = x;
                    % approx
                    Base = chebpoly_base(S.nA+1, S.d_A*(A_next_u(k) - S.extmin_A) - 1);
                    V_next_u(k) = sum(coeff(x_next,:).*Base,2);
                
                    % value function
                    Value(k) = u(k) + G.beta * ((prob_work*V_next(k))+(1-prob_work)*V_next_u(k));
                    Value_u(k) = u_u(k) + G.beta * ((prob_work*V_next(k))+(1-prob_work)*V_next_u(k));
                end

                % validate
                Value(A_next < min(A_wide)) = NaN;
                Value(A_next > max(A_wide)) = NaN;
                Value_u(A_next_u < min(A_wide)) = NaN;
                Value_u(A_next_u > max(A_wide)) = NaN;

                % save optimal consumption
                [V_max, index] = max(Value);
                c_max = c_vector(index);
                [V_max_u, index_u] = max(Value_u);
                c_max_u = cu_vector(index_u);
                c_star_aux = [c_max,c_max_u];
                
                % save optimal value and labor choice
                [V_star, index_l] = max([V_max,V_max_u]);

                % policy functions
                c_star(j,i,x,t) = c_star_aux(index_l);
                v_star(j,i,x,t) = V_star;
                l_star(j,i,x,t) = index_l;

                % save outside
                A_out(j,i,x,t) = sum(A_next < A_min) + sum(A_next > A_max);
                A_out_u(j,i,x,t) = sum(A_next_u < A_min) + sum(A_next_u > A_max);
            end
        end
    end
    
    % Integrate
    for x = 1:1:G.n_exp
        W(:,x,t) = pi^(-1/2)*v_star(:,:,x,t)*S.weight;
    end
end

end