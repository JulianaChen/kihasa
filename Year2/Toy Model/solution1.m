function [c_star, V_star] = solution1(G,P,S)

% expanded assets
A_min=S.extmin_A;
A_max=S.extmax_A;
A_wide = S.SS_A;
A_wide(1) = A_min; %min(S.SS_A);
A_wide(G.n_assets)= A_max; %max(S.SS_A);

% fixed (constant) income
wage = log(2000);

% Terminal Value Function:
%TVF = G.lambda*(1-exp(-S.SS_A));
TVF = P.lambda*(S.SS_A).^(1-P.sigma)/(1-P.sigma);

tic
% loop for time (20):
for t = G.n_period-1:-1:1
    t
    toc
    
    % Coefficients for Polynomial Approximation
    if t==G.n_period-1
        Emax = TVF;
        Num = Emax*S.T_A;
        Den = S.T2_A;
        coeff = Num./Den';
    else
        Emax = W(:,t+1);
        Num = Emax'*S.T_A;
        Den = S.T2_A;
        coeff = Num./Den';
    end
    
    % loop for shocks (3):
    for i = 1:1:G.n_shocks
        i
        
        % shocks
        shock = S.shocks(i);
        
        % wages
        wage_i = wage + shock;
        
        % loop over assets (15):
        for j = 1:1:G.n_assets
            j
            
            % household assets
            A_j = S.SS_A(j);
            
            % consumption vector
            c_min = 0.5*exp(wage_i); %this needs to be a parameter/function = state transfer to poor/unemployed
            c_max = exp(wage_i) + A_j; %check if min consumption matters ()
            %c_vector = linspace(c_min,c_max,G.n_cons); %try this non linear consumption
            
            % non-linear consumption vector
            c1 = linspace(c_min,0.75*c_max,G.n_cons-5);
            c2 = linspace(0.75*c_max,c_max,6);
            c_vector = [c1, c2(2:6)];
            
            % loop for consumtion (30):
            for k = 1:1:G.n_cons
                k
                
                % consumption
                c = c_vector(k);
                
                % utility
                u(k) = (c^(1-P.sigma))/(1-P.sigma);
                
                % assets
                A_next(k) = (1+G.r) * (A_j + exp(wage_i) - c);
                
                % approx
                Base = chebpoly_base(S.nA+1, S.d_A*(A_next(k) - S.extmin_A) - 1);
                V_next = sum(coeff.*Base,2);
                
                % value function
                Value(k) = u(k) + G.beta * V_next;
            end
            
            % validate
            Value(A_next < min(A_wide)) = NaN;
            Value(A_next > max(A_wide)) = NaN;
            
            % save optimal
            [V_max, index] = max(Value);
            c_max = c_vector(index);
            
            % save choice
            c_star(j,i,t) = c_max;
            V_star(j,i,t) = V_max;
        
            % save outside
            A_out(j,i,t) = sum(A_next < A_min) + sum(A_next > A_max);
        end
    end
    % Integrate
    W(:,t) = pi^(-1/2)*V_star(:,:,t)*S.weight;
end

end
