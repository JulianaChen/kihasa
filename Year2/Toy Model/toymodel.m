clear all; clc;

%% Parameters

sigma = 0.8; 
lambda = 1;
beta = 0.8;
var = 0.2;
r=0.07;

Ne = 3;
n_shocks = 3;
n_period = 20;
n_assets = 15;
n_cons = 30;

G = struct('Ne',Ne,'sigma',sigma,'lambda',lambda,'beta',beta,'r',r,...
    'n_assets',n_assets,'n_period',n_period,'n_shocks',n_shocks,'n_cons',n_cons);

%% State Space

% shocks
[e, wt] = GaussHermite(Ne);
eps = sqrt(2)*e*var; 
shocks = eps;
weight = wt;
% basis for shocks
zeps = 2*(eps-eps(1))/(eps(Ne,1)-eps(1))-1; 
Teps=chebpoly_base(Ne-1,zeps);
T2eps = diag(Teps'*Teps);

% assets
assets_lb  = 0;
assets_int = 8000;
assets_ub  = 40000;
% values for polynomial approximation
[assets,nA,extmin_A,extmax_A,d_A,T_A,T2_A] = cheby_values(G.n_assets,assets_ub,assets_lb);
% SS_A = assets;
% asset vector
assets1 = linspace(assets_lb,assets_int,n_assets-3);
assets2 = linspace(assets_int,assets_ub,4);
SS_A = [assets1,assets2([2:4])];

S = struct('Teps',Teps,'T2eps',T2eps,'shocks',shocks,'weight',weight,'nA',nA,...
    'extmin_A',extmin_A,'extmax_A',extmax_A,'d_A',d_A,'T_A',T_A,'T2_A',T2_A,...
    'SS_A',SS_A);

%% Policy Function

% expanded assets
A_min=S.extmin_A;
A_max=S.extmax_A;
A_wide = S.SS_A;
A_wide(1) = A_min; %min(S.SS_A);
A_wide(G.n_assets)= A_max; %max(S.SS_A);

% fixed (constant) income
wage = log(2000);

% Terminal Value Function
%TVF = G.lambda*(1-exp(-S.SS_A));
TVF = G.lambda*(S.SS_A).^(1-G.sigma)/(1-G.sigma);

tic
% loop for time (20):
for t = G.n_period-1:-1:1
    t
    toc
    
    % Coefficients for Polynomial Approximation
    if t==G.n_period-1
%         Emax_linear = TVF;
        Emax = TVF;
        Num = Emax*S.T_A;
        Den = S.T2_A;
        coeff = Num./Den';
    else
%         Emax_linear = W_linear(:,t+1);
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
            c1 = linspace(c_min,0.75*c_max,n_cons-5);
            c2 = linspace(0.75*c_max,c_max,6);
            c_vector = [c1, c2(2:6)];
            
            % loop for consumtion (30):
            for k = 1:1:G.n_cons
                k
                
                % consumption
                c = c_vector(k);
                
                % utility
                u(k) = (c^(1-G.sigma))/(1-G.sigma);
                
                % assets
                A_next(k) = (1+G.r) * (A_j + exp(wage_i) - c);
                
                % approx
                Base = chebpoly_base(S.nA+1, S.d_A*(A_next(k) - S.extmin_A) - 1);
                V_next = sum(coeff.*Base,2);
                V_next_aux(k) = V_next;
%                 V_next_linear = interpn(A_wide,Emax_linear,A_next(k));
%                 V_next_aux_linear(k) = V_next_linear;
                % value function
                Value(k) = u(k) + G.beta * V_next;
%                 Value_linear(k) = u(k) + G.beta * V_next_linear;
            end
            
            % validate
            Value(A_next < min(A_wide)) = NaN;
            Value(A_next > max(A_wide)) = NaN;
%             Value_linear(A_next < min(A_wide)) = NaN;
%             Value_linear(A_next > max(A_wide)) = NaN;
            
            % save optimal
            [V_max, index] = max(Value);
            c_max = c_vector(index);
%             [V_max_linear, index_linear] = max(Value_linear);
%             c_max_linear = c_vector(index_linear);
            
            % save choice
            c_star(j,i,t) = c_max;
            V_star(j,i,t) = V_max;
%             c_star_linear(j,i,t) = c_max_linear;
%             V_star_linear(j,i,t) = V_max_linear;
            
            % save outside
            A_out(j,i,t) = sum(A_next < A_min) + sum(A_next > A_max);
        end
    end
    % Integrate
    W(:,t) = pi^(-1/2)*V_star(:,:,t)*S.weight;
%     W_linear(:,t) = pi^(-1/2)*V_star_linear(:,:,t)*S.weight;
end

cons=squeeze(c_star(:,1,:))
%cons_lin=squeeze(c_star_linear(:,1,:))

tau=1:1:19

plot(tau,cons(1,:),tau,cons(15,:))
plot(tau,cons(2,:),tau,cons(15,:))
