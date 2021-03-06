function [alpC,alpR,alpN,alpU,alpM]=polfunc_approx_small(C,R,N,U,M,S,G)

%T_sim=kron(S.T_A, kron(S.T_H,kron(S.T_K,kron(S.Teps_r,kron(S.Teps_n,S.Teps_i)))));
%Den= kron(S.T2_A, kron(S.T2_H,kron(S.T2_K,kron(S.T2eps_r,kron(S.T2eps_n,S.T2eps_i)))));
T_sim=kron(S.T_A, kron(S.T_H,kron(S.T_K,kron(S.Teps_r,S.Teps_n))));
Den= kron(S.T2_A, kron(S.T2_H,kron(S.T2_K,kron(S.T2eps_r,S.T2eps_n))));

%%% Polynomial Bases and Derivatives %%%% 

 % square of T multiplied
for t=1:1:G.n_period-1
	for x = 1:1:(G.n_matstat*G.n_wrkexp)
	NumC(:,x,t) = C(:,x,t)'*T_sim;  
	alpC(:,x,t) = NumC(:,x,t)./Den;
	NumR(:,x,t) = R(:,x,t)'*T_sim; 
	alpR(:,x,t) = NumR(:,x,t)./Den;
	NumN(:,x,t) = N(:,x,t)'*T_sim; 
	alpN(:,x,t) = NumN(:,x,t)./Den;
	NumU(:,x,t) = U(:,x,t)'*T_sim; 
	alpU(:,x,t) = NumU(:,x,t)./Den;
	NumM(:,x,t) = M(:,x,t)'*T_sim; 
	alpM(:,x,t) = NumM(:,x,t)./Den;
	end
end
