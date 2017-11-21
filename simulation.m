function [c_s,r_s,n_s,u_s,m_s,a_s,k_s,wr_s,wn_s] = simulation(params,alpC,alpR,alpN,alpU,alpM,G,S,type,wh)

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


% Initial Conditions
wh_s=wh;
a_s = zeros(G.n_pop,G.n_period); 
k_s= zeros(G.n_pop,G.n_period); 
exp_s= zeros(G.n_pop,G.n_period); 
m_s= zeros(G.n_pop,G.n_period); 

for t=1:1:G.nper
    for n=1:1:G.npop
     
    epssim_r(n,t)=sqrt(2)*G.Eps(1,n,t)'*sigma_r;
    epssim_n(n,t)=sqrt(2)*G.Eps(2,n,t)'*sigma_n;
    epssim_i(n,t)=sqrt(2)*G.Eps(3,n,t)'*sigma_i;
        
    T_eps_r=chebpoly_base(G.Ne-1, 2*(epssim_r(n,t) - G.eps_r(1))/(G.eps_r(G.Ne)-G.eps_r(1)) - 1);
    T_eps_n=chebpoly_base(G.Ne-1, 2*(epssim_n(n,t) - G.eps_n(1))/(G.eps_n(G.Ne)-G.eps_n(1)) - 1);
    T_eps_i=chebpoly_base(G.Ne-1, 2*(epssim_i(n,t) - G.eps_i(1))/(G.eps_i(G.Ne)-G.eps_i(1)) - 1);
    T_a=chebpoly_base(S.nA+1, S.dA*(a_s(n,t) - S.extmin_A) - 1);
    T_h=chebpoly_base(S.nH+1, S.dH*(wh_s(n) - S.extmin_H) - 1);
    T_k=chebpoly_base(S.nK+1, S.dK*(k_s(n,t) - S.extmin_K) - 1);
    T_s=kron(T_a,kron(T_h,kron(T_k,kron(Teps_r,kron(Teps_n,Teps_i)))));
    
    % Locate in the experience/marriage vector    
    x=(1-m_s(n,t))*10 + exp_s(n,t)+1;
    
    % Calculate wages
    wr_s(n,t) = exp(alpha01_r + alpha02_r*(abi==2) + alpha11_r*(edu==2) + alpha12_r*(edu==3) + alpha2_r*log(1+exp_s(n,t)) + epssim_r(n,t)); 
    wn_s(n,t) = exp(alpha01_n + alpha02_n*(abi==2) + alpha11_n*(edu==2) + alpha12_n*(edu==3) + alpha2_n*log(1+exp_s(n,t)) + epssim_n(n,t));
    
    % Optimal Choices
    c_s(n,t)=sum(alpC(:,x,t,type).*T_s,2);
    r(n,t)=sum(alpR(:,x,t,type).*T_s,2);
    n(n,t)=sum(alpN(:,x,t,type).*T_s,2);
    u(n,t)=sum(alpU(:,x,t,type).*T_s,2);
    [v, Ind] = max([r_s(n,t), n_s(n,t), u_s(n,t)])
    if Ind==1  
     r_s(n,t)=1;
     n_s(n,t)=0;
     u_s(n,t)=0;
     exp_s(n,t+1)=exp_s(n,t)+1;
    elseif Ind==2
     r_s(n,t)=0;
     n_s(n,t)=1;
     u_s(n,t)=0;
     exp_s(n,t+1)=exp_s(n,t)+1;
    else
     r_s(n,t)=0;
     n_s(n,t)=0;
     u_s(n,t)=1;
     exp_s(n,t+1)=exp_s(n,t);
    end
    
    if m_s(n,t)==0
       marr(n,t)=sum(alpM(:,x,t,type).*T_s,2);
       if marr(n,t)<0.5 
          m_s(n,t+1)=0;
       else
          m_s(n,t+1)=1;
       end
    else
        m_s(n,t+1)=m_s(n,t);
    end
   
    % Transition functions
    a_s(n,t+1)= (1+P.r)*(a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*wh_s(n,t) + epssim_i(n,t) - m_s(n,t)*G.Inv - c_s(n,t));
    k_s(n,t+1)=(gamma1*k_s(n,t)^phi + (1-gamma1)*G.Inv^phi)^(1/phi);
    wh_s(n,t+1)=wh_s(n,t); 
    
    n
    
    end
    t
end



end

