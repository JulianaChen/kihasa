function [c_s,r_s,n_s,u_s,m_s,a_s,k_s,wr_s,wn_s] = simulation(params,C,R,N,U,M,G,S,abi,edu,type,wh)

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

for t=1:1:G.n_period-1
    for n=1:1:G.n_pop
     
    epssim_r(n,t)=sqrt(2)*G.Eps(1,n,t)'*sigma_r;
    epssim_n(n,t)=sqrt(2)*G.Eps(2,n,t)'*sigma_n;
%     T_eps_r=chebpoly_base(G.Ne-1, 2*(epssim_r(n,t) - S.eps_r(1))/(S.eps_r(G.Ne)-S.eps_r(1)) - 1);
%     T_eps_n=chebpoly_base(G.Ne-1, 2*(epssim_n(n,t) - S.eps_n(1))/(S.eps_n(G.Ne)-S.eps_n(1)) - 1);
%     T_a=chebpoly_base(S.nA+1, S.d_A*(a_s(n,t) - S.extmin_A) - 1);
%     T_h=chebpoly_base(S.nH+1, S.d_H*(wh_s(n) - S.extmin_H) - 1);
%     T_k=chebpoly_base(S.nK+1, S.d_K*(k_s(n,t) - S.extmin_K) - 1);
%     T_s=kron(T_a,kron(T_h,kron(T_k,kron(T_eps_r,T_eps_n))));
    
    % Locate in the experience/marriage vector    
    x=(1-m_s(n,t))*10 + exp_s(n,t)+1;
    
    % Calculate wages
    wr_s(n,t) = exp(alpha01_r + alpha02_r*(abi(n)==2) + alpha11_r*(edu(n)==2) + alpha12_r*(edu(n)==3) + alpha2_r*log(1+exp_s(n,t)) + epssim_r(n,t)); 
    wn_s(n,t) = exp(alpha01_n + alpha02_n*(abi(n)==2) + alpha11_n*(edu(n)==2) + alpha12_n*(edu(n)==3) + alpha2_n*log(1+exp_s(n,t)) + epssim_n(n,t));
    
    % Optimal Choices
    c_s(n,t)=sum(alpC(:,x,t,type).*T_s',1);
    rr(n,t)=sum(alpR(:,x,t,type).*T_s',1);
    nn(n,t)=sum(alpN(:,x,t,type).*T_s',1);
    uu(n,t)=sum(alpU(:,x,t,type).*T_s',1);
    [v, Ind] = max([rr(n,t), nn(n,t), uu(n,t)])
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
       marr(n,t)=sum(alpM(:,x,t,type).*T_s',1);
       if marr(n,t)<0.5 
          m_s(n,t+1)=0;
       else
          m_s(n,t+1)=1;
       end
    else
        m_s(n,t+1)=m_s(n,t);
    end
   
    % Transition functions
    a_s(n,t+1)= (1+G.r)*(a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*wh_s(n,t) - m_s(n,t)*G.Inv - c_s(n,t));
    k_s(n,t+1)=(gamma1*k_s(n,t)^phi + (1-gamma1)*G.Inv^phi)^(1/phi);
    wh_s(n,t+1)=wh_s(n,t); 
    
    n
    
    end
    t
end



end

