function [c_slin,r_slin,n_slin,u_slin,m_slin,a_slin,k_slin,wr_slin,wn_slin] = simulation(params,alpC,alpR,alpN,alpU,alpM,G,S,abi,edu,type,wh_s,C_lin,M_lin,R_lin,N_lin,U_lin)

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
wh_s=wh_s;
a_s = zeros(G.n_pop,G.n_period); 
a_slin = zeros(G.n_pop,G.n_period); 
k_s= zeros(G.n_pop,G.n_period); 
k_slin= zeros(G.n_pop,G.n_period); 
exp_s= zeros(G.n_pop,G.n_period); 
exp_slin= zeros(G.n_pop,G.n_period); 
m_s= zeros(G.n_pop,G.n_period); 
m_slin= zeros(G.n_pop,G.n_period); 
c_s = zeros(G.n_pop,G.n_period); 
c_slin = zeros(G.n_pop,G.n_period); 

for t=1:1:G.n_period-1
    for n=1:1:G.n_pop
     
    epssim_r(n,t)=sqrt(2)*G.Eps(1,n,t)'*sigma_r;
    epssim_n(n,t)=sqrt(2)*G.Eps(2,n,t)'*sigma_n;
    
    if epssim_r(n,t)<S.eps_r(1) || epssim_r(n,t)>S.eps_r(3)
        eps_rg=epssim_r(n,t)*[-1;0;1];
    else 
        eps_rg=S.eps_r;
    end
    if epssim_n(n,t)<S.eps_n(1) || epssim_n(n,t)>S.eps_n(3)
        eps_ng=epssim_n(n,t)*[-1;0;1];
    else 
        eps_ng=S.eps_n;
    end
    
    T_eps_r=chebpoly_base(G.Ne-1, 2*(epssim_r(n,t) - eps_rg(1))/(eps_rg(G.Ne)-eps_rg(1)) - 1);
    T_eps_n=chebpoly_base(G.Ne-1, 2*(epssim_n(n,t) - eps_ng(1))/(eps_ng(G.Ne)-eps_ng(1)) - 1);
    T_a=chebpoly_base(S.nA+1, 2*(a_s(n,t) - S.assets(1))/(S.assets(G.n_assets)-S.assets(1)) - 1);
    T_h=chebpoly_base(S.nH+1, S.d_H*(wh_s(n) - S.extmin_H) - 1);
    T_k=chebpoly_base(S.nK+1, S.d_K*(k_s(n,t) - S.extmin_K) - 1);
    T_s=kron(T_a,kron(T_h,kron(T_k,kron(T_eps_r,T_eps_n))));
    
    % Locate in the experience/marriage vector    
    x=min(20,(1-m_s(n,t))*10 + exp_s(n,t)+1);
    
    % Calculate wages
    wr_s(n,t) = exp(alpha01_r + alpha02_r*(abi(n)==2) + alpha11_r*(edu(n)==2) + alpha12_r*(edu(n)==3) + alpha2_r*log(1+exp_s(n,t)) + epssim_r(n,t)); 
    wn_s(n,t) = exp(alpha01_n + alpha02_n*(abi(n)==2) + alpha11_n*(edu(n)==2) + alpha12_n*(edu(n)==3) + alpha2_n*log(1+exp_s(n,t)) + epssim_n(n,t));

    wr_slin(n,t) = exp(alpha01_r + alpha02_r*(abi(n)==2) + alpha11_r*(edu(n)==2) + alpha12_r*(edu(n)==3) + alpha2_r*log(1+exp_slin(n,t)) + epssim_r(n,t)); 
    wn_slin(n,t) = exp(alpha01_n + alpha02_n*(abi(n)==2) + alpha11_n*(edu(n)==2) + alpha12_n*(edu(n)==3) + alpha2_n*log(1+exp_slin(n,t)) + epssim_n(n,t));
    
    wr_slin_nshock(n,t) = exp(alpha01_r + alpha02_r*(abi(n)==2) + alpha11_r*(edu(n)==2) + alpha12_r*(edu(n)==3) + alpha2_r*log(1+exp_slin(n,t))); 
    wn_slin_nshock(n,t) = exp(alpha01_n + alpha02_n*(abi(n)==2) + alpha11_n*(edu(n)==2) + alpha12_n*(edu(n)==3) + alpha2_n*log(1+exp_slin(n,t)));
   
    asset_g=S.assets;
%     if a_slin(n,t)>S.assets(G.n_assets)
%         asset_g(G.n_assets)=a_slin(n,t);
%     end
       
    % Optimal Choices
    cc_s(n,t)=sum(alpC(:,x,t,type(n)).*T_s',1);
    cc_slin(n,t)= interpn(S.childK,asset_g,S.hwages,eps_rg,eps_ng,C_lin(:,:,:,:,:,x,t,type(n)),k_slin(n,t),a_slin(n,t),wh_s(n,t),epssim_r(n,t),epssim_n(n,t));
    rr(n,t)=sum(alpR(:,x,t,type(n)).*T_s',1);
    rr_lin(n,t)=interpn(S.childK,asset_g,S.hwages,eps_rg,eps_ng,R_lin(:,:,:,:,:,x,t,type(n)),k_slin(n,t),a_slin(n,t),wh_s(n,t),epssim_r(n,t),epssim_n(n,t));
    nn(n,t)=sum(alpN(:,x,t,type(n)).*T_s',1);
    nn_lin(n,t)=interpn(S.childK,asset_g,S.hwages,eps_rg,eps_ng,N_lin(:,:,:,:,:,x,t,type(n)),k_slin(n,t),a_slin(n,t),wh_s(n,t),epssim_r(n,t),epssim_n(n,t));
    uu(n,t)=sum(alpU(:,x,t,type(n)).*T_s',1);
    uu_lin(n,t)=interpn(S.childK,asset_g,S.hwages,eps_rg,eps_ng,U_lin(:,:,:,:,:,x,t,type(n)),k_slin(n,t),a_slin(n,t),wh_s(n,t),epssim_r(n,t),epssim_n(n,t));
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
    
    [v_lin, Ind_lin] = max([rr_lin(n,t), nn_lin(n,t), uu_lin(n,t)]);
    if Ind_lin==1  
     r_slin(n,t)=1;
     n_slin(n,t)=0;
     u_slin(n,t)=0;
     exp_slin(n,t+1)=exp_slin(n,t)+1;
    elseif Ind_lin==2
     r_slin(n,t)=0;
     n_slin(n,t)=1;
     u_slin(n,t)=0;
     exp_slin(n,t+1)=exp_slin(n,t)+1;
    else
     r_slin(n,t)=0;
     n_slin(n,t)=0;
     u_slin(n,t)=1;
     exp_slin(n,t+1)=exp_slin(n,t);
    end
    
    if m_s(n,t)==0
       marr(n,t)=sum(alpM(:,x,t,type(n)).*T_s',1);
       if marr(n,t)<0.5 
          m_s(n,t+1)=0;
       else
          m_s(n,t+1)=1;
       end
    else
        m_s(n,t+1)=m_s(n,t);
    end
    
    if m_slin(n,t)==0
       marr_lin(n,t)=interpn(S.childK,asset_g,S.hwages,eps_rg,eps_ng,M_lin(:,:,:,:,:,x,t,type(n)),k_slin(n,t),a_slin(n,t),wh_s(n,t),epssim_r(n,t),epssim_n(n,t));
       if marr_lin(n,t)<0.5 
          m_slin(n,t+1)=0;
       else
          m_slin(n,t+1)=1;
       end
    else
        m_slin(n,t+1)=m_slin(n,t);
    end
   
    % Transition functions
    if (1+G.r)*(a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*wh_s(n) - m_s(n,t)*G.Inv-cc_s(n,t))<S.assets(1)
        c_s(n,t)=a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*wh_s(n) - m_s(n,t)*G.Inv - S.assets(1)/(1+G.r);    
    else
        c_s(n,t)=cc_s(n,t);
    end
    a_s(n,t+1)=(1+G.r)*(a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*wh_s(n) - m_s(n,t)*G.Inv-c_s(n,t));
    k_s(n,t+1)=(gamma1*k_s(n,t)^phi + (1-gamma1)*G.Inv^phi)^(1/phi);
    
    % Transition functions linear approx
    if (1+G.r)*(a_slin(n,t) + r_slin(n,t)*wr_slin(n,t) + n_slin(n,t)*wn_slin(n,t) + m_slin(n,t)*wh_s(n) - m_slin(n,t)*G.Inv-cc_slin(n,t))<S.assets(1)
        c_slin(n,t)=a_slin(n,t) + r_slin(n,t)*wr_slin(n,t) + n_slin(n,t)*wn_slin(n,t) + m_slin(n,t)*wh_s(n) - m_slin(n,t)*G.Inv - S.assets(1)/(1+G.r);    
    elseif (1+G.r)*(a_slin(n,t) + r_slin(n,t)*wr_slin(n,t) + n_slin(n,t)*wn_slin(n,t) + m_slin(n,t)*wh_s(n) - m_slin(n,t)*G.Inv-cc_slin(n,t))>S.assets(G.n_assets)
        c_slin(n,t)=a_slin(n,t) + r_slin(n,t)*wr_slin(n,t) + n_slin(n,t)*wn_slin(n,t) + m_slin(n,t)*wh_s(n) - m_slin(n,t)*G.Inv - S.assets(G.n_assets)/(1+G.r)+0.1;    
    else       
        c_slin(n,t)=cc_slin(n,t);
    end
    a_slin(n,t+1)= (1+G.r)*(a_slin(n,t) + r_slin(n,t)*wr_slin(n,t) + n_slin(n,t)*wn_slin(n,t) + m_slin(n,t)*wh_s(n) - m_slin(n,t)*G.Inv - c_slin(n,t));
    k_slin(n,t+1)=(gamma1*k_slin(n,t)^phi + (1-gamma1)*G.Inv^phi)^(1/phi);
        
    wh_s(n,t+1)=wh_s(n,t);   
    
    n
    
    end
    t
end



end

