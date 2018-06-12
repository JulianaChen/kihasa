function [c_s,r_s,n_s,u_s,m_s,a_s,wh_s,inv_s,wr_s,wn_s] = simulation(params,G,S,abi,edu,type,C,M,R,N,U)

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
      a_s(:,1)= %assign an initial value based on education and ability%
exp_s= zeros(G.n_pop,G.n_period); 
m_s= zeros(G.n_pop,G.n_period); 
c_s = zeros(G.n_pop,G.n_period); 


for t=1:1:G.n_period-1
    for n=1:1:G.n_pop
     
    epssim_r(n,t)=sqrt(2)*G.Eps(1,n,t)'*sigma_r;
    epssim_n(n,t)=sqrt(2)*G.Eps(2,n,t)'*sigma_n;
       
    % Calculate wages
    wr_s(n,t) = exp(alpha01_r + alpha02_r*(abi(n)==2) + alpha11_r*(edu(n)==2) + alpha12_r*(edu(n)==3) + alpha2_r*log(1+exp_s(n,t)) + epssim_r(n,t)); 
    wn_s(n,t) = exp(alpha01_n + alpha02_n*(abi(n)==2) + alpha11_n*(edu(n)==2) + alpha12_n*(edu(n)==3) + alpha2_n*log(1+exp_s(n,t)) + epssim_n(n,t));
   
    asset_g=S.assets;
%     if a_s(n,t)>S.assets(G.n_assets)
%         asset_g(G.n_assets)=a_slin(n,t);
%     end

% Locate in the experience/marriage vector    
    x=min(20,(1-m_s(n,t))*10 + exp_s(n,t)+1);        
                   
% Optimal Choices
    cc_s(n,t)= interpn(asset_g, S.eps_r, S.eps_n, C(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t));
    rr_s(n,t)= interpn(asset_g, S.eps_r, S.eps_n, R(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t));
    nn_s(n,t)= interpn(asset_g, S.eps_r, S.eps_n, N(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t));
    uu_s(n,t)= interpn(asset_g, S.eps_r, S.eps_n, U(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t));
      
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
       marr(n,t)=interpn(asset_g,S.eps_r, S.eps_n,M_(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t));
       if marr(n,t)<0.5 
          m_s(n,t+1)=0;
       else
          m_s(n,t+1)=1;
       end
    else
        m_s(n,t+1)=m_s(n,t);
    end


% Draw Husband wage wh_s(n,t)


% Draw child investment inv_s(n,t)



% Optimal simulated consumption with validations

    if (1+G.r)*(a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*wh_s(n,t) - m_s(n,t)*inv_s(n,t)-cc_s(n,t))<S.assets(1)
       c_s(n,t)=a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*wh_s(n,t) - m_s(n,t)*inv_s(n,t)- S.assets(1)/(1+G.r);    
    elseif (1+G.r)*(a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*wh_s(n,t) - m_s(n,t)*inv_s(n,t)-cc_s(n,t))>S.assets(G.n_assets)
        c_s(n,t)=a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*wh_s(n,t) - m_s(n,t)*inv_s(n,t) - S.assets(G.n_assets)/(1+G.r);    
    else       
        c_s(n,t)=cc_s(n,t);
    end

%  Transition for assets (Budget Constraint)    
    
    a_s(n,t+1)= (1+G.r)*(a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*wh_s(n,t) - m_s(n,t)*inv_s(n,t) - c_s(n,t));

    
    n
    
    end
    t
end



end

