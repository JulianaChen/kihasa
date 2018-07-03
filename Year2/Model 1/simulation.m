function [c_s,r_s,n_s,u_s,m_s,a_s,wh_s,inv_s,wr_s,wn_s] = simulation(params,G,S,abi,edu,type,C,M,R,N,U)

% Parameters CHECK THIS
alpha01_r=params(9); % wage return for the ability type, regular
alpha01_n=params(10); % wage return for the ability type, non-regular
alpha02_r=params(11); % additional return for high ability type, regular
alpha02_n=params(12); % additional return for high ability type, non-regular
alpha11_r=params(13); % wage return to 2yr college, regular
alpha11_n=params(14); % wage return to 2yr college, non-regular
alpha12_r=params(15); % wage return to 4yr college, regular
alpha12_n=params(16); % wage return to 4yr college, non-regular
alpha2_r=params(17); % wage return to general work experience, regular
alpha2_n=params(18); % wage return to general work experience, non-regular
sigma_r = params(19); % shock, regular
sigma_n = params(20); % shock, non-regular
sigma_i = params(21); % shock, unemployed

%% ADD THE PHIS FOR PROBABILITY OF HAVING KIDS
%% BRING ETAS AND IOTAS


% Initial Conditions\
a_s = zeros(G.n_pop,G.n_period); 
    a_s(:,1)= %Run regressions of mean and sd of assets at age 18-20 on education and ability, and then draw from that distribution. 
exp_s= zeros(G.n_pop,G.n_period); 
m_s= zeros(G.n_pop,G.n_period); 
n_s= zeros(G.n_pop,G.n_period); 

for t=1:1:G.n_period-1
    for n=1:1:G.n_pop
     
    epssim_r(n,t)=sqrt(2)*G.Eps(1,n,t)'*sigma_r;
    epssim_n(n,t)=sqrt(2)*G.Eps(2,n,t)'*sigma_n;
       
    % Calculate wages
    wr_s(n,t) = exp(alpha01_r + alpha02_r*(abi(n)==2) + alpha11_r*(edu(n)==2) + alpha12_r*(edu(n)==3) + alpha2_r*log(1+exp_s(n,t)) + epssim_r(n,t)); 
    wn_s(n,t) = exp(alpha01_n + alpha02_n*(abi(n)==2) + alpha11_n*(edu(n)==2) + alpha12_n*(edu(n)==3) + alpha2_n*log(1+exp_s(n,t)) + epssim_n(n,t));

  
% Locate in the experience/marriage/children vector   
if m_s(n,t)==0 
    x=min(30, 20 + exp_s(n,t)+1);  
elseif m_s(n,t)==1
	if n_s(n,t)==1
	   x=min(10, exp_s(n,t)+1); 
	elseif n_s(n,t)==2
	   x=min(20, 10 + exp_s(n,t)+1);  
	end
end
                            
% Optimal Choices
    cc_s(n,t)= interpn(S.assets, S.eps_r, S.eps_n, C(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t));
    rr_s(n,t)= interpn(S.assets, S.eps_r, S.eps_n, R(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t));
    nn_s(n,t)= interpn(S.assets, S.eps_r, S.eps_n, N(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t));
    uu_s(n,t)= interpn(S.assets, S.eps_r, S.eps_n, U(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t));
      
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
    
    prob_2kids_r = normcdf(phi10 + phi11*(edu==2) + phi12*(edu==3) + phi13*exp_s(n,t));
    prob_2kids_n = normcdf(phi20 + phi21*(edu==2) + phi22*(edu==3) + phi23*exp_s(n,t));
    prob_2kids_u = normcdf(phi30 + phi31*(edu==2) + phi32*(edu==3) + phi33*exp_s(n,t));
       
    if m_s(n,t)==0
       marr(n,t)=interpn(S.assets,S.eps_r, S.eps_n,M_(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t));
       if marr(n,t)<0.5 
          m_s(n,t+1)=0;
       else
          m_s(n,t+1)=1;
          n_s(n,t+1)=1;
       end
    else
        m_s(n,t+1)=m_s(n,t);
        if n_s(n,t)==1 && r_s(n,t)==1;
           if prob_2kids_r<0.5
               n_s(n,t+1)=1
           else
               n_s(n,t+1)=2
        elseif n_s(n,t)==1 && n_s(n,t)==1;
           if prob_2kids_n<0.5
               n_s(n,t+1)=1
           else
               n_s(n,t+1)=2
        elseif n_s(n,t)==1 && u_s(n,t)==1;
           if prob_2kids_u<0.5
               n_s(n,t+1)=1
           else
               n_s(n,t+1)=2                                        
        elseif n_s(n,t)==2;
           n_s(n,t+1)==2;       
    end


% Draw husband wages
	age = 18*(edu==1) + 20*(edu==2) + 22*(edu==3) + t;
	wh_mean = eta01 + eta02*(abi==2) + eta11*(edu==2) + eta12*(edu==3) + eta2*age;
	wh_sd = eta03 + eta04*(abi==2) + eta21*(edu==2) + eta22*(edu==3) + eta3*age;
	wh_s(n,t) = normrnd(wh_mean,wh_sd);

% Draw child investments
	Inv_mean = iota01 + iota02*(abi==2) + iota11*(edu==2) + iota12*(edu==3) + iota2*age;
	Inv_sd = iota03 + iota04*(abi==2) + iota21*(edu==2) + iota22*(edu==3) + iota3*age;
	inv_s(n,t) = normrnd(Inv_mean,Inv_sd);

% Optimal simulated consumption with validations

    a_tmw(n,t)= (1+G.r)*(a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*wh_s(n,t) - n_s(n,t)*inv_s(n,t)-cc_s(n,t))
    if a_tmw(n,t)<S.assets(1)
       c_s(n,t)=max(0,a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*wh_s(n,t) - n_s(n,t)*inv_s(n,t)- S.assets(1)/(1+G.r));  
      
    elseif a_tmw(n,t)>S.assets(G.n_assets)
       c_s(n,t)=max(0,a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*wh_s(n,t) - n_s(n,t)*inv_s(n,t) - S.assets(G.n_assets)/(1+G.r));    
        
    else       
       c_s(n,t)=cc_s(n,t);
    end

%  Transition for assets (Budget Constraint)    
    
    a_s(n,t+1)= (1+G.r)*(a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*wh_s(n,t) - n_s(n,t)*inv_s(n,t) - c_s(n,t));
   
    n    
    end
    t
end

end

