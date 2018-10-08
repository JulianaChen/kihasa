function [c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,exp_s] = simulation(eparams,G,S,P,abi,edu,type,C,M,R,N,U)

%% Parameters

alpha01_r=eparams(13);
alpha01_n=eparams(14);
alpha02_r=eparams(15);
alpha02_n=eparams(16);
alpha11_r=eparams(17);
alpha11_n=eparams(18);
alpha12_r=eparams(19);
alpha12_n=eparams(20);
alpha2_r=eparams(21);
alpha2_n=eparams(22);

sigma_r=eparams(23);
sigma_n=eparams(24);
sigma_i=eparams(25);

phi10=P.phi10;
phi11=P.phi11;
phi12=P.phi12;
phi13=P.phi13;
phi20=P.phi20;
phi21=P.phi21;
phi22=P.phi22;
phi23=P.phi23;
phi30=P.phi30;
phi31=P.phi31;
phi32=P.phi32;
phi33=P.phi33;

eta01=P.eta01;
eta02=P.eta02;
eta11=P.eta11;
eta12=P.eta12;
eta2=P.eta2;
iota01=P.iota01;
iota02=P.iota02;
iota11=P.iota11;
iota12=P.iota12;
iota2=P.iota2;

rho01 = 5112; % assets at age 18-20 (mean)
rho02 = 2464; % assets at age 18-20 (mean)
rho11 = 2133; % assets at age 18-20 (mean)
rho12 = 8671; % assets at age 18-20 (mean)
rho03 = 2628; % assets at age 18-20 (sd)
rho04 = -1479; % assets at age 18-20 (sd)
rho21 = 158.9; % assets at age 18-20 (sd)
rho22 = -16.16; % assets at age 18-20 (sd)

%% Initial Conditions

exp_s = zeros(G.n_pop,G.n_period); 
m_s = zeros(G.n_pop,G.n_period); 
ch_s = zeros(G.n_pop,G.n_period); 
a_s = zeros(G.n_pop,G.n_period); 

%Run regressions of mean and sd of assets at age 18-20 on education/ability
a_mean = rho01 + rho02*(abi==2) + rho11*(edu==2) + rho12*(edu==3);
a_sd = rho03 + rho04*(abi==2) + rho21*(edu==2) + rho22*(edu==3);
a_s(:,1) = 0; %normrnd(a_mean, a_sd); % draw from distribution

%Ensure initial assets are within the bounds
for n=1:G.n_pop
    if a_s(n,1)<min(S.SS_A)
        a_s(n,1)= min(S.SS_A);
    end
    if a_s(:,1)>max(S.SS_A)
        a_s(:,1)= max(S.SS_A);
    end
end

%Wider Assets Vector
A_wide = S.SS_A;
A_wide(1) = min(S.SS_A)-50000;
A_wide(G.n_assets)= max(S.SS_A)+50000;

% Use the basis for assets and shocks from s_space.m
T_sim=kron(S.T_A,kron(S.Teps_i,kron(S.Teps_r,S.Teps_n)));
Den=kron(S.T2_A,kron(S.T2eps_i,kron(S.T2eps_n,S.T2eps_r)));

% Reshape the policy functions
for z = 1:1:G.n_incond
    for t = 1:1:G.n_period-1
        for x = 1:1:(G.n_matstat*G.n_wrkexp)
            C_rsp(:,x,t,z) = reshape(C(:,:,:,:,x,t,z),[],1);
            R_rsp(:,x,t,z) = reshape(R(:,:,:,:,x,t,z),[],1);
            N_rsp(:,x,t,z) = reshape(N(:,:,:,:,x,t,z),[],1);
            U_rsp(:,x,t,z) = reshape(U(:,:,:,:,x,t,z),[],1);
            M_rsp(:,x,t,z) = reshape(M(:,:,:,:,x,t,z),[],1);
            CH_rsp(:,x,t,z)= reshape(CH(:,:,:,:,x,t,z),[],1);
        end
    end
end

%%% Polynomial Bases and Derivatives %%%% 

for z = 1:1:G.n_incond
    for t = 1:1:G.n_period-1
        for x = 1:1:(G.n_matstat*G.n_wrkexp)
            NumC(:,x,t,z) = C_rsp(:,x,t,z)'*T_sim;  
            alpC(:,x,t,z) = NumC(:,x,t,z)./Den;
            NumR(:,x,t,z) = R_rsp(:,x,t,z)'*T_sim; 
            alpR(:,x,t,z) = NumR(:,x,t,z)./Den;
            NumN(:,x,t,z) = N_rsp(:,x,t,z)'*T_sim; 
            alpN(:,x,t,z) = NumN(:,x,t,z)./Den;
            NumU(:,x,t,z) = U_rsp(:,x,t,z)'*T_sim; 
            alpU(:,x,t,z) = NumU(:,x,t,z)./Den;       
            NumM(:,x,t,z) = M_rsp(:,x,t,z)'*T_sim; 
            alpM(:,x,t,z) = NumM(:,x,t,z)./Den;
            NumCH(:,x,t,z)= CH_rsp(:,x,t,z)'*T_sim;
            alpCH(:,x,t,z)= NumCH(:,x,t,z)./Den;
        end
	end
end

%% Loop over all periods/individuals

for n=1:G.n_pop
for t=1:1:G.n_period-1
    
    %for n=1:n%1:G.n_pop
   
%     % Obtain sector shocks
%     epssim_r(n,t)=sqrt(2)*G.Eps(1,n,t)'*sigma_r;
%     epssim_n(n,t)=sqrt(2)*G.Eps(2,n,t)'*sigma_n;
%     
%     if epssim_r(n,t)<S.eps_r(1) || epssim_r(n,t)>S.eps_r(3)
%         eps_rg=epssim_r(n,t)*[-1;0;1];
%     else 
%         eps_rg=S.eps_r;
%     end
%     if epssim_n(n,t)<S.eps_n(1) || epssim_n(n,t)>S.eps_n(3)
%         eps_ng=epssim_n(n,t)*[-1;0;1];
%     else 
%         eps_ng=S.eps_n;
%     end

    % draw the shocks and adjust the grid if necessary
     
    epssim_i(n,t)=sqrt(2)*G.Eps(1,n,t)'*sigma_i;
    
    if epssim_i(n,t)<S.eps_i(1) || epssim_i(n,t)>S.eps_i(3)
        eps_ig=epssim_i(n,t)*[-1;0;1];
    else 
        eps_ig=S.eps_i;
    end
    
    epssim_r(n,t)=sqrt(2)*G.Eps(1,n,t)'*sigma_r;
    
    if epssim_r(n,t)<S.eps_r(1) || epssim_r(n,t)>S.eps_r(3)
        eps_rg=epssim_r(n,t)*[-1;0;1];
    else 
        eps_rg=S.eps_r;
    end
    
    epssim_n(n,t)=sqrt(2)*G.Eps(1,n,t)'*sigma_n;
    
    if epssim_n(n,t)<S.eps_n(1) || epssim_n(n,t)>S.eps_n(3)
        eps_ng=epssim_n(n,t)*[-1;0;1];
    else 
        eps_ng=S.eps_n;
    end
    
    % Calculate wages
    wr_s(n,t) = exp(alpha01_r + alpha02_r*(abi(n)==2) + alpha11_r*(edu(n)==2) + alpha12_r*(edu(n)==3) + alpha2_r*log(1+exp_s(n,t)) + epssim_r(n,t)); 
    wn_s(n,t) = exp(alpha01_n + alpha02_n*(abi(n)==2) + alpha11_n*(edu(n)==2) + alpha12_n*(edu(n)==3) + alpha2_n*log(1+exp_s(n,t)) + epssim_n(n,t));

    % Locate in the experience/marriage/children vector   
    if m_s(n,t)==0 
        x=min(30, 20 + exp_s(n,t)+1);  
    elseif m_s(n,t)==1
        if ch_s(n,t)==1
           x=min(10, exp_s(n,t)+1); 
        elseif ch_s(n,t)==2
           x=min(20, 10 + exp_s(n,t)+1);  
        end
    end
                            
%     % Optimal Choices (using linear interpolation)
%     cc_s(n,t)= interpn(A_wide, eps_rg, eps_ng, C(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t));
%     rr_s(n,t)= interpn(A_wide, eps_rg, eps_ng, R(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t)); 
%     nn_s(n,t)= interpn(A_wide, eps_rg, eps_ng, N(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t)); 
%     uu_s(n,t)= interpn(A_wide, eps_rg, eps_ng, U(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t)); 
 
    % Evaluate the basis in the specific values of assets and wage shocks
    % for individual n at period t
    
    T_eps_i=chebpoly_base(G.Ne-1, 2*(epssim_i(n,t) - eps_ig(1))/(eps_ig(G.Ne)-eps_ig(1)) - 1);
    T_eps_r=chebpoly_base(G.Ne-1, 2*(epssim_r(n,t) - eps_rg(1))/(eps_rg(G.Ne)-eps_rg(1)) - 1);
    T_eps_n=chebpoly_base(G.Ne-1, 2*(epssim_n(n,t) - eps_ng(1))/(eps_ng(G.Ne)-eps_ng(1)) - 1);  
    T_a=chebpoly_base(S.nA+1, 2*(a_s(n,t) - S.SS_A(1))/(S.SS_A(G.n_assets)-S.SS_A(1)) - 1);
    
    % New basis
    T_s=kron(T_a,kron(T_eps_i,kron(T_eps_r,T_eps_n)));  
    
    % Interpolate policy functions
    cc_s(n,t)=sum(alpC(:,x,t,type(n)).*T_s',1);
    rr_s(n,t)=sum(alpR(:,x,t,type(n)).*T_s',1);
    nn_s(n,t)=sum(alpN(:,x,t,type(n)).*T_s',1); 
    uu_s(n,t)=sum(alpU(:,x,t,type(n)).*T_s',1);     
    mm_s(n,t)=sum(alpM(:,x,t,type(n)).*T_s',1);    
    chh_s(n,t)=sum(alpCH(:,x,t,type(n)).*T_s',1);
    
    % add unemployed 

    [v, Ind] = max([rr_s(n,t), nn_s(n,t), uu_s(n,t)]);
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
    
    % Probability of Children
%     prob_2kids_r = normcdf(phi10 + phi11*(edu(n)==2) + phi12*(edu(n)==3) + phi13*exp_s(n,t));
%     prob_2kids_n = normcdf(phi20 + phi21*(edu(n)==2) + phi22*(edu(n)==3) + phi23*exp_s(n,t));
%     prob_2kids_u = normcdf(phi30 + phi31*(edu(n)==2) + phi32*(edu(n)==3) + phi33*exp_s(n,t));
       
    if m_s(n,t)==0 
       %marr(n,t)=interpn(A_wide, eps_rg, eps_ng,M(:,:,:,x,t,type(n)),a_s(n,t),epssim_r(n,t),epssim_n(n,t));
       if mm_s(n,t)<=0.5 
          m_s(n,t+1)=0;
          ch_s(n,t+1)=0;
       else
          m_s(n,t+1)=1;
          ch_s(n,t+1)=1;
       end
    else
        m_s(n,t+1)=m_s(n,t);
        if chh_s(n,t) < 0.5
            ch_s(n,t+1)=ch_s(n,t);
        else
            ch_s(n,t+1)=ch_s(n,t)+1;
        end
%         if ch_s(n,t)==1 && r_s(n,t)==1
%            if prob_2kids_r<0.5
%                ch_s(n,t+1)=1;
%            else
%                ch_s(n,t+1)=2;
%            end
%         elseif ch_s(n,t)==1 && n_s(n,t)==1
%            if prob_2kids_n<0.5
%                ch_s(n,t+1)=1;
%            else
%                ch_s(n,t+1)=2;
%            end
%         elseif ch_s(n,t)==1 && u_s(n,t)==1
%            if prob_2kids_u<0.5
%                ch_s(n,t+1)=1;
%            else
%                ch_s(n,t+1)=2;
%            end
%         elseif ch_s(n,t)==2
%            ch_s(n,t+1)=2;       
%         end
    end
         
    % Find age of woman
    age = 18*(edu(n)==1) + 20*(edu(n)==2) + 22*(edu(n)==3) + t;

    % Draw husband wages
	wh_mean = eta01 + eta02*(abi(n)==2) + eta11*(edu(n)==2) + eta12*(edu(n)==3) + eta2*age;
	wh_sd = 0.7022257; %eta03 + eta04*(abi(n)==2) + eta21*(edu(n)==2) + eta22*(edu(n)==3) + eta3*age;
	wh_s(n,t) = normrnd(wh_mean,wh_sd);

    % Draw child investments
	Inv_mean = iota01 + iota02*(abi(n)==2) + iota11*(edu(n)==2) + iota12*(edu(n)==3) + iota2*age;
	Inv_sd = 0.9270494; %iota03 + iota04*(abi(n)==2) + iota21*(edu(n)==2) + iota22*(edu(n)==3) + iota3*age;
	inv_s(n,t) = normrnd(Inv_mean,Inv_sd);

    % Optimal simulated consumption with validations
    a_tmw(n,t)= (1+G.r)*(a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*exp(wh_s(n,t)) - ch_s(n,t)*exp(inv_s(n,t)) - cc_s(n,t));
    
    if a_tmw(n,t)<S.SS_A(1) % or make into -10
       c_s(n,t)=max(0,a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*exp(wh_s(n,t)) - ch_s(n,t)*exp(inv_s(n,t))- S.SS_A(1)/(1+G.r));  
      
    elseif a_tmw(n,t)>S.SS_A(G.n_assets)
       c_s(n,t)=max(0,a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*exp(wh_s(n,t)) - ch_s(n,t)*exp(inv_s(n,t)) - S.SS_A(G.n_assets)/(1+G.r));    
        
    else       
       c_s(n,t)=cc_s(n,t);
    end

    % Transition for assets (Budget Constraint)    
    a_s(n,t+1)= (1+G.r)*(a_s(n,t) + r_s(n,t)*wr_s(n,t) + n_s(n,t)*wn_s(n,t) + m_s(n,t)*exp(wh_s(n,t)) - ch_s(n,t)*exp(inv_s(n,t)) - c_s(n,t));    
end
end

% bound children to 2
ch_s(ch_s > 2) = 2;

end