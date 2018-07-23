function JTN=mdf(params,G,S,abi,edu,type,types,data)  %JTN Objective Function to be minimized

params_all = 
S = sspace(params0,G);

parfor z=1:1:G.n_incond
    z
    [C(:,:,:,:,:,z),M(:,:,:,:,:,z),R(:,:,:,:,:,z),N(:,:,:,:,:,z),U(:,:,:,:,:,z)]= solution(G,types(z,1),types(z,2),S,params);     toc

end

[c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s] = simulation(params,G,S,abi,edu,type,C,M,R,N,U);

Sbeta=moments(c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s);

% Data Moments

Dbeta=data(:,1);
Dse=srqt(data(:,2));

JTN=nansum(((Dbeta-Sbeta)./Dse).^2);


end