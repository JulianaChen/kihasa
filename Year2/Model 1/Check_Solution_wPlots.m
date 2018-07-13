%% Plots

% 1.1 consumption (y) v assets (x) - should be increasing
% 1.2 wages (y) v experience (x) - high and low education (2 lines)
% 1.3 wages (y) v experience (x) - high and low ability (2 lines)
% 1.4 participation (y) v shock to wages (x) - for each sector (3 lines)
% 1.5 marriage (y) v husband earnings (x)
% 1.6 number of children (from state X) vs husband earnings (x)

%% reshape policy functions
for t=1:1:G.n_period-1
    for x=1:1:G.n_matstat*G.n_wrkexp
        C_rsp(:,:,x,t) = reshape(C(:,:,:,x,t),[G.n_assets,9,1,1]);
        M_rsp(:,:,x,t) = reshape(M(:,:,:,x,t),[G.n_assets,9,1,1]);
        R_rsp(:,:,x,t) = reshape(R(:,:,:,x,t),[G.n_assets,9,1,1]);
        N_rsp(:,:,x,t) = reshape(N(:,:,:,x,t),[G.n_assets,9,1,1]);
        U_rsp(:,:,x,t) = reshape(U(:,:,:,x,t),[G.n_assets,9,1,1]);
    end
end
for t=1:1:G.n_period-1
    t
    [S.SS_A', C_rsp(:,:,1,t)]
end

% husband wages
wh_mean = eta01 + eta02*(abi==2) + eta11*(edu==2) + eta12*(edu==3) + eta2*age;
wh_sd = eta03 + eta04*(abi==2) + eta21*(edu==2) + eta22*(edu==3) + eta3*age;
wh = normrnd(wh_mean,wh_sd); 

plot(age,wh)
title('Husband Wages by age');
xlabel('Age');
saveas(gcf,'wh_age.png');

