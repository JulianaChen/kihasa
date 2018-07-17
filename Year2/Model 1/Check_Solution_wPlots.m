%% Plots

% 1.1 consumption (y) v assets (x) - should be increasing
% 1.2 wages (y) v experience (x) - high and low education (2 lines)
% 1.3 wages (y) v experience (x) - high and low ability (2 lines)
% 1.4 participation (y) v shock to wages (x) - for each sector (3 lines)
% 1.5 marriage (y) v husband earnings (x)
% 1.6 number of children (from state X) vs husband earnings (x)

%% pick type

C1 = C(:,:,:,:,:,1);
M1 = M(:,:,:,:,:,1);
R1 = R(:,:,:,:,:,1);
N1 = N(:,:,:,:,:,1);
U1 = U(:,:,:,:,:,1);

C4 = C(:,:,:,:,:,4);
M4 = M(:,:,:,:,:,4);
R4 = R(:,:,:,:,:,4);
N4 = N(:,:,:,:,:,4);
U4 = U(:,:,:,:,:,4);

%% reshape policy functions
for t=1:1:G.n_period-1
    for x=1:1:G.n_matstat*G.n_wrkexp
        C_rsp(:,:,x,t) = reshape(C1(:,:,:,x,t),[G.n_assets,9,1,1]);
        M_rsp(:,:,x,t) = reshape(M1(:,:,:,x,t),[G.n_assets,9,1,1]);
        R_rsp(:,:,x,t) = reshape(R1(:,:,:,x,t),[G.n_assets,9,1,1]);
        N_rsp(:,:,x,t) = reshape(N1(:,:,:,x,t),[G.n_assets,9,1,1]);
        U_rsp(:,:,x,t) = reshape(U1(:,:,:,x,t),[G.n_assets,9,1,1]);
    end
end

%find married women
for t=1:6:G.n_period-1
    for x=1:4:G.n_matstat*G.n_wrkexp
        t
        x
        sum(M_rsp(:,:,x,t))
    end
end

for t=1:6:G.n_period-1
    for x=1:4:G.n_matstat*G.n_wrkexp
        t
        x
        [S.SS_A', C_rsp(:,:,x,t)]
        [S.SS_A', M_rsp(:,:,x,t)]
        [S.SS_A', R_rsp(:,:,x,t)]
        [S.SS_A', N_rsp(:,:,x,t)]
        [S.SS_A', U_rsp(:,:,x,t)]
    end
end

%print simulated assets and marriage
xlswrite('check_July17.xls',[type,a_s(:,1:6:G.n_period-1),m_s(:,1:6:G.n_period-1)],'simulated')
%one type (instead of all rows :, pick only rows where (type==1))
type1 = [a_s(type==1,1:6:G.n_period-1),m_s(type==1,1:6:G.n_period-1)];


%print married women
for t=1:6:G.n_period-1
    for x=1:4:G.n_matstat*G.n_wrkexp
        t
        x
        sheetname = strcat('marr_',num2str(x),'_',num2str(t));
        xlswrite('check.xls',[t,x,sum(M_rsp(:,:,x,t))],sheetname)
    end
end

x=5;
t=5;
        [S.SS_A', C_rsp(:,:,x,t)]
        [S.SS_A', M_rsp(:,:,x,t)]
        [S.SS_A', R_rsp(:,:,x,t)]
        [S.SS_A', N_rsp(:,:,x,t)]
        [S.SS_A', U_rsp(:,:,x,t)]

% xlswrite('check.xls',[S.SS_A', C_rsp(:,:,x,t)],strcat('cons',x,t))
% xlswrite('check.xls',[S.SS_A', M_rsp(:,:,x,t)],strcat('marr',x,t))
% xlswrite('check.xls',[S.SS_A', R_rsp(:,:,x,t)],strcat('wreg',x,t))
% xlswrite('check.xls',[S.SS_A', N_rsp(:,:,x,t)],strcat('nreg',x,t))
% xlswrite('check.xls',[S.SS_A', U_rsp(:,:,x,t)],strcat('unem',x,t))
