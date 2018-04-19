clear all; clc;

%% choose aproximation
approx = 'linear';

%% choose output version
version = '_wVstar';

% without V*
%datafile = strcat('solution_',approx,'.mat')

% with V*
datafile = strcat('solution_',approx,'2.mat')


%% load data solution
load(datafile)

%% policy functions as functions of assets
filename = strcat(approx,'_assets',version,'.xls')

% fix other variables
childK = 1;
Hwages = 1;
shockr = 1;
shockn = 1;
type = 1;

% period 1
period = 1

for i = 1:20 % experience & marriage
    V_assets_1(:,i) = V_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
    C_assets_1(:,i) = C_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
    M_assets_1(:,i) = M_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
    R_assets_1(:,i) = R_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
    N_assets_1(:,i) = N_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
    U_assets_1(:,i) = U_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
end

A1 = [unique(S.SS_A),V_assets_1,C_assets_1,M_assets_1,R_assets_1,N_assets_1,U_assets_1];
xlswrite(filename,A1,'period1')

% period 9
period = 9;

for i = 1:20 % experience & marriage
    V_assets_9(:,i) = V_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
    C_assets_9(:,i) = C_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
    M_assets_9(:,i) = M_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
    R_assets_9(:,i) = R_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
    N_assets_9(:,i) = N_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
    U_assets_9(:,i) = U_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
end

A9 = [unique(S.SS_A),V_assets_9,C_assets_9,M_assets_9,R_assets_9,N_assets_9,U_assets_9];
xlswrite(filename,A9,'period9')

% period 19
period = 19;

for i = 1:20 % experience & marriage
    V_assets_19(:,i) = V_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
    C_assets_19(:,i) = C_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
    M_assets_19(:,i) = M_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
    R_assets_19(:,i) = R_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
    N_assets_19(:,i) = N_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
    U_assets_19(:,i) = U_lin(childK,:,Hwages,shockr,shockn,i,period,type)';
end

A19 = [unique(S.SS_A),V_assets_19,C_assets_19,M_assets_19,R_assets_19,N_assets_19,U_assets_19];
xlswrite(filename,A19,'period19')

%% policy functions as functions of husband's earnings
filename = strcat(approx,'_hwages',version,'.xls')

% fix other variables
assets = 5;
childK = 1;
shockr = 1;
shockn = 1;
type = 1;

% period 1
period = 1;

for i = 1:20 % experience & marriage
    V_hwages_1(:,i) = V_lin(childK,assets,:,shockr,shockn,i,period,type);
    C_hwages_1(:,i) = C_lin(childK,assets,:,shockr,shockn,i,period,type);
    M_hwages_1(:,i) = M_lin(childK,assets,:,shockr,shockn,i,period,type);
    R_hwages_1(:,i) = R_lin(childK,assets,:,shockr,shockn,i,period,type);
    N_hwages_1(:,i) = N_lin(childK,assets,:,shockr,shockn,i,period,type);
    U_hwages_1(:,i) = U_lin(childK,assets,:,shockr,shockn,i,period,type);
end

H1 = [unique(S.SS_H),V_hwages_1,C_hwages_1,M_hwages_1,R_hwages_1,N_hwages_1,U_hwages_1];
xlswrite(filename,H1,'period1')

% period 9
period = 9;

for i = 1:20 % experience & marriage
    V_hwages_9(:,i) = V_lin(childK,assets,:,shockr,shockn,i,period,type);
    C_hwages_9(:,i) = C_lin(childK,assets,:,shockr,shockn,i,period,type);
    M_hwages_9(:,i) = M_lin(childK,assets,:,shockr,shockn,i,period,type);
    R_hwages_9(:,i) = R_lin(childK,assets,:,shockr,shockn,i,period,type);
    N_hwages_9(:,i) = N_lin(childK,assets,:,shockr,shockn,i,period,type);
    U_hwages_9(:,i) = U_lin(childK,assets,:,shockr,shockn,i,period,type);
end

H9 = [unique(S.SS_H),V_hwages_9,C_hwages_9,M_hwages_9,R_hwages_9,N_hwages_9,U_hwages_9];
xlswrite(filename,H9,'period9')

% period 19
period = 19;

for i = 1:20 % experience & marriage
    V_hwages_19(:,i) = V_lin(childK,assets,:,shockr,shockn,i,period,type);
    C_hwages_19(:,i) = C_lin(childK,assets,:,shockr,shockn,i,period,type);
    M_hwages_19(:,i) = M_lin(childK,assets,:,shockr,shockn,i,period,type);
    R_hwages_19(:,i) = R_lin(childK,assets,:,shockr,shockn,i,period,type);
    N_hwages_19(:,i) = N_lin(childK,assets,:,shockr,shockn,i,period,type);
    U_hwages_19(:,i) = U_lin(childK,assets,:,shockr,shockn,i,period,type);
end

H19 = [unique(S.SS_H),V_hwages_19,C_hwages_19,M_hwages_19,R_hwages_19,N_hwages_19,U_hwages_19];
xlswrite(filename,H19,'period19')

%% consumption as a function of child's human capital
filename = strcat(approx,'_childK',version,'.xls')

% fix other variables
assets = 5;
Hwages = 1;
shockr = 1;
shockn = 1;
type = 1;

% period 1
period = 1;

for i = 1:20
    V_childK_1(:,i) = V_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
    C_childK_1(:,i) = C_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
    M_childK_1(:,i) = M_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
    R_childK_1(:,i) = R_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
    N_childK_1(:,i) = N_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
    U_childK_1(:,i) = U_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
end

K1 = [unique(S.SS_K),V_childK_1,C_childK_1,M_childK_1,R_childK_1,N_childK_1,U_childK_1];
xlswrite(filename,K1,'period1')

% period 9
period = 9;

for i = 1:20
    V_childK_9(:,i) = V_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
    C_childK_9(:,i) = C_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
    M_childK_9(:,i) = M_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
    R_childK_9(:,i) = R_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
    N_childK_9(:,i) = N_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
    U_childK_9(:,i) = U_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
end

K9 = [unique(S.SS_K),V_childK_9,C_childK_9,M_childK_9,R_childK_9,N_childK_9,U_childK_9];
xlswrite(filename,K9,'period9')

% period 19
period = 19;

for i = 1:20
    V_childK_19(:,i) = V_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
    C_childK_19(:,i) = C_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
    M_childK_19(:,i) = M_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
    R_childK_19(:,i) = R_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
    N_childK_19(:,i) = N_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
    U_childK_19(:,i) = U_lin(:,assets,Hwages,shockr,shockn,i,period,type)';
end

K19 = [unique(S.SS_K),V_childK_19,C_childK_19,M_childK_19,R_childK_19,N_childK_19,U_childK_19];
xlswrite(filename,K19,'period19')