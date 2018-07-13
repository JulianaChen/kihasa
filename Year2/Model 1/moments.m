function [SBeta] = moments(c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s)

t=[22:40];
n_hs=sum(edu==1);
n_col2=sum(edu==2);
n_col4=sum(edu==3);

%age groups(18-22(t=1-6),23-28(t=7-12),29-35(t=13-19))

%% Work participation by age and sector (2 x 4 age groups)

r_s2=r_s;
n_s2=n_s;
u_s2=1-r_s2-n_s2;

% by age
ls=r_s2+n_s2;
work_part=sum(ls,1)./G.n_pop;

% By age and work sector
work_reg=sum(r_s2,1)./G.n_pop;
work_nreg=sum(n_s2,1)./G.n_pop;
unem=sum(u_s2,1)./G.n_pop;

%% Work participation by age and married (1 x 4 age groups)

m_s2=m_s(:,1:19);
work_marr=sum((m_s2==1).*ls,1)./sum(m_s2==1);
work_marr(isnan(work_marr))=0;
work_single=sum((m_s2==0).*ls,1)./sum(m_s2==0);
work_single(work_single>1)=1;

%% Work participation by age and number of children (2 x 4 age groups)

ch_s2=ch_s(:,1:19);
work_ch0=sum((ch_s2==0).*ls,1)./sum(ch_s2==0);
work_ch1=sum((ch_s2==1).*ls,1)./sum(ch_s2==1);
work_ch2=sum((ch_s2==2).*ls,1)./sum(ch_s2==2);

%% Work participation by age and education (2 x 4 age groups)

work_hs=sum(ls(edu==1,1:19),1)./n_hs;
work_col2=sum(ls(edu==2,1:19),1)./n_col2;
work_col4=sum(ls(edu==3,1:19),1)./n_col4;

%% Fraction married by age (4 age groups)

marriage=sum(m_s2(:,[1:19]),1)./G.n_pop;

%% Fraction married by age and education (2x4 age groups)

marr_hs=sum(m_s2(edu==1,1:19),1)./n_hs;
marr_col2=sum(m_s2(edu==2,1:19),1)./n_col2;
marr_col4=sum(m_s2(edu==3,1:19),1)./n_col4;

%% Fraction married by age and sector (2x4 age groups)

marr_reg=sum(m_s2.*r_s,1)./sum(r_s>0);
marr_nreg=sum(m_s2.*n_s,1)./sum(n_s>0);
marr_unem=sum(m_s2.*u_s,1)./sum(u_s>0);
marr_unem(isnan(marr_unem))=0;

%% Assets tertiles by age (2x4 age groups) %% how to do tertiles?

asset_age=(sum(a_s,1)./G.n_pop)*1000000/10000;

%% Asset tertiles by age and marital status (2x1x4 age groups) %% need to check

n_marr=sum(m_s2,1);
n_single=G.n_pop-n_marr;

asset_marr=(sum(a_s(m_s2==1),1)./n_marr)*1000000/10000;
asset_single=(sum(a_s(m_s2==0),1)./n_single)*1000000/10000;

%% Asset tertiles by age and education (2x2x4 age groups)

asset_hs=(sum(a_s(edu==1,1:19),1)./n_hs)*1000000/10000;
asset_col2=(sum(a_s(edu==2,1:19),1)./n_col2)*1000000/10000;
asset_col4=(sum(a_s(edu==3,1:19),1)./n_col4)*1000000/10000;

%% Experience by age (1x4 age groups)
%% Experience by marital status and education (2x3)
%% Transitions by general experience (3x3x5 bins 0,1,2-4,5-9,10+)
%% Transitions between sectors by marital status (3x3x2)
%% Transitions Between sectors by education (3x3x3)
%% Husband earnings

% Mean and variance of earnings by wife?s age (2x4 age groups) %% variance?

wh_s2=wh_s;
lhw= sum(wh_s2,1) ./sum(wh_s2~=0,1); 

% Mean earnings by wife?s education (2 x 3)

wh_hs=wh_s2.*(edu==1);
wh_col2=wh_s2.*(edu==2);
wh_col4=wh_s2.*(edu==3);

wh_hssel= sum(wh_hs,1) ./ sum(wh_hs~=0,1);
wh_col2sel= sum(wh_col2,1) ./ sum(wh_col2~=0,1);
wh_col4sel= sum(wh_col4,1) ./ sum(wh_col4~=0,1);

%% Child Investments

% Mean and variance of investment by women?s age (2x4 age groups) %% variance?

inv_s2=inv_s;
inv_sel= sum(inv_s2,1) ./sum(inv_s2~=0,1);

% Mean investment by wife?s education (2 x 3)

inv_hs=inv_s2.*(edu==1);
inv_col2=inv_s2.*(edu==2);
inv_col4=inv_s2.*(edu==3);

inv_hssel= sum(inv_hs,1) ./ sum(inv_hs~=0,1);
inv_col2sel= sum(inv_col2,1) ./ sum(inv_col2~=0,1);
inv_col4sel= sum(inv_col4,1) ./ sum(inv_col4~=0,1);

%% Wages 

% Mean wage by sector, education and ability type (2x3x2)

% By work sector
wr=wr_s.*r_s;
wn=wn_s.*n_s;
wr_sel= sum(wr,1) ./ sum(wr~=0,1);
wn_sel= sum(wn,1) ./ sum(wn~=0,1);

% Mean wage by sector and ability type (2x2)

abi_l=(abi==1);
abi_h=(abi==2);

wr_l=wr_s.*r_s.*abi_l;
wr_h=wr_s.*r_s.*abi_h;
wn_l=wn_s.*n_s.*abi_l;
wn_h=wn_s.*n_s.*abi_h;

wr_lsel= sum(wr_l,1) ./ sum(wr_l~=0,1);
wr_hsel= sum(wr_h,1) ./ sum(wr_h~=0,1);
wn_lsel= sum(wn_l,1) ./ sum(wn_l~=0,1);
wn_hsel= sum(wn_h,1) ./ sum(wn_h~=0,1);

% Mean wage by sector and education (2x3)

wr_hs=wr_s.*r_s.*(edu==1);
wr_col2=wr_s.*r_s.*(edu==2);
wr_col4=wr_s.*r_s.*(edu==3);

wn_hs=wn_s.*n_s.*(edu==1);
wn_col2=wn_s.*n_s.*(edu==2);
wn_col4=wn_s.*n_s.*(edu==3);

wr_hssel= sum(wr_hs,1) ./ sum(wr_hs~=0,1);
wr_col2sel= sum(wr_col2,1) ./ sum(wr_col2~=0,1);
wr_col4sel= sum(wr_col4,1) ./ sum(wr_col4~=0,1);

wn_hssel= sum(wn_hs,1) ./ sum(wn_hs~=0,1);
wn_col2sel= sum(wn_col2,1) ./ sum(wn_col2~=0,1);
wn_col4sel= sum(wn_col4,1) ./ sum(wn_col4~=0,1);

% Mean wage by sector and general experience (2x 5 bins 0,1,2-4,5-9,10+)

% Variance of observed earnings (2-sector-specific)


