function [P]=plots(c_s,r_s,ch_s,u_s,m_s,a_s,wh_s,inv_s,wr_s,wn_s,abi,edu,type)

t=[22:40];
n_hs=sum(edu==1);
n_col2=sum(edu==2);
n_col4=sum(edu==3);

%% Work Participation

r_s2=r_s;
r_s2(1:2)=r_s2(1:2)*0.4;
n_s2=n_s;
u_s2=1-r_s2-n_s2;
% Overall
ls=r_s2+n_s2;
work_part=sum(ls,1)./G.n_pop;
work_part2=work_part;
work_part2(1)=work_part(1);
work_part2(2)=work_part(2);
work_partf=smooth(work_part2); % To use 'smooth', the following product must be licensed, installed, and enabled: Curve Fitting Toolbox
plot(t,work_partf) % this will not work bc "smooth" AND vector lengths?
axis([22 40 0 1.2])
title('Work particiption by age');
xlabel('Age');
saveas(gcf,'work_part.png');

% By work sector
work_reg=smooth(sum(r_s2,1)./G.n_pop);
work_nreg=smooth(sum(n_s2,1)./G.n_pop);
unem=smooth(sum(u_s2,1)./G.n_pop);
plot(t,work_reg)
hold on
plot(t,work_nreg)
hold off
hold on
plot(t,unem)
hold off
axis([22 40 0 1.2])
title('Sector participation by age');
xlabel('Age');
hleg6 = legend('Regular','Non-regular','Unemployment');
saveas(gcf,'sector_part.png');

% By marital status
m_s2=m_s(:,1:19);
work_marr=smooth(sum((m_s2==1).*ls,1)./sum(m_s2==1));
work_marr(isnan(work_marr))=0;
work_single=smooth(sum((m_s2==0).*ls,1)./sum(m_s2==0));
work_single(work_single>1)=1;

plot(t,work_marr)
hold on
plot(t,work_single)
hold off
axis([22 40 0 1.2])
title('Work participation by marital status');
xlabel('Age');
hleg6 = legend('Married','Single');
saveas(gcf,'work_marr.png');

% By education
work_hs=smooth(sum(ls(edu==1,1:19),1)./n_hs);
work_col2=smooth(sum(ls(edu==2,1:19),1)./n_col2);
work_col4=smooth(sum(ls(edu==3,1:19),1)./n_col4);
plot(t,work_hs)
hold on
plot(t,work_col2)
hold off
hold on
plot(t,work_col4)
hold off
axis([22 40 0 1.2])
title('Work participation by education');
xlabel('Age');
hleg6 = legend('College 2 yrs','College 4 yrs','High School');
saveas(gcf,'work_edu.png');


%% Marriage

% Overall
marriage=sum(m_s2(:,[1:19]),1)./G.n_pop;
plot(t,marriage)
axis([22 40 0 1.2])
title('Fraction married by age');
xlabel('Age');
saveas(gcf,'fmarried.png');

% By education
marr_hs=sum(m_s2(edu==1,1:19),1)./n_hs;
marr_col2=sum(m_s2(edu==2,1:19),1)./n_col2;
marr_col4=sum(m_s2(edu==3,1:19),1)./n_col4;
plot(t,marr_hs)
hold on
plot(t,marr_col2)
hold off
hold on
plot(t,marr_col4)
hold off
axis([22 40 0 1.2])
title('Fraction married by education');
xlabel('Age');
hleg6 = legend('College 4 yrs','College 2 yrs','High School');
saveas(gcf,'marr_edu.png');

% By work sector
marr_reg=smooth(sum(m_s2.*r_s,1)./sum(r_s>0));
marr_nreg=smooth(sum(m_s2.*n_s,1)./sum(n_s>0));
marr_unem=smooth(sum(m_s2.*u_s,1)./sum(u_s>0));
marr_unem(isnan(marr_unem))=0;
plot(t,marr_reg)
hold on
plot(t,marr_nreg)
hold off
hold on
plot(t,marr_unem)
hold off
axis([22 40 0 1.2])
title('Fraction married by employment');
xlabel('Age');
hleg6 = legend('Regular','Non-regular','Unemployed');
saveas(gcf,'marr_work.png');

%% Assets

% By age

asset_age=smooth((sum(a_s,1)./G.n_pop)*1000000/10000);
plot(t,asset_age(1:19))
axis([22 40 0 50000])
title('Assets by age');
xlabel('Age');
saveas(gcf,'asset_age.png');

% By education
asset_hs=smooth((sum(a_s(edu==1,1:19),1)./n_hs)*1000000/10000);
asset_col2=smooth((sum(a_s(edu==2,1:19),1)./n_col2)*1000000/10000);
asset_col4=smooth((sum(a_s(edu==3,1:19),1)./n_col4)*1000000/10000);
plot(t,asset_hs)
hold on
plot(t,asset_col2)
hold off
hold on
plot(t,asset_col4)
hold off
axis([22 40 0 50000])
title('Assets by education');
xlabel('Age');
hleg6 = legend('College 2 yrs','College 4 yrs','High School');
saveas(gcf,'asset_edu.png');

%% Wages

% By work sector
wr=wr_s.*r_s;
wn=wn_s.*n_s;
wr_sel= sum(wr,1) ./ sum(wr~=0,1);
wn_sel= sum(wn,1) ./ sum(wn~=0,1);
lw_r=smooth(log(wr_sel));
lw_n=smooth(log(wn_sel));
%wages=[mean(wr_s); mean(wn_s);wr_sel;wn_sel;lw_r;lw_n]; 

plot(t,lw_r)
hold on
plot(t,lw_n)
hold off
axis([22 40 3 6])
title('Wages by sector');
xlabel('Age');
hleg6 = legend('Regular','Non-regular');
saveas(gcf,'lwages_sector.png');

% By work sector and ability
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
lw_rl=smooth(log(wr_lsel));
lw_rh=smooth(log(wr_hsel));
lw_nl=smooth(log(wn_lsel));
lw_nh=smooth(log(wn_hsel));

plot(t,lw_rl)
hold on
plot(t,lw_rh)
hold off
axis([22 40 4 6])
title('Log wages regular sector by ability');
xlabel('Age');
hleg6 = legend('Low ability','High ability');
saveas(gcf,'wagereg_ability.png');

plot(t,lw_nl)
hold on
plot(t,lw_nh)
hold off
axis([22 40 4 6])
title('Log wages non-regular sector by ability');
xlabel('Age');
hleg6 = legend('Low ability','High ability');
saveas(gcf,'wagenreg_ability.png');

% By work sector and education
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

lw_rhs=smooth(log(wr_hssel));
lw_rcol2=smooth(log(wr_col2sel));
lw_rcol4=smooth(log(wr_col4sel));
lw_nhs=smooth(log(wn_hssel));
lw_ncol2=smooth(log(wn_col2sel));
lw_ncol4=smooth(log(wn_col4sel));

plot(t,lw_rhs)
hold on
plot(t,lw_rcol2)
hold off
hold on
plot(t,lw_rcol4)
hold off
axis([22 40 4 7])
title('Log wages regular sector by education');
xlabel('Age');
hleg6 = legend('High School','College 2 yrs','College 4 yrs');
saveas(gcf,'wagereg_edu.png');

plot(t,lw_nhs)
hold on
plot(t,lw_ncol2)
hold off
hold on
plot(t,lw_ncol4)
hold off
axis([22 40 4 7])
title('Log wages non-regular sector by education');
xlabel('Age');
hleg6 = legend('High School','College 2 yrs','College 4 yrs');
saveas(gcf,'wagenreg_edu.png');

lwage_edu=[lw_rhs;lw_rcol2;lw_rcol4;lw_nhs;lw_ncol2;lw_ncol4];

%% Husband Wages

% Overall
wh_s2=wh_s;
wh_sel= sum(wh_s2,1) ./sum(wh_s2~=0,1);
lhw= smooth(log(wh_sel)); 

plot(t,lhw)
axis([22 40 3 6])
title('Husband Wages by sector');
xlabel('Age');
saveas(gcf,'hw_lwages_sector.png');

% By education
wh_hs=wh_s2.*(edu==1);
wh_col2=wh_s2.*(edu==2);
wh_col4=wh_s2.*(edu==3);

wh_hssel= sum(wh_hs,1) ./ sum(wh_hs~=0,1);
wh_col2sel= sum(wh_col2,1) ./ sum(wh_col2~=0,1);
wh_col4sel= sum(wh_col4,1) ./ sum(wh_col4~=0,1);

lhw_hs=smooth(log(wh_hssel));
lhw_col2=smooth(log(wh_col2sel));
lhw_col4=smooth(log(wh_col4sel));

plot(t,lhw_hs)
hold on
plot(t,lhw_col2)
hold off
hold on
plot(t,lhw_col4)
hold off
axis([22 40 4 7])
title('Log husband wages by womens education');
xlabel('Age');
hleg6 = legend('High School','College 2 yrs','College 4 yrs');
saveas(gcf,'hwage_edu.png');

hwage_edu=[lhw_hs;lhw_col2;lhw_col4];

%% Child Investment

% Overall
inv_s2=inv_s;
inv_sel= sum(inv_s2,1) ./sum(inv_s2~=0,1);
inv= smooth(log(inv_sel)); 

plot(t,inv)
axis([22 40 3 6])
title('Child investment by womens age');
xlabel('Age');
saveas(gcf,'child_inv_age.png');

% By education
inv_hs=inv_s2.*(edu==1);
inv_col2=inv_s2.*(edu==2);
inv_col4=inv_s2.*(edu==3);

inv_hssel= sum(inv_hs,1) ./ sum(inv_hs~=0,1);
inv_col2sel= sum(inv_col2,1) ./ sum(inv_col2~=0,1);
inv_col4sel= sum(inv_col4,1) ./ sum(inv_col4~=0,1);

linv_hs=smooth(log(inv_hssel));
linv_col2=smooth(log(inv_col2sel));
linv_col4=smooth(log(inv_col4sel));

plot(t,linv_hs)
hold on
plot(t,linv_col2)
hold off
hold on
plot(t,linv_col4)
hold off
axis([22 40 4 7])
title('Log child investment by womens education');
xlabel('Age');
hleg6 = legend('High School','College 2 yrs','College 4 yrs');
saveas(gcf,'child_inv_edu.png');

inv_edu=[linv_hs;linv_col2;linv_col4];

end