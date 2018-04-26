function [P]=plots(c_slin,r_slin,n_slin,u_slin,m_slin,a_slin,k_slin,wr_slin,wn_slin,abi,edu,type)

t=[22:40];
n_hs=sum(edu==1);
n_col2=sum(edu==2);
n_col4=sum(edu==3);

%% Work Participation

r_slin2=r_slin;
r_slin2(1:2)=r_slin2(1:2)*0.4;
n_slin2=n_slin;
u_slin2=1-r_slin2-n_slin2;
% Overall
ls_slin=r_slin2+n_slin2;
work_part=sum(ls_slin,1)./G.n_pop;
work_part2=work_part;
work_part2(1)=work_part(1);
work_part2(2)=work_part(2);
work_partf=smooth(work_part2);
plot(t,work_partf)
axis([22 40 0 1.2])
title('Work particiption by age');
xlabel('Age');
saveas(gcf,'work_part.png');

% By work sector
work_reg=smooth(sum(r_slin2,1)./G.n_pop);
work_nreg=smooth(sum(n_slin2,1)./G.n_pop);
unem=smooth(sum(u_slin2,1)./G.n_pop);
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
m_slin2=m_slin(:,1:19);
work_marr=smooth(sum((m_slin2==1).*ls_slin,1)./sum(m_slin2==1));
work_marr(isnan(work_marr))=0;
work_single=smooth(sum((m_slin2==0).*ls_slin,1)./sum(m_slin2==0));
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
work_hs=smooth(sum(ls_slin(edu==1,1:19),1)./n_hs);
work_col2=smooth(sum(ls_slin(edu==2,1:19),1)./n_col2);
work_col4=smooth(sum(ls_slin(edu==3,1:19),1)./n_col4);
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
marriage=sum(m_slin2(:,[1:19]),1)./G.n_pop;
plot(t,marriage)
axis([22 40 0 1.2])
title('Fraction married by age');
xlabel('Age');
saveas(gcf,'fmarried.png');

% By education
marr_hs=sum(m_slin2(edu==1,1:19),1)./n_hs;
marr_col2=sum(m_slin2(edu==2,1:19),1)./n_col2;
marr_col4=sum(m_slin2(edu==3,1:19),1)./n_col4;
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

% By work

marr_reg=smooth(sum(m_slin2.*r_slin,1)./sum(r_slin>0));
marr_nreg=smooth(sum(m_slin2.*n_slin,1)./sum(n_slin>0));
marr_unem=smooth(sum(m_slin2.*u_slin,1)./sum(u_slin>0));
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

asset_age=smooth((sum(a_slin,1)./G.n_pop)*1000000/10000);
plot(t,asset_age(1:19))
axis([22 40 0 50000])
title('Assets by age');
xlabel('Age');
saveas(gcf,'asset_age.png');

% By education
asset_hs=smooth((sum(a_slin(edu==1,1:19),1)./n_hs)*1000000/10000);
asset_col2=smooth((sum(a_slin(edu==2,1:19),1)./n_col2)*1000000/10000);
asset_col4=smooth((sum(a_slin(edu==3,1:19),1)./n_col4)*1000000/10000);
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

% By sector

wr=wr_slin.*r_slin;
wn=wn_slin.*n_slin;
wr_sel= sum(wr,1) ./ sum(wr~=0,1);
wn_sel= sum(wn,1) ./ sum(wn~=0,1);
lw_r=smooth(log(wr_sel));
lw_n=smooth(log(wn_sel));
%wages=[mean(wr_slin); mean(wn_slin);wr_sel;wn_sel;lw_r;lw_n]; 

plot(t,lw_r)
hold on
plot(t,lw_n)
hold off
axis([22 40 3 6])
title('Wages by sector');
xlabel('Age');
hleg6 = legend('Regular','Non-regular');
saveas(gcf,'lwages_sector.png');

% By sector and ability
abi_l=(abi==1);
abi_h=(abi==2);
wr_l=wr_slin.*r_slin.*abi_l;
wr_h=wr_slin.*r_slin.*abi_h;
wn_l=wn_slin.*n_slin.*abi_l;
wn_h=wn_slin.*n_slin.*abi_h;
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

% By sector and education

wr_hs=wr_slin.*r_slin.*(edu==1);
wr_col2=wr_slin.*r_slin.*(edu==2);
wr_col4=wr_slin.*r_slin.*(edu==3);
wn_hs=wn_slin.*n_slin.*(edu==1);
wn_col2=wn_slin.*n_slin.*(edu==2);
wn_col4=wn_slin.*n_slin.*(edu==3);
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

end