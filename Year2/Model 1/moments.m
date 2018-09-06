function [SBeta] = moments(c_s,r_s,n_s,u_s,m_s,ch_s,a_s,wh_s,inv_s,wr_s,wn_s,exp_s)

t=[1:19];
n_hs=sum(edu==1);
n_col2=sum(edu==2);
n_col4=sum(edu==3);

%% age groups(18-22(t=1-6),23-28(t=7-12),29-35(t=13-19))
% age = 18*(edu==1) + 20*(edu==2) + 22*(edu==3) + t - 1;
%  
% age2=age;
% 
% for t=1:1:G.n_period-1;
%    for n=1:1:G.n_pop;
%         
%     if age2(n,t) < 23;
%         age2(n,t) = 1;
%     elseif age2(n,t) > 22 & age(n,t) < 29;
%         age2(n,t) = 2; 
%     elseif age2(n,t) > 30 & age(n,t) < 36;
%         age2(n,t) = 3;
%     else age2(n,t) > 35;
%         age2(n,t) = 4; 
%     end
%     end
% end

% save age2.mat
load('age2.mat')

%% experience bins
% exp_s2=exp_s

% for t=1:1:G.n_period-1
%    for n=1:1:G.n_pop
%         
%     if exp_s2(n,t)==0;
%         exp_s2(n,t)=0;
%     elseif exp_s2(n,t)==1; 
%         exp_s2(n,t)=1;
%     elseif exp_s2(n,t) > 1 & exp_s2(n,t) < 5;
%         exp_s2(n,t)=2;
%     elseif exp_s2(n,t) > 4 & exp_s2(n,t) < 10;
%         exp_s2(n,t)=3; 
%     else exp_s2(n,t)>9; 
%         exp_s2(n,t)=4;
%     end
%     end
% end

% save exp_s2.mat
load('exp_s2.mat')

%%Asset tertiles
tertile1=prctile(a_s,33)
tertile2=prctile(a_s,63)
tertile3=prctile(a_s,93)

for t=1:1:G.n_period-1;
   for n=1:1:G.n_pop;
        
    if a_s(n,t) < tertile1(:,t);
        a_s3(n,t) = 1;
    elseif a_s(n,t) >= tertile1(:,t) & a_s(n,t) < tertile2(:,t);
        a_s3(n,t) = 2; 
    else a_s(n,t) >= tertile2(:,t);
        a_s3(n,t) = 3; 
    end
    end
end

%% Work participation by age and sector (2 x 4 age groups)

r_s2=r_s;
n_s2=n_s;
u_s2=1-r_s2-n_s2;

% By age 
ls=r_s2+n_s2;
work_part=sum(ls,1)./G.n_pop;

work_part_age1=sum(sum(ls.*(age2==1)))/sum(sum(ls,1));
work_part_age2=sum(sum(ls.*(age2==2)))/sum(sum(ls,1));
work_part_age3=sum(sum(ls.*(age2==3)))/sum(sum(ls,1));
work_part_age4=sum(sum(ls.*(age2==4)))/sum(sum(ls,1));

work_part_age=[work_part_age1;work_part_age2;work_part_age3;work_part_age4]

% By age and sector

work_reg=sum(r_s2,1)./G.n_pop;
work_nreg=sum(n_s2,1)./G.n_pop;
unem=sum(u_s2,1)./G.n_pop;

%By age group and sector

work_reg_age1=sum(r_s2.*(age2==1))/sum(r_s2,1);
work_reg_age2=sum(r_s2.*(age2==2))/sum(r_s2,1);
work_reg_age3=sum(r_s2.*(age2==3))/sum(r_s2,1);
work_reg_age4=sum(r_s2.*(age2==4))/sum(r_s2,1);

work_nreg_age1=sum(n_s2.*(age2==1))/sum(n_s2,1);
work_nreg_age2=sum(n_s2.*(age2==2))/sum(n_s2,1);
work_nreg_age3=sum(n_s2.*(age2==3))/sum(n_s2,1);
work_nreg_age4=sum(n_s2.*(age2==4))/sum(n_s2,1);

work_unem_age1=sum(u_s2.*(age2==1))/sum(u_s2,1);
work_unem_age2=sum(u_s2.*(age2==2))/sum(u_s2,1);
work_unem_age3=sum(u_s2.*(age2==3))/sum(u_s2,1);
work_unem_age4=sum(u_s2.*(age2==4))/sum(u_s2,1);

work_reg_age=[work_reg_age1;work_reg_age2;work_reg_age3;work_reg_age4];
work_nreg_age=[work_nreg_age1;work_nreg_age2;work_nreg_age3;work_nreg_age4];
work_unem_age=[work_unem_age1;work_unem_age2;work_unem_age3;work_unem_age1];

%% Work participation by age and married (1 x 4 age groups)

%By age and married

m_s2=m_s(:,1:19);
work_marr=sum((m_s2==1).*ls,1)./sum(m_s2==1);
work_marr(isnan(work_marr))=0;
work_single=sum((m_s2==0).*ls,1)./sum(m_s2==0);
work_single(work_single>1)=1;

%By age group and married

work_marr_age1=sum((m_s2==1).*(age2==1).*ls,1)/sum((m_s2==1).*(age2==1));
work_marr_age2=sum((m_s2==1).*(age2==2).*ls,1)/sum((m_s2==1).*(age2==2));
work_marr_age3=sum((m_s2==1).*(age2==3).*ls,1)/sum((m_s2==1).*(age2==3));
work_marr_age4=sum((m_s2==1).*(age2==4).*ls,1)/sum((m_s2==1).*(age2==4));

work_single_age1=sum((m_s2==0).*(age2==1).*ls,1)/sum((m_s2==0).*(age2==1));
work_single_age2=sum((m_s2==0).*(age2==2).*ls,1)/sum((m_s2==0).*(age2==2));
work_single_age3=sum((m_s2==0).*(age2==3).*ls,1)/sum((m_s2==0).*(age2==3));
work_single_age4=sum((m_s2==0).*(age2==4).*ls,1)/sum((m_s2==0).*(age2==4));

work_marr_age=[work_marr_age1;work_marr_age2;work_marr_age3;work_marr_age4];
work_single_age=[work_single_age1;work_single_age2;work_single_age3;work_single_age4];

%% Work participation by age and number of children (2 x 4 age groups)

%By age and number of children

ch_s2=ch_s(:,1:19);
work_ch0=sum((ch_s2==0).*ls,1)./sum(ch_s2==0);
work_ch1=sum((ch_s2==1).*ls,1)./sum(ch_s2==1);
work_ch2=sum((ch_s2==2).*ls,1)./sum(ch_s2==2);

%By age group and number of children

work_ch0_age1=sum((ch_s2==0).*(age2==1).*ls,1)/sum((ch_s2==0).*(age2==1));
work_ch0_age2=sum((ch_s2==0).*(age2==2).*ls,1)/sum((ch_s2==0).*(age2==2));
work_ch0_age3=sum((ch_s2==0).*(age2==3).*ls,1)/sum((ch_s2==0).*(age2==3));
work_ch0_age4=sum((ch_s2==0).*(age2==4).*ls,1)/sum((ch_s2==0).*(age2==4));

work_ch1_age1=sum((ch_s2==1).*(age2==1).*ls,1)/sum((ch_s2==1).*(age2==1));
work_ch1_age2=sum((ch_s2==1).*(age2==2).*ls,1)/sum((ch_s2==1).*(age2==2));
work_ch1_age3=sum((ch_s2==1).*(age2==3).*ls,1)/sum((ch_s2==1).*(age2==3));
work_ch1_age4=sum((ch_s2==1).*(age2==4).*ls,1)/sum((ch_s2==1).*(age2==4));

work_ch2_age1=sum((ch_s2==2).*(age2==1).*ls,1)/sum((ch_s2==2).*(age2==1));
work_ch2_age2=sum((ch_s2==2).*(age2==2).*ls,1)/sum((ch_s2==2).*(age2==2));
work_ch2_age3=sum((ch_s2==2).*(age2==3).*ls,1)/sum((ch_s2==2).*(age2==3));
work_ch2_age4=sum((ch_s2==2).*(age2==4).*ls,1)/sum((ch_s2==2).*(age2==4));

work_ch0_age=[work_ch0_age1;work_ch0_age2;work_ch0_age3;work_ch0_age4];
work_ch1_age=[work_ch1_age1;work_ch1_age2;work_ch1_age3;work_ch1_age4];
work_ch2_age=[work_ch2_age1;work_ch2_age2;work_ch2_age3;work_ch2_age4];

%% Work participation by age and education (2 x 4 age groups)

%By age and education

work_hs=sum(ls(edu==1,1:19),1)./n_hs;
work_col2=sum(ls(edu==2,1:19),1)./n_col2;
work_col4=sum(ls(edu==3,1:19),1)./n_col4;

%By age group and education

work_hs_age1=sum((edu==1).*(age2==1).*ls)/sum((edu==1).*(age2==1));
work_hs_age2=sum((edu==1).*(age2==2).*ls)/sum((edu==1).*(age2==2));
work_hs_age3=sum((edu==1).*(age2==3).*ls)/sum((edu==1).*(age2==3));
work_hs_age4=sum((edu==1).*(age2==4).*ls)/sum((edu==1).*(age2==4));

work_col2_age1=sum((edu==2).*(age2==1).*ls)/sum((edu==2).*(age2==1));
work_col2_age2=sum((edu==2).*(age2==2).*ls)/sum((edu==2).*(age2==2));
work_col2_age3=sum((edu==2).*(age2==3).*ls)/sum((edu==2).*(age2==3));
work_col2_age4=sum((edu==2).*(age2==4).*ls)/sum((edu==2).*(age2==4));

work_col4_age1=sum((edu==3).*(age2==1).*ls)/sum((edu==3).*(age2==1));
work_col4_age2=sum((edu==3).*(age2==2).*ls)/sum((edu==3).*(age2==2));
work_col4_age3=sum((edu==3).*(age2==3).*ls)/sum((edu==3).*(age2==3));
work_col4_age4=sum((edu==3).*(age2==4).*ls)/sum((edu==3).*(age2==4));

work_hs_age=[work_hs_age1;work_hs_age2;work_hs_age3;work_hs_age4];
work_col2_age=[work_col2_age1;work_col2_age2;work_col2_age3;work_col2_age4];
work_col4_age=[work_col4_age1;work_col4_age2;work_col4_age3;work_col4_age4];

%% Fraction married by age (4 age groups)

%By age

marriage=sum(m_s2(:,[1:19]),1)./G.n_pop;

%By age group

marriage_age1=sum(sum(m_s2(:,[1:19]).*(age2==1),1))/G.n_pop;
marriage_age2=sum(sum(m_s2(:,[1:19]).*(age2==2),1))/G.n_pop;
marriage_age3=sum(sum(m_s2(:,[1:19]).*(age2==3),1))/G.n_pop;
marriage_age4=sum(sum(m_s2(:,[1:19]).*(age2==4),1))/G.n_pop;

marriage_age=[marriage_age1;marriage_age2;marriage_age3;marriage_age4]

%% Fraction married by age and education (2x4 age groups)

%By age and education

marr_hs=sum(m_s2(edu==1,1:19),1)./n_hs;
marr_col2=sum(m_s2(edu==2,1:19),1)./n_col2;
marr_col4=sum(m_s2(edu==3,1:19),1)./n_col4;

%By age group and education

marr_hs_age1=sum(m_s2.*(edu==1).*(age2==1))/sum((edu==2).*(age2==1));
marr_hs_age2=sum(m_s2.*(edu==1).*(age2==2))/sum((edu==2).*(age2==2));
marr_hs_age3=sum(m_s2.*(edu==1).*(age2==3))/sum((edu==2).*(age2==3));
marr_hs_age4=sum(m_s2.*(edu==1).*(age2==4))/sum((edu==2).*(age2==4));

marr_col2_age1=sum(m_s2.*(edu==2).*(age2==1))/sum((edu==2).*(age2==1));
marr_col2_age2=sum(m_s2.*(edu==2).*(age2==2))/sum((edu==2).*(age2==2));
marr_col2_age3=sum(m_s2.*(edu==2).*(age2==3))/sum((edu==2).*(age2==3));
marr_col2_age4=sum(m_s2.*(edu==2).*(age2==4))/sum((edu==2).*(age2==4));

marr_col4_age1=sum(m_s2.*(edu==3).*(age2==1))/sum((edu==2).*(age2==1));
marr_col4_age2=sum(m_s2.*(edu==3).*(age2==2))/sum((edu==2).*(age2==1));
marr_col4_age3=sum(m_s2.*(edu==3).*(age2==3))/sum((edu==2).*(age2==1));
marr_col4_age4=sum(m_s2.*(edu==3).*(age2==4))/sum((edu==2).*(age2==1));

marr_hs_age=[marr_hs_age1;marr_hs_age2;marr_hs_age3;marr_hs_age4];
marr_col2_age=[marr_col2_age1;marr_col2_age2;marr_col2_age3;marr_col2_age4];
marr_col4_age=[marr_col4_age1;marr_col4_age2;marr_col4_age3;marr_col4_age4];

%% Fraction married by age and sector (2x4 age groups)

%By age and sector

marr_reg=sum(m_s2.*r_s,1)./sum(r_s>0);
marr_nreg=sum(m_s2.*n_s,1)./sum(n_s>0);
marr_unem=sum(m_s2.*u_s,1)./sum(u_s>0);
marr_unem(isnan(marr_unem))=0;

% By age group and sector

marr_reg_age1=sum(sum(m_s2.*r_s.*(age2==1),1))./sum(sum((r_s>0).*(age2==1)));
marr_reg_age2=sum(sum(m_s2.*r_s.*(age2==2),1))./sum(sum((r_s>0).*(age2==2)));
marr_reg_age3=sum(sum(m_s2.*r_s.*(age2==3),1))./sum(sum((r_s>0).*(age2==3)));
marr_reg_age4=sum(sum(m_s2.*r_s.*(age2==4),1))./sum(sum((r_s>0).*(age2==4)));

marr_nreg_age1=sum(sum(m_s2.*n_s.*(age2==1),1))./sum(sum((n_s>0).*(age2==1)));
marr_nreg_age2=sum(sum(m_s2.*n_s.*(age2==2),1))./sum(sum((n_s>0).*(age2==2)));
marr_nreg_age3=sum(sum(m_s2.*n_s.*(age2==3),1))./sum(sum((n_s>0).*(age2==3)));
marr_nreg_age4=sum(sum(m_s2.*n_s.*(age2==4),1))./sum(sum((n_s>0).*(age2==4)));

marr_unem_age1=sum(sum(m_s2.*u_s.*(age2==1),1))./sum(sum((u_s>0).*(age2==1)));
marr_unem_age2=sum(sum(m_s2.*u_s.*(age2==2),1))./sum(sum((u_s>0).*(age2==2)));
marr_unem_age3=sum(sum(m_s2.*u_s.*(age2==3),1))./sum(sum((u_s>0).*(age2==3)));
marr_unem_age4=sum(sum(m_s2.*u_s.*(age2==4),1))./sum(sum((u_s>0).*(age2==4)));

marr_reg_age=[marr_reg_age1;marr_reg_age2;marr_reg_age3;marr_reg_age4];
marr_nreg_age=[marr_nreg_age1;marr_nreg_age2;marr_nreg_age3;marr_nreg_age4];
marr_unem_age=[marr_unem_age1;marr_unem_age2;marr_unem_age3;marr_unem_age4];

%% Assets by age (2x4 age groups) %% how to do tertiles? use mean/median

a_s2=a_s(:,1:19);

asset_age1_1=sum((age2==1).*(a_s3==1));
asset_age1_2=sum((age2==2).*(a_s3==1));
asset_age1_3=sum((age2==3).*(a_s3==1));
asset_age1_4=sum((age2==4).*(a_s3==1));

asset_age2_1=sum((age2==1).*(a_s3==2));
asset_age2_2=sum((age2==2).*(a_s3==2));
asset_age2_3=sum((age2==3).*(a_s3==2));
asset_age2_4=sum((age2==4).*(a_s3==2));

asset_age3_1=sum((age2==1).*(a_s3==3));
asset_age3_2=sum((age2==2).*(a_s3==3));
asset_age3_3=sum((age2==3).*(a_s3==3));
asset_age3_4=sum((age2==4).*(a_s3==3));

%By age

asset_age=(sum(a_s2,1)./G.n_pop);

%By age group, 1st tertile

asset_agegr1_1=nanmean(nansum(a_s2.*(age2==1).*(a_s3==1))./asset_age1_1);
asset_agegr1_2=nanmean(nansum(a_s2.*(age2==1).*(a_s3==1))./asset_age1_2);
asset_agegr1_3=nanmean(nansum(a_s2.*(age2==1).*(a_s3==1))./asset_age1_3);
asset_agegr1_4=nanmean(nansum(a_s2.*(age2==1).*(a_s3==1))./asset_age1_4);

%By age group, 2nd tertile

asset_agegr2_1=nanmean(nansum(a_s2.*(age2==1).*(a_s3==2))./asset_age2_1);
asset_agegr2_2=nanmean(nansum(a_s2.*(age2==1).*(a_s3==2))./asset_age2_2);
asset_agegr2_3=nanmean(nansum(a_s2.*(age2==1).*(a_s3==2))./asset_age2_3);
asset_agegr2_4=nanmean(nansum(a_s2.*(age2==1).*(a_s3==2))./asset_age2_4);

%By age group, 3rd tertile

asset_agegr3_1=nanmean(nansum(a_s2.*(age2==1).*(a_s3==3))./asset_age3_1);
asset_agegr3_2=nanmean(nansum(a_s2.*(age2==1).*(a_s3==3))./asset_age3_2);
asset_agegr3_3=nanmean(nansum(a_s2.*(age2==1).*(a_s3==3))./asset_age3_3);
asset_agegr3_4=nanmean(nansum(a_s2.*(age2==1).*(a_s3==3))./asset_age3_4);

asset_agegr1=[asset_agegr1_1;asset_agegr1_2;asset_agegr1_3;asset_agegr1_4];
asset_agegr2=[asset_agegr2_1;asset_agegr2_2;asset_agegr2_3;asset_agegr2_4];
asset_agegr3=[asset_agegr3_1;asset_agegr3_2;asset_agegr3_3;asset_agegr3_4];

%% Asset by age and marital status (2x1x4 age groups) %% need to check

%By age and married

% n_marr=sum(m_s2==1);
% n_single=G.n_pop-n_marr;
% 
% asset_marr=sum((a_s2.*(m_s2==1)),1)./n_marr;
% asset_single=sum((a_s2.*(m_s2==0)),1)./n_single;

%By age group and married, 1st tertile 

n_marr_age1_1=sum(sum(m_s2.*(age2==1).*(a_s3==1)));
n_marr_age1_2=sum(sum(m_s2.*(age2==2).*(a_s3==1)));
n_marr_age1_3=sum(sum(m_s2.*(age2==3).*(a_s3==1)));
n_marr_age1_4=sum(sum(m_s2.*(age2==4).*(a_s3==1)));

asset_marr_age1_1=sum(sum(a_s2.*(m_s2==1).*(age2==1).*(a_s3==1)))/n_marr_age1_1;
asset_marr_age1_2=sum(sum(a_s2.*(m_s2==1).*(age2==2).*(a_s3==1)))/n_marr_age1_1;
asset_marr_age1_3=sum(sum(a_s2.*(m_s2==1).*(age2==3).*(a_s3==1)))/n_marr_age1_1;
asset_marr_age1_4=sum(sum(a_s2.*(m_s2==1).*(age2==4).*(a_s3==1)))/n_marr_age1_1;

n_single_age1_1=G.n_pop-n_marr_age1_1;
n_single_age1_2=G.n_pop-n_marr_age1_2;
n_single_age1_3=G.n_pop-n_marr_age1_3;
n_single_age1_4=G.n_pop-n_marr_age1_4;

asset_single_age1_1=sum(sum(a_s2.*(m_s2==0).*(age2==1).*(a_s3==1)))/n_single_age1_1;
asset_single_age1_2=sum(sum(a_s2.*(m_s2==0).*(age2==2).*(a_s3==1)))/n_single_age1_1;
asset_single_age1_3=sum(sum(a_s2.*(m_s2==0).*(age2==3).*(a_s3==1)))/n_single_age1_1;
asset_single_age1_4=sum(sum(a_s2.*(m_s2==0).*(age2==4).*(a_s3==1)))/n_single_age1_1;

asset_marr_age1=[asset_marr_age1_1;asset_marr_age1_2;asset_marr_age1_3;asset_marr_age1_4];
asset_single_age1=[asset_single_age1_1;asset_single_age1_2;asset_single_age1_3;asset_single_age1_4];

%By age group and married, 2nd tertile 

n_marr_age2_1=sum(sum(m_s2.*(age2==1).*(a_s3==2)));
n_marr_age2_2=sum(sum(m_s2.*(age2==2).*(a_s3==2)));
n_marr_age2_3=sum(sum(m_s2.*(age2==3).*(a_s3==2)));
n_marr_age2_4=sum(sum(m_s2.*(age2==4).*(a_s3==2)));

asset_marr_age2_1=sum(sum(a_s2.*(m_s2==1).*(age2==1).*(a_s3==2)))/n_marr_age2_1;
asset_marr_age2_2=sum(sum(a_s2.*(m_s2==1).*(age2==2).*(a_s3==2)))/n_marr_age2_1;
asset_marr_age2_3=sum(sum(a_s2.*(m_s2==1).*(age2==3).*(a_s3==2)))/n_marr_age2_1;
asset_marr_age2_4=sum(sum(a_s2.*(m_s2==1).*(age2==4).*(a_s3==2)))/n_marr_age2_1;

n_single_age2_1=G.n_pop-n_marr_age2_1;
n_single_age2_2=G.n_pop-n_marr_age2_2;
n_single_age2_3=G.n_pop-n_marr_age2_3;
n_single_age2_4=G.n_pop-n_marr_age2_4;

asset_single_age2_1=sum(sum(a_s2.*(m_s2==0).*(age2==1).*(a_s3==2)))/n_single_age2_1;
asset_single_age2_2=sum(sum(a_s2.*(m_s2==0).*(age2==2).*(a_s3==2)))/n_single_age2_1;
asset_single_age2_3=sum(sum(a_s2.*(m_s2==0).*(age2==3).*(a_s3==2)))/n_single_age2_1;
asset_single_age2_4=sum(sum(a_s2.*(m_s2==0).*(age2==4).*(a_s3==2)))/n_single_age2_1;

asset_marr_age2=[asset_marr_age2_1;asset_marr_age2_2;asset_marr_age2_3;asset_marr_age2_4];
asset_single_age2=[asset_single_age2_1;asset_single_age2_2;asset_single_age2_3;asset_single_age2_4];

%By age group and married, 3rd tertile 

n_marr_age3_1=sum(sum(m_s2.*(age2==1).*(a_s3==3)));
n_marr_age3_2=sum(sum(m_s2.*(age2==2).*(a_s3==3)));
n_marr_age3_3=sum(sum(m_s2.*(age2==3).*(a_s3==3)));
n_marr_age3_4=sum(sum(m_s2.*(age2==4).*(a_s3==3)));

asset_marr_age3_1=sum(sum(a_s2.*(m_s2==1).*(age2==1).*(a_s3==3)))/n_marr_age3_1;
asset_marr_age3_2=sum(sum(a_s2.*(m_s2==1).*(age2==2).*(a_s3==3)))/n_marr_age3_1;
asset_marr_age3_3=sum(sum(a_s2.*(m_s2==1).*(age2==3).*(a_s3==3)))/n_marr_age3_1;
asset_marr_age3_4=sum(sum(a_s2.*(m_s2==1).*(age2==4).*(a_s3==3)))/n_marr_age3_1;

n_single_age3_1=G.n_pop-n_marr_age3_1;
n_single_age3_2=G.n_pop-n_marr_age3_2;
n_single_age3_3=G.n_pop-n_marr_age3_3;
n_single_age3_4=G.n_pop-n_marr_age3_4;

asset_single_age3_1=sum(sum(a_s2.*(m_s2==0).*(age2==1).*(a_s3==3)))/n_single_age3_1;
asset_single_age3_2=sum(sum(a_s2.*(m_s2==0).*(age2==2).*(a_s3==3)))/n_single_age3_1;
asset_single_age3_3=sum(sum(a_s2.*(m_s2==0).*(age2==3).*(a_s3==3)))/n_single_age3_1;
asset_single_age3_4=sum(sum(a_s2.*(m_s2==0).*(age2==4).*(a_s3==3)))/n_single_age3_1;

asset_marr_age3=[asset_marr_age3_1;asset_marr_age3_2;asset_marr_age3_3;asset_marr_age3_4];
asset_single_age3=[asset_single_age3_1;asset_single_age3_2;asset_single_age3_3;asset_single_age3_4];


%% Asset by age and education (2x2x4 age groups)

%By age and education

% asset_hs=(sum(a_s(edu==1,1:19),1)./n_hs);
% asset_col2=(sum(a_s(edu==2,1:19),1)./n_col2);
% asset_col4=(sum(a_s(edu==3,1:19),1)./n_col4);

%By age group and education, 1st tertile

n_hs_age1_1=sum(sum((age2==1).*(edu==1).*(a_s3==1)));
n_hs_age1_2=sum(sum((age2==2).*(edu==1).*(a_s3==1)));
n_hs_age1_3=sum(sum((age2==3).*(edu==1).*(a_s3==1)));
n_hs_age1_4=sum(sum((age2==4).*(edu==1).*(a_s3==1)));

n_col2_age1_1=sum(sum((age2==1).*(edu==2).*(a_s3==1)));
n_col2_age1_2=sum(sum((age2==2).*(edu==2).*(a_s3==1)));
n_col2_age1_3=sum(sum((age2==3).*(edu==2).*(a_s3==1)));
n_col2_age1_4=sum(sum((age2==4).*(edu==2).*(a_s3==1)));

n_col4_age1_1=sum(sum((age2==1).*(edu==3).*(a_s3==1)));
n_col4_age1_2=sum(sum((age2==2).*(edu==3).*(a_s3==1)));
n_col4_age1_3=sum(sum((age2==3).*(edu==3).*(a_s3==1)));
n_col4_age1_4=sum(sum((age2==4).*(edu==3).*(a_s3==1)));

asset_hs_age1_1=nansum(sum(a_s2.*(edu==1).*(age2==1).*(a_s3==1)))./n_hs_age1_1;
asset_hs_age1_2=nansum(sum(a_s2.*(edu==1).*(age2==2).*(a_s3==1)))./n_hs_age1_2;
asset_hs_age1_3=nansum(sum(a_s2.*(edu==1).*(age2==3).*(a_s3==1)))./n_hs_age1_3;
asset_hs_age1_4=nansum(sum(a_s2.*(edu==1).*(age2==4).*(a_s3==1)))./n_hs_age1_4;

asset_col2_age1_1=nansum(sum(a_s2.*(edu==2).*(age2==1).*(a_s3==1)))./n_col2_age1_1;
asset_col2_age1_2=nansum(sum(a_s2.*(edu==2).*(age2==2).*(a_s3==1)))./n_col2_age1_2;
asset_col2_age1_3=nansum(sum(a_s2.*(edu==2).*(age2==3).*(a_s3==1)))./n_col2_age1_3;
asset_col2_age1_4=nansum(sum(a_s2.*(edu==2).*(age2==4).*(a_s3==1)))./n_col2_age1_4;

asset_col4_age1_1=nansum(sum(a_s2.*(edu==3).*(age2==1).*(a_s3==1)))./n_col4_age1_1;
asset_col4_age1_2=nansum(sum(a_s2.*(edu==3).*(age2==2).*(a_s3==1)))./n_col4_age1_2;
asset_col4_age1_3=nansum(sum(a_s2.*(edu==3).*(age2==3).*(a_s3==1)))./n_col4_age1_3;
asset_col4_age1_4=nansum(sum(a_s2.*(edu==3).*(age2==4).*(a_s3==1)))./n_col4_age1_4;

asset_hs_age1=[asset_hs_age1_1;asset_hs_age1_2;asset_hs_age1_3;asset_hs_age1_4];
asset_col2_age1=[asset_col2_age1_1;asset_col2_age1_2;asset_col2_age1_3;asset_col2_age1_4];
asset_col4_age1=[asset_col4_age1_1;asset_col4_age1_2;asset_col4_age1_3;asset_col4_age1_4];

%By age group and education, 2nd tertile

n_hs_age2_1=sum(sum((age2==1).*(edu==1).*(a_s3==2)));
n_hs_age2_2=sum(sum((age2==2).*(edu==1).*(a_s3==2)));
n_hs_age2_3=sum(sum((age2==3).*(edu==1).*(a_s3==2)));
n_hs_age2_4=sum(sum((age2==4).*(edu==1).*(a_s3==2)));

n_col2_age2_1=sum(sum((age2==1).*(edu==2).*(a_s3==2)));
n_col2_age2_2=sum(sum((age2==2).*(edu==2).*(a_s3==2)));
n_col2_age2_3=sum(sum((age2==3).*(edu==2).*(a_s3==2)));
n_col2_age2_4=sum(sum((age2==4).*(edu==2).*(a_s3==2)));

n_col4_age2_1=sum(sum((age2==1).*(edu==3).*(a_s3==2)));
n_col4_age2_2=sum(sum((age2==2).*(edu==3).*(a_s3==2)));
n_col4_age2_3=sum(sum((age2==3).*(edu==3).*(a_s3==2)));
n_col4_age2_4=sum(sum((age2==4).*(edu==3).*(a_s3==2)));

asset_hs_age2_1=nansum(sum(a_s2.*(edu==1).*(age2==1).*(a_s3==2)))./n_hs_age2_1;
asset_hs_age2_2=nansum(sum(a_s2.*(edu==1).*(age2==2).*(a_s3==2)))./n_hs_age2_2;
asset_hs_age2_3=nansum(sum(a_s2.*(edu==1).*(age2==3).*(a_s3==2)))./n_hs_age2_3;
asset_hs_age2_4=nansum(sum(a_s2.*(edu==1).*(age2==4).*(a_s3==2)))./n_hs_age2_4;

asset_col2_age2_1=nansum(sum(a_s2.*(edu==2).*(age2==1).*(a_s3==2)))./n_col2_age2_1;
asset_col2_age2_2=nansum(sum(a_s2.*(edu==2).*(age2==2).*(a_s3==2)))./n_col2_age2_2;
asset_col2_age2_3=nansum(sum(a_s2.*(edu==2).*(age2==3).*(a_s3==2)))./n_col2_age2_3;
asset_col2_age2_4=nansum(sum(a_s2.*(edu==2).*(age2==4).*(a_s3==2)))./n_col2_age2_4;

asset_col4_age2_1=nansum(sum(a_s2.*(edu==3).*(age2==1).*(a_s3==2)))./n_col4_age2_1;
asset_col4_age2_2=nansum(sum(a_s2.*(edu==3).*(age2==2).*(a_s3==2)))./n_col4_age2_2;
asset_col4_age2_3=nansum(sum(a_s2.*(edu==3).*(age2==3).*(a_s3==2)))./n_col4_age2_3;
asset_col4_age2_4=nansum(sum(a_s2.*(edu==3).*(age2==4).*(a_s3==2)))./n_col4_age2_4;

asset_hs_age2=[asset_hs_age2_1;asset_hs_age2_2;asset_hs_age2_3;asset_hs_age2_4];
asset_col2_age2=[asset_col2_age2_1;asset_col2_age2_2;asset_col2_age2_3;asset_col2_age2_4];
asset_col4_age2=[asset_col4_age2_1;asset_col4_age2_2;asset_col4_age2_3;asset_col4_age2_4];

%By age group and education, 3rd tertile

n_hs_age3_1=sum(sum((age2==1).*(edu==1).*(a_s3==3)));
n_hs_age3_2=sum(sum((age2==2).*(edu==1).*(a_s3==3)));
n_hs_age3_3=sum(sum((age2==3).*(edu==1).*(a_s3==3)));
n_hs_age3_4=sum(sum((age2==4).*(edu==1).*(a_s3==3)));

n_col2_age3_1=sum(sum((age2==1).*(edu==2).*(a_s3==3)));
n_col2_age3_2=sum(sum((age2==2).*(edu==2).*(a_s3==3)));
n_col2_age3_3=sum(sum((age2==3).*(edu==2).*(a_s3==3)));
n_col2_age3_4=sum(sum((age2==4).*(edu==2).*(a_s3==3)));

n_col4_age3_1=sum(sum((age2==1).*(edu==3).*(a_s3==3)));
n_col4_age3_2=sum(sum((age2==2).*(edu==3).*(a_s3==3)));
n_col4_age3_3=sum(sum((age2==3).*(edu==3).*(a_s3==3)));
n_col4_age3_4=sum(sum((age2==4).*(edu==3).*(a_s3==3)));

asset_hs_age3_1=nansum(sum(a_s2.*(edu==1).*(age2==1).*(a_s3==3)))./n_hs_age3_1;
asset_hs_age3_2=nansum(sum(a_s2.*(edu==1).*(age2==2).*(a_s3==3)))./n_hs_age3_2;
asset_hs_age3_3=nansum(sum(a_s2.*(edu==1).*(age2==3).*(a_s3==3)))./n_hs_age3_3;
asset_hs_age3_4=nansum(sum(a_s2.*(edu==1).*(age2==4).*(a_s3==3)))./n_hs_age3_4;

asset_col2_age3_1=nansum(sum(a_s2.*(edu==2).*(age2==1).*(a_s3==3)))./n_col2_age3_1;
asset_col2_age3_2=nansum(sum(a_s2.*(edu==2).*(age2==2).*(a_s3==3)))./n_col2_age3_2;
asset_col2_age3_3=nansum(sum(a_s2.*(edu==2).*(age2==3).*(a_s3==3)))./n_col2_age3_3;
asset_col2_age3_4=nansum(sum(a_s2.*(edu==2).*(age2==4).*(a_s3==3)))./n_col2_age3_4;

asset_col4_age3_1=nansum(sum(a_s2.*(edu==3).*(age2==1).*(a_s3==3)))./n_col4_age3_1;
asset_col4_age3_2=nansum(sum(a_s2.*(edu==3).*(age2==2).*(a_s3==3)))./n_col4_age3_2;
asset_col4_age3_3=nansum(sum(a_s2.*(edu==3).*(age2==3).*(a_s3==3)))./n_col4_age3_3;
asset_col4_age3_4=nansum(sum(a_s2.*(edu==3).*(age2==4).*(a_s3==3)))./n_col4_age3_4;

asset_hs_age3=[asset_hs_age3_1;asset_hs_age3_2;asset_hs_age3_3;asset_hs_age3_4];
asset_col2_age3=[asset_col2_age3_1;asset_col2_age3_2;asset_col2_age3_3;asset_col2_age3_4];
asset_col4_age3=[asset_col4_age3_1;asset_col4_age3_2;asset_col4_age3_3;asset_col4_age3_4];

asset_output1=[asset_agegr1;asset_agegr2;asset_agegr3];
asset_output2=[asset_marr_age1;asset_marr_age2;asset_marr_age3;asset_single_age1;asset_single_age2;asset_single_age3];
asset_output3=[asset_hs_age1;asset_col2_age1;asset_col4_age1;asset_hs_age2;asset_col2_age2;asset_col4_age2;asset_hs_age3;asset_col2_age3;asset_col4_age3];
%% Experience by age (1x4 age groups)

exp_s3=exp_s2(:,1:19);

exp_age1=sum(exp_s3(age2==1),1)/sum(sum(age2==1));
exp_age2=sum(exp_s3(age2==2),1)/sum(sum(age2==2));
exp_age3=sum(exp_s3(age2==3),1)/sum(sum(age2==3));
exp_age4=sum(exp_s3(age2==4),1)/sum(sum(age2==4));

%% Experience by marital status and education (2x3)

% By marital status 

exp_marr=sum(sum(exp_s3.*(m_s2==1),1))./sum(n_marr);
exp_single=sum(sum(exp_s3.*(m_s2==0),1))./sum(n_single);

% By education

exp_hs=sum(sum(exp_s3(edu==1,1:19),1))./n_hs;
exp_col2=sum(sum(exp_s3(edu==2,1:19),1))./n_col2;
exp_col4=sum(sum(exp_s3(edu==3,1:19),1))./n_col4;

% By marital status and education 

n_marr_hs=sum((m_s2==1).*(edu==1));
n_marr_col2=sum((m_s2==1).*(edu==2));
n_marr_col4=sum((m_s2==1).*(edu==3));

exp_marr_hs=nansum(sum(exp_s3.*(m_s2==1).*(edu==1)))./sum(n_marr_hs);
exp_marr_col2=nansum(sum(exp_s3.*(m_s2==1).*(edu==2)))./sum(n_marr_col2);
exp_marr_col4=nansum(sum(exp_s3.*(m_s2==1).*(edu==3)))./sum(n_marr_col4);

n_single_hs=sum((m_s2==0).*(edu==1));
n_single_col2=sum((m_s2==0).*(edu==2));
n_single_col4=sum((m_s2==0).*(edu==3));

exp_single_hs=nansum(sum(exp_s3.*(m_s2==0).*(edu==1)))./sum(n_single_hs);
exp_single_col2=nansum(sum(exp_s3.*(m_s2==0).*(edu==2)))./sum(n_single_col2);
exp_single_col4=nansum(sum(exp_s3.*(m_s2==0).*(edu==3)))./sum(n_single_col4);

exp_marr_edu=[exp_marr_hs;exp_marr_col2;exp_marr_col4];
exp_single_edu=[exp_single_hs;exp_single_col2;exp_single_col4];

%% Transitions by general experience (3x3x5 bins 0,1,2-4,5-9,10+)
    %r_s, n_s, u_s -> rr_s, rn_s, ru_s, nr_s, nn_s, nu_s, ur_s, un_s, uu_s
    
%% Transitions between sectors by marital status (3x3x2)
    %r_s, n_s, u_s -> rr_s, rn_s, ru_s, nr_s, nn_s, nu_s, ur_s, un_s, uu_s
    
%% Transitions Between sectors by education (3x3x3)
    %r_s, n_s, u_s -> rr_s, rn_s, ru_s, nr_s, nn_s, nu_s, ur_s, un_s, uu_s

%% Husband earnings

% Mean and variance of earnings by wife's age (2x4 age groups) 

wh_s2=wh_s;
lhw= sum(wh_s2,1) ./ sum(wh_s2~=0,1); 

n_age1=sum(age2==1);
n_age2=sum(age2==2);
n_age3=sum(age2==3);
n_age4=sum(age2==4);

wh_age1=nansum(sum(wh_s2.*(age2==1))) ./ sum(n_age1);
wh_age2=nansum(sum(wh_s2.*(age2==2))) ./ sum(n_age2);
wh_age3=nansum(sum(wh_s2.*(age2==3))) ./ sum(n_age3);
wh_age4=nansum(sum(wh_s2.*(age2==4))) ./ sum(n_age4);

wh_age=[wh_age1;wh_age2;wh_age3;wh_age4];

var_wh1=var(wh_s(age2==1));
var_wh2=var(wh_s(age2==2));
var_wh3=var(wh_s(age2==3));
var_wh4=var(wh_s(age2==4));

wh_var=[var_wh1;var_wh2;var_wh3;var_wh4];

% Mean earnings by wife's education (2x3)

wh_hs=wh_s2.*(edu==1);
wh_col2=wh_s2.*(edu==2);
wh_col4=wh_s2.*(edu==3);

wh_hssel= sum(sum(wh_hs,1) ./ sum(wh_hs~=0,1))/t;
wh_col2sel= sum(sum(wh_col2,1) ./ sum(wh_col2~=0,1))/t;
wh_col4sel= sum(sum(wh_col4,1) ./ sum(wh_col4~=0,1))/t;

wh_edu=[wh_hssel;wh_col2sel;wh_col4sel];

%% Child Investments

% Mean and variance of investment by women's age (2x4 age groups) 

inv_s2=inv_s;
inv_sel= sum(inv_s2,1) ./sum(inv_s2~=0,1);

inv_age1=nansum(sum(inv_s2.*(age2==1))) ./ sum(n_age1);
inv_age2=nansum(sum(inv_s2.*(age2==2))) ./ sum(n_age2);
inv_age3=nansum(sum(inv_s2.*(age2==3))) ./ sum(n_age3);
inv_age4=nansum(sum(inv_s2.*(age2==4))) ./ sum(n_age4);

inv_age=[inv_age1;inv_age2;inv_age3;inv_age4];

var_inv1=var(inv_s2(age2==1));
var_inv2=var(inv_s2(age2==2));
var_inv3=var(inv_s2(age2==3));
var_inv4=var(inv_s2(age2==4));

inv_var=[var_inv1;var_inv2;var_inv3;var_inv4];

% Mean investment by wife's education (2x3)
inv_hs=inv_s2.*(edu==1);
inv_col2=inv_s2.*(edu==2);
inv_col4=inv_s2.*(edu==3);

inv_hssel= sum(sum(inv_hs,1) ./ sum(inv_hs~=0,1))/t;
inv_col2sel= sum(sum(inv_col2,1) ./ sum(inv_col2~=0,1))/t;
inv_col4sel= sum(sum(inv_col4,1) ./ sum(inv_col4~=0,1))/t;

inv_edu=[inv_hssel;inv_col2sel;inv_col4sel];

%% Wages 

% Mean wage by work sector

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

wr_lsel= nansum(sum(wr_l,1) ./ sum(wr_l~=0,1))/t;
wr_hsel= nansum(sum(wr_h,1) ./ sum(wr_h~=0,1))/t;
wn_lsel= nansum(sum(wn_l,1) ./ sum(wn_l~=0,1))/t;
wn_hsel= nansum(sum(wn_h,1) ./ sum(wn_h~=0,1))/t;

wg_output1=[wr_hsel;wn_hsel;wr_lsel;wn_lsel];

% Mean wage by sector and education (2x3)

wr_hs=wr_s.*r_s.*(edu==1);
wr_col2=wr_s.*r_s.*(edu==2);
wr_col4=wr_s.*r_s.*(edu==3);

wn_hs=wn_s.*n_s.*(edu==1);
wn_col2=wn_s.*n_s.*(edu==2);
wn_col4=wn_s.*n_s.*(edu==3);

wr_hssel= nansum(sum(wr_hs,1) ./ sum(wr_hs~=0,1))/t;
wr_col2sel= nansum(sum(wr_col2,1) ./ sum(wr_col2~=0,1))/t;
wr_col4sel= nansum(sum(wr_col4,1) ./ sum(wr_col4~=0,1))/t;

wn_hssel= nansum(sum(wn_hs,1) ./ sum(wn_hs~=0,1))/t;
wn_col2sel= nansum(sum(wn_col2,1) ./ sum(wn_col2~=0,1))/t;
wn_col4sel= nansum(sum(wn_col4,1) ./ sum(wn_col4~=0,1))/t;

wg_output2=[wr_hssel;wr_col2sel;wr_col4sel;wn_hssel;wn_col2sel;wn_col4sel];

% % Mean wage by sector, education and ability type (2x3x2)
% 
% n_reg_hs_abi1=sum((r_s2==1).*(edu==1).*(abi==1));
% n_reg_hs_abi0=sum((r_s2==1).*(edu==1).*(abi==0));
% 
% wr_hs_abi1=sum(wr_h(edu==1,1:19)) ./ n_reg_hs_abi1 *1000000/10000;
% wr_hs_abi0=sum(wr_l(edu==1,1:19)) ./ n_reg_hs_abi0 *1000000/10000;
% 
% n_nreg_hs_abi1=sum((n_s2==1).*(edu==1).*(abi==1));
% n_nreg_hs_abi0=sum((n_s2==1).*(edu==1).*(abi==0));
% 
% wn_hs_abi1=sum(wn_h(edu==1,1:19)) ./ n_nreg_hs_abi1 *1000000/10000;
% wn_hs_abi0=sum(wn_l(edu==1,1:19)) ./ n_nreg_hs_abi0 *1000000/10000;
% 
% wr_hs_abi=[wr_hs_abi1,wr_hs_abi0]
% wn_hs_abi=[wn_hs_abi1,wn_hs_abi0]
% 
% n_reg_col2_abi1=sum((r_s2==1).*(edu==2).*(abi==1));
% n_reg_col2_abi0=sum((r_s2==1).*(edu==2).*(abi==0));
% 
% wr_col2_abi1=sum(wr_h(edu==2,1:19)) ./ n_reg_col2_abi1 *1000000/10000;
% wr_col2_abi0=sum(wr_l(edu==2,1:19)) ./ n_reg_col2_abi0 *1000000/10000;
% 
% n_nreg_col2_abi1=sum((n_s2==1).*(edu==2).*(abi==1));
% n_nreg_col2_abi0=sum((n_s2==1).*(edu==2).*(abi==0));
% 
% wn_col2_abi1=sum(wn_h(edu==2,1:19)) ./ n_nreg_col2_abi1 *1000000/10000;
% wn_col2_abi0=sum(wn_l(edu==2,1:19)) ./ n_nreg_col2_abi0 *1000000/10000;
% 
% wr_col2_abi=[wr_col2_abi1,wr_col2_abi0]
% wn_col2_abi=[wn_col2_abi1,wn_col2_abi0]
% 
% n_reg_col4_abi1=sum((r_s2==1).*(edu==3).*(abi==1));
% n_reg_col4_abi0=sum((r_s2==1).*(edu==3).*(abi==0));
% 
% wr_col4_abi1=sum(wr_h(edu==3,1:19)) ./ n_reg_col4_abi1 *1000000/10000;
% wr_col4_abi0=sum(wr_l(edu==3,1:19)) ./ n_reg_col4_abi0 *1000000/10000;
% 
% n_nreg_col4_abi1=sum((n_s2==1).*(edu==3).*(abi==1));
% n_nreg_col4_abi0=sum((n_s2==1).*(edu==3).*(abi==0));
% 
% wn_col4_abi1=sum(wn_h(edu==3,1:19)) ./ n_nreg_col4_abi1 *1000000/10000;
% wn_col4_abi0=sum(wn_l(edu==3,1:19)) ./ n_nreg_col4_abi0 *1000000/10000;
% 
% wr_col4_abi=[wr_col4_abi1,wr_col4_abi0]
% wn_col4_abi=[wn_col4_abi1,wn_col4_abi0]

% Mean wage by sector and general experience (2x5 bins 0,1,2-4,5-9,10+)

wr_exp0=sum(sum(wr_s.*r_s.*(exp_s3==0)))/sum(r_s(exp_s3==0));
wr_exp1=sum(sum(wr_s.*r_s.*(exp_s3==1)))/sum(r_s(exp_s3==1));
wr_exp2=sum(sum(wr_s.*r_s.*(exp_s3==2)))/sum(r_s(exp_s3==2));
wr_exp3=sum(sum(wr_s.*r_s.*(exp_s3==3)))/sum(r_s(exp_s3==3));
wr_exp4=sum(sum(wr_s.*r_s.*(exp_s3==4)))/sum(r_s(exp_s3==4));

wn_exp0=sum(sum(wn_s.*n_s.*(exp_s3==0)))/sum(n_s(exp_s3==0));
wn_exp1=sum(sum(wn_s.*n_s.*(exp_s3==1)))/sum(n_s(exp_s3==1));
wn_exp2=sum(sum(wn_s.*n_s.*(exp_s3==2)))/sum(n_s(exp_s3==2));
wn_exp3=sum(sum(wn_s.*n_s.*(exp_s3==3)))/sum(n_s(exp_s3==3));
wn_exp4=sum(sum(wn_s.*n_s.*(exp_s3==4)))/sum(n_s(exp_s3==4));

wr_exp=[wr_exp0;wr_exp1;wr_exp2;wr_exp3;wr_exp4];
wn_exp=[wn_exp0;wn_exp1;wn_exp2;wn_exp3;wn_exp4];

wg_output3=[wr_exp;wn_exp];
    
% Variance of observed earnings (2-sector-specific)

var_rs=sum(var(wr_s))/t;
var_ns=sum(var(wn_s))/t;

wg_var=[var_rs;var_ns];

%% Output vector
work_part_output=[work_reg_age;work_nreg_age;work_marr_age;work_single_age;work_ch0_age;work_ch1_age;work_ch2_age;work_hs_age;work_col2_age;work_col4_age];
married_output=[marriage_age;marr_hs_age;marr_col2_age;marr_col4_age;marr_reg_age;marr_nreg_age;marr_unem_age];
asset_output=[asset_output1;asset_output2;asset_output3];
exp_output=[exp_age1;exp_age2;exp_age3;exp_age4;exp_marr_edu;exp_single_edu];
hw_inv_output=[wh_age;wh_var;wh_edu;inv_age;inv_var;inv_edu];
wg_output=[wg_output1;wg_output2;wg_output3;wg_var];

sim_moments=[work_part_output;married_output;asset_output;exp_output;hw_inv_output;wg_output]

