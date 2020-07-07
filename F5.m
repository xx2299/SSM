% % SSM_n for models and the observation as a function of mean SSM for observation over different time scales
% Remove regions with SSM less than 0.1
avg_mdif1_try1(find(mean_SM<0.1))=nan;
avg_mdif2_try1(find(mean_SM<0.1))=nan;
avg_mdif3_try1(find(mean_SM<0.1))=nan;
cv_nsm1(find(mean_SM<0.1))=nan;
cv_nsm2(find(mean_SM<0.1))=nan;
cv_nsm3(find(mean_SM<0.1))=nan;
mean_SM(find(mean_SM<0.1))=nan;
Flattened_mean_SM=mean_SM(:)';
MappedF_mean_SM=mapminmax(Flattened_mean_SM,0,1);
Mapped_mean_SM=reshape(MappedF_mean_SM,size(mean_SM));
% Separate SSM into 20 bins of equal size
dlevels=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1];
account=length(dlevels)-1;
idx=cell(1,account);
% Calculate averaged biases and coefficient of variation of SSM_n located in each bin
spdif_comp1=cell(1,account);
spdif_slope_avg1=cell(1,account);
spdif_slope_std1=cell(1,account);
for k=1:account
    idx{k}=find(Mapped_mean_SM>=dlevels(k) & Mapped_mean_SM<=dlevels(k+1));
    spdif_comp1{k}=avg_mdif1_try1(idx{k});
    spdif_slope_avg1{k}=mean(spdif_comp1{k},'omitnan');
    spdif_slope_std1{k}=std(spdif_comp1{k},'omitnan');
end
spdif_comp2=cell(1,account);
spdif_slope_avg2=cell(1,account);
spdif_slope_std2=cell(1,account);
for k=1:account
    spdif_comp2{k}=avg_mdif2_try1(idx{k});
    spdif_slope_avg2{k}=mean(spdif_comp2{k},'omitnan');
    spdif_slope_std2{k}=std(spdif_comp2{k},'omitnan');
end
spdif_comp3=cell(1,account);
spdif_slope_avg3=cell(1,account);
spdif_slope_std3=cell(1,account);
for k=1:account
    spdif_comp3{k}=avg_mdif3_try1(idx{k});
    spdif_slope_avg3{k}=mean(spdif_comp3{k},'omitnan');
    spdif_slope_std3{k}=std(spdif_comp3{k},'omitnan');
end
spdif_slope_Avg1=cell2mat(spdif_slope_avg1);
spdif_slope_Std1=cell2mat(spdif_slope_std1);
spdif_slope_Avg2=cell2mat(spdif_slope_avg2);
spdif_slope_Std2=cell2mat(spdif_slope_std2);
spdif_slope_Avg3=cell2mat(spdif_slope_avg3);
spdif_slope_Std3=cell2mat(spdif_slope_std3);
hpdif_comp1=cell(1,account);
hpdif_slope_avg1=cell(1,account);
hpdif_slope_std1=cell(1,account);
for k=1:account
    hpdif_comp1{k}=cv_nsm1(idx{k});
    hpdif_slope_avg1{k}=mean(hpdif_comp1{k},'omitnan');
    hpdif_slope_std1{k}=std(hpdif_comp1{k},'omitnan');
end
hpdif_comp2=cell(1,account);
hpdif_slope_avg2=cell(1,account);
hpdif_slope_std2=cell(1,account);
for k=1:account
    hpdif_comp2{k}=cv_nsm2(idx{k});
    hpdif_slope_avg2{k}=mean(hpdif_comp2{k},'omitnan');
    hpdif_slope_std2{k}=std(hpdif_comp2{k},'omitnan');
end
hpdif_comp3=cell(1,account);
hpdif_slope_avg3=cell(1,account);
hpdif_slope_std3=cell(1,account);
for k=1:account
    hpdif_comp3{k}=cv_nsm3(idx{k});
    hpdif_slope_avg3{k}=mean(hpdif_comp3{k},'omitnan');
    hpdif_slope_std3{k}=std(hpdif_comp3{k},'omitnan');
end
hpdif_slope_Avg1=cell2mat(hpdif_slope_avg1);
hpdif_slope_Std1=cell2mat(hpdif_slope_std1);
hpdif_slope_Avg2=cell2mat(hpdif_slope_avg2);
hpdif_slope_Std2=cell2mat(hpdif_slope_std2);
hpdif_slope_Avg3=cell2mat(hpdif_slope_avg3);
hpdif_slope_Std3=cell2mat(hpdif_slope_std3);
% Plot the diagram (Figure 5)
x=[0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,0.825,0.875,0.925,0.975];
plot1=subplot(2,3,1);
shadedErrorBar(x,spdif_slope_Avg1,spdif_slope_Std1,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
set(gca,'YTick',[-0.8,-0.4,0,0.4,0.8])
xlabel('(a)')
ylabel('SSM_n bias AVG')
title('1/30<f<1/7 day^{-1}')
plot2=subplot(2,3,2);
shadedErrorBar(x,spdif_slope_Avg2,spdif_slope_Std2,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
set(gca,'YTick',[-0.8,-0.4,0,0.4,0.8])
xlabel('(b)')
title('1/90<f<1/30 day^{-1}')
plot3=subplot(2,3,3);
shadedErrorBar(x,spdif_slope_Avg3,spdif_slope_Std3,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
set(gca,'YTick',[-0.8,-0.4,0,0.4,0.8])
xlabel('(c)')
title('1/365<f<1/90 day^{-1}')
plot4=subplot(2,3,4);
shadedErrorBar(x,hpdif_slope_Avg1,hpdif_slope_Std1,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
xlabel('(d)')
ylabel('SSM_n CV')
plot5=subplot(2,3,5);
shadedErrorBar(x,hpdif_slope_Avg2,hpdif_slope_Std2,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
xlabel('(e)')
plot6=subplot(2,3,6);
shadedErrorBar(x,hpdif_slope_Avg3,hpdif_slope_Std3,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
xlabel('(f)')
hAxis=axes('visible','off');
h=text(-0.05,0.5,'Normalized global mean SSM');
set(h,'fontsize',11,'VerticalAlignment','middle')
axis([plot1 plot2 plot3],[0 1 -0.8 0.8])
axis([plot4 plot5],[0 1 0 1])
axis(plot6,[0 1 0 0.4])