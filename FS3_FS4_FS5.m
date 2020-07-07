% % Global distribution of sand content and clay content of SM and their relationship with model biases
% Read sand content of SM for the first two layers from observation (GSDE)
sand_lonData=ncread('SAND5min.nc','lon');
sand_latData=ncread('SAND5min.nc','lat');
sand=ncread('SAND5min.nc','SAND');
sand1=sand(:,:,1);
sand1=single(sand1);
sand1(sand1==-999)=nan;
sand1=sand1/100;
sand1rot=rot90(sand1,3);
sand1rf=fliplr(sand1rot);
sand2=sand(:,:,2);
sand2=single(sand2);
sand2(sand2==-999)=nan;
sand2=sand2/100;
sand2rot=rot90(sand2,3);
sand2rf=fliplr(sand2rot);
% Read clay content of SM for the first two layers from observation (GSDE)
clay_lonData=ncread('CLAY5min.nc','lon');
clay_latData=ncread('CLAY5min.nc','lat');
clay=ncread('CLAY5min.nc','CLAY');
clay1=clay(:,:,1);
clay1=single(clay1);
clay1(clay1==-999)=nan;
clay1=clay1/100;
clay1rot=rot90(clay1,3);
clay1rf=fliplr(clay1rot);
clay2=clay(:,:,2);
clay2=single(clay2);
clay2(clay2==-999)=nan;
clay2=clay2/100;
clay2rot=rot90(clay2,3);
clay2rf=fliplr(clay2rot);
% Global map of sand content and clay content for the first two layers (Figure S3)
detalgx=-180:180;
detalgy=-90:90;
subplot(2,2,1)
map1=imagesc(detalgx,detalgy,sand1rf);
set(map1,'alphadata',~isnan(sand1rf));
colorbar;
caxis([0,1]);
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
title('GSDE soil sand content percentage layer 1 (4.5cm)')
subplot(2,2,2)
map4=imagesc(detalgx,detalgy,clay1rf);
set(map4,'alphadata',~isnan(clay1rf));
colorbar;
caxis([0,1]);
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
title('GSDE soil clay content percentage layer 1 (4.5cm)')
subplot(2,2,3)
map4=imagesc(detalgx,detalgy,sand2rf);
set(map4,'alphadata',~isnan(sand2rf));
colorbar;
caxis([0,1]);
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
title('GSDE soil sand content percentage layer 2 (9.1cm)')
subplot(2,2,4)
map4=imagesc(detalgx,detalgy,clay2rf);
set(map4,'alphadata',~isnan(clay2rf));
colorbar;
caxis([0,1]);
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
title('GSDE soil clay content percentage layer 2 (9.1cm)')
sgtitle('GSDE soil sand and clay content')
axes('position',[0.41,0.14,0.15,0.76])
axis off
colorbar('eastoutside')
caxis([0,1])
hold off
% Change the spatial resolution of sand content to the same as SMAP
sand_avg=(sand1rf+sand2rf)/2;
[v1,v2]=size(sand_avg);
sand_avg_try1=imresize(sand_avg,[v1/2,v2/2],'bilinear');
[w1,w2]=size(sand_avg_try1);
sand_londata=imresize(sand_lonData,[v2/2,1],'bilinear');
sand_latdata=imresize(sand_latData,[v1/2,1],'bilinear');
lon=sand_londata.';
count=1;
for i=1:w1
    sand_lat(count:count+w2-1,1)=sand_latdata(i,1);
    sand_lon(count:count+w2-1,1)=lon;
    count=count+w2;
end
sand_lat=single(sand_lat);
sand_avg_used=sand_avg_try1.';
sand1_f=reshape(sand_avg_used,[w1*w2,1]);
sand_content2=[sand_lat,sand_lon,sand1_f];
[m1,n1]=size(sand_content2);
[m2,n2]=size(SMAP1);
dif(m1,1)=0;
sand_new(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(sand_content2(i,1)-SMAP1(j,1))+abs(sand_content2(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    sand_new(j,1)=sand_content2(idx,3);
    disp(j)
end
new=reshape(sand_new,[964,406]);
sand_content_new=new.';
normalized_sand=sand_content_new./(max(max(sand_content_new)));
% Change the spatial resolution of clay content to the same as SMAP
clay_avg=(clay1rf+clay2rf)/2;
[v1,v2]=size(clay_avg);
clay_avg_try1=imresize(clay_avg,[v1/2,v2/2],'bilinear');
[w1,w2]=size(clay_avg_try1);
clay_londata=imresize(clay_lonData,[v2/2,1],'bilinear');
clay_latdata=imresize(clay_latData,[v1/2,1],'bilinear');
lon=clay_londata.';
count=1;
for i=1:w1
    clay_lat(count:count+w2-1,1)=clay_latdata(i,1);
    clay_lon(count:count+w2-1,1)=lon;
    count=count+w2;
end
clay_lat=single(clay_lat);
clay_avg_used=clay_avg_try1.';
clay1_f=reshape(clay_avg_used,[w1*w2,1]);
clay_content2=[clay_lat,clay_lon,clay1_f];
[m1,n1]=size(clay_content2);
[m2,n2]=size(SMAP1);
dif(m1,1)=0;
clay_new(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(clay_content2(i,1)-SMAP1(j,1))+abs(clay_content2(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    clay_new(j,1)=clay_content2(idx,3);
    disp(j)
end
new=reshape(clay_new,[964,406]);
clay_content_new=new.';
normalized_clay=clay_content_new./(max(max(clay_content_new)));
% Remove regions with SSM less than 0.1
normalized_sand(find(mean_SM<0.1))=nan;
normalized_clay(find(mean_SM<0.1))=nan;
avg_mdif1_try1(find(mean_SM<0.1))=nan;
avg_mdif2_try1(find(mean_SM<0.1))=nan;
avg_mdif3_try1(find(mean_SM<0.1))=nan;
avg_mhpdif1(find(mean_SM<0.1))=nan;
avg_mhpdif2(find(mean_SM<0.1))=nan;
avg_mhpdif3(find(mean_SM<0.1))=nan;
avg_hhpdif1(find(mean_SM<0.1))=nan;
avg_hhpdif2(find(mean_SM<0.1))=nan;
avg_hhpdif3(find(mean_SM<0.1))=nan;
% Separate sand content into 20 bins of equal size
x=[0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,0.825,0.875,0.925,0.975];
dlevels=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1];
% Calculate averaged biases of SSM_n, ET_n, P_n located in each bin
account=length(dlevels)-1;
sand_idx1=cell(1,account);
sand_comp1=cell(1,account);
sand_spdif_slope_avg1=cell(1,account);
sand_spdif_slope_std1=cell(1,account);
for k=1:account
    sand_idx1{k}=find(normalized_sand>=dlevels(k) & normalized_sand<=dlevels(k+1));
    sand_comp1{k}=avg_mdif1_try1(sand_idx1{k});
    sand_spdif_slope_avg1{k}=mean(sand_comp1{k},'omitnan');
    sand_spdif_slope_std1{k}=std(sand_comp1{k},'omitnan');
end
sand_idx2=cell(1,account);
sand_comp2=cell(1,account);
sand_spdif_slope_avg2=cell(1,account);
sand_spdif_slope_std2=cell(1,account);
for k=1:account
    sand_idx2{k}=find(normalized_sand>=dlevels(k) & normalized_sand<=dlevels(k+1));
    sand_comp2{k}=avg_mdif2_try1(sand_idx2{k});
    sand_spdif_slope_avg2{k}=mean(sand_comp2{k},'omitnan');
    sand_spdif_slope_std2{k}=std(sand_comp2{k},'omitnan');
end
sand_idx3=cell(1,account);
sand_comp3=cell(1,account);
sand_spdif_slope_avg3=cell(1,account);
sand_spdif_slope_std3=cell(1,account);
for k=1:account
    sand_idx3{k}=find(normalized_sand>=dlevels(k) & normalized_sand<=dlevels(k+1));
    sand_comp3{k}=avg_mdif3_try1(sand_idx3{k});
    sand_spdif_slope_avg3{k}=mean(sand_comp3{k},'omitnan');
    sand_spdif_slope_std3{k}=std(sand_comp3{k},'omitnan');
end
sand_spdif_slope_Avg1=cell2mat(sand_spdif_slope_avg1);
sand_spdif_slope_Std1=cell2mat(sand_spdif_slope_std1);
sand_spdif_slope_Avg2=cell2mat(sand_spdif_slope_avg2);
sand_spdif_slope_Std2=cell2mat(sand_spdif_slope_std2);
sand_spdif_slope_Avg3=cell2mat(sand_spdif_slope_avg3);
sand_spdif_slope_Std3=cell2mat(sand_spdif_slope_std3);
sand_idx1=cell(1,account);
sand_comp1=cell(1,account);
sand_hpdif_slope_avg1=cell(1,account);
sand_hpdif_slope_std1=cell(1,account);
for k=1:account
    sand_idx1{k}=find(normalized_sand>=dlevels(k) & normalized_sand<=dlevels(k+1));
    sand_comp1{k}=avg_mhpdif1(sand_idx1{k});
    sand_hpdif_slope_avg1{k}=mean(sand_comp1{k},'omitnan');
    sand_hpdif_slope_std1{k}=std(sand_comp1{k},'omitnan');
end
sand_idx2=cell(1,account);
sand_comp2=cell(1,account);
sand_hpdif_slope_avg2=cell(1,account);
sand_hpdif_slope_std2=cell(1,account);
for k=1:account
    sand_idx2{k}=find(normalized_sand>=dlevels(k) & normalized_sand<=dlevels(k+1));
    sand_comp2{k}=avg_mhpdif2(sand_idx2{k});
    sand_hpdif_slope_avg2{k}=mean(sand_comp2{k},'omitnan');
    sand_hpdif_slope_std2{k}=std(sand_comp2{k},'omitnan');
end
sand_idx3=cell(1,account);
sand_comp3=cell(1,account);
sand_hpdif_slope_avg3=cell(1,account);
sand_hpdif_slope_std3=cell(1,account);
for k=1:account
    sand_idx3{k}=find(normalized_sand>=dlevels(k) & normalized_sand<=dlevels(k+1));
    sand_comp3{k}=avg_mhpdif3(sand_idx3{k});
    sand_hpdif_slope_avg3{k}=mean(sand_comp3{k},'omitnan');
    sand_hpdif_slope_std3{k}=std(sand_comp3{k},'omitnan');
end
sand_hpdif_slope_Avg1=cell2mat(sand_hpdif_slope_avg1);
sand_hpdif_slope_Std1=cell2mat(sand_hpdif_slope_std1);
sand_hpdif_slope_Avg2=cell2mat(sand_hpdif_slope_avg2);
sand_hpdif_slope_Std2=cell2mat(sand_hpdif_slope_std2);
sand_hpdif_slope_Avg3=cell2mat(sand_hpdif_slope_avg3);
sand_hpdif_slope_Std3=cell2mat(sand_hpdif_slope_std3);
sand_idx1=cell(1,account);
sand_comp1=cell(1,account);
sand_ppdif_slope_avg1=cell(1,account);
sand_ppdif_slope_std1=cell(1,account);
for k=1:account
    sand_idx1{k}=find(normalized_sand>=dlevels(k) & normalized_sand<=dlevels(k+1));
    sand_comp1{k}=avg_hhpdif1(sand_idx1{k});
    sand_ppdif_slope_avg1{k}=mean(sand_comp1{k},'omitnan');
    sand_ppdif_slope_std1{k}=std(sand_comp1{k},'omitnan');
end
sand_idx2=cell(1,account);
sand_comp2=cell(1,account);
sand_ppdif_slope_avg2=cell(1,account);
sand_ppdif_slope_std2=cell(1,account);
for k=1:account
    sand_idx2{k}=find(normalized_sand>=dlevels(k) & normalized_sand<=dlevels(k+1));
    sand_comp2{k}=avg_hhpdif2(sand_idx2{k});
    sand_ppdif_slope_avg2{k}=mean(sand_comp2{k},'omitnan');
    sand_ppdif_slope_std2{k}=std(sand_comp2{k},'omitnan');
end
sand_idx3=cell(1,account);
sand_comp3=cell(1,account);
sand_ppdif_slope_avg3=cell(1,account);
sand_ppdif_slope_std3=cell(1,account);
for k=1:account
    sand_idx3{k}=find(normalized_sand>=dlevels(k) & normalized_sand<=dlevels(k+1));
    sand_comp3{k}=avg_hhpdif3(sand_idx3{k});
    sand_ppdif_slope_avg3{k}=mean(sand_comp3{k},'omitnan');
    sand_ppdif_slope_std3{k}=std(sand_comp3{k},'omitnan');
end
sand_ppdif_slope_Avg1=cell2mat(sand_ppdif_slope_avg1);
sand_ppdif_slope_Std1=cell2mat(sand_ppdif_slope_std1);
sand_ppdif_slope_Avg2=cell2mat(sand_ppdif_slope_avg2);
sand_ppdif_slope_Std2=cell2mat(sand_ppdif_slope_std2);
sand_ppdif_slope_Avg3=cell2mat(sand_ppdif_slope_avg3);
sand_ppdif_slope_Std3=cell2mat(sand_ppdif_slope_std3);
% Plot the diagram (Figure S4)
plot1=subplot(3,3,1);
shadedErrorBar(x,sand_spdif_slope_Avg1,sand_spdif_slope_Std1,'lineprops',{'-ro','MarkerFaceColor','r'});
grid on
xlabel('(a)')
ylabel('SSM_n diff')
title('1/30<f<1/7 day^{-1}')
plot2=subplot(3,3,2);
shadedErrorBar(x,sand_spdif_slope_Avg2,sand_spdif_slope_Std2,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
xlabel('(b)')
title('1/90<f<1/30 day^{-1}')
plot3=subplot(3,3,3);
shadedErrorBar(x,sand_spdif_slope_Avg3,sand_spdif_slope_Std3,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
xlabel('(c)')
title('1/365<f<1/90 day^{-1}')
plot4=subplot(3,3,4);
shadedErrorBar(x,sand_hpdif_slope_Avg1,sand_hpdif_slope_Std1,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
set(gca,'YTick',[-3,-1.5,0,1.5,3])
xlabel('(d)')
ylabel('H_S_E_P_n diff')
plot5=subplot(3,3,5);
shadedErrorBar(x,sand_hpdif_slope_Avg2,sand_hpdif_slope_Std2,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
set(gca,'YTick',[-3,-1.5,0,1.5,3])
xlabel('(e)')
plot6=subplot(3,3,6);
shadedErrorBar(x,sand_hpdif_slope_Avg3,sand_hpdif_slope_Std3,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
set(gca,'YTick',[-3,-1.5,0,1.5,3])
xlabel('(f)')
plot7=subplot(3,3,7);
shadedErrorBar(x,sand_ppdif_slope_Avg1,sand_ppdif_slope_Std1,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
xlabel('(g)')
ylabel('H_E_E_P_n diff')
plot8=subplot(3,3,8);
shadedErrorBar(x,sand_ppdif_slope_Avg2,sand_ppdif_slope_Std2,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
xlabel('(h)')
plot9=subplot(3,3,9);
shadedErrorBar(x,sand_ppdif_slope_Avg3,sand_ppdif_slope_Std3,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
xlabel('(i)')
axis([plot1 plot2 plot3],[0 1 -1 1])
axis([plot4 plot5 plot6],[0 1 -3 3])
axis([plot7 plot8 plot9],[0 1 -1 1])
hAxis=axes('visible','off');
h=text(-0.05,0.15,'Soil sand content');
set(h,'fontsize',11,'VerticalAlignment','middle')
% Separate clay content into 20 bins of equal size
x=[0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,0.825,0.875,0.925,0.975];
dlevels=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1];
% Calculate averaged biases of SSM_n, ET_n, P_n located in each bin
account=length(dlevels)-1;
clay_idx1=cell(1,account);
clay_comp1=cell(1,account);
clay_spdif_slope_avg1=cell(1,account);
clay_spdif_slope_std1=cell(1,account);
for k=1:account
    clay_idx1{k}=find(normalized_clay>=dlevels(k) & normalized_clay<=dlevels(k+1));
    clay_comp1{k}=avg_mdif1_try1(clay_idx1{k});
    clay_spdif_slope_avg1{k}=mean(clay_comp1{k},'omitnan');
    clay_spdif_slope_std1{k}=std(clay_comp1{k},'omitnan');
end
clay_idx2=cell(1,account);
clay_comp2=cell(1,account);
clay_spdif_slope_avg2=cell(1,account);
clay_spdif_slope_std2=cell(1,account);
for k=1:account
    clay_idx2{k}=find(normalized_clay>=dlevels(k) & normalized_clay<=dlevels(k+1));
    clay_comp2{k}=avg_mdif2_try1(clay_idx2{k});
    clay_spdif_slope_avg2{k}=mean(clay_comp2{k},'omitnan');
    clay_spdif_slope_std2{k}=std(clay_comp2{k},'omitnan');
end
clay_idx3=cell(1,account);
clay_comp3=cell(1,account);
clay_spdif_slope_avg3=cell(1,account);
clay_spdif_slope_std3=cell(1,account);
for k=1:account
    clay_idx3{k}=find(normalized_clay>=dlevels(k) & normalized_clay<=dlevels(k+1));
    clay_comp3{k}=avg_mdif3_try1(clay_idx3{k});
    clay_spdif_slope_avg3{k}=mean(clay_comp3{k},'omitnan');
    clay_spdif_slope_std3{k}=std(clay_comp3{k},'omitnan');
end
clay_spdif_slope_Avg1=cell2mat(clay_spdif_slope_avg1);
clay_spdif_slope_Std1=cell2mat(clay_spdif_slope_std1);
clay_spdif_slope_Avg2=cell2mat(clay_spdif_slope_avg2);
clay_spdif_slope_Std2=cell2mat(clay_spdif_slope_std2);
clay_spdif_slope_Avg3=cell2mat(clay_spdif_slope_avg3);
clay_spdif_slope_Std3=cell2mat(clay_spdif_slope_std3);
clay_idx1=cell(1,account);
clay_comp1=cell(1,account);
clay_hpdif_slope_avg1=cell(1,account);
clay_hpdif_slope_std1=cell(1,account);
for k=1:account
    clay_idx1{k}=find(normalized_clay>=dlevels(k) & normalized_clay<=dlevels(k+1));
    clay_comp1{k}=avg_mhpdif1(clay_idx1{k});
    clay_hpdif_slope_avg1{k}=mean(clay_comp1{k},'omitnan');
    clay_hpdif_slope_std1{k}=std(clay_comp1{k},'omitnan');
end
clay_idx2=cell(1,account);
clay_comp2=cell(1,account);
clay_hpdif_slope_avg2=cell(1,account);
clay_hpdif_slope_std2=cell(1,account);
for k=1:account
    clay_idx2{k}=find(normalized_clay>=dlevels(k) & normalized_clay<=dlevels(k+1));
    clay_comp2{k}=avg_mhpdif2(clay_idx2{k});
    clay_hpdif_slope_avg2{k}=mean(clay_comp2{k},'omitnan');
    clay_hpdif_slope_std2{k}=std(clay_comp2{k},'omitnan');
end
clay_idx3=cell(1,account);
clay_comp3=cell(1,account);
clay_hpdif_slope_avg3=cell(1,account);
clay_hpdif_slope_std3=cell(1,account);
for k=1:account
    clay_idx3{k}=find(normalized_clay>=dlevels(k) & normalized_clay<=dlevels(k+1));
    clay_comp3{k}=avg_mhpdif3(clay_idx3{k});
    clay_hpdif_slope_avg3{k}=mean(clay_comp3{k},'omitnan');
    clay_hpdif_slope_std3{k}=std(clay_comp3{k},'omitnan');
end
clay_hpdif_slope_Avg1=cell2mat(clay_hpdif_slope_avg1);
clay_hpdif_slope_Std1=cell2mat(clay_hpdif_slope_std1);
clay_hpdif_slope_Avg2=cell2mat(clay_hpdif_slope_avg2);
clay_hpdif_slope_Std2=cell2mat(clay_hpdif_slope_std2);
clay_hpdif_slope_Avg3=cell2mat(clay_hpdif_slope_avg3);
clay_hpdif_slope_Std3=cell2mat(clay_hpdif_slope_std3);
clay_idx1=cell(1,account);
clay_comp1=cell(1,account);
clay_ppdif_slope_avg1=cell(1,account);
clay_ppdif_slope_std1=cell(1,account);
for k=1:account
    clay_idx1{k}=find(normalized_clay>=dlevels(k) & normalized_clay<=dlevels(k+1));
    clay_comp1{k}=avg_hhpdif1(clay_idx1{k});
    clay_ppdif_slope_avg1{k}=mean(clay_comp1{k},'omitnan');
    clay_ppdif_slope_std1{k}=std(clay_comp1{k},'omitnan');
end
clay_idx2=cell(1,account);
clay_comp2=cell(1,account);
clay_ppdif_slope_avg2=cell(1,account);
clay_ppdif_slope_std2=cell(1,account);
for k=1:account
    clay_idx2{k}=find(normalized_clay>=dlevels(k) & normalized_clay<=dlevels(k+1));
    clay_comp2{k}=avg_hhpdif2(clay_idx2{k});
    clay_ppdif_slope_avg2{k}=mean(clay_comp2{k},'omitnan');
    clay_ppdif_slope_std2{k}=std(clay_comp2{k},'omitnan');
end
clay_idx3=cell(1,account);
clay_comp3=cell(1,account);
clay_ppdif_slope_avg3=cell(1,account);
clay_ppdif_slope_std3=cell(1,account);
for k=1:account
    clay_idx3{k}=find(normalized_clay>=dlevels(k) & normalized_clay<=dlevels(k+1));
    clay_comp3{k}=avg_hhpdif3(clay_idx3{k});
    clay_ppdif_slope_avg3{k}=mean(clay_comp3{k},'omitnan');
    clay_ppdif_slope_std3{k}=std(clay_comp3{k},'omitnan');
end
clay_ppdif_slope_Avg1=cell2mat(clay_ppdif_slope_avg1);
clay_ppdif_slope_Std1=cell2mat(clay_ppdif_slope_std1);
clay_ppdif_slope_Avg2=cell2mat(clay_ppdif_slope_avg2);
clay_ppdif_slope_Std2=cell2mat(clay_ppdif_slope_std2);
clay_ppdif_slope_Avg3=cell2mat(clay_ppdif_slope_avg3);
clay_ppdif_slope_Std3=cell2mat(clay_ppdif_slope_std3);
% Plot the diagram (Figure S5)
plot1=subplot(3,3,1);
shadedErrorBar(x,clay_spdif_slope_Avg1,clay_spdif_slope_Std1,'lineprops',{'-ro','MarkerFaceColor','r'});
grid on
xlabel('(a)')
ylabel('SSM_n diff')
title('1/30<f<1/7 day^{-1}')
plot2=subplot(3,3,2);
shadedErrorBar(x,clay_spdif_slope_Avg2,clay_spdif_slope_Std2,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
xlabel('(b)')
title('1/90<f<1/30 day^{-1}')
plot3=subplot(3,3,3);
shadedErrorBar(x,clay_spdif_slope_Avg3,clay_spdif_slope_Std3,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
xlabel('(c)')
title('1/365<f<1/90 day^{-1}')
plot4=subplot(3,3,4);
shadedErrorBar(x,clay_hpdif_slope_Avg1,clay_hpdif_slope_Std1,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
set(gca,'YTick',[-3,-1.5,0,1.5,3])
xlabel('(d)')
ylabel('H_S_E_P_n diff')
plot5=subplot(3,3,5);
shadedErrorBar(x,clay_hpdif_slope_Avg2,clay_hpdif_slope_Std2,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
set(gca,'YTick',[-3,-1.5,0,1.5,3])
xlabel('(e)')
plot6=subplot(3,3,6);
shadedErrorBar(x,clay_hpdif_slope_Avg3,clay_hpdif_slope_Std3,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
set(gca,'YTick',[-3,-1.5,0,1.5,3])
xlabel('(f)')
plot7=subplot(3,3,7);
shadedErrorBar(x,clay_ppdif_slope_Avg1,clay_ppdif_slope_Std1,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
xlabel('(g)')
ylabel('H_E_E_P_n diff')
plot8=subplot(3,3,8);
shadedErrorBar(x,clay_ppdif_slope_Avg2,clay_ppdif_slope_Std2,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
xlabel('(h)')
plot9=subplot(3,3,9);
shadedErrorBar(x,clay_ppdif_slope_Avg3,clay_ppdif_slope_Std3,'lineprops',{'-ro','MarkerFaceColor','r'})
grid on
xlabel('(i)')
axis([plot1 plot2 plot3],[0 1 -1 1])
axis([plot4 plot5 plot6],[0 1 -3 3])
axis([plot7 plot8 plot9],[0 1 -1 1])
hAxis=axes('visible','off');
h=text(-0.05,0.15,'Soil clay content');
set(h,'fontsize',11,'VerticalAlignment','middle')