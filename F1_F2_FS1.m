% % Global SSM_n over different time scales
% Read SSM observation data (SMAP)
cd('E:\SMAP_L3_data')
latitude=h5read('SMAP_L3_SM_P_20190122_R16022_001.h5','/Soil_Moisture_Retrieval_Data_AM/latitude');
longitude=h5read('SMAP_L3_SM_P_20190122_R16022_001.h5','/Soil_Moisture_Retrieval_Data_AM/longitude');
lon=longitude(:,1);
lonSMAP=lon.';
[~,n]=size(latitude);
lat(1,n)=0;
count=1;
for i=1:n
    latidx=find(latitude(:,i)~=-9999);
    idx=latidx(1);
    lat(1,count)=latitude(idx,i);
    count=count+1;
end
latSMAP=lat.';
datadir='F:\SMAP_L3_data\';
filelist=dir([datadir,'*.h5']);
k=length(filelist);
for m=1:k
    filename=[datadir,filelist(m).name];
    SM_AM{m}=h5read(filename,'/Soil_Moisture_Retrieval_Data_AM/soil_moisture');
end
A=[];
for i=1:k
    eval(['A',num2str(i),'=','SM_AM{1,i}']);
    eval(['A',num2str(i),'=','reshape(A',num2str(i),',[1,391384])']);
    eval(['A=[A;A',num2str(i),'];']);
end
A(A<0)=nan;
% Fill SMAP missing data
PA=double(PA);
[m,n]=size(PA);
Z(m,n)=0;
count=1;
for i=1:n
    data=fillmissing(PA(:,i),'movmedian',10);
    Z(:,count)=fillgaps(data,80,40);
    count=count+1;
    disp(i);
end
Z=single(Z);
Z(Z<0.02)=0.02;
Z(Z>0.5)=0.5;
% FFT of SSM observation
Y=fft(Z);
Ayy=abs(Y).^2;
Ayy=Ayy(1:547,:);
N=1093;
Ayy=Ayy/(N/2);
Fs=1;
F=((1:N)-1)*Fs/N;
x=F;
a=find(abs(x-1/7)<=0.0004);
b=find(abs(x-1/30)<=0.0004);
c=find(abs(x-1/90)<=0.0004);
d=find(abs(x-1/365)<=0.0002);
e=2;
SMAP_ybar1=nansum(Ayy(b:a,:));
YBAR1=reshape(SMAP_ybar1,[964,406]);
Fybar1=YBAR1.';
SMAP_ybar2=nansum(Ayy(c:b,:));
YBAR2=reshape(SMAP_ybar2,[964,406]);
Fybar2=YBAR2.';
SMAP_ybar3=nansum(Ayy(d:c,:));
YBAR3=reshape(SMAP_ybar3,[964,406]);
Fybar3=YBAR3.';
smf1=SMAP_ybar1./(SMAP_ybar1+SMAP_ybar2+SMAP_ybar3);
smf2=SMAP_ybar2./(SMAP_ybar1+SMAP_ybar2+SMAP_ybar3);
smf3=SMAP_ybar3./(SMAP_ybar1+SMAP_ybar2+SMAP_ybar3);
SMF1=Fybar1./(Fybar1+Fybar2+Fybar3);
SMF2=Fybar2./(Fybar1+Fybar2+Fybar3);
SMF3=Fybar3./(Fybar1+Fybar2+Fybar3);

% % Global ET_n over different time scales
% Read ET observation data (GLEAM)
N=1097;
Fs=1;
EA(N,1036800)=0;
count=1;
for q=1:720
    for p=1:1440
        C=EData(q,p,:);
        EA(:,count)=C(:);
        count=count+1;
    end
end
EA=single(EA);
% FFT of ET observation
Y=fft(EA);
Y=Y(1:549,:);
Ayy=abs(Y).^2;
Ayy=Ayy/((N+1)/2);
F=((1:N)-1)*Fs/N;
x=F(1:(N+1)/2);
a=find(abs(x-1/7)<=0.0003);
b=find(abs(x-1/30)<=0.0004);
c=find(abs(x-1/90)<=0.0002);
d=find(abs(x-1/365)<=0.00001);
ybar1=nansum(Ayy(b:a,:));
YBAR1=reshape(ybar1,[1440,720]);
Fybar1=YBAR1.';
ybar2=nansum(Ayy(c:b,:));
YBAR2=reshape(ybar2,[1440,720]);
Fybar2=YBAR2.';
ybar3=nansum(Ayy(d:c,:));
YBAR3=reshape(ybar3,[1440,720]);
Fybar3=YBAR3.';
% Change the spatial resolution of ET_n to the same as SSM_n
ef1=ybar1./(ybar1+ybar2+ybar3);
ef2=ybar2./(ybar1+ybar2+ybar3);
ef3=ybar3./(ybar1+ybar2+ybar3);
count=1;
for i=1:720
    GLEAM_lat(count:count+1439,1)=latGLEAM(i,1);
    GLEAM_lon(count:count+1439,1)=lonGLEAM;
    count=count+1440;
end
GLEAM_lat=single(GLEAM_lat);
GLEAM_lon=single(GLEAM_lon);
ef1=ef1.';
ef2=ef2.';
ef3=ef3.';
GLEAM1=[GLEAM_lat,GLEAM_lon,ef1];
GLEAM2=[GLEAM_lat,GLEAM_lon,ef2];
GLEAM3=[GLEAM_lat,GLEAM_lon,ef3];
lon=lonSMAP.';
count=1;
for i=1:406
    SMAP_lat(count:count+963,1)=latSMAP(i,1);
    SMAP_lon(count:count+963,1)=lon;
    count=count+964;
end
SMAP_lat=single(SMAP_lat);
smf1=smf1.';
smf2=smf2.';
smf3=smf3.';
SMAP1=[SMAP_lat,SMAP_lon,smf1];
SMAP2=[SMAP_lat,SMAP_lon,smf2];
SMAP3=[SMAP_lat,SMAP_lon,smf3];
[m1,~]=size(GLEAM1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
GLEAM_new_f1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GLEAM1(i,1)-SMAP1(j,1))+abs(GLEAM1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GLEAM_new_f1(j,1)=GLEAM1(idx,3);
    disp(j)
end
new_f1=reshape(GLEAM_new_f1,[964,406]);
GLEAM_F1=new_f1.';
[m1,~]=size(GLEAM2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
GLEAM_new_f2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GLEAM2(i,1)-SMAP2(j,1))+abs(GLEAM2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GLEAM_new_f2(j,1)=GLEAM2(idx,3);
    disp(j)
end
new_f2=reshape(GLEAM_new_f2,[964,406]);
GLEAM_F2=new_f2.';
[m1,~]=size(GLEAM3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
GLEAM_new_f3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GLEAM3(i,1)-SMAP3(j,1))+abs(GLEAM3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GLEAM_new_f3(j,1)=GLEAM3(idx,3);
    disp(j)
end
new_f3=reshape(GLEAM_new_f3,[964,406]);
GLEAM_F3=new_f3.';

% % Global P_n over different time scales 
% Read P observation data (ERA5)
tp=ncread('ERA5_P.nc','tp');
tp=single(tp);
PA(t,6483600)=0;
count=1;
for q=1:1801
    for p=1:3600
        C=tp(p,q,:);
        PA(:,count)=C(:);
        count=count+1;
    end
end
PA=single(PA);
[m,n]=size(PA);
temp=isnan(PA);
temp=single(temp);
temp_sum=sum(temp);
Z(m,n)=0;
count=1;
for i=1:n
    if temp_sum(1,i)==0
        Z(:,count)=PA(:,i);
    else
        data=fillmissing(PA(:,i),'movmedian',10);
        data=double(data);
        Z(:,count)=fillgaps(data,80,40);
    end
    count=count+1;
    disp(i);
end
% FFT of P observation
N=t;
Fs=8;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
Z=single(Z);
Y=fft(Z);
Ayy=abs(Y).^2;
Ayy=Ayy/(N/2);
a=find(abs(x-1/7)<=0.0005);
b=find(abs(x-1/30)<=0.00045);
c=find(abs(x-1/90)<=0.0002);
d=find(abs(x-1/365)<=0.0002);
ybar1=nansum(Ayy(b:a,:));
ybar2=nansum(Ayy(c:b,:));
ybar3=nansum(Ayy(d:c,:));
% Change the spatial resolution of P_n to the same as SSM_n
temp=isnan(PA);
temp=single(temp);
temp_sum=sum(temp);
idx=find(temp_sum==120);
ybar1_tot(1,6483600)=0;
ybar1_tot(1,idx)=nan;
ybar1_tot=single(ybar1_tot);
nan_temp=isnan(ybar1_tot);
nan_temp=single(nan_temp);
idx1=find(nan_temp==0);
ybar1_tot(1,idx1)=ybar1;
ybar2_tot(1,6483600)=0;
ybar2_tot(1,idx)=nan;
ybar2_tot=single(ybar2_tot);
nan_temp=isnan(ybar2_tot);
nan_temp=single(nan_temp);
idx2=find(nan_temp==0);
ybar2_tot(1,idx2)=ybar2;
ybar3_tot(1,6483600)=0;
ybar3_tot(1,idx)=nan;
ybar3_tot=single(ybar3_tot);
nan_temp=isnan(ybar3_tot);
nan_temp=single(nan_temp);
idx3=find(nan_temp==0);
ybar3_tot(1,idx3)=ybar3;
tp1_1=reshape(ybar1_tot,[3600 1801]);
tp1_2=flipud(tp1_1);
tp1_F=rot90(tp1_2,-1);
a1=tp1_F(:,1:1800);
a2=tp1_F(:,1801:3600);
a3=[a2,a1];
tp2_1=reshape(ybar2_tot,[3600 1801]);
tp2_2=flipud(tp2_1);
tp2_F=rot90(tp2_2,-1);
b1=tp2_F(:,1:1800);
b2=tp2_F(:,1801:3600);
b3=[b2,b1];
tp3_1=reshape(ybar3_tot,[3600 1801]);
tp3_2=flipud(tp3_1);
tp3_F=rot90(tp3_2,-1);
c1=tp3_F(:,1:1800);
c2=tp3_F(:,1801:3600);
c3=[c2,c1];
PF1=a3./(a3+b3+c3);
PF2=b3./(a3+b3+c3);
PF3=c3./(a3+b3+c3);
lon=linspace(-179.95,179.95,3600);
lon=lon.';
count=1;
for i=1:1801
    ERA5_lat(count:count+3599,1)=latERA5(i,1);
    ERA5_lon(count:count+3599,1)=lon;
    count=count+3600;
end
ERA5_lat=single(ERA5_lat);
ERA5_lon=single(ERA5_lon);
ERA5_tp1_F=PF1.';
ERA5_PF1=reshape(ERA5_tp1_F,[6483600,1]);
ERA5_tp2_F=PF2.';
ERA5_PF2=reshape(ERA5_tp2_F,[6483600,1]);
ERA5_tp3_F=PF3.';
ERA5_PF3=reshape(ERA5_tp3_F,[6483600,1]);
ERAF1=[ERA5_lat,ERA5_lon,ERA5_PF1];
ERAF2=[ERA5_lat,ERA5_lon,ERA5_PF2];
ERAF3=[ERA5_lat,ERA5_lon,ERA5_PF3];
lon=lonSMAP.';
count=1;
for i=1:406
    SMAP_lat(count:count+963,1)=latSMAP(i,1);
    SMAP_lon(count:count+963,1)=lon;
    count=count+964;
end
SMAP_lat=single(SMAP_lat);
smf1=smf1.';
smf2=smf2.';
smf3=smf3.';
SMAP1=[SMAP_lat,SMAP_lon,smf1];
SMAP2=[SMAP_lat,SMAP_lon,smf2];
SMAP3=[SMAP_lat,SMAP_lon,smf3];
[m1,~]=size(ERAF1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
ERA5_new_f1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(ERAF1(i,1)-SMAP1(j,1))+abs(ERAF1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    ERA5_new_f1(j,1)=ERAF1(idx,3);
    disp(j)
end
new_f1=reshape(ERA5_new_f1,[964,406]);
ERA5_F1=new_f1.';
[m1,~]=size(ERAF2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
ERA5_new_f2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(ERAF2(i,1)-SMAP2(j,1))+abs(ERAF2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    ERA5_new_f2(j,1)=ERAF2(idx,3);
    disp(j)
end
new_f2=reshape(ERA5_new_f2,[964,406]);
ERA5_F2=new_f2.';
[m1,~]=size(ERAF3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
ERA5_new_f3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(ERAF3(i,1)-SMAP3(j,1))+abs(ERAF3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    ERA5_new_f3(j,1)=ERAF3(idx,3);
    disp(j)
end
new_f3=reshape(ERA5_new_f3,[964,406]);
ERA5_F3=new_f3.';

% Mean SSM for observation after spatiotemporal normalization (Figure S1)
meanA=nanmean(A);
meanS=(meanA-min(meanA))/(max(meanA)-min(meanA));
meanS1=reshape(meanS,[964,406]);
mean_SM=meanS1.';
detalgx=lonSMAP;
detalgy=latSMAP;
map=pcolor(detalgx,detalgy,mean_SM);
set(map,'alphadata',~isnan(mean_SM));
shading interp
colorbar
caxis([0,1])
xlabel('Longitude')
ylabel('Latitude')
set(gca,'Ydir','normal')
set(gca,'YTick',[-84,-60,-30,0,30,60,84])
% SSM_n,ET_n,and P_n over different time scales£¨Figure 1£©
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
GLEAM_F1_new=GLEAM_F1./temp;
GLEAM_F2_new=GLEAM_F2./temp;
GLEAM_F3_new=GLEAM_F3./temp;
ERA5_F1_new=ERA5_F1./temp;
ERA5_F2_new=ERA5_F2./temp;
ERA5_F3_new=ERA5_F3./temp;
SMF1(find(mean_SM<0.1))=0;
SMF2(find(mean_SM<0.1))=0;
SMF3(find(mean_SM<0.1))=0;
detalgx=lonSMAP;
detalgy=latSMAP;
subplot(3,3,1)
map1=pcolor(detalgx,detalgy,SMF1);
set(map1,'alphadata',~isnan(SMF1))
shading interp
colorbar
caxis([0,1])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(a)')
ylabel('SSM_n')
title('1/30<f<1/7 day^-1')
subplot(3,3,2)
map2=pcolor(detalgx,detalgy,SMF2);
set(map2,'alphadata',~isnan(SMF2))
shading interp
colorbar
caxis([0,1])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(b)')
title('1/90<f<1/30 day^-1')
subplot(3,3,3)
map3=pcolor(detalgx,detalgy,SMF3);
set(map3,'alphadata',~isnan(SMF3))
shading interp
colorbar
caxis([0,1])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(c)')
title('1/365<f<1/90 day^-1')
subplot(3,3,4)
map4=pcolor(detalgx,detalgy,GLEAM_F1_new);
set(map4,'alphadata',~isnan(GLEAM_F1_new))
shading interp
colorbar
caxis([0,1])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(d)')
ylabel('ET_n')
subplot(3,3,5)
map5=pcolor(detalgx,detalgy,GLEAM_F2_new);
set(map5,'alphadata',~isnan(GLEAM_F2_new))
shading interp
colorbar
caxis([0,1])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(e)')
subplot(3,3,6)
map6=pcolor(detalgx,detalgy,GLEAM_F3_new);
set(map6,'alphadata',~isnan(GLEAM_F3_new))
shading interp
colorbar
caxis([0,1])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(f)')
subplot(3,3,7)
map7=pcolor(detalgx,detalgy,ERA5_F1_new);
set(map7,'alphadata',~isnan(ERA5_F1_new))
shading interp
colorbar
caxis([0,1])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(g)')
ylabel('P_n')
subplot(3,3,8)
map8=pcolor(detalgx,detalgy,ERA5_F2_new);
set(map8,'alphadata',~isnan(ERA5_F2_new))
shading interp
colorbar
caxis([0,1])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(h)')
subplot(3,3,9)
map9=pcolor(detalgx,detalgy,ERA5_F3_new);
set(map9,'alphadata',~isnan(ERA5_F3_new))
shading interp
colorbar
caxis([0,1])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(i)')
axes('position',[0.81,0.11,0.15,0.75])
axis off
c=colorbar('eastoutside');
c.Label.FontSize=5;
% H_SEP_n and H_EEP_n over different time scales (Figure 2)
SEPR1=SMF1./(GLEAM_F1_new+ERA5_F1_new);
SEPR2=SMF2./(GLEAM_F2_new+ERA5_F2_new);
SEPR3=SMF3./(GLEAM_F3_new+ERA5_F3_new);
EEPR1=GLEAM_F1_new./(GLEAM_F1_new+ERA5_F1_new);
EEPR2=GLEAM_F2_new./(GLEAM_F2_new+ERA5_F2_new);
EEPR3=GLEAM_F3_new./(GLEAM_F3_new+ERA5_F3_new);
nSEPR1=SEPR1./(SEPR1+SEPR2+SEPR3);
nSEPR2=SEPR2./(SEPR1+SEPR2+SEPR3);
nSEPR3=SEPR3./(SEPR1+SEPR2+SEPR3);
nSEPR1(find(mean_SM<0.1))=0;
nSEPR2(find(mean_SM<0.1))=0;
nSEPR3(find(mean_SM<0.1))=0;
nEEPR1=EEPR1./(EEPR1+EEPR2+EEPR3);
nEEPR2=EEPR2./(EEPR1+EEPR2+EEPR3);
nEEPR3=EEPR3./(EEPR1+EEPR2+EEPR3);
nEEPR1(find(mean_SM<0.1))=0;
nEEPR2(find(mean_SM<0.1))=0;
nEEPR3(find(mean_SM<0.1))=0;
MAXn1=[SEPR1,SEPR2,SEPR3];
MAX1=max(max(MAXn1));
MAXn2=[EEPR1,EEPR2,EEPR3];
MAX2=max(max(MAXn2));
detalgx=lonSMAP;
detalgy=latSMAP;
[LON,LAT]=meshgrid(detalgx,detalgy);
subplot(2,3,1)
map1=pcolor(detalgx,detalgy,nSEPR1);
set(map1,'alphadata',~isnan(nSEPR1))
colorbar
caxis([0,1])
shading interp
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(a)')
ylabel('SSM_n/(ET_n+P_n)')
title('1/30<f<1/7 day^{-1}')
subplot(2,3,2)
map2=pcolor(detalgx,detalgy,nSEPR2);
set(map2,'alphadata',~isnan(nSEPR2))
colorbar
caxis([0,1])
shading interp
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(b)')
title('1/90<f<1/30 day^{-1}')
subplot(2,3,3)
map3=pcolor(detalgx,detalgy,nSEPR3);
set(map3,'alphadata',~isnan(nSEPR3))
colorbar
caxis([0,1])
shading interp
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(c)')
title('1/365<f<1/90 day^{-1}');
subplot(2,3,4)
map4=pcolor(detalgx,detalgy,nEEPR1);
set(map4,'alphadata',~isnan(nEEPR1))
colorbar
caxis([0,1])
shading interp
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(d)')
ylabel('ET_n/(ET_n+P_n)')
subplot(2,3,5)
map5=pcolor(detalgx,detalgy,nEEPR2);
set(map5,'alphadata',~isnan(nEEPR2))
colorbar
caxis([0,1])
shading interp
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(e)')
subplot(2,3,6)
map3=pcolor(detalgx,detalgy,nEEPR3);
set(map3,'alphadata',~isnan(nEEPR3))
colorbar
caxis([0,1])
shading interp
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(f)')
axes('position',[0.41,0.14,0.15,0.34])
axis off
colorbar('eastoutside')
caxis([0,1])
hold on
axes('position',[0.41,0.24,0.15,0.34])
axis off
colorbar('eastoutside')
caxis([0,1])
hold off