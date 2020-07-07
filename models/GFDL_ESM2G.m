% % SSM_n of GFDL_ESM2G (mrsos)
% Get SSM_n over different time scales
N=52925;
Fs=1;
A(52925,12960)=0;
count=1;
for q=1:90
    for p=1:144
        C=mrsosData(p,q,:);
        A(:,count)=C(:);
        count=count+1;
    end
end
A(A==0)=nan;
A=single(A);
[row,col] = size(A);
dA(row,col)=0;
for i=1:col
    dA(:,i)=detrend(A(:,i));
end
dA=single(dA);
A1=dA(:,1:3240);
Y1=fft(A1);
Y1=Y1(1:26463,:);
Ayy1=abs(Y1).^2;
A2=dA(:,3241:6480);
Y2=fft(A2);
Y2=Y2(1:26463,:);
Ayy2=abs(Y2).^2;
A3=dA(:,6481:9720);
Y3=fft(A3);
Y3=Y3(1:26463,:);
Ayy3=abs(Y3).^2;
A4=dA(:,9721:12960);
Y4=fft(A4);
Y4=Y4(1:26463,:);
Ayy4=abs(Y4).^2;
GFDL_ESM2G_Ayy=[Ayy1,Ayy2,Ayy3,Ayy4];
GFDL_ESM2G_Ayy=GFDL_ESM2G_Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:(N+1)/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
GFDL_ESM2G_ybar1=sum(GFDL_ESM2G_Ayy(b:a,:));
YBAR1=reshape(GFDL_ESM2G_ybar1,[144,90]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
GFDL_ESM2G_ybar2=sum(GFDL_ESM2G_Ayy(c:b,:));
YBAR2=reshape(GFDL_ESM2G_ybar2,[144,90]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
GFDL_ESM2G_ybar3=sum(GFDL_ESM2G_Ayy(d:c,:));
YBAR3=reshape(GFDL_ESM2G_ybar3,[144,90]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
GFDL_ESM2G_mrsosf1=GFDL_ESM2G_ybar1./(GFDL_ESM2G_ybar1+GFDL_ESM2G_ybar2+GFDL_ESM2G_ybar3);
GFDL_ESM2G_mrsosf2=GFDL_ESM2G_ybar2./(GFDL_ESM2G_ybar1+GFDL_ESM2G_ybar2+GFDL_ESM2G_ybar3);
GFDL_ESM2G_mrsosf3=GFDL_ESM2G_ybar3./(GFDL_ESM2G_ybar1+GFDL_ESM2G_ybar2+GFDL_ESM2G_ybar3);
GFDL_ESM2G_mrsosF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
GFDL_ESM2G_mrsosF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
GFDL_ESM2G_mrsosF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(GFDL_ESM2G_lonData);
[q,~]=size(GFDL_ESM2G_latData);
m=GFDL_ESM2G_lonData(p);
GFDL_ESM2G_lon=linspace(-m/2,m/2,p);
GFDL_ESM2G_lon=GFDL_ESM2G_lon.';
BNU_lon1=GFDL_ESM2G_lon(1:p/2,1);
BNU_lon2=GFDL_ESM2G_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    GFDL_ESM2G_lat_D(count:count+(p-1),1)=GFDL_ESM2G_latData(i,1);
    GFDL_ESM2G_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
GFDL_ESM2G_mrsosf1=GFDL_ESM2G_mrsosf1.';
GFDL_ESM2G_mrsosf2=GFDL_ESM2G_mrsosf2.';
GFDL_ESM2G_mrsosf3=GFDL_ESM2G_mrsosf3.';
GFDL_ESM2G_alpha1=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_mrsosf1];
GFDL_ESM2G_alpha2=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_mrsosf2];
GFDL_ESM2G_alpha3=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_mrsosf3];
[m1,~]=size(GFDL_ESM2G_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
GFDL_ESM2G_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha1(i,1)-SMAP1(j,1))+abs(GFDL_ESM2G_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_P1(j,1)=GFDL_ESM2G_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(GFDL_ESM2G_new_P1,[964,406]);
GFDL_ESM2G_new_mrsosf1=new_P1.';
[m1,~]=size(GFDL_ESM2G_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
GFDL_ESM2G_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha2(i,1)-SMAP2(j,1))+abs(GFDL_ESM2G_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_P2(j,1)=GFDL_ESM2G_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(GFDL_ESM2G_new_P2,[964,406]);
GFDL_ESM2G_new_mrsosf2=new_mf2.';
[m1,~]=size(GFDL_ESM2G_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BCC_CSM_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha3(i,1)-SMAP3(j,1))+abs(GFDL_ESM2G_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P3(j,1)=GFDL_ESM2G_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(BCC_CSM_new_P3,[964,406]);
GFDL_ESM2G_new_mrsosf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
GFDL_ESM2G_mrsosf1_new=GFDL_ESM2G_new_mrsosf1./temp;
GFDL_ESM2G_mrsosf2_new=GFDL_ESM2G_new_mrsosf2./temp;
GFDL_ESM2G_mrsosf3_new=GFDL_ESM2G_new_mrsosf3./temp;

% % ET_n of GFDL_ESM2G (hfls)
% Get ET_n over different time scales
N=52925;
Fs=1;
A(52925,12960)=0;
count=1;
for q=1:90
    for p=1:144
        C=hflsData(p,q,:);
        A(:,count)=C(:);
        count=count+1;
    end
end
A=single(A);
[row,col] = size(A);
dA(row,col)=0;
for i=1:col
    dA(:,i)=detrend(A(:,i));
end
dA=single(dA);
A1=dA(:,1:3240);
Y1=fft(A1);
Y1=Y1(1:26463,:);
Ayy1=abs(Y1).^2;
A2=dA(:,3241:6480);
Y2=fft(A2);
Y2=Y2(1:26463,:);
Ayy2=abs(Y2).^2;
A3=dA(:,6481:9720);
Y3=fft(A3);
Y3=Y3(1:26463,:);
Ayy3=abs(Y3).^2;
A4=dA(:,9721:12960);
Y4=fft(A4);
Y4=Y4(1:26463,:);
Ayy4=abs(Y4).^2;
GFDL_ESM2G_Ayy=[Ayy1,Ayy2,Ayy3,Ayy4];
GFDL_ESM2G_Ayy=GFDL_ESM2G_Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:(N+1)/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
GFDL_ESM2G_ybar1=sum(GFDL_ESM2G_Ayy(b:a,:));
YBAR1=reshape(GFDL_ESM2G_ybar1,[144,90]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
GFDL_ESM2G_ybar2=sum(GFDL_ESM2G_Ayy(c:b,:));
YBAR2=reshape(GFDL_ESM2G_ybar2,[144,90]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
GFDL_ESM2G_ybar3=sum(GFDL_ESM2G_Ayy(d:c,:));
YBAR3=reshape(GFDL_ESM2G_ybar3,[144,90]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
GFDL_ESM2G_hflsf1=GFDL_ESM2G_ybar1./(GFDL_ESM2G_ybar1+GFDL_ESM2G_ybar2+GFDL_ESM2G_ybar3);
GFDL_ESM2G_hflsf2=GFDL_ESM2G_ybar2./(GFDL_ESM2G_ybar1+GFDL_ESM2G_ybar2+GFDL_ESM2G_ybar3);
GFDL_ESM2G_hflsf3=GFDL_ESM2G_ybar3./(GFDL_ESM2G_ybar1+GFDL_ESM2G_ybar2+GFDL_ESM2G_ybar3);
GFDL_ESM2G_hflsF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
GFDL_ESM2G_hflsF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
GFDL_ESM2G_hflsF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's ET_n to the same as SMAP
[p,~]=size(GFDL_ESM2G_lonData);
[q,~]=size(GFDL_ESM2G_latData);
m=GFDL_ESM2G_lonData(p);
GFDL_ESM2G_lon=linspace(-m/2,m/2,p);
GFDL_ESM2G_lon=GFDL_ESM2G_lon.';
BNU_lon1=GFDL_ESM2G_lon(1:p/2,1);
BNU_lon2=GFDL_ESM2G_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    GFDL_ESM2G_lat_D(count:count+(p-1),1)=GFDL_ESM2G_latData(i,1);
    GFDL_ESM2G_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
GFDL_ESM2G_hflsf1=GFDL_ESM2G_hflsf1.';
GFDL_ESM2G_hflsf2=GFDL_ESM2G_hflsf2.';
GFDL_ESM2G_hflsf3=GFDL_ESM2G_hflsf3.';
GFDL_ESM2G_alpha1=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_hflsf1];
GFDL_ESM2G_alpha2=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_hflsf2];
GFDL_ESM2G_alpha3=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_hflsf3];
[m1,~]=size(GFDL_ESM2G_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
GFDL_ESM2G_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha1(i,1)-SMAP1(j,1))+abs(GFDL_ESM2G_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_P1(j,1)=GFDL_ESM2G_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(GFDL_ESM2G_new_P1,[964,406]);
GFDL_ESM2G_new_hflsf1=new_P1.';
[m1,~]=size(GFDL_ESM2G_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
GFDL_ESM2G_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha2(i,1)-SMAP2(j,1))+abs(GFDL_ESM2G_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_P2(j,1)=GFDL_ESM2G_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(GFDL_ESM2G_new_P2,[964,406]);
GFDL_ESM2G_new_hflsf2=new_mf2.';
[m1,~]=size(GFDL_ESM2G_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BCC_CSM_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha3(i,1)-SMAP3(j,1))+abs(GFDL_ESM2G_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P3(j,1)=GFDL_ESM2G_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(BCC_CSM_new_P3,[964,406]);
GFDL_ESM2G_new_hflsf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
GFDL_ESM2G_hflsf1_new=GFDL_ESM2G_new_hflsf1./temp;
GFDL_ESM2G_hflsf2_new=GFDL_ESM2G_new_hflsf2./temp;
GFDL_ESM2G_hflsf3_new=GFDL_ESM2G_new_hflsf3./temp;

% % P_n of GFDL_ESM2G (pr)
% Get P_n over different time scales
N=52925;
Fs=1;
A(52925,12960)=0;
count=1;
for q=1:90
    for p=1:144
        C=prData(p,q,:);
        A(:,count)=C(:);
        count=count+1;
    end
    disp(q)
end
A=single(A);
[row,col] = size(A);
dA(row,col)=0;
for i=1:col
    dA(:,i)=detrend(A(:,i));
end
dA=single(dA);
A1=dA(:,1:3240);
Y1=fft(A1);
Y1=Y1(1:26463,:);
Ayy1=abs(Y1).^2;
A2=dA(:,3241:6480);
Y2=fft(A2);
Y2=Y2(1:26463,:);
Ayy2=abs(Y2).^2;
A3=dA(:,6481:9720);
Y3=fft(A3);
Y3=Y3(1:26463,:);
Ayy3=abs(Y3).^2;
A4=dA(:,9721:12960);
Y4=fft(A4);
Y4=Y4(1:26463,:);
Ayy4=abs(Y4).^2;
Ayy=[Ayy1,Ayy2,Ayy3,Ayy4];
Ayy=Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:(N+1)/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
GFDL_ESM2G_ybar1=sum(GFDL_ESM2G_Ayy(b:a,:));
YBAR1=reshape(GFDL_ESM2G_ybar1,[144,90]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
GFDL_ESM2G_ybar2=sum(GFDL_ESM2G_Ayy(c:b,:));
YBAR2=reshape(GFDL_ESM2G_ybar2,[144,90]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
GFDL_ESM2G_ybar3=sum(GFDL_ESM2G_Ayy(d:c,:));
YBAR3=reshape(GFDL_ESM2G_ybar3,[144,90]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
GFDL_ESM2G_prf1=GFDL_ESM2G_ybar1./(GFDL_ESM2G_ybar1+GFDL_ESM2G_ybar2+GFDL_ESM2G_ybar3);
GFDL_ESM2G_prf2=GFDL_ESM2G_ybar2./(GFDL_ESM2G_ybar1+GFDL_ESM2G_ybar2+GFDL_ESM2G_ybar3);
GFDL_ESM2G_prf3=GFDL_ESM2G_ybar3./(GFDL_ESM2G_ybar1+GFDL_ESM2G_ybar2+GFDL_ESM2G_ybar3);
GFDL_ESM2G_prF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
GFDL_ESM2G_prF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
GFDL_ESM2G_prF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's P_n to the same as SMAP
[p,~]=size(GFDL_ESM2G_lonData);
[q,~]=size(GFDL_ESM2G_latData);
m=GFDL_ESM2G_lonData(p);
GFDL_ESM2G_lon=linspace(-m/2,m/2,p);
GFDL_ESM2G_lon=GFDL_ESM2G_lon.';
BNU_lon1=GFDL_ESM2G_lon(1:p/2,1);
BNU_lon2=GFDL_ESM2G_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    GFDL_ESM2G_lat_D(count:count+(p-1),1)=GFDL_ESM2G_latData(i,1);
    GFDL_ESM2G_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
GFDL_ESM2G_prf1=GFDL_ESM2G_prf1.';
GFDL_ESM2G_prf2=GFDL_ESM2G_prf2.';
GFDL_ESM2G_prf3=GFDL_ESM2G_prf3.';
GFDL_ESM2G_alpha1=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_prf1];
GFDL_ESM2G_alpha2=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_prf2];
GFDL_ESM2G_alpha3=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_prf3];
[m1,~]=size(GFDL_ESM2G_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
GFDL_ESM2G_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha1(i,1)-SMAP1(j,1))+abs(GFDL_ESM2G_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_P1(j,1)=GFDL_ESM2G_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(GFDL_ESM2G_new_P1,[964,406]);
GFDL_ESM2G_new_prf1=new_P1.';
[m1,~]=size(GFDL_ESM2G_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
GFDL_ESM2G_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha2(i,1)-SMAP2(j,1))+abs(GFDL_ESM2G_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_P2(j,1)=GFDL_ESM2G_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(GFDL_ESM2G_new_P2,[964,406]);
GFDL_ESM2G_new_prf2=new_mf2.';
[m1,~]=size(GFDL_ESM2G_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
GFDL_ESM2G_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha3(i,1)-SMAP3(j,1))+abs(GFDL_ESM2G_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_P3(j,1)=GFDL_ESM2G_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(GFDL_ESM2G_new_P3,[964,406]);
GFDL_ESM2G_new_prf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
GFDL_ESM2G_prf1_new=GFDL_ESM2G_new_prf1./temp;
GFDL_ESM2G_prf2_new=GFDL_ESM2G_new_prf2./temp;
GFDL_ESM2G_prf3_new=GFDL_ESM2G_new_prf3./temp;

% Differences of SSM_n between models and observations over different time scales
GFDL_ESM2G_mdif1=GFDL_ESM2G_mF1-SMF1;
GFDL_ESM2G_mdif2=GFDL_ESM2G_mF2-SMF2;
GFDL_ESM2G_mdif3=GFDL_ESM2G_mF3-SMF3;
[m,n]=size(SMF1);
GFDL_ESM2G_mdif1_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF1(:,i);
    temp2=GFDL_ESM2G_mdif1(:,i);
    idx1=isnan(temp1);
    idx1=double(idx1);
    idx2=isnan(temp2);
    idx2=double(idx2);
    n_nan1=numel(find(isnan(temp1)));
    m_nan1=m-n_nan1;
    n_nan2=numel(find(isnan(temp2)));
    m_nan2=m-n_nan2;
    if m_nan1==0 || m_nan2==0
        idx1=nan;
    else
        use=temp2(idx2==0);
        new=imresize(use,[m_nan1,1],'bilinear');
        idx1(idx1==1)=nan;
        idx1(idx1==0)=new;
    end
    GFDL_ESM2G_mdif1_try1(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF2);
GFDL_ESM2G_mdif2_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF2(:,i);
    temp2=GFDL_ESM2G_mdif2(:,i);
    idx1=isnan(temp1);
    idx1=double(idx1);
    idx2=isnan(temp2);
    idx2=double(idx2);
    n_nan1=numel(find(isnan(temp1)));
    m_nan1=m-n_nan1;
    n_nan2=numel(find(isnan(temp2)));
    m_nan2=m-n_nan2;
    if m_nan1==0 || m_nan2==0
        idx1=nan;
    else
        use=temp2(idx2==0);
        new=imresize(use,[m_nan1,1],'bilinear');
        idx1(idx1==1)=nan;
        idx1(idx1==0)=new;
    end
    GFDL_ESM2G_mdif2_try1(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF3);
GFDL_ESM2G_mdif3_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF3(:,i);
    temp2=GFDL_ESM2G_mdif3(:,i);
    idx1=isnan(temp1);
    idx1=double(idx1);
    idx2=isnan(temp2);
    idx2=double(idx2);
    n_nan1=numel(find(isnan(temp1)));
    m_nan1=m-n_nan1;
    n_nan2=numel(find(isnan(temp2)));
    m_nan2=m-n_nan2;
    if m_nan1==0 || m_nan2==0
        idx1=nan;
    else
        use=temp2(idx2==0);
        new=imresize(use,[m_nan1,1],'bilinear');
        idx1(idx1==1)=nan;
        idx1(idx1==0)=new;
    end
    GFDL_ESM2G_mdif3_try1(:,count)=idx1;
    count=count+1;
end

% % SSM_kw of GFDL_ESM2G (mrsos)
% Get SSM_kw over different time scales
N=52925;
F=((1:N)-1)*Fs/N;
x=F(1:(N+1)/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
[m1,n1]=size(GFDL_ESM2G_mrsosA);
[m2,n2]=size(GFDL_ESM2G_mrsosF1);
GFDL_ESM2G_mrsosp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=GFDL_ESM2G_mrsosA(:,i);
    temp2=isnan(GFDL_ESM2G_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        GFDL_ESM2G_mrsosp1(1,count)=nan;
    else
        A=GFDL_ESM2G_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m1+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        GFDL_ESM2G_mrsosp1(1,count)=p1(1);
    end
    count=count+1;
end
mrsosp11=reshape(GFDL_ESM2G_mrsosp1,[n2,m2]);
p11=rot90(mrsosp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
GFDL_ESM2G_mrsosP1=[P2 P1];
GFDL_ESM2G_mrsosp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=GFDL_ESM2G_mrsosA(:,i);
    temp2=isnan(GFDL_ESM2G_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        GFDL_ESM2G_mrsosp2(1,count)=nan;
    else
        A=GFDL_ESM2G_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m1+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        GFDL_ESM2G_mrsosp2(1,count)=p2(1);
    end
    count=count+1;
end
mrsosp22=reshape(GFDL_ESM2G_mrsosp2,[n2,m2]);
p22=rot90(mrsosp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
GFDL_ESM2G_mrsosP2=[P2 P1];
GFDL_ESM2G_mrsosp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=GFDL_ESM2G_mrsosA(:,i);
    temp2=isnan(GFDL_ESM2G_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        GFDL_ESM2G_mrsosp3(1,count)=nan;
    else
        A=GFDL_ESM2G_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m1+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        GFDL_ESM2G_mrsosp3(1,count)=p3(1);
    end
    count=count+1;
end
mrsosp33=reshape(GFDL_ESM2G_mrsosp3,[n2,m2]);
p33=rot90(mrsosp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
GFDL_ESM2G_mrsosP3=[P2 P1];
% Change the spatial resolution of the model's SSM_kw to the same as SMAP
[p,~]=size(GFDL_ESM2G_lonData);
[q,~]=size(GFDL_ESM2G_latData);
m=GFDL_ESM2G_lonData(p);
GFDL_ESM2G_lon=linspace(-m/2,m/2,p);
GFDL_ESM2G_lon=GFDL_ESM2G_lon.';
BNU_lon1=GFDL_ESM2G_lon(1:p/2,1);
BNU_lon2=GFDL_ESM2G_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    GFDL_ESM2G_lat_D(count:count+(p-1),1)=GFDL_ESM2G_latData(i,1);
    GFDL_ESM2G_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
GFDL_ESM2G_mrsosp1=GFDL_ESM2G_mrsosp1.';
GFDL_ESM2G_mrsosp2=GFDL_ESM2G_mrsosp2.';
GFDL_ESM2G_mrsosp3=GFDL_ESM2G_mrsosp3.';
GFDL_ESM2G_alpha1=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_mrsosp1];
GFDL_ESM2G_alpha2=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_mrsosp2];
GFDL_ESM2G_alpha3=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_mrsosp3];
[m1,~]=size(GFDL_ESM2G_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
GFDL_ESM2G_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha1(i,1)-SMAP1(j,1))+abs(GFDL_ESM2G_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_P1(j,1)=GFDL_ESM2G_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(GFDL_ESM2G_new_P1,[964,406]);
GFDL_ESM2G_new_mrsosP1=new_P1.';
[m1,~]=size(GFDL_ESM2G_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
GFDL_ESM2G_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha2(i,1)-SMAP2(j,1))+abs(GFDL_ESM2G_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_P2(j,1)=GFDL_ESM2G_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(GFDL_ESM2G_new_P2,[964,406]);
GFDL_ESM2G_new_mrsosP2=new_mf2.';
[m1,~]=size(GFDL_ESM2G_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BCC_CSM_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha3(i,1)-SMAP3(j,1))+abs(GFDL_ESM2G_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P3(j,1)=GFDL_ESM2G_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(BCC_CSM_new_P3,[964,406]);
GFDL_ESM2G_new_mrsosP3=new_mf3.';
[m,n]=size(SMF1);
GFDL_ESM2G_mrsosP1_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF1(:,i);
    new_temp2=GFDL_ESM2G_new_mrsosP1(:,i);
    idx1=isnan(new_temp1);
    idx1=double(idx1);
    idx2=isnan(new_temp2);
    idx2=double(idx2);
    n_nan1=numel(find(isnan(new_temp1)));
    m_nan1=m-n_nan1;
    n_nan2=numel(find(isnan(new_temp2)));
    m_nan2=m-n_nan2;
    if m_nan1==0 || m_nan2==0
        idx1=nan;
    else
        use=new_temp2(idx2==0);
        new=imresize(use,[m_nan1,1],'bilinear');
        idx1(idx1==1)=nan;
        idx1(idx1==0)=new;
    end
    GFDL_ESM2G_mrsosP1_new(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF2);
GFDL_ESM2G_mrsosP2_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF2(:,i);
    new_temp2=GFDL_ESM2G_new_mrsosP2(:,i);
    idx1=isnan(new_temp1);
    idx1=double(idx1);
    idx2=isnan(new_temp2);
    idx2=double(idx2);
    n_nan1=numel(find(isnan(new_temp1)));
    m_nan1=m-n_nan1;
    n_nan2=numel(find(isnan(new_temp2)));
    m_nan2=m-n_nan2;
    if m_nan1==0 || m_nan2==0
        idx1=nan;
    else
        use=new_temp2(idx2==0);
        new=imresize(use,[m_nan1,1],'bilinear');
        idx1(idx1==1)=nan;
        idx1(idx1==0)=new;
    end
    GFDL_ESM2G_mrsosP2_new(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF3);
GFDL_ESM2G_mrsosP3_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF3(:,i);
    new_temp2=GFDL_ESM2G_new_mrsosP3(:,i);
    idx1=isnan(new_temp1);
    idx1=double(idx1);
    idx2=isnan(new_temp2);
    idx2=double(idx2);
    n_nan1=numel(find(isnan(new_temp1)));
    m_nan1=m-n_nan1;
    n_nan2=numel(find(isnan(new_temp2)));
    m_nan2=m-n_nan2;
    if m_nan1==0 || m_nan2==0
        idx1=nan;
    else
        use=new_temp2(idx2==0);
        new=imresize(use,[m_nan1,1],'bilinear');
        idx1(idx1==1)=nan;
        idx1(idx1==0)=new;
    end
    GFDL_ESM2G_mrsosP3_new(:,count)=idx1;
    count=count+1;
end
% Differences of SSM_kw between models and observations over different time scales
GFDL_ESM2G_SMPdif1=GFDL_ESM2G_mrsosP1_new-SMP1;
GFDL_ESM2G_SMPdif2=GFDL_ESM2G_mrsosP2_new-SMP2;
GFDL_ESM2G_SMPdif3=GFDL_ESM2G_mrsosP3_new-SMP3;

% % ET_kw of GFDL_ESM2G (hfls)
% Get ET_kw over different time scales
N=52925;
F=((1:N)-1)*Fs/N;
x=F(1:(N+1)/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
[m1,n1]=size(GFDL_ESM2G_hflsA);
[m2,n2]=size(GFDL_ESM2G_mrsosF1);
GFDL_ESM2G_hflsp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=GFDL_ESM2G_hflsA(:,i);
    temp2=isnan(GFDL_ESM2G_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        GFDL_ESM2G_hflsp1(1,count)=nan;
    else
        A=GFDL_ESM2G_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m1+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        GFDL_ESM2G_hflsp1(1,count)=p1(1);
    end
    count=count+1;
end
hflsp11=reshape(GFDL_ESM2G_hflsp1,[n2,m2]);
p11=rot90(hflsp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
GFDL_ESM2G_hflsP1=[P2 P1];
GFDL_ESM2G_hflsp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=GFDL_ESM2G_hflsA(:,i);
    temp2=isnan(GFDL_ESM2G_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        GFDL_ESM2G_hflsp2(1,count)=nan;
    else
        A=GFDL_ESM2G_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m1+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        GFDL_ESM2G_hflsp2(1,count)=p2(1);
    end
    count=count+1;
end
hflsp22=reshape(GFDL_ESM2G_hflsp2,[n2,m2]);
p22=rot90(hflsp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
GFDL_ESM2G_hflsP2=[P2 P1];
GFDL_ESM2G_hflsp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=GFDL_ESM2G_hflsA(:,i);
    temp2=isnan(GFDL_ESM2G_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        GFDL_ESM2G_hflsp3(1,count)=nan;
    else
        A=GFDL_ESM2G_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m1+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        GFDL_ESM2G_hflsp3(1,count)=p3(1);
    end
    count=count+1;
end
hflsp33=reshape(GFDL_ESM2G_hflsp3,[n2,m2]);
p33=rot90(hflsp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
GFDL_ESM2G_hflsP3=[P2 P1];
% Change the spatial resolution of the model's ET_kw to the same as SMAP
[p,~]=size(GFDL_ESM2G_lonData);
[q,~]=size(GFDL_ESM2G_latData);
m=GFDL_ESM2G_lonData(p);
GFDL_ESM2G_lon=linspace(-m/2,m/2,p);
GFDL_ESM2G_lon=GFDL_ESM2G_lon.';
BNU_lon1=GFDL_ESM2G_lon(1:p/2,1);
BNU_lon2=GFDL_ESM2G_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    GFDL_ESM2G_lat_D(count:count+(p-1),1)=GFDL_ESM2G_latData(i,1);
    GFDL_ESM2G_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
GFDL_ESM2G_hflsp1=GFDL_ESM2G_hflsp1.';
GFDL_ESM2G_hflsp2=GFDL_ESM2G_hflsp2.';
GFDL_ESM2G_hflsp3=GFDL_ESM2G_hflsp3.';
GFDL_ESM2G_alpha1=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_hflsp1];
GFDL_ESM2G_alpha2=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_hflsp2];
GFDL_ESM2G_alpha3=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_hflsp3];
[m1,~]=size(GFDL_ESM2G_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
GFDL_ESM2G_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha1(i,1)-SMAP1(j,1))+abs(GFDL_ESM2G_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_P1(j,1)=GFDL_ESM2G_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(GFDL_ESM2G_new_P1,[964,406]);
GFDL_ESM2G_new_hflsP1=new_P1.';
[m1,~]=size(GFDL_ESM2G_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
GFDL_ESM2G_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha2(i,1)-SMAP2(j,1))+abs(GFDL_ESM2G_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_P2(j,1)=GFDL_ESM2G_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(GFDL_ESM2G_new_P2,[964,406]);
GFDL_ESM2G_new_hflsP2=new_mf2.';
[m1,~]=size(GFDL_ESM2G_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
GFDL_ESM2G_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha3(i,1)-SMAP3(j,1))+abs(GFDL_ESM2G_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_P3(j,1)=GFDL_ESM2G_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(GFDL_ESM2G_new_P3,[964,406]);
GFDL_ESM2G_new_hflsP3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
GFDL_ESM2G_hflsP1_new=GFDL_ESM2G_new_hflsP1./temp;
GFDL_ESM2G_hflsP2_new=GFDL_ESM2G_new_hflsP2./temp;
GFDL_ESM2G_hflsP3_new=GFDL_ESM2G_new_hflsP3./temp;
% Differences of ET_kw between models and observations over different time scales
GFDL_ESM2G_EPdif1=GFDL_ESM2G_hflsP1_new-GLEAM_P1;
GFDL_ESM2G_EPdif2=GFDL_ESM2G_hflsP2_new-GLEAM_P2;
GFDL_ESM2G_EPdif3=GFDL_ESM2G_hflsP3_new-GLEAM_P3;

% % P_kw of GFDL_ESM2G (pr)
% Get P_kw over different time scales
N=52925;
F=((1:N)-1)*Fs/N;
x=F(1:(N+1)/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
[m1,n1]=size(GFDL_ESM2G_prA);
[m2,n2]=size(GFDL_ESM2G_mrsosF1);
GFDL_ESM2G_prp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=GFDL_ESM2G_prA(:,i);
    temp2=isnan(GFDL_ESM2G_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        GFDL_ESM2G_prp1(1,count)=nan;
    else
        A=GFDL_ESM2G_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m1+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        GFDL_ESM2G_prp1(1,count)=p1(1);
    end
    count=count+1;
end
prp11=reshape(GFDL_ESM2G_prp1,[n2,m2]);
p11=rot90(prp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
GFDL_ESM2G_prP1=[P2 P1];
GFDL_ESM2G_prp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=GFDL_ESM2G_prA(:,i);
    temp2=isnan(GFDL_ESM2G_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        GFDL_ESM2G_prp2(1,count)=nan;
    else
        A=GFDL_ESM2G_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m1+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        GFDL_ESM2G_prp2(1,count)=p2(1);
    end
    count=count+1;
end
prp22=reshape(GFDL_ESM2G_prp2,[n2,m2]);
p22=rot90(prp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
GFDL_ESM2G_prP2=[P2 P1];
GFDL_ESM2G_prp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=GFDL_ESM2G_prA(:,i);
    temp2=isnan(GFDL_ESM2G_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        GFDL_ESM2G_prp3(1,count)=nan;
    else
        A=GFDL_ESM2G_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m1+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        GFDL_ESM2G_prp3(1,count)=p3(1);
    end
    count=count+1;
end
prp33=reshape(GFDL_ESM2G_prp3,[n2,m2]);
p33=rot90(prp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
GFDL_ESM2G_prP3=[P2 P1];
% Change the spatial resolution of the model's P_kw to the same as SMAP
[p,~]=size(GFDL_ESM2G_lonData);
[q,~]=size(GFDL_ESM2G_latData);
m=GFDL_ESM2G_lonData(p);
GFDL_ESM2G_lon=linspace(-m/2,m/2,p);
GFDL_ESM2G_lon=GFDL_ESM2G_lon.';
BNU_lon1=GFDL_ESM2G_lon(1:p/2,1);
BNU_lon2=GFDL_ESM2G_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    GFDL_ESM2G_lat_D(count:count+(p-1),1)=GFDL_ESM2G_latData(i,1);
    GFDL_ESM2G_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
GFDL_ESM2G_prp1=GFDL_ESM2G_prp1.';
GFDL_ESM2G_prp2=GFDL_ESM2G_prp2.';
GFDL_ESM2G_prp3=GFDL_ESM2G_prp3.';
GFDL_ESM2G_alpha1=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_prp1];
GFDL_ESM2G_alpha2=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_prp2];
GFDL_ESM2G_alpha3=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_prp3];
[m1,~]=size(GFDL_ESM2G_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
GFDL_ESM2G_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha1(i,1)-SMAP1(j,1))+abs(GFDL_ESM2G_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_P1(j,1)=GFDL_ESM2G_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(GFDL_ESM2G_new_P1,[964,406]);
GFDL_ESM2G_new_prP1=new_P1.';
[m1,~]=size(GFDL_ESM2G_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
GFDL_ESM2G_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha2(i,1)-SMAP2(j,1))+abs(GFDL_ESM2G_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_P2(j,1)=GFDL_ESM2G_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(GFDL_ESM2G_new_P2,[964,406]);
GFDL_ESM2G_new_prP2=new_mf2.';
[m1,~]=size(GFDL_ESM2G_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
GFDL_ESM2G_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_alpha3(i,1)-SMAP3(j,1))+abs(GFDL_ESM2G_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_P3(j,1)=GFDL_ESM2G_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(GFDL_ESM2G_new_P3,[964,406]);
GFDL_ESM2G_new_prP3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
GFDL_ESM2G_prP1_new=GFDL_ESM2G_new_prP1./temp;
GFDL_ESM2G_prP2_new=GFDL_ESM2G_new_prP2./temp;
GFDL_ESM2G_prP3_new=GFDL_ESM2G_new_prP3./temp;
% Differences of P_kw between models and observations over different time scales
GFDL_ESM2G_PPdif1=GFDL_ESM2G_prP1_new-ERA5_P1;
GFDL_ESM2G_PPdif2=GFDL_ESM2G_prP2_new-ERA5_P2;
GFDL_ESM2G_PPdif3=GFDL_ESM2G_prP3_new-ERA5_P3;