% % SSM_n of CSIRO_Mk3.6 (mrsos)
% Get SSM_n over different time scales
N=56940;
Fs=1;
A1(29200,18432)=0;
count=1;
for q=1:96
    for p=1:192
        C=mrsosData1234(p,q,:);
        A1(:,count)=C(:);
        count=count+1;
    end
end
A1(A1==3.395080566406250e-04)=nan;
A1=single(A1);
[row,col] = size(A1);
dA1(row,col)=0;
for i=1:col
    dA1(:,i)=detrend(A1(:,i));
end
dA1=single(dA1);
A2(27740,18432)=0;
for q=1:96
    for p=1:192
        C=mrsosData5678(p,q,:);
        A2(:,count)=C(:);
        count=count+1;
    end
end
A2(A2==3.395080566406250e-04)=nan;
A2=single(A2);
[row,col] = size(A2);
dA2(row,col)=0;
for i=1:col
    dA2(:,i)=detrend(A2(:,i));
end
dA2=single(dA2);
A11=dA1(:,1:4608);
A12=dA1(:,4609:9216);
A13=dA1(:,9217:13824);
A14=dA1(:,13825:18432);
A21=dA2(:,1:4608);
A22=dA2(:,4609:9216);
A23=dA2(:,9217:13824);
A24=dA2(:,13825:18432);
A1121=[A11;A21];
Y1=fft(A1121);
Y1=Y1(1:28470,:);
Ayy1=abs(Y1).^2;
A1222=[A12;A22];
Y2=fft(A1222);
Y2=Y2(1:28470,:);
Ayy2=abs(Y2).^2;
A1323=[A13;A23];
Y3=fft(A1323);
Y3=Y3(1:28470,:);
Ayy3=abs(Y3).^2;
A1424=[A14;A24];
Y4=fft(A1424);
Y4=Y4(1:28470,:);
Ayy4=abs(Y4).^2;
CSIRO_Mk_Ayy=[Ayy1,Ayy2,Ayy3,Ayy4];
CSIRO_Mk_Ayy=CSIRO_Mk_Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
CSIRO_Mk_ybar1=sum(CSIRO_Mk_Ayy(b:a,:));
YBAR1=reshape(CSIRO_Mk_ybar1,[192,96]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
CSIRO_Mk_ybar2=sum(CSIRO_Mk_Ayy(c:b,:));
YBAR2=reshape(CSIRO_Mk_ybar2,[192,96]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
CSIRO_Mk_ybar3=sum(CSIRO_Mk_Ayy(d:c,:));
YBAR3=reshape(CSIRO_Mk_ybar3,[192,96]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
CSIRO_Mk_mrsosf1=CSIRO_Mk_ybar1./(CSIRO_Mk_ybar1+CSIRO_Mk_ybar2+CSIRO_Mk_ybar3);
CSIRO_Mk_mrsosf2=CSIRO_Mk_ybar2./(CSIRO_Mk_ybar1+CSIRO_Mk_ybar2+CSIRO_Mk_ybar3);
CSIRO_Mk_mrsosf3=CSIRO_Mk_ybar3./(CSIRO_Mk_ybar1+CSIRO_Mk_ybar2+CSIRO_Mk_ybar3);
CSIRO_Mk_mrsosF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
CSIRO_Mk_mrsosF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
CSIRO_Mk_mrsosF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(CSIRO_Mk_lonData);
[q,~]=size(CSIRO_Mk_latData);
m=CSIRO_Mk_lonData(p);
CSIRO_Mk_lon=linspace(-m/2,m/2,p);
CSIRO_Mk_lon=CSIRO_Mk_lon.';
BNU_lon1=CSIRO_Mk_lon(1:p/2,1);
BNU_lon2=CSIRO_Mk_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    CSIRO_Mk_lat_D(count:count+(p-1),1)=CSIRO_Mk_latData(i,1);
    CSIRO_Mk_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
CSIRO_Mk_mrsosf1=CSIRO_Mk_mrsosf1.';
CSIRO_Mk_mrsosf2=CSIRO_Mk_mrsosf2.';
CSIRO_Mk_mrsosf3=CSIRO_Mk_mrsosf3.';
CSIRO_Mk_alpha1=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_mrsosf1];
CSIRO_Mk_alpha2=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_mrsosf2];
CSIRO_Mk_alpha3=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_mrsosf3];
[m1,~]=size(CSIRO_Mk_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
CSIRO_Mk_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha1(i,1)-SMAP1(j,1))+abs(CSIRO_Mk_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P1(j,1)=CSIRO_Mk_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(CSIRO_Mk_new_P1,[964,406]);
CSIRO_Mk_new_mrsosf1=new_P1.';
[m1,~]=size(CSIRO_Mk_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
CSIRO_Mk_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha2(i,1)-SMAP2(j,1))+abs(CSIRO_Mk_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P2(j,1)=CSIRO_Mk_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(CSIRO_Mk_new_P2,[964,406]);
CSIRO_Mk_new_mrsosf2=new_mf2.';
[m1,~]=size(CSIRO_Mk_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
CSIRO_Mk_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha3(i,1)-SMAP3(j,1))+abs(CSIRO_Mk_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P3(j,1)=CSIRO_Mk_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(CSIRO_Mk_new_P3,[964,406]);
CSIRO_Mk_new_mrsosf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
CSIRO_Mk_mrsosf1_new=CSIRO_Mk_new_mrsosf1./temp;
CSIRO_Mk_mrsosf2_new=CSIRO_Mk_new_mrsosf2./temp;
CSIRO_Mk_mrsosf3_new=CSIRO_Mk_new_mrsosf3./temp;

% % ET_n of CSIRO_Mk3.6 (hfls)
% Get ET_n over different time scales
N=56940;
Fs=1;
A1(29200,18432)=0;
count=1;
for q=1:96
    for p=1:192
        C=hflsData1234(p,q,:);
        A1(:,count)=C(:);
        count=count+1;
    end
end
A1=single(A1);
[row,col] = size(A1);
dA1(row,col)=0;
for i=1:col
    dA1(:,i)=detrend(A1(:,i));
end
dA1=single(dA1);
A2(27740,18432)=0;
for q=1:96
    for p=1:192
        C=hflsData5678(p,q,:);
        A2(:,count)=C(:);
        count=count+1;
    end
end
A2=single(A2);
[row,col] = size(A2);
dA2(row,col)=0;
for i=1:col
    dA2(:,i)=detrend(A2(:,i));
end
dA2=single(dA2);
A11=dA1(:,1:4608);
A12=dA1(:,4609:9216);
A13=dA1(:,9217:13824);
A14=dA1(:,13825:18432);
A21=dA2(:,1:4608);
A22=dA2(:,4609:9216);
A23=dA2(:,9217:13824);
A24=dA2(:,13825:18432);
A1121=[A11;A21];
Y1=fft(A1121);
Y1=Y1(1:28470,:);
Ayy1=abs(Y1).^2;
A1222=[A12;A22];
Y2=fft(A1222);
Y2=Y2(1:28470,:);
Ayy2=abs(Y2).^2;
A1323=[A13;A23];
Y3=fft(A1323);
Y3=Y3(1:28470,:);
Ayy3=abs(Y3).^2;
A1424=[A14;A24];
Y4=fft(A1424);
Y4=Y4(1:28470,:);
Ayy4=abs(Y4).^2;
Ayy=[Ayy1,Ayy2,Ayy3,Ayy4];
Ayy=Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
CSIRO_Mk_ybar1=sum(Ayy(b:a,:));
YBAR1=reshape(CSIRO_Mk_ybar1,[192,96]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
CSIRO_Mk_ybar2=sum(Ayy(c:b,:));
YBAR2=reshape(CSIRO_Mk_ybar2,[192,96]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
CSIRO_Mk_ybar3=sum(Ayy(d:c,:));
YBAR3=reshape(CSIRO_Mk_ybar3,[192,96]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
CSIRO_Mk_hflsf1=CSIRO_Mk_ybar1./(CSIRO_Mk_ybar1+CSIRO_Mk_ybar2+CSIRO_Mk_ybar3);
CSIRO_Mk_hflsf2=CSIRO_Mk_ybar2./(CSIRO_Mk_ybar1+CSIRO_Mk_ybar2+CSIRO_Mk_ybar3);
CSIRO_Mk_hflsf3=CSIRO_Mk_ybar3./(CSIRO_Mk_ybar1+CSIRO_Mk_ybar2+CSIRO_Mk_ybar3);
CSIRO_Mk_hflsF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
CSIRO_Mk_hflsF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
CSIRO_Mk_hflsF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's ET_n to the same as SMAP
[p,~]=size(CSIRO_Mk_lonData);
[q,~]=size(CSIRO_Mk_latData);
m=CSIRO_Mk_lonData(p);
CSIRO_Mk_lon=linspace(-m/2,m/2,p);
CSIRO_Mk_lon=CSIRO_Mk_lon.';
BNU_lon1=CSIRO_Mk_lon(1:p/2,1);
BNU_lon2=CSIRO_Mk_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    CSIRO_Mk_lat_D(count:count+(p-1),1)=CSIRO_Mk_latData(i,1);
    CSIRO_Mk_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
CSIRO_Mk_hflsf1=CSIRO_Mk_hflsf1.';
CSIRO_Mk_hflsf2=CSIRO_Mk_hflsf2.';
CSIRO_Mk_hflsf3=CSIRO_Mk_hflsf3.';
CSIRO_Mk_alpha1=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_hflsf1];
CSIRO_Mk_alpha2=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_hflsf2];
CSIRO_Mk_alpha3=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_hflsf3];
[m1,~]=size(CSIRO_Mk_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
CSIRO_Mk_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha1(i,1)-SMAP1(j,1))+abs(CSIRO_Mk_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P1(j,1)=CSIRO_Mk_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(CSIRO_Mk_new_P1,[964,406]);
CSIRO_Mk_new_hflsf1=new_P1.';
[m1,~]=size(CSIRO_Mk_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
CSIRO_Mk_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha2(i,1)-SMAP2(j,1))+abs(CSIRO_Mk_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P2(j,1)=CSIRO_Mk_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(CSIRO_Mk_new_P2,[964,406]);
CSIRO_Mk_new_hflsf2=new_mf2.';
[m1,~]=size(CSIRO_Mk_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
CSIRO_Mk_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha3(i,1)-SMAP3(j,1))+abs(CSIRO_Mk_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P3(j,1)=CSIRO_Mk_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(CSIRO_Mk_new_P3,[964,406]);
CSIRO_Mk_new_hflsf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
CSIRO_Mk_hflsf1_new=CSIRO_Mk_new_hflsf1./temp;
CSIRO_Mk_hflsf2_new=CSIRO_Mk_new_hflsf2./temp;
CSIRO_Mk_hflsf3_new=CSIRO_Mk_new_hflsf3./temp;

% % P_n of CSIRO_Mk3.6 (pr)
% Get P_n over different time scales
N=56940;
Fs=1;
A1(29200,18432)=0;
count=1;
for q=1:96
    for p=1:192
        C=prData1234(p,q,:);
        A1(:,count)=C(:);
        count=count+1;
    end
    disp(q)
end
A1=single(A1);
[row,col] = size(A1);
dA1(row,col)=0;
for i=1:col
    dA1(:,i)=detrend(A1(:,i));
end
dA1=single(dA1);
count=1;
A2(27740,18432)=0;
for q=1:96
    for p=1:192
        C=prData5678(p,q,:);
        A2(:,count)=C(:);
        count=count+1;
    end
    disp(q)
end
A2=single(A2);
A=[A1;A2];
[row,col] = size(A2);
dA2(row,col)=0;
for i=1:col
    dA2(:,i)=detrend(A2(:,i));
end
dA2=single(dA2);
A11=dA1(:,1:4608);
A12=dA1(:,4609:9216);
A13=dA1(:,9217:13824);
A14=dA1(:,13825:18432);
A21=dA2(:,1:4608);
A22=dA2(:,4609:9216);
A23=dA2(:,9217:13824);
A24=dA2(:,13825:18432);
A1121=[A11;A21];
Y1=fft(A1121);
Y1=Y1(1:28470,:);
Ayy1=abs(Y1).^2;
A1222=[A12;A22];
Y2=fft(A1222);
Y2=Y2(1:28470,:);
Ayy2=abs(Y2).^2;
A1323=[A13;A23];
Y3=fft(A1323);
Y3=Y3(1:28470,:);
Ayy3=abs(Y3).^2;
A1424=[A14;A24];
Y4=fft(A1424);
Y4=Y4(1:28470,:);
Ayy4=abs(Y4).^2;
Ayy=[Ayy1,Ayy2,Ayy3,Ayy4];
Ayy=Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
CSIRO_Mk_ybar1=sum(CSIRO_Mk_Ayy(b:a,:));
YBAR1=reshape(CSIRO_Mk_ybar1,[192,96]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
CSIRO_Mk_ybar2=sum(CSIRO_Mk_Ayy(c:b,:));
YBAR2=reshape(CSIRO_Mk_ybar2,[192,96]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
CSIRO_Mk_ybar3=sum(CSIRO_Mk_Ayy(d:c,:));
YBAR3=reshape(CSIRO_Mk_ybar3,[192,96]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
CSIRO_Mk_prf1=CSIRO_Mk_ybar1./(CSIRO_Mk_ybar1+CSIRO_Mk_ybar2+CSIRO_Mk_ybar3);
CSIRO_Mk_prf2=CSIRO_Mk_ybar2./(CSIRO_Mk_ybar1+CSIRO_Mk_ybar2+CSIRO_Mk_ybar3);
CSIRO_Mk_prf3=CSIRO_Mk_ybar3./(CSIRO_Mk_ybar1+CSIRO_Mk_ybar2+CSIRO_Mk_ybar3);
CSIRO_Mk_prF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
CSIRO_Mk_prF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
CSIRO_Mk_prF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's P_n to the same as SMAP
[p,~]=size(CSIRO_Mk_lonData);
[q,~]=size(CSIRO_Mk_latData);
m=CSIRO_Mk_lonData(p);
CSIRO_Mk_lon=linspace(-m/2,m/2,p);
CSIRO_Mk_lon=CSIRO_Mk_lon.';
BNU_lon1=CSIRO_Mk_lon(1:p/2,1);
BNU_lon2=CSIRO_Mk_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    CSIRO_Mk_lat_D(count:count+(p-1),1)=CSIRO_Mk_latData(i,1);
    CSIRO_Mk_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
CSIRO_Mk_prf1=CSIRO_Mk_prf1.';
CSIRO_Mk_prf2=CSIRO_Mk_prf2.';
CSIRO_Mk_prf3=CSIRO_Mk_prf3.';
CSIRO_Mk_alpha1=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_prf1];
CSIRO_Mk_alpha2=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_prf2];
CSIRO_Mk_alpha3=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_prf3];
[m1,~]=size(CSIRO_Mk_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
CSIRO_Mk_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha1(i,1)-SMAP1(j,1))+abs(CSIRO_Mk_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P1(j,1)=CSIRO_Mk_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(CSIRO_Mk_new_P1,[964,406]);
CSIRO_Mk_new_prf1=new_P1.';
[m1,~]=size(CSIRO_Mk_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
CSIRO_Mk_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha2(i,1)-SMAP2(j,1))+abs(CSIRO_Mk_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P2(j,1)=CSIRO_Mk_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(CSIRO_Mk_new_P2,[964,406]);
CSIRO_Mk_new_prf2=new_mf2.';
[m1,~]=size(CSIRO_Mk_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
CSIRO_Mk_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha3(i,1)-SMAP3(j,1))+abs(CSIRO_Mk_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P3(j,1)=CSIRO_Mk_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(CSIRO_Mk_new_P3,[964,406]);
CSIRO_Mk_new_prf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
CSIRO_Mk_prf1_new=CSIRO_Mk_new_prf1./temp;
CSIRO_Mk_prf2_new=CSIRO_Mk_new_prf2./temp;
CSIRO_Mk_prf3_new=CSIRO_Mk_new_prf3./temp;

% Differences of SSM_n between models and observations over different time scales
CSIRO_Mk_mdif1=CSIRO_Mk_mF1-SMF1;
CSIRO_Mk_mdif2=CSIRO_Mk_mF2-SMF2;
CSIRO_Mk_mdif3=CSIRO_Mk_mF3-SMF3;
[m,n]=size(SMF1);
CSIRO_Mk_mdif1_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF1(:,i);
    temp2=CSIRO_Mk_mdif1(:,i);
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
    CSIRO_Mk_mdif1_try1(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF2);
CSIRO_Mk_mdif2_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF2(:,i);
    temp2=CSIRO_Mk_mdif2(:,i);
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
    CSIRO_Mk_mdif2_try1(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF3);
CSIRO_Mk_mdif3_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF3(:,i);
    temp2=CSIRO_Mk_mdif3(:,i);
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
    CSIRO_Mk_mdif3_try1(:,count)=idx1;
    count=count+1;
end

% % SSM_kw of CSIRO_Mk3.6 (mrsos)
% Get SSM_kw over different time scales
[m1,n1]=size(CSIRO_Mk_mrsosA);
[m2,n2]=size(CSIRO_Mk_mrsosF1);
CSIRO_Mk_mrsosp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=CSIRO_Mk_mrsosA(:,i);
    temp2=isnan(CSIRO_Mk_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CSIRO_Mk_mrsosp1(1,count)=nan;
    else
        A=CSIRO_Mk_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        CSIRO_Mk_mrsosp1(1,count)=p1(1);
    end
    count=count+1;
end
mrsosp11=reshape(CSIRO_Mk_mrsosp1,[n2,m2]);
p11=rot90(mrsosp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
CSIRO_Mk_mrsosP1=[P2 P1];
CSIRO_Mk_mrsosp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=CSIRO_Mk_mrsosA(:,i);
    temp2=isnan(CSIRO_Mk_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CSIRO_Mk_mrsosp2(1,count)=nan;
    else
        A=CSIRO_Mk_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        CSIRO_Mk_mrsosp2(1,count)=p2(1);
    end
    count=count+1;
end
mrsosp22=reshape(CSIRO_Mk_mrsosp2,[n2,m2]);
p22=rot90(mrsosp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
CSIRO_Mk_mrsosP2=[P2 P1];
CSIRO_Mk_mrsosp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=CSIRO_Mk_mrsosA(:,i);
    temp2=isnan(CSIRO_Mk_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CSIRO_Mk_mrsosp3(1,count)=nan;
    else
        A=CSIRO_Mk_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        CSIRO_Mk_mrsosp3(1,count)=p3(1);
    end
    count=count+1;
end
mrsosp33=reshape(CSIRO_Mk_mrsosp3,[n2,m2]);
p33=rot90(mrsosp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
CSIRO_Mk_mrsosP3=[P2 P1];
% Change the spatial resolution of the model's SSM_kw to the same as SMAP
[p,~]=size(CSIRO_Mk_lonData);
[q,~]=size(CSIRO_Mk_latData);
m=CSIRO_Mk_lonData(p);
CSIRO_Mk_lon=linspace(-m/2,m/2,p);
CSIRO_Mk_lon=CSIRO_Mk_lon.';
BNU_lon1=CSIRO_Mk_lon(1:p/2,1);
BNU_lon2=CSIRO_Mk_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    CSIRO_Mk_lat_D(count:count+(p-1),1)=CSIRO_Mk_latData(i,1);
    CSIRO_Mk_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
CSIRO_Mk_mrsosp1=CSIRO_Mk_mrsosp1.';
CSIRO_Mk_mrsosp2=CSIRO_Mk_mrsosp2.';
CSIRO_Mk_mrsosp3=CSIRO_Mk_mrsosp3.';
CSIRO_Mk_alpha1=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_mrsosp1];
CSIRO_Mk_alpha2=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_mrsosp2];
CSIRO_Mk_alpha3=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_mrsosp3];
[m1,~]=size(CSIRO_Mk_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
CSIRO_Mk_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha1(i,1)-SMAP1(j,1))+abs(CSIRO_Mk_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P1(j,1)=CSIRO_Mk_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(CSIRO_Mk_new_P1,[964,406]);
CSIRO_Mk_new_mrsosP1=new_P1.';
[m1,~]=size(CSIRO_Mk_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
CSIRO_Mk_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha2(i,1)-SMAP2(j,1))+abs(CSIRO_Mk_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P2(j,1)=CSIRO_Mk_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(CSIRO_Mk_new_P2,[964,406]);
CSIRO_Mk_new_mrsosP2=new_mf2.';
[m1,~]=size(CSIRO_Mk_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
CSIRO_Mk_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha3(i,1)-SMAP3(j,1))+abs(CSIRO_Mk_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P3(j,1)=CSIRO_Mk_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(CSIRO_Mk_new_P3,[964,406]);
CSIRO_Mk_new_mrsosP3=new_mf3.';
[m,n]=size(SMF1);
CSIRO_Mk_mrsosP1_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF1(:,i);
    new_temp2=CSIRO_Mk_new_mrsosP1(:,i);
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
    CSIRO_Mk_mrsosP1_new(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF2);
CSIRO_Mk_mrsosP2_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF2(:,i);
    new_temp2=CSIRO_Mk_new_mrsosP2(:,i);
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
    CSIRO_Mk_mrsosP2_new(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF3);
CSIRO_Mk_mrsosP3_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF3(:,i);
    new_temp2=CSIRO_Mk_new_mrsosP3(:,i);
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
    CSIRO_Mk_mrsosP3_new(:,count)=idx1;
    count=count+1;
end
% Differences of SSM_kw between models and observations over different time scales
CSIRO_Mk_SMPdif1=CSIRO_Mk_mrsosP1_new-SMP1;
CSIRO_Mk_SMPdif2=CSIRO_Mk_mrsosP2_new-SMP2;
CSIRO_Mk_SMPdif3=CSIRO_Mk_mrsosP3_new-SMP3;

% % ET_kw of CSIRO_Mk3.6 (hfls)
% Get ET_kw over different time scales
N=56940;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
[m1,n1]=size(CSIRO_Mk_hflsA);
[m2,n2]=size(CSIRO_Mk_mrsosF1);
CSIRO_Mk_hflsp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=CSIRO_Mk_hflsA(:,i);
    temp2=isnan(CSIRO_Mk_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CSIRO_Mk_hflsp1(1,count)=nan;
    else
        A=CSIRO_Mk_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        CSIRO_Mk_hflsp1(1,count)=p1(1);
    end
    count=count+1;
end
hflsp11=reshape(CSIRO_Mk_hflsp1,[n2,m2]);
p11=rot90(hflsp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
CSIRO_Mk_hflsP1=[P2 P1];
CSIRO_Mk_hflsp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=CSIRO_Mk_hflsA(:,i);
    temp2=isnan(CSIRO_Mk_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CSIRO_Mk_hflsp2(1,count)=nan;
    else
        A=CSIRO_Mk_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        CSIRO_Mk_hflsp2(1,count)=p2(1);
    end
    count=count+1;
end
hflsp22=reshape(CSIRO_Mk_hflsp2,[n2,m2]);
p22=rot90(hflsp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
CSIRO_Mk_hflsP2=[P2 P1];
CSIRO_Mk_hflsp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=CSIRO_Mk_hflsA(:,i);
    temp2=isnan(CSIRO_Mk_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CSIRO_Mk_hflsp3(1,count)=nan;
    else
        A=CSIRO_Mk_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        CSIRO_Mk_hflsp3(1,count)=p3(1);
    end
    count=count+1;
end
hflsp33=reshape(CSIRO_Mk_hflsp3,[n2,m2]);
p33=rot90(hflsp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
CSIRO_Mk_hflsP3=[P2 P1];
% Change the spatial resolution of the model's ET_kw to the same as SMAP
[p,~]=size(CSIRO_Mk_lonData);
[q,~]=size(CSIRO_Mk_latData);
m=CSIRO_Mk_lonData(p);
CSIRO_Mk_lon=linspace(-m/2,m/2,p);
CSIRO_Mk_lon=CSIRO_Mk_lon.';
BNU_lon1=CSIRO_Mk_lon(1:p/2,1);
BNU_lon2=CSIRO_Mk_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    CSIRO_Mk_lat_D(count:count+(p-1),1)=CSIRO_Mk_latData(i,1);
    CSIRO_Mk_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
CSIRO_Mk_hflsp1=CSIRO_Mk_hflsp1.';
CSIRO_Mk_hflsp2=CSIRO_Mk_hflsp2.';
CSIRO_Mk_hflsp3=CSIRO_Mk_hflsp3.';
CSIRO_Mk_alpha1=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_hflsp1];
CSIRO_Mk_alpha2=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_hflsp2];
CSIRO_Mk_alpha3=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_hflsp3];
[m1,~]=size(CSIRO_Mk_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
CSIRO_Mk_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha1(i,1)-SMAP1(j,1))+abs(CSIRO_Mk_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P1(j,1)=CSIRO_Mk_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(CSIRO_Mk_new_P1,[964,406]);
CSIRO_Mk_new_hflsP1=new_P1.';
[m1,~]=size(CSIRO_Mk_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
CSIRO_Mk_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha2(i,1)-SMAP2(j,1))+abs(CSIRO_Mk_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P2(j,1)=CSIRO_Mk_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(CSIRO_Mk_new_P2,[964,406]);
CSIRO_Mk_new_hflsP2=new_mf2.';
[m1,~]=size(CSIRO_Mk_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
CSIRO_Mk_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha3(i,1)-SMAP3(j,1))+abs(CSIRO_Mk_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P3(j,1)=CSIRO_Mk_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(CSIRO_Mk_new_P3,[964,406]);
CSIRO_Mk_new_hflsP3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
CSIRO_Mk_hflsP1_new=CSIRO_Mk_new_hflsP1./temp;
CSIRO_Mk_hflsP2_new=CSIRO_Mk_new_hflsP2./temp;
CSIRO_Mk_hflsP3_new=CSIRO_Mk_new_hflsP3./temp;
% Differences of ET_kw between models and observations over different time scales
CSIRO_Mk_EPdif1=CSIRO_Mk_hflsP1_new-GLEAM_P1;
CSIRO_Mk_EPdif2=CSIRO_Mk_hflsP2_new-GLEAM_P2;
CSIRO_Mk_EPdif3=CSIRO_Mk_hflsP3_new-GLEAM_P3;

% % P_kw of CSIRO_Mk3.6 (pr)
% Get P_kw over different time scales
N=56940;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
e=find(abs(x-1/1458)<=0.00001);
[m1,n1]=size(CSIRO_Mk_prA);
[m2,n2]=size(CSIRO_Mk_mrsosF1);
CSIRO_Mk_prp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=CSIRO_Mk_prA(:,i);
    temp2=isnan(CSIRO_Mk_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CSIRO_Mk_prp1(1,count)=nan;
    else
        A=CSIRO_Mk_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        CSIRO_Mk_prp1(1,count)=p1(1);
    end
    count=count+1;
end
prp11=reshape(CSIRO_Mk_prp1,[n2,m2]);
p11=rot90(prp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
CSIRO_Mk_prP1=[P2 P1];
CSIRO_Mk_prp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=CSIRO_Mk_prA(:,i);
    temp2=isnan(CSIRO_Mk_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CSIRO_Mk_prp2(1,count)=nan;
    else
        A=CSIRO_Mk_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        CSIRO_Mk_prp2(1,count)=p2(1);
    end
    count=count+1;
end
prp22=reshape(CSIRO_Mk_prp2,[n2,m2]);
p22=rot90(prp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
CSIRO_Mk_prP2=[P2 P1];
CSIRO_Mk_prp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=CSIRO_Mk_prA(:,i);
    temp2=isnan(CSIRO_Mk_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CSIRO_Mk_prp3(1,count)=nan;
    else
        A=CSIRO_Mk_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        CSIRO_Mk_prp3(1,count)=p3(1);
    end
    count=count+1;
end
prp33=reshape(CSIRO_Mk_prp3,[n2,m2]);
p33=rot90(prp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
CSIRO_Mk_prP3=[P2 P1];
% Change the spatial resolution of the model's P_kw to the same as SMAP
[p,~]=size(CSIRO_Mk_lonData);
[q,~]=size(CSIRO_Mk_latData);
m=CSIRO_Mk_lonData(p);
CSIRO_Mk_lon=linspace(-m/2,m/2,p);
CSIRO_Mk_lon=CSIRO_Mk_lon.';
BNU_lon1=CSIRO_Mk_lon(1:p/2,1);
BNU_lon2=CSIRO_Mk_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    CSIRO_Mk_lat_D(count:count+(p-1),1)=CSIRO_Mk_latData(i,1);
    CSIRO_Mk_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
CSIRO_Mk_prp1=CSIRO_Mk_prp1.';
CSIRO_Mk_prp2=CSIRO_Mk_prp2.';
CSIRO_Mk_prp3=CSIRO_Mk_prp3.';
CSIRO_Mk_alpha1=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_prp1];
CSIRO_Mk_alpha2=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_prp2];
CSIRO_Mk_alpha3=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_prp3];
[m1,~]=size(CSIRO_Mk_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
CSIRO_Mk_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha1(i,1)-SMAP1(j,1))+abs(CSIRO_Mk_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P1(j,1)=CSIRO_Mk_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(CSIRO_Mk_new_P1,[964,406]);
CSIRO_Mk_new_prP1=new_P1.';
[m1,~]=size(CSIRO_Mk_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
CSIRO_Mk_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha2(i,1)-SMAP2(j,1))+abs(CSIRO_Mk_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P2(j,1)=CSIRO_Mk_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(CSIRO_Mk_new_P2,[964,406]);
CSIRO_Mk_new_prP2=new_mf2.';
[m1,~]=size(CSIRO_Mk_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
CSIRO_Mk_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_alpha3(i,1)-SMAP3(j,1))+abs(CSIRO_Mk_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_P3(j,1)=CSIRO_Mk_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(CSIRO_Mk_new_P3,[964,406]);
CSIRO_Mk_new_prP3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
CSIRO_Mk_prP1_new=CSIRO_Mk_new_prP1./temp;
CSIRO_Mk_prP2_new=CSIRO_Mk_new_prP2./temp;
CSIRO_Mk_prP3_new=CSIRO_Mk_new_prP3./temp;
% Differences of P_kw between models and observations over different time scales
CSIRO_Mk_PPdif1=CSIRO_Mk_prP1_new-ERA5_P1;
CSIRO_Mk_PPdif2=CSIRO_Mk_prP2_new-ERA5_P2;
CSIRO_Mk_PPdif3=CSIRO_Mk_prP3_new-ERA5_P3;