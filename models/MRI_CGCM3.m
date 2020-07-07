% % SSM_n of MRI_CGCM3 (mrsos)
% Get SSM_n over different time scales
N=20454;
Fs=1;
A1(10957,51200)=0;
count=1;
for q=1:160
    for p=1:320
        C=mrsosData123(p,q,:);
        A1(:,count)=C(:);
        count=count+1;
    end
end
A1(A1==0)=nan;
A1=single(A1);
[row,col] = size(A1);
dA1(row,col)=0;
for i=1:col
    dA1(:,i)=detrend(A1(:,i));
end
dA1=single(dA1);
A2(9497,51200)=0;
for q=1:160
    for p=1:320
        C=mrsosData456(p,q,:);
        A2(:,count)=C(:);
        count=count+1;
    end
end
A2(A2==0)=nan;
A2=single(A2);
[row,col] = size(A2);
dA2(row,col)=0;
for i=1:col
    dA2(:,i)=detrend(A2(:,i));
end
dA2=single(dA2);
A11=dA1(:,1:12800);
A12=dA1(:,12801:25600);
A13=dA1(:,25601:38400);
A14=dA1(:,38401:51200);
A21=dA2(:,1:12800);
A22=dA2(:,12801:25600);
A23=dA2(:,25601:38400);
A24=dA2(:,38401:51200);
A1121=[A11;A21];
Y1=fft(A1121);
Y1=Y1(1:10227,:);
Ayy1=abs(Y1).^2;
A1222=[A12;A22];
Y2=fft(A1222);
Y2=Y2(1:10227,:);
Ayy2=abs(Y2).^2;
A1323=[A13;A23];
Y3=fft(A1323);
Y3=Y3(1:10227,:);
Ayy3=abs(Y3).^2;
A1424=[A14;A24];
Y4=fft(A1424);
Y4=Y4(1:10227,:);
Ayy4=abs(Y4).^2;
MRI_CGCM3_Ayy=[Ayy1,Ayy2,Ayy3,Ayy4];
MRI_CGCM3_Ayy=MRI_CGCM3_Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00002);
d=find(abs(x-1/365)<=0.00001);
MRI_CGCM3_ybar1=sum(MRI_CGCM3_Ayy(b:a,:));
YBAR1=reshape(MRI_CGCM3_ybar1,[320,160]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
MRI_CGCM3_ybar2=sum(MRI_CGCM3_Ayy(c:b,:));
YBAR2=reshape(MRI_CGCM3_ybar2,[320,160]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
MRI_CGCM3_ybar3=sum(MRI_CGCM3_Ayy(d:c,:));
YBAR3=reshape(MRI_CGCM3_ybar3,[320,160]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
MRI_CGCM3_mrsosf1=MRI_CGCM3_ybar1./(MRI_CGCM3_ybar1+MRI_CGCM3_ybar2+MRI_CGCM3_ybar3);
MRI_CGCM3_mrsosf2=MRI_CGCM3_ybar2./(MRI_CGCM3_ybar1+MRI_CGCM3_ybar2+MRI_CGCM3_ybar3);
MRI_CGCM3_mrsosf3=MRI_CGCM3_ybar3./(MRI_CGCM3_ybar1+MRI_CGCM3_ybar2+MRI_CGCM3_ybar3);
MRI_CGCM3_mrsosF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
MRI_CGCM3_mrsosF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
MRI_CGCM3_mrsosF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(MRI_CGCM3_lonData);
[q,~]=size(MRI_CGCM3_latData);
m=MRI_CGCM3_lonData(p);
MRI_CGCM3_lon=linspace(-m/2,m/2,p);
MRI_CGCM3_lon=MRI_CGCM3_lon.';
BNU_lon1=MRI_CGCM3_lon(1:p/2,1);
BNU_lon2=MRI_CGCM3_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    MRI_CGCM3_lat_D(count:count+(p-1),1)=MRI_CGCM3_latData(i,1);
    MRI_CGCM3_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
MRI_CGCM3_mrsosf1=MRI_CGCM3_mrsosf1.';
MRI_CGCM3_mrsosf2=MRI_CGCM3_mrsosf2.';
MRI_CGCM3_mrsosf3=MRI_CGCM3_mrsosf3.';
MRI_CGCM3_alpha1=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_mrsosf1];
MRI_CGCM3_alpha2=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_mrsosf2];
MRI_CGCM3_alpha3=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_mrsosf3];
[m1,~]=size(MRI_CGCM3_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
MRI_CGCM3_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha1(i,1)-SMAP1(j,1))+abs(MRI_CGCM3_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P1(j,1)=MRI_CGCM3_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(MRI_CGCM3_new_P1,[964,406]);
MRI_CGCM3_new_mrsosf1=new_P1.';
[m1,~]=size(MRI_CGCM3_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
MRI_CGCM3_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha2(i,1)-SMAP2(j,1))+abs(MRI_CGCM3_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P2(j,1)=MRI_CGCM3_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(MRI_CGCM3_new_P2,[964,406]);
MRI_CGCM3_new_mrsosf2=new_mf2.';
[m1,~]=size(MRI_CGCM3_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
MRI_CGCM3_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha3(i,1)-SMAP3(j,1))+abs(MRI_CGCM3_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P3(j,1)=MRI_CGCM3_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(MRI_CGCM3_new_P3,[964,406]);
MRI_CGCM3_new_mrsosf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
MRI_CGCM3_mrsosf1_new=MRI_CGCM3_new_mrsosf1./temp;
MRI_CGCM3_mrsosf2_new=MRI_CGCM3_new_mrsosf2./temp;
MRI_CGCM3_mrsosf3_new=MRI_CGCM3_new_mrsosf3./temp;

% % ET_n of MRI_CGCM3 (hfls)
% Get ET_n over different time scales
N=20454;
Fs=1;
A1(10957,51200)=0;
count=1;
for q=1:160
    for p=1:320
        C=hflsData123(p,q,:);
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
A2(9497,51200)=0;
for q=1:160
    for p=1:320
        C=hflsData456(p,q,:);
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
A11=dA1(:,1:12800);
A12=dA1(:,12801:25600);
A13=dA1(:,25601:38400);
A14=dA1(:,38401:51200);
A21=dA2(:,1:12800);
A22=dA2(:,12801:25600);
A23=dA2(:,25601:38400);
A24=dA2(:,38401:51200);
A1121=[A11;A21];
Y1=fft(A1121);
Y1=Y1(1:10227,:);
Ayy1=abs(Y1).^2;
A1222=[A12;A22];
Y2=fft(A1222);
Y2=Y2(1:10227,:);
Ayy2=abs(Y2).^2;
A1323=[A13;A23];
Y3=fft(A1323);
Y3=Y3(1:10227,:);
Ayy3=abs(Y3).^2;
A1424=[A14;A24];
Y4=fft(A1424);
Y4=Y4(1:10227,:);
Ayy4=abs(Y4).^2;
MRI_CGCM3_Ayy=[Ayy1,Ayy2,Ayy3,Ayy4];
MRI_CGCM3_Ayy=MRI_CGCM3_Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00002);
d=find(abs(x-1/365)<=0.00001);
MRI_CGCM3_ybar1=sum(MRI_CGCM3_Ayy(b:a,:));
YBAR1=reshape(MRI_CGCM3_ybar1,[320,160]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
MRI_CGCM3_ybar2=sum(MRI_CGCM3_Ayy(c:b,:));
YBAR2=reshape(MRI_CGCM3_ybar2,[320,160]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
MRI_CGCM3_ybar3=sum(MRI_CGCM3_Ayy(d:c,:));
YBAR3=reshape(MRI_CGCM3_ybar3,[320,160]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
MRI_CGCM3_hflsf1=MRI_CGCM3_ybar1./(MRI_CGCM3_ybar1+MRI_CGCM3_ybar2+MRI_CGCM3_ybar3);
MRI_CGCM3_hflsf2=MRI_CGCM3_ybar2./(MRI_CGCM3_ybar1+MRI_CGCM3_ybar2+MRI_CGCM3_ybar3);
MRI_CGCM3_hflsf3=MRI_CGCM3_ybar3./(MRI_CGCM3_ybar1+MRI_CGCM3_ybar2+MRI_CGCM3_ybar3);
MRI_CGCM3_hflsF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
MRI_CGCM3_hflsF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
MRI_CGCM3_hflsF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's ET_n to the same as SMAP
[p,~]=size(MRI_CGCM3_lonData);
[q,~]=size(MRI_CGCM3_latData);
m=MRI_CGCM3_lonData(p);
MRI_CGCM3_lon=linspace(-m/2,m/2,p);
MRI_CGCM3_lon=MRI_CGCM3_lon.';
BNU_lon1=MRI_CGCM3_lon(1:p/2,1);
BNU_lon2=MRI_CGCM3_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    MRI_CGCM3_lat_D(count:count+(p-1),1)=MRI_CGCM3_latData(i,1);
    MRI_CGCM3_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
MRI_CGCM3_hflsf1=MRI_CGCM3_hflsf1.';
MRI_CGCM3_hflsf2=MRI_CGCM3_hflsf2.';
MRI_CGCM3_hflsf3=MRI_CGCM3_hflsf3.';
MRI_CGCM3_alpha1=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_hflsf1];
MRI_CGCM3_alpha2=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_hflsf2];
MRI_CGCM3_alpha3=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_hflsf3];
[m1,~]=size(MRI_CGCM3_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
MRI_CGCM3_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha1(i,1)-SMAP1(j,1))+abs(MRI_CGCM3_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P1(j,1)=MRI_CGCM3_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(MRI_CGCM3_new_P1,[964,406]);
MRI_CGCM3_new_hflsf1=new_P1.';
[m1,~]=size(MRI_CGCM3_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
MRI_CGCM3_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha2(i,1)-SMAP2(j,1))+abs(MRI_CGCM3_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P2(j,1)=MRI_CGCM3_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(MRI_CGCM3_new_P2,[964,406]);
MRI_CGCM3_new_hflsf2=new_mf2.';
[m1,~]=size(MRI_CGCM3_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
MRI_CGCM3_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha3(i,1)-SMAP3(j,1))+abs(MRI_CGCM3_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P3(j,1)=MRI_CGCM3_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(MRI_CGCM3_new_P3,[964,406]);
MRI_CGCM3_new_hflsf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
MRI_CGCM3_hflsf1_new=MRI_CGCM3_new_hflsf1./temp;
MRI_CGCM3_hflsf2_new=MRI_CGCM3_new_hflsf2./temp;
MRI_CGCM3_hflsf3_new=MRI_CGCM3_new_hflsf3./temp;

% % P_n of MRI_CGCM3 (pr)
% Get P_n over different time scales
N=20454;
Fs=1;
A1(10957,51200)=0;
count=1;
for q=1:160
    for p=1:320
        C=prData123(p,q,:);
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
A2(9497,51200)=0;
count=1;
for q=1:160
    for p=1:320
        C=prData456(p,q,:);
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
A11=dA1(:,1:12800);
A12=dA1(:,12801:25600);
A13=dA1(:,25601:38400);
A14=dA1(:,38401:51200);
A21=dA2(:,1:12800);
A22=dA2(:,12801:25600);
A23=dA2(:,25601:38400);
A24=dA2(:,38401:51200);
A1121=[A11;A21];
Y1=fft(A1121);
Y1=Y1(1:10227,:);
Ayy1=abs(Y1).^2;
A1222=[A12;A22];
Y2=fft(A1222);
Y2=Y2(1:10227,:);
Ayy2=abs(Y2).^2;
A1323=[A13;A23];
Y3=fft(A1323);
Y3=Y3(1:10227,:);
Ayy3=abs(Y3).^2;
A1424=[A14;A24];
Y4=fft(A1424);
Y4=Y4(1:10227,:);
Ayy4=abs(Y4).^2;
MRI_CGCM3_Ayy=[Ayy1,Ayy2,Ayy3,Ayy4];
MRI_CGCM3_Ayy=MRI_CGCM3_Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00002);
d=find(abs(x-1/365)<=0.00001);
MRI_CGCM3_ybar1=sum(MRI_CGCM3_Ayy(b:a,:));
YBAR1=reshape(MRI_CGCM3_ybar1,[320,160]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
MRI_CGCM3_ybar2=sum(MRI_CGCM3_Ayy(c:b,:));
YBAR2=reshape(MRI_CGCM3_ybar2,[320,160]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
MRI_CGCM3_ybar3=sum(MRI_CGCM3_Ayy(d:c,:));
YBAR3=reshape(MRI_CGCM3_ybar3,[320,160]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
MRI_CGCM3_prf1=MRI_CGCM3_ybar1./(MRI_CGCM3_ybar1+MRI_CGCM3_ybar2+MRI_CGCM3_ybar3);
MRI_CGCM3_prf2=MRI_CGCM3_ybar2./(MRI_CGCM3_ybar1+MRI_CGCM3_ybar2+MRI_CGCM3_ybar3);
MRI_CGCM3_prf3=MRI_CGCM3_ybar3./(MRI_CGCM3_ybar1+MRI_CGCM3_ybar2+MRI_CGCM3_ybar3);
MRI_CGCM3_prF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
MRI_CGCM3_prF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
MRI_CGCM3_prF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's P_n to the same as SMAP
[p,~]=size(MRI_CGCM3_lonData);
[q,~]=size(MRI_CGCM3_latData);
m=MRI_CGCM3_lonData(p);
MRI_CGCM3_lon=linspace(-m/2,m/2,p);
MRI_CGCM3_lon=MRI_CGCM3_lon.';
BNU_lon1=MRI_CGCM3_lon(1:p/2,1);
BNU_lon2=MRI_CGCM3_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    MRI_CGCM3_lat_D(count:count+(p-1),1)=MRI_CGCM3_latData(i,1);
    MRI_CGCM3_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
MRI_CGCM3_prf1=MRI_CGCM3_prf1.';
MRI_CGCM3_prf2=MRI_CGCM3_prf2.';
MRI_CGCM3_prf3=MRI_CGCM3_prf3.';
MRI_CGCM3_alpha1=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_prf1];
MRI_CGCM3_alpha2=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_prf2];
MRI_CGCM3_alpha3=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_prf3];
[m1,~]=size(MRI_CGCM3_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
MRI_CGCM3_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha1(i,1)-SMAP1(j,1))+abs(MRI_CGCM3_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P1(j,1)=MRI_CGCM3_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(MRI_CGCM3_new_P1,[964,406]);
MRI_CGCM3_new_prf1=new_P1.';
[m1,~]=size(MRI_CGCM3_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
MRI_CGCM3_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha2(i,1)-SMAP2(j,1))+abs(MRI_CGCM3_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P2(j,1)=MRI_CGCM3_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(MRI_CGCM3_new_P2,[964,406]);
MRI_CGCM3_new_prf2=new_mf2.';
[m1,~]=size(MRI_CGCM3_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
MRI_CGCM3_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha3(i,1)-SMAP3(j,1))+abs(MRI_CGCM3_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P3(j,1)=MRI_CGCM3_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(MRI_CGCM3_new_P3,[964,406]);
MRI_CGCM3_new_prf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
MRI_CGCM3_prf1_new=MRI_CGCM3_new_prf1./temp;
MRI_CGCM3_prf2_new=MRI_CGCM3_new_prf2./temp;
MRI_CGCM3_prf3_new=MRI_CGCM3_new_prf3./temp;

% Differences of SSM_n between models and observations over different time scales
MRI_CGCM3_mdif1=MRI_CGCM3_mF1-SMF1;
MRI_CGCM3_mdif2=MRI_CGCM3_mF2-SMF2;
MRI_CGCM3_mdif3=MRI_CGCM3_mF3-SMF3;
[m,n]=size(SMF1);
MRI_CGCM3_mdif1_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF1(:,i);
    temp2=MRI_CGCM3_mdif1(:,i);
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
    MRI_CGCM3_mdif1_try1(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF2);
MRI_CGCM3_mdif2_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF2(:,i);
    temp2=MRI_CGCM3_mdif2(:,i);
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
    MRI_CGCM3_mdif2_try1(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF3);
MRI_CGCM3_mdif3_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF3(:,i);
    temp2=MRI_CGCM3_mdif3(:,i);
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
    MRI_CGCM3_mdif3_try1(:,count)=idx1;
    count=count+1;
end

% % SSM_kw of MRI_CGCM3 (mrsos)
% Get SSM_kw over different time scales
N=20454;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00002);
d=find(abs(x-1/365)<=0.00001);
[m1,n1]=size(MRI_CGCM3_mrsosA);
[m2,n2]=size(MRI_CGCM3_mrsosF1);
MRI_CGCM3_mrsosp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=MRI_CGCM3_mrsosA(:,i);
    temp2=isnan(MRI_CGCM3_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        MRI_CGCM3_mrsosp1(1,count)=nan;
    else
        A=MRI_CGCM3_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        MRI_CGCM3_mrsosp1(1,count)=p1(1);
    end
    count=count+1;
end
mrsosp11=reshape(MRI_CGCM3_mrsosp1,[n2,m2]);
p11=rot90(mrsosp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
MRI_CGCM3_mrsosP1=[P2 P1];
MRI_CGCM3_mrsosp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=MRI_CGCM3_mrsosA(:,i);
    temp2=isnan(MRI_CGCM3_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        MRI_CGCM3_mrsosp2(1,count)=nan;
    else
        A=MRI_CGCM3_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        MRI_CGCM3_mrsosp2(1,count)=p2(1);
    end
    count=count+1;
end
mrsosp22=reshape(MRI_CGCM3_mrsosp2,[n2,m2]);
p22=rot90(mrsosp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
MRI_CGCM3_mrsosP2=[P2 P1];
MRI_CGCM3_mrsosp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=MRI_CGCM3_mrsosA(:,i);
    temp2=isnan(MRI_CGCM3_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        MRI_CGCM3_mrsosp3(1,count)=nan;
    else
        A=MRI_CGCM3_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        MRI_CGCM3_mrsosp3(1,count)=p3(1);
    end
    count=count+1;
end
mrsosp33=reshape(MRI_CGCM3_mrsosp3,[n2,m2]);
p33=rot90(mrsosp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
MRI_CGCM3_mrsosP3=[P2 P1];
% Change the spatial resolution of the model's SSM_kw to the same as SMAP
[p,~]=size(MRI_CGCM3_lonData);
[q,~]=size(MRI_CGCM3_latData);
m=MRI_CGCM3_lonData(p);
MRI_CGCM3_lon=linspace(-m/2,m/2,p);
MRI_CGCM3_lon=MRI_CGCM3_lon.';
BNU_lon1=MRI_CGCM3_lon(1:p/2,1);
BNU_lon2=MRI_CGCM3_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    MRI_CGCM3_lat_D(count:count+(p-1),1)=MRI_CGCM3_latData(i,1);
    MRI_CGCM3_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
MRI_CGCM3_mrsosp1=MRI_CGCM3_mrsosp1.';
MRI_CGCM3_mrsosp2=MRI_CGCM3_mrsosp2.';
MRI_CGCM3_mrsosp3=MRI_CGCM3_mrsosp3.';
MRI_CGCM3_alpha1=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_mrsosp1];
MRI_CGCM3_alpha2=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_mrsosp2];
MRI_CGCM3_alpha3=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_mrsosp3];
[m1,~]=size(MRI_CGCM3_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
MRI_CGCM3_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha1(i,1)-SMAP1(j,1))+abs(MRI_CGCM3_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P1(j,1)=MRI_CGCM3_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(MRI_CGCM3_new_P1,[964,406]);
MRI_CGCM3_new_mrsosP1=new_P1.';
[m1,~]=size(MRI_CGCM3_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
MRI_CGCM3_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha2(i,1)-SMAP2(j,1))+abs(MRI_CGCM3_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P2(j,1)=MRI_CGCM3_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(MRI_CGCM3_new_P2,[964,406]);
MRI_CGCM3_new_mrsosP2=new_mf2.';
[m1,~]=size(MRI_CGCM3_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
MRI_CGCM3_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha3(i,1)-SMAP3(j,1))+abs(MRI_CGCM3_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P3(j,1)=MRI_CGCM3_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(MRI_CGCM3_new_P3,[964,406]);
MRI_CGCM3_new_mrsosP3=new_mf3.';
[m,n]=size(SMF1);
MRI_CGCM3_mrsosP1_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF1(:,i);
    new_temp2=MRI_CGCM3_new_mrsosP1(:,i);
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
    MRI_CGCM3_mrsosP1_new(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF2);
MRI_CGCM3_mrsosP2_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF2(:,i);
    new_temp2=MRI_CGCM3_new_mrsosP2(:,i);
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
    MRI_CGCM3_mrsosP2_new(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF3);
MRI_CGCM3_mrsosP3_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF3(:,i);
    new_temp2=MRI_CGCM3_new_mrsosP3(:,i);
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
    MRI_CGCM3_mrsosP3_new(:,count)=idx1;
    count=count+1;
end
% Differences of SSM_kw between models and observations over different time scales
MRI_CGCM3_SMPdif1=MRI_CGCM3_mrsosP1_new-SMP1;
MRI_CGCM3_SMPdif2=MRI_CGCM3_mrsosP2_new-SMP2;
MRI_CGCM3_SMPdif3=MRI_CGCM3_mrsosP3_new-SMP3;

% % ET_kw of MRI_CGCM3 (hfls)
% Get ET_kw over different time scales
N=20454;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00002);
d=find(abs(x-1/365)<=0.00001);
[m1,n1]=size(MRI_CGCM3_hflsA);
[m2,n2]=size(MRI_CGCM3_mrsosF1);
MRI_CGCM3_hflsp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=MRI_CGCM3_hflsA(:,i);
    temp2=isnan(MRI_CGCM3_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        MRI_CGCM3_hflsp1(1,count)=nan;
    else
        A=MRI_CGCM3_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        MRI_CGCM3_hflsp1(1,count)=p1(1);
    end
    count=count+1;
end
hflsp11=reshape(MRI_CGCM3_hflsp1,[n2,m2]);
p11=rot90(hflsp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
MRI_CGCM3_hflsP1=[P2 P1];
MRI_CGCM3_hflsp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=MRI_CGCM3_hflsA(:,i);
    temp2=isnan(MRI_CGCM3_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        MRI_CGCM3_hflsp2(1,count)=nan;
    else
        A=MRI_CGCM3_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        MRI_CGCM3_hflsp2(1,count)=p2(1);
    end
    count=count+1;
end
hflsp22=reshape(MRI_CGCM3_hflsp2,[n2,m2]);
p22=rot90(hflsp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
MRI_CGCM3_hflsP2=[P2 P1];
MRI_CGCM3_hflsp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=MRI_CGCM3_hflsA(:,i);
    temp2=isnan(MRI_CGCM3_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        MRI_CGCM3_hflsp3(1,count)=nan;
    else
        A=MRI_CGCM3_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        MRI_CGCM3_hflsp3(1,count)=p3(1);
    end
    count=count+1;
end
hflsp33=reshape(MRI_CGCM3_hflsp3,[n2,m2]);
p33=rot90(hflsp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
MRI_CGCM3_hflsP3=[P2 P1];
% Change the spatial resolution of the model's ET_kw to the same as SMAP
[p,~]=size(MRI_CGCM3_lonData);
[q,~]=size(MRI_CGCM3_latData);
m=MRI_CGCM3_lonData(p);
MRI_CGCM3_lon=linspace(-m/2,m/2,p);
MRI_CGCM3_lon=MRI_CGCM3_lon.';
BNU_lon1=MRI_CGCM3_lon(1:p/2,1);
BNU_lon2=MRI_CGCM3_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    MRI_CGCM3_lat_D(count:count+(p-1),1)=MRI_CGCM3_latData(i,1);
    MRI_CGCM3_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
MRI_CGCM3_hflsp1=MRI_CGCM3_hflsp1.';
MRI_CGCM3_hflsp2=MRI_CGCM3_hflsp2.';
MRI_CGCM3_hflsp3=MRI_CGCM3_hflsp3.';
MRI_CGCM3_alpha1=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_hflsp1];
MRI_CGCM3_alpha2=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_hflsp2];
MRI_CGCM3_alpha3=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_hflsp3];
[m1,~]=size(MRI_CGCM3_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
MRI_CGCM3_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha1(i,1)-SMAP1(j,1))+abs(MRI_CGCM3_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P1(j,1)=MRI_CGCM3_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(MRI_CGCM3_new_P1,[964,406]);
MRI_CGCM3_new_hflsP1=new_P1.';
[m1,~]=size(MRI_CGCM3_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
MRI_CGCM3_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha2(i,1)-SMAP2(j,1))+abs(MRI_CGCM3_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P2(j,1)=MRI_CGCM3_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(MRI_CGCM3_new_P2,[964,406]);
MRI_CGCM3_new_hflsP2=new_mf2.';
[m1,~]=size(MRI_CGCM3_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
MRI_CGCM3_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha3(i,1)-SMAP3(j,1))+abs(MRI_CGCM3_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P3(j,1)=MRI_CGCM3_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(MRI_CGCM3_new_P3,[964,406]);
MRI_CGCM3_new_hflsP3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
MRI_CGCM3_hflsP1_new=MRI_CGCM3_new_hflsP1./temp;
MRI_CGCM3_hflsP2_new=MRI_CGCM3_new_hflsP2./temp;
MRI_CGCM3_hflsP3_new=MRI_CGCM3_new_hflsP3./temp;
% Differences of ET_kw between models and observations over different time scales
MRI_CGCM3_EPdif1=MRI_CGCM3_hflsP1_new-GLEAM_P1;
MRI_CGCM3_EPdif2=MRI_CGCM3_hflsP2_new-GLEAM_P2;
MRI_CGCM3_EPdif3=MRI_CGCM3_hflsP3_new-GLEAM_P3;

% % P_kw of MRI_CGCM3 (pr)
% Get P_kw over different time scales
N=20454;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00002);
d=find(abs(x-1/365)<=0.00001);
[m1,n1]=size(MRI_CGCM3_prA);
[m2,n2]=size(MRI_CGCM3_mrsosF1);
MRI_CGCM3_prp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=MRI_CGCM3_prA(:,i);
    temp2=isnan(MRI_CGCM3_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        MRI_CGCM3_prp1(1,count)=nan;
    else
        A=MRI_CGCM3_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        MRI_CGCM3_prp1(1,count)=p1(1);
    end
    count=count+1;
end
prp11=reshape(MRI_CGCM3_prp1,[n2,m2]);
p11=rot90(prp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
MRI_CGCM3_prP1=[P2 P1];
MRI_CGCM3_prp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=MRI_CGCM3_prA(:,i);
    temp2=isnan(MRI_CGCM3_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        MRI_CGCM3_prp2(1,count)=nan;
    else
        A=MRI_CGCM3_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        MRI_CGCM3_prp2(1,count)=p2(1);
    end
    count=count+1;
end
prp22=reshape(MRI_CGCM3_prp2,[n2,m2]);
p22=rot90(prp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
MRI_CGCM3_prP2=[P2 P1];
MRI_CGCM3_prp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=MRI_CGCM3_prA(:,i);
    temp2=isnan(MRI_CGCM3_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        MRI_CGCM3_prp3(1,count)=nan;
    else
        A=MRI_CGCM3_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        MRI_CGCM3_prp3(1,count)=p3(1);
    end
    count=count+1;
end
prp33=reshape(MRI_CGCM3_prp3,[n2,m2]);
p33=rot90(prp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
MRI_CGCM3_prP3=[P2 P1];
% Change the spatial resolution of the model's P_kw to the same as SMAP
[p,~]=size(MRI_CGCM3_lonData);
[q,~]=size(MRI_CGCM3_latData);
m=MRI_CGCM3_lonData(p);
MRI_CGCM3_lon=linspace(-m/2,m/2,p);
MRI_CGCM3_lon=MRI_CGCM3_lon.';
BNU_lon1=MRI_CGCM3_lon(1:p/2,1);
BNU_lon2=MRI_CGCM3_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    MRI_CGCM3_lat_D(count:count+(p-1),1)=MRI_CGCM3_latData(i,1);
    MRI_CGCM3_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
MRI_CGCM3_prp1=MRI_CGCM3_prp1.';
MRI_CGCM3_prp2=MRI_CGCM3_prp2.';
MRI_CGCM3_prp3=MRI_CGCM3_prp3.';
MRI_CGCM3_alpha1=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_prp1];
MRI_CGCM3_alpha2=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_prp2];
MRI_CGCM3_alpha3=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_prp3];
[m1,~]=size(MRI_CGCM3_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
MRI_CGCM3_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha1(i,1)-SMAP1(j,1))+abs(MRI_CGCM3_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P1(j,1)=MRI_CGCM3_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(MRI_CGCM3_new_P1,[964,406]);
MRI_CGCM3_new_prP1=new_P1.';
[m1,~]=size(MRI_CGCM3_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
MRI_CGCM3_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha2(i,1)-SMAP2(j,1))+abs(MRI_CGCM3_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P2(j,1)=MRI_CGCM3_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(MRI_CGCM3_new_P2,[964,406]);
MRI_CGCM3_new_prP2=new_mf2.';
[m1,~]=size(MRI_CGCM3_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
MRI_CGCM3_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_alpha3(i,1)-SMAP3(j,1))+abs(MRI_CGCM3_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_P3(j,1)=MRI_CGCM3_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(MRI_CGCM3_new_P3,[964,406]);
MRI_CGCM3_new_prP3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
MRI_CGCM3_prP1_new=MRI_CGCM3_new_prP1./temp;
MRI_CGCM3_prP2_new=MRI_CGCM3_new_prP2./temp;
MRI_CGCM3_prP3_new=MRI_CGCM3_new_prP3./temp;
% Differences of P_kw between models and observations over different time scales
MRI_CGCM3_PPdif1=MRI_CGCM3_prP1_new-ERA5_P1;
MRI_CGCM3_PPdif2=MRI_CGCM3_prP2_new-ERA5_P2;
MRI_CGCM3_PPdif3=MRI_CGCM3_prP3_new-ERA5_P3;