% % SSM_n of CanESM2 (mrsos)
% Get SSM_n over different time scales
N=56940;
Fs=1;
A(56940,8192)=0;
count=1;
for q=1:64
    for p=1:128
        C=mrsosData(p,q,:);
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
A1=dA(:,1:2048);
Y1=fft(A1);
Y1=Y1(1:28470,:);
Ayy1=abs(Y1).^2;
A2=dA(:,2049:4096);
Y2=fft(A2);
Y2=Y2(1:28470,:);
Ayy2=abs(Y2).^2;
A3=dA(:,4097:6144);
Y3=fft(A3);
Y3=Y3(1:28470,:);
Ayy3=abs(Y3).^2;
A4=dA(:,6145:8192);
Y4=fft(A4);
Y4=Y4(1:28470,:);
Ayy4=abs(Y4).^2;
CanESM2_Ayy=[Ayy1,Ayy2,Ayy3,Ayy4];
CanESM2_Ayy=CanESM2_Ayy/(N/2);
CanESM2_Ayy(Ayy==0)=nan;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
CanESM2_ybar1=sum(CanESM2_Ayy(b:a,:));
YBAR1=reshape(CanESM2_ybar1,[128,64]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
CanESM2_ybar2=sum(CanESM2_Ayy(c:b,:));
YBAR2=reshape(CanESM2_ybar2,[128,64]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
CanESM2_ybar3=sum(CanESM2_Ayy(d:c,:));
YBAR3=reshape(CanESM2_ybar3,[128,64]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
CanESM2_mrsosf1=CanESM2_ybar1./(CanESM2_ybar1+CanESM2_ybar2+CanESM2_ybar3);
CanESM2_mrsosf2=CanESM2_ybar2./(CanESM2_ybar1+CanESM2_ybar2+CanESM2_ybar3);
CanESM2_mrsosf3=CanESM2_ybar3./(CanESM2_ybar1+CanESM2_ybar2+CanESM2_ybar3);
CanESM2_mrsosF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
CanESM2_mrsosF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
CanESM2_mrsosF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(CanESM2_lonData);
[q,~]=size(CanESM2_latData);
m=CanESM2_lonData(p);
CanESM2_lon=linspace(-m/2,m/2,p);
CanESM2_lon=CanESM2_lon.';
BNU_lon1=CanESM2_lon(1:p/2,1);
BNU_lon2=CanESM2_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    CanESM2_lat_D(count:count+(p-1),1)=CanESM2_latData(i,1);
    CanESM2_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
CanESM2_mrsosf1=CanESM2_mrsosf1.';
CanESM2_mrsosf2=CanESM2_mrsosf2.';
CanESM2_mrsosf3=CanESM2_mrsosf3.';
CanESM2_alpha1=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_mrsosf1];
CanESM2_alpha2=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_mrsosf2];
CanESM2_alpha3=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_mrsosf3];
[m1,~]=size(CanESM2_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
CanESM2_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha1(i,1)-SMAP1(j,1))+abs(CanESM2_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P1(j,1)=CanESM2_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(CanESM2_new_P1,[964,406]);
CanESM2_new_mrsosf1=new_P1.';
[m1,~]=size(CanESM2_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
CanESM2_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha2(i,1)-SMAP2(j,1))+abs(CanESM2_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P2(j,1)=CanESM2_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(CanESM2_new_P2,[964,406]);
CanESM2_new_mrsosf2=new_mf2.';
[m1,~]=size(CanESM2_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
CanESM2_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha3(i,1)-SMAP3(j,1))+abs(CanESM2_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P3(j,1)=CanESM2_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(CanESM2_new_P3,[964,406]);
CanESM2_new_mrsosf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
CanESM2_mrsosf1_new=CanESM2_new_mrsosf1./temp;
CanESM2_mrsosf2_new=CanESM2_new_mrsosf2./temp;
CanESM2_mrsosf3_new=CanESM2_new_mrsosf3./temp;

% % ET_n of CanESM2 (hfls)
% Get ET_n over different time scales
N=56940;
Fs=1;
A(56940,8192)=0;
count=1;
for q=1:64
    for p=1:128
        C=hflsData(p,q,:);
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
A1=dA(:,1:2048);
Y1=fft(A1);
Y1=Y1(1:28470,:);
Ayy1=abs(Y1).^2;
A2=dA(:,2049:4096);
Y2=fft(A2);
Y2=Y2(1:28470,:);
Ayy2=abs(Y2).^2;
A3=dA(:,4097:6144);
Y3=fft(A3);
Y3=Y3(1:28470,:);
Ayy3=abs(Y3).^2;
A4=dA(:,6145:8192);
Y4=fft(A4);
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
CanESM2_ybar1=sum(Ayy(b:a,:));
YBAR1=reshape(CanESM2_ybar1,[128,64]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
CanESM2_ybar2=sum(Ayy(c:b,:));
YBAR2=reshape(CanESM2_ybar2,[128,64]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
CanESM2_ybar3=sum(Ayy(d:c,:));
YBAR3=reshape(CanESM2_ybar3,[128,64]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
CanESM2_hflsf1=CanESM2_ybar1./(CanESM2_ybar1+CanESM2_ybar2+CanESM2_ybar3);
CanESM2_hflsf2=CanESM2_ybar2./(CanESM2_ybar1+CanESM2_ybar2+CanESM2_ybar3);
CanESM2_hflsf3=CanESM2_ybar3./(CanESM2_ybar1+CanESM2_ybar2+CanESM2_ybar3);
CanESM2_hflsF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
CanESM2_hflsF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
CanESM2_hflsF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's ET_n to the same as SMAP
[p,~]=size(CanESM2_lonData);
[q,~]=size(CanESM2_latData);
m=CanESM2_lonData(p);
CanESM2_lon=linspace(-m/2,m/2,p);
CanESM2_lon=CanESM2_lon.';
BNU_lon1=CanESM2_lon(1:p/2,1);
BNU_lon2=CanESM2_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    CanESM2_lat_D(count:count+(p-1),1)=CanESM2_latData(i,1);
    CanESM2_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
CanESM2_hflsf1=CanESM2_hflsf1.';
CanESM2_hflsf2=CanESM2_hflsf2.';
CanESM2_hflsf3=CanESM2_hflsf3.';
CanESM2_alpha1=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_hflsf1];
CanESM2_alpha2=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_hflsf2];
CanESM2_alpha3=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_hflsf3];
[m1,~]=size(CanESM2_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
CanESM2_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha1(i,1)-SMAP1(j,1))+abs(CanESM2_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P1(j,1)=CanESM2_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(CanESM2_new_P1,[964,406]);
CanESM2_new_hflsf1=new_P1.';
[m1,~]=size(CanESM2_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
CanESM2_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha2(i,1)-SMAP2(j,1))+abs(CanESM2_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P2(j,1)=CanESM2_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(CanESM2_new_P2,[964,406]);
CanESM2_new_hflsf2=new_mf2.';
[m1,~]=size(CanESM2_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
CanESM2_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha3(i,1)-SMAP3(j,1))+abs(CanESM2_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P3(j,1)=CanESM2_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(CanESM2_new_P3,[964,406]);
CanESM2_new_hflsf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
CanESM2_hflsf1_new=CanESM2_new_hflsf1./temp;
CanESM2_hflsf2_new=CanESM2_new_hflsf2./temp;
CanESM2_hflsf3_new=CanESM2_new_hflsf3./temp;

% % P_n of CanESM2 (pr)
% Get P_n over different time scales
N=56940;
Fs=1;
A(56940,8192)=0;
count=1;
for q=1:64
    for p=1:128
        C=prData(p,q,:);
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
A1=dA(:,1:2048);
Y1=fft(A1);
Y1=Y1(1:28470,:);
Ayy1=abs(Y1).^2;
A2=dA(:,2049:4096);
Y2=fft(A2);
Y2=Y2(1:28470,:);
Ayy2=abs(Y2).^2;
A3=dA(:,4097:6144);
Y3=fft(A3);
Y3=Y3(1:28470,:);
Ayy3=abs(Y3).^2;
A4=dA(:,6145:8192);
Y4=fft(A4);
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
CanESM2_ybar1=sum(CanESM2_Ayy(b:a,:));
YBAR1=reshape(CanESM2_ybar1,[128,64]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
CanESM2_ybar2=sum(CanESM2_Ayy(c:b,:));
YBAR2=reshape(CanESM2_ybar2,[128,64]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
CanESM2_ybar3=sum(CanESM2_Ayy(d:c,:));
YBAR3=reshape(CanESM2_ybar3,[128,64]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
CanESM2_prf1=CanESM2_ybar1./(CanESM2_ybar1+CanESM2_ybar2+CanESM2_ybar3);
CanESM2_prf2=CanESM2_ybar2./(CanESM2_ybar1+CanESM2_ybar2+CanESM2_ybar3);
CanESM2_prf3=CanESM2_ybar3./(CanESM2_ybar1+CanESM2_ybar2+CanESM2_ybar3);
CanESM2_prF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
CanESM2_prF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
CanESM2_prF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's P_n to the same as SMAP
[p,~]=size(CanESM2_lonData);
[q,~]=size(CanESM2_latData);
m=CanESM2_lonData(p);
CanESM2_lon=linspace(-m/2,m/2,p);
CanESM2_lon=CanESM2_lon.';
BNU_lon1=CanESM2_lon(1:p/2,1);
BNU_lon2=CanESM2_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    CanESM2_lat_D(count:count+(p-1),1)=CanESM2_latData(i,1);
    CanESM2_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
CanESM2_prf1=CanESM2_prf1.';
CanESM2_prf2=CanESM2_prf2.';
CanESM2_prf3=CanESM2_prf3.';
CanESM2_alpha1=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_prf1];
CanESM2_alpha2=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_prf2];
CanESM2_alpha3=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_prf3];
[m1,~]=size(CanESM2_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
CanESM2_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha1(i,1)-SMAP1(j,1))+abs(CanESM2_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P1(j,1)=CanESM2_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(CanESM2_new_P1,[964,406]);
CanESM2_new_prf1=new_P1.';
[m1,~]=size(CanESM2_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
CanESM2_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha2(i,1)-SMAP2(j,1))+abs(CanESM2_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P2(j,1)=CanESM2_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(CanESM2_new_P2,[964,406]);
CanESM2_new_prf2=new_mf2.';
[m1,~]=size(CanESM2_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
CanESM2_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha3(i,1)-SMAP3(j,1))+abs(CanESM2_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P3(j,1)=CanESM2_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(CanESM2_new_P3,[964,406]);
CanESM2_new_prf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
CanESM2_prf1_new=CanESM2_new_prf1./temp;
CanESM2_prf2_new=CanESM2_new_prf2./temp;
CanESM2_prf3_new=CanESM2_new_prf3./temp;

% Differences of SSM_n between models and observations over different time scales
CanESM2_mdif1=CanESM2_mF1-SMF1;
CanESM2_mdif2=CanESM2_mF2-SMF2;
CanESM2_mdif3=CanESM2_mF3-SMF3;
[m,n]=size(SMF1);
CanESM2_mdif1_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF1(:,i);
    temp2=CanESM2_mdif1(:,i);
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
    CanESM2_mdif1_try1(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF2);
CanESM2_mdif2_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF2(:,i);
    temp2=CanESM2_mdif2(:,i);
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
    CanESM2_mdif2_try1(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF3);
CanESM2_mdif3_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF3(:,i);
    temp2=CanESM2_mdif3(:,i);
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
    CanESM2_mdif3_try1(:,count)=idx1;
    count=count+1;
end

% % SSM_kw of CanESM2 (mrsos)
% Get SSM_kw over different time scales
N=56940;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
[m1,n1]=size(CanESM2_mrsosA);
[m2,n2]=size(CanESM2_mrsosF1);
CanESM2_mrsosp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=CanESM2_mrsosA(:,i);
    temp2=isnan(CanESM2_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CanESM2_mrsosp1(1,count)=nan;
    else
        A=CanESM2_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        CanESM2_mrsosp1(1,count)=p1(1);
    end
    count=count+1;
end
mrsosp11=reshape(CanESM2_mrsosp1,[n2,m2]);
p11=rot90(mrsosp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
CanESM2_mrsosP1=[P2 P1];
CanESM2_mrsosp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=CanESM2_mrsosA(:,i);
    temp2=isnan(CanESM2_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CanESM2_mrsosp2(1,count)=nan;
    else
        A=CanESM2_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        CanESM2_mrsosp2(1,count)=p2(1);
    end
    count=count+1;
end
mrsosp22=reshape(CanESM2_mrsosp2,[n2,m2]);
p22=rot90(mrsosp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
CanESM2_mrsosP2=[P2 P1];
CanESM2_mrsosp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=CanESM2_mrsosA(:,i);
    temp2=isnan(CanESM2_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CanESM2_mrsosp3(1,count)=nan;
    else
        A=CanESM2_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        CanESM2_mrsosp3(1,count)=p3(1);
    end
    count=count+1;
end
mrsosp33=reshape(CanESM2_mrsosp3,[n2,m2]);
p33=rot90(mrsosp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
CanESM2_mrsosP3=[P2 P1];
% Change the spatial resolution of the model's SSM_kw to the same as SMAP
[p,~]=size(CanESM2_lonData);
[q,~]=size(CanESM2_latData);
m=CanESM2_lonData(p);
CanESM2_lon=linspace(-m/2,m/2,p);
CanESM2_lon=CanESM2_lon.';
BNU_lon1=CanESM2_lon(1:p/2,1);
BNU_lon2=CanESM2_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    CanESM2_lat_D(count:count+(p-1),1)=CanESM2_latData(i,1);
    CanESM2_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
CanESM2_mrsosp1=CanESM2_mrsosp1.';
CanESM2_mrsosp2=CanESM2_mrsosp2.';
CanESM2_mrsosp3=CanESM2_mrsosp3.';
CanESM2_alpha1=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_mrsosp1];
CanESM2_alpha2=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_mrsosp2];
CanESM2_alpha3=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_mrsosp3];
[m1,~]=size(CanESM2_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
CanESM2_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha1(i,1)-SMAP1(j,1))+abs(CanESM2_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P1(j,1)=CanESM2_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(CanESM2_new_P1,[964,406]);
CanESM2_new_mrsosP1=new_P1.';
[m1,~]=size(CanESM2_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
CanESM2_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha2(i,1)-SMAP2(j,1))+abs(CanESM2_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P2(j,1)=CanESM2_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(CanESM2_new_P2,[964,406]);
CanESM2_new_mrsosP2=new_mf2.';
[m1,~]=size(CanESM2_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
CanESM2_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha3(i,1)-SMAP3(j,1))+abs(CanESM2_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P3(j,1)=CanESM2_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(CanESM2_new_P3,[964,406]);
CanESM2_new_mrsosP3=new_mf3.';
[m,n]=size(SMF1);
CanESM2_mrsosP1_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF1(:,i);
    new_temp2=CanESM2_new_mrsosP1(:,i);
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
    CanESM2_mrsosP1_new(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF2);
CanESM2_mrsosP2_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF2(:,i);
    new_temp2=CanESM2_new_mrsosP2(:,i);
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
    CanESM2_mrsosP2_new(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF3);
CanESM2_mrsosP3_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF3(:,i);
    new_temp2=CanESM2_new_mrsosP3(:,i);
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
    CanESM2_mrsosP3_new(:,count)=idx1;
    count=count+1;
end
% Differences of SSM_kw between models and observations over different time scales
CanESM2_SMPdif1=CanESM2_mrsosP1_new-SMP1;
CanESM2_SMPdif2=CanESM2_mrsosP2_new-SMP2;
CanESM2_SMPdif3=CanESM2_mrsosP3_new-SMP3;

% % ET_kw of CanESM2 (hfls)
% Get ET_kw over different time scales
N=56940;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
[m1,n1]=size(CanESM2_hflsA);
[m2,n2]=size(CanESM2_mrsosF1);
CanESM2_hflsp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=CanESM2_hflsA(:,i);
    temp2=isnan(CanESM2_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CanESM2_hflsp1(1,count)=nan;
    else
        A=CanESM2_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        CanESM2_hflsp1(1,count)=p1(1);
    end
    count=count+1;
end
hflsp11=reshape(CanESM2_hflsp1,[n2,m2]);
p11=rot90(hflsp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
CanESM2_hflsP1=[P2 P1];
CanESM2_hflsp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=CanESM2_hflsA(:,i);
    temp2=isnan(CanESM2_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CanESM2_hflsp2(1,count)=nan;
    else
        A=CanESM2_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        CanESM2_hflsp2(1,count)=p2(1);
    end
    count=count+1;
end
hflsp22=reshape(CanESM2_hflsp2,[n2,m2]);
p22=rot90(hflsp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
CanESM2_hflsP2=[P2 P1];
CanESM2_hflsp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=CanESM2_hflsA(:,i);
    temp2=isnan(CanESM2_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CanESM2_hflsp3(1,count)=nan;
    else
        A=CanESM2_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        CanESM2_hflsp3(1,count)=p3(1);
    end
    count=count+1;
end
hflsp33=reshape(CanESM2_hflsp3,[n2,m2]);
p33=rot90(hflsp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
CanESM2_hflsP3=[P2 P1];
% Change the spatial resolution of the model's ET_kw to the same as SMAP
[p,~]=size(CanESM2_lonData);
[q,~]=size(CanESM2_latData);
m=CanESM2_lonData(p);
CanESM2_lon=linspace(-m/2,m/2,p);
CanESM2_lon=CanESM2_lon.';
BNU_lon1=CanESM2_lon(1:p/2,1);
BNU_lon2=CanESM2_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    CanESM2_lat_D(count:count+(p-1),1)=CanESM2_latData(i,1);
    CanESM2_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
CanESM2_hflsp1=CanESM2_hflsp1.';
CanESM2_hflsp2=CanESM2_hflsp2.';
CanESM2_hflsp3=CanESM2_hflsp3.';
CanESM2_alpha1=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_hflsp1];
CanESM2_alpha2=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_hflsp2];
CanESM2_alpha3=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_hflsp3];
[m1,~]=size(CanESM2_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
CanESM2_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha1(i,1)-SMAP1(j,1))+abs(CanESM2_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P1(j,1)=CanESM2_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(CanESM2_new_P1,[964,406]);
CanESM2_new_hflsP1=new_P1.';
[m1,~]=size(CanESM2_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
CanESM2_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha2(i,1)-SMAP2(j,1))+abs(CanESM2_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P2(j,1)=CanESM2_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(CanESM2_new_P2,[964,406]);
CanESM2_new_hflsP2=new_mf2.';
[m1,~]=size(CanESM2_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
CanESM2_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha3(i,1)-SMAP3(j,1))+abs(CanESM2_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P3(j,1)=CanESM2_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(CanESM2_new_P3,[964,406]);
CanESM2_new_hflsP3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
CanESM2_hflsP1_new=CanESM2_new_hflsP1./temp;
CanESM2_hflsP2_new=CanESM2_new_hflsP2./temp;
CanESM2_hflsP3_new=CanESM2_new_hflsP3./temp;
% Differences of ET_kw between models and observations over different time scales
CanESM2_EPdif1=CanESM2_hflsP1_new-GLEAM_P1;
CanESM2_EPdif2=CanESM2_hflsP2_new-GLEAM_P2;
CanESM2_EPdif3=CanESM2_hflsP3_new-GLEAM_P3;

% % P_kw of CanESM2 (pr)
% Get P_kw over different time scales
N=56940;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
[m1,n1]=size(CanESM2_prA);
[m2,n2]=size(CanESM2_mrsosF1);
CanESM2_prp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=CanESM2_prA(:,i);
    temp2=isnan(CanESM2_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CanESM2_prp1(1,count)=nan;
    else
        A=CanESM2_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        CanESM2_prp1(1,count)=p1(1);
    end
    count=count+1;
end
prp11=reshape(CanESM2_prp1,[n2,m2]);
p11=rot90(prp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
CanESM2_prP1=[P2 P1];
CanESM2_prp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=CanESM2_prA(:,i);
    temp2=isnan(CanESM2_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CanESM2_prp2(1,count)=nan;
    else
        A=CanESM2_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        CanESM2_prp2(1,count)=p2(1);
    end
    count=count+1;
end
prp22=reshape(CanESM2_prp2,[n2,m2]);
p22=rot90(prp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
CanESM2_prP2=[P2 P1];
CanESM2_prp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=CanESM2_prA(:,i);
    temp2=isnan(CanESM2_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        CanESM2_prp3(1,count)=nan;
    else
        A=CanESM2_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        CanESM2_prp3(1,count)=p3(1);
    end
    count=count+1;
end
prp33=reshape(CanESM2_prp3,[n2,m2]);
p33=rot90(prp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
CanESM2_prP3=[P2 P1];
% Change the spatial resolution of the model's P_kw to the same as SMAP
[p,~]=size(CanESM2_lonData);
[q,~]=size(CanESM2_latData);
m=CanESM2_lonData(p);
CanESM2_lon=linspace(-m/2,m/2,p);
CanESM2_lon=CanESM2_lon.';
BNU_lon1=CanESM2_lon(1:p/2,1);
BNU_lon2=CanESM2_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    CanESM2_lat_D(count:count+(p-1),1)=CanESM2_latData(i,1);
    CanESM2_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
CanESM2_prp1=CanESM2_prp1.';
CanESM2_prp2=CanESM2_prp2.';
CanESM2_prp3=CanESM2_prp3.';
CanESM2_alpha1=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_prp1];
CanESM2_alpha2=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_prp2];
CanESM2_alpha3=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_prp3];
[m1,~]=size(CanESM2_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
CanESM2_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha1(i,1)-SMAP1(j,1))+abs(CanESM2_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P1(j,1)=CanESM2_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(CanESM2_new_P1,[964,406]);
CanESM2_new_prP1=new_P1.';
[m1,~]=size(CanESM2_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
CanESM2_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha2(i,1)-SMAP2(j,1))+abs(CanESM2_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P2(j,1)=CanESM2_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(CanESM2_new_P2,[964,406]);
CanESM2_new_prP2=new_mf2.';
[m1,~]=size(CanESM2_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
CanESM2_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_alpha3(i,1)-SMAP3(j,1))+abs(CanESM2_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_P3(j,1)=CanESM2_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(CanESM2_new_P3,[964,406]);
CanESM2_new_prP3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
CanESM2_prP1_new=CanESM2_new_prP1./temp;
CanESM2_prP2_new=CanESM2_new_prP2./temp;
CanESM2_prP3_new=CanESM2_new_prP3./temp;
% Differences of P_kw between models and observations over different time scales
CanESM2_PPdif1=CanESM2_prP1_new-ERA5_P1;
CanESM2_PPdif2=CanESM2_prP2_new-ERA5_P2;
CanESM2_PPdif3=CanESM2_prP3_new-ERA5_P3;