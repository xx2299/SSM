% % SSM_n of NorESM1_M (mrsos)
% Get SSM_n over different time scales
N=20440;
Fs=1;
A(20440,13824)=0;
count=1;
for q=1:96
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
A1=dA(:,1:4608);
Y1=fft(A1);
Y1=Y1(1:10220,:);
A2=dA(:,4609:9216);
Y2=fft(A2);
Y2=Y2(1:10220,:);
A3=dA(:,9217:13824);
Y3=fft(A3);
Y3=Y3(1:10220,:);
Y=[Y1,Y2,Y3];
NorESM1_M_Ayy=abs(Y).^2;
NorESM1_M_Ayy=NorESM1_M_Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00002);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
NorESM1_M_ybar1=sum(NorESM1_M_Ayy(b:a,:));
YBAR1=reshape(NorESM1_M_ybar1,[144,96]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
NorESM1_M_ybar2=sum(NorESM1_M_Ayy(c:b,:));
YBAR2=reshape(NorESM1_M_ybar2,[144,96]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
NorESM1_M_ybar3=sum(NorESM1_M_Ayy(d:c,:));
YBAR3=reshape(NorESM1_M_ybar3,[144,96]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
NorESM1_M_mrsosf1=NorESM1_M_ybar1./(NorESM1_M_ybar1+NorESM1_M_ybar2+NorESM1_M_ybar3);
NorESM1_M_mrsosf2=NorESM1_M_ybar2./(NorESM1_M_ybar1+NorESM1_M_ybar2+NorESM1_M_ybar3);
NorESM1_M_mrsosf3=NorESM1_M_ybar3./(NorESM1_M_ybar1+NorESM1_M_ybar2+NorESM1_M_ybar3);
NorESM1_M_mrsosF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
NorESM1_M_mrsosF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
NorESM1_M_mrsosF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(NorESM1_M_lonData);
[q,~]=size(NorESM1_M_latData);
m=NorESM1_M_lonData(p);
NorESM1_M_lon=linspace(-m/2,m/2,p);
NorESM1_M_lon=NorESM1_M_lon.';
BNU_lon1=NorESM1_M_lon(1:p/2,1);
BNU_lon2=NorESM1_M_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    NorESM1_M_lat_D(count:count+(p-1),1)=NorESM1_M_latData(i,1);
    NorESM1_M_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
NorESM1_M_mrsosf1=NorESM1_M_mrsosf1.';
NorESM1_M_mrsosf2=NorESM1_M_mrsosf2.';
NorESM1_M_mrsosf3=NorESM1_M_mrsosf3.';
NorESM1_M_alpha1=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_mrsosf1];
NorESM1_M_alpha2=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_mrsosf2];
NorESM1_M_alpha3=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_mrsosf3];
[m1,~]=size(NorESM1_M_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
NorESM1_M_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha1(i,1)-SMAP1(j,1))+abs(NorESM1_M_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P1(j,1)=NorESM1_M_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(NorESM1_M_new_P1,[964,406]);
NorESM1_M_new_mrsosf1=new_P1.';
[m1,~]=size(NorESM1_M_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
NorESM1_M_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha2(i,1)-SMAP2(j,1))+abs(NorESM1_M_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P2(j,1)=NorESM1_M_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(NorESM1_M_new_P2,[964,406]);
NorESM1_M_new_mrsosf2=new_mf2.';
[m1,~]=size(NorESM1_M_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
NorESM1_M_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha3(i,1)-SMAP3(j,1))+abs(NorESM1_M_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P3(j,1)=NorESM1_M_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(NorESM1_M_new_P3,[964,406]);
NorESM1_M_new_mrsosf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
NorESM1_M_mrsosf1_new=NorESM1_M_new_mrsosf1./temp;
NorESM1_M_mrsosf2_new=NorESM1_M_new_mrsosf2./temp;
NorESM1_M_mrsosf3_new=NorESM1_M_new_mrsosf3./temp;

% % ET_n of NorESM1_M (hfls)
% Get ET_n over different time scales
N=20440;
Fs=1;
A(20440,13824)=0;
count=1;
for q=1:96
    for p=1:144
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
A1=dA(:,1:4608);
Y1=fft(A1);
Y1=Y1(1:10220,:);
A2=dA(:,4609:9216);
Y2=fft(A2);
Y2=Y2(1:10220,:);
A3=dA(:,9217:13824);
Y3=fft(A3);
Y3=Y3(1:10220,:);
Y=[Y1,Y2,Y3];
NorESM1_M_Ayy=abs(Y).^2;
NorESM1_M_Ayy=NorESM1_M_Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00002);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
NorESM1_M_ybar1=sum(NorESM1_M_Ayy(b:a,:));
YBAR1=reshape(NorESM1_M_ybar1,[144,96]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
NorESM1_M_ybar2=sum(NorESM1_M_Ayy(c:b,:));
YBAR2=reshape(NorESM1_M_ybar2,[144,96]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
NorESM1_M_ybar3=sum(NorESM1_M_Ayy(d:c,:));
YBAR3=reshape(NorESM1_M_ybar3,[144,96]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
NorESM1_M_hflsf1=NorESM1_M_ybar1./(NorESM1_M_ybar1+NorESM1_M_ybar2+NorESM1_M_ybar3);
NorESM1_M_hflsf2=NorESM1_M_ybar2./(NorESM1_M_ybar1+NorESM1_M_ybar2+NorESM1_M_ybar3);
NorESM1_M_hflsf3=NorESM1_M_ybar3./(NorESM1_M_ybar1+NorESM1_M_ybar2+NorESM1_M_ybar3);
NorESM1_M_hflsF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
NorESM1_M_hflsF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
NorESM1_M_hflsF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's ET_n to the same as SMAP
[p,~]=size(NorESM1_M_lonData);
[q,~]=size(NorESM1_M_latData);
m=NorESM1_M_lonData(p);
NorESM1_M_lon=linspace(-m/2,m/2,p);
NorESM1_M_lon=NorESM1_M_lon.';
BNU_lon1=NorESM1_M_lon(1:p/2,1);
BNU_lon2=NorESM1_M_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    NorESM1_M_lat_D(count:count+(p-1),1)=NorESM1_M_latData(i,1);
    NorESM1_M_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
NorESM1_M_hflsf1=NorESM1_M_hflsf1.';
NorESM1_M_hflsf2=NorESM1_M_hflsf2.';
NorESM1_M_hflsf3=NorESM1_M_hflsf3.';
NorESM1_M_alpha1=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_hflsf1];
NorESM1_M_alpha2=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_hflsf2];
NorESM1_M_alpha3=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_hflsf3];
[m1,~]=size(NorESM1_M_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
NorESM1_M_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha1(i,1)-SMAP1(j,1))+abs(NorESM1_M_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P1(j,1)=NorESM1_M_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(NorESM1_M_new_P1,[964,406]);
NorESM1_M_new_hflsf1=new_P1.';
[m1,~]=size(NorESM1_M_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
NorESM1_M_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha2(i,1)-SMAP2(j,1))+abs(NorESM1_M_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P2(j,1)=NorESM1_M_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(NorESM1_M_new_P2,[964,406]);
NorESM1_M_new_hflsf2=new_mf2.';
[m1,~]=size(NorESM1_M_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
NorESM1_M_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha3(i,1)-SMAP3(j,1))+abs(NorESM1_M_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P3(j,1)=NorESM1_M_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(NorESM1_M_new_P3,[964,406]);
NorESM1_M_new_hflsf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
NorESM1_M_hflsf1_new=NorESM1_M_new_hflsf1./temp;
NorESM1_M_hflsf2_new=NorESM1_M_new_hflsf2./temp;
NorESM1_M_hflsf3_new=NorESM1_M_new_hflsf3./temp;

% % P_n of NorESM1_M (pr)
% Get P_n over different time scales
N=20440;
Fs=1;
A(20440,13824)=0;
count=1;
for q=1:96
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
A1=dA(:,1:4608);
Y1=fft(A1);
Y1=Y1(1:10220,:);
A2=dA(:,4609:9216);
Y2=fft(A2);
Y2=Y2(1:10220,:);
A3=dA(:,9217:13824);
Y3=fft(A3);
Y3=Y3(1:10220,:);
Y=[Y1,Y2,Y3];
NorESM1_M_Ayy=abs(Y).^2;
NorESM1_M_Ayy=NorESM1_M_Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00002);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
NorESM1_M_ybar1=sum(NorESM1_M_Ayy(b:a,:));
YBAR1=reshape(NorESM1_M_ybar1,[144,96]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
NorESM1_M_ybar2=sum(NorESM1_M_Ayy(c:b,:));
YBAR2=reshape(NorESM1_M_ybar2,[144,96]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
NorESM1_M_ybar3=sum(NorESM1_M_Ayy(d:c,:));
YBAR3=reshape(NorESM1_M_ybar3,[144,96]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
NorESM1_M_prf1=NorESM1_M_ybar1./(NorESM1_M_ybar1+NorESM1_M_ybar2+NorESM1_M_ybar3);
NorESM1_M_prf2=NorESM1_M_ybar2./(NorESM1_M_ybar1+NorESM1_M_ybar2+NorESM1_M_ybar3);
NorESM1_M_prf3=NorESM1_M_ybar3./(NorESM1_M_ybar1+NorESM1_M_ybar2+NorESM1_M_ybar3);
NorESM1_M_prF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
NorESM1_M_prF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
NorESM1_M_prF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's P_n to the same as SMAP
[p,~]=size(NorESM1_M_lonData);
[q,~]=size(NorESM1_M_latData);
m=NorESM1_M_lonData(p);
NorESM1_M_lon=linspace(-m/2,m/2,p);
NorESM1_M_lon=NorESM1_M_lon.';
BNU_lon1=NorESM1_M_lon(1:p/2,1);
BNU_lon2=NorESM1_M_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    NorESM1_M_lat_D(count:count+(p-1),1)=NorESM1_M_latData(i,1);
    NorESM1_M_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
NorESM1_M_prf1=NorESM1_M_prf1.';
NorESM1_M_prf2=NorESM1_M_prf2.';
NorESM1_M_prf3=NorESM1_M_prf3.';
NorESM1_M_alpha1=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_prf1];
NorESM1_M_alpha2=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_prf2];
NorESM1_M_alpha3=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_prf3];
[m1,~]=size(NorESM1_M_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
NorESM1_M_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha1(i,1)-SMAP1(j,1))+abs(NorESM1_M_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P1(j,1)=NorESM1_M_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(NorESM1_M_new_P1,[964,406]);
NorESM1_M_new_prf1=new_P1.';
[m1,~]=size(NorESM1_M_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
NorESM1_M_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha2(i,1)-SMAP2(j,1))+abs(NorESM1_M_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P2(j,1)=NorESM1_M_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(NorESM1_M_new_P2,[964,406]);
NorESM1_M_new_prf2=new_mf2.';
[m1,~]=size(NorESM1_M_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
NorESM1_M_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha3(i,1)-SMAP3(j,1))+abs(NorESM1_M_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P3(j,1)=NorESM1_M_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(NorESM1_M_new_P3,[964,406]);
NorESM1_M_new_prf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
NorESM1_M_prf1_new=NorESM1_M_new_prf1./temp;
NorESM1_M_prf2_new=NorESM1_M_new_prf2./temp;
NorESM1_M_prf3_new=NorESM1_M_new_prf3./temp;

% Differences of SSM_n between models and observations over different time scales
NorESM1_M_mdif1=NorESM1_M_mF1-SMF1;
NorESM1_M_mdif2=NorESM1_M_mF2-SMF2;
NorESM1_M_mdif3=NorESM1_M_mF3-SMF3;
[m,n]=size(SMF1);
NorESM1_M_mdif1_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF1(:,i);
    temp2=NorESM1_M_mdif1(:,i);
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
    NorESM1_M_mdif1_try1(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF2);
NorESM1_M_mdif2_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF2(:,i);
    temp2=NorESM1_M_mdif2(:,i);
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
    NorESM1_M_mdif2_try1(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF3);
NorESM1_M_mdif3_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF3(:,i);
    temp2=NorESM1_M_mdif3(:,i);
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
    NorESM1_M_mdif3_try1(:,count)=idx1;
    count=count+1;
end

% % SSM_kw of NorESM1_M (mrsos)
% Get SSM_kw over different time scales
N=20440;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00002);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
[m1,n1]=size(NorESM1_M_mrsosA);
[m2,n2]=size(NorESM1_M_mrsosF1);
NorESM1_M_mrsosp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=NorESM1_M_mrsosA(:,i);
    temp2=isnan(NorESM1_M_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        NorESM1_M_mrsosp1(1,count)=nan;
    else
        A=NorESM1_M_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        NorESM1_M_mrsosp1(1,count)=p1(1);
    end
    count=count+1;
end
mrsosp11=reshape(NorESM1_M_mrsosp1,[n2,m2]);
p11=rot90(mrsosp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
NorESM1_M_mrsosP1=[P2 P1];
NorESM1_M_mrsosp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=NorESM1_M_mrsosA(:,i);
    temp2=isnan(NorESM1_M_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        NorESM1_M_mrsosp2(1,count)=nan;
    else
        A=NorESM1_M_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        NorESM1_M_mrsosp2(1,count)=p2(1);
    end
    count=count+1;
end
mrsosp22=reshape(NorESM1_M_mrsosp2,[n2,m2]);
p22=rot90(mrsosp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
NorESM1_M_mrsosP2=[P2 P1];
NorESM1_M_mrsosp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=NorESM1_M_mrsosA(:,i);
    temp2=isnan(NorESM1_M_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        NorESM1_M_mrsosp3(1,count)=nan;
    else
        A=NorESM1_M_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        NorESM1_M_mrsosp3(1,count)=p3(1);
    end
    count=count+1;
end
mrsosp33=reshape(NorESM1_M_mrsosp3,[n2,m2]);
p33=rot90(mrsosp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
NorESM1_M_mrsosP3=[P2 P1];
% Change the spatial resolution of the model's SSM_kw to the same as SMAP
[p,~]=size(NorESM1_M_lonData);
[q,~]=size(NorESM1_M_latData);
m=NorESM1_M_lonData(p);
NorESM1_M_lon=linspace(-m/2,m/2,p);
NorESM1_M_lon=NorESM1_M_lon.';
BNU_lon1=NorESM1_M_lon(1:p/2,1);
BNU_lon2=NorESM1_M_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    NorESM1_M_lat_D(count:count+(p-1),1)=NorESM1_M_latData(i,1);
    NorESM1_M_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
NorESM1_M_mrsosp1=NorESM1_M_mrsosp1.';
NorESM1_M_mrsosp2=NorESM1_M_mrsosp2.';
NorESM1_M_mrsosp3=NorESM1_M_mrsosp3.';
NorESM1_M_alpha1=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_mrsosp1];
NorESM1_M_alpha2=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_mrsosp2];
NorESM1_M_alpha3=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_mrsosp3];
[m1,~]=size(NorESM1_M_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
NorESM1_M_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha1(i,1)-SMAP1(j,1))+abs(NorESM1_M_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P1(j,1)=NorESM1_M_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(NorESM1_M_new_P1,[964,406]);
NorESM1_M_new_mrsosP1=new_P1.';
[m1,~]=size(NorESM1_M_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
NorESM1_M_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha2(i,1)-SMAP2(j,1))+abs(NorESM1_M_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P2(j,1)=NorESM1_M_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(NorESM1_M_new_P2,[964,406]);
NorESM1_M_new_mrsosP2=new_mf2.';
[m1,~]=size(NorESM1_M_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
NorESM1_M_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha3(i,1)-SMAP3(j,1))+abs(NorESM1_M_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P3(j,1)=NorESM1_M_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(NorESM1_M_new_P3,[964,406]);
NorESM1_M_new_mrsosP3=new_mf3.';
[m,n]=size(SMF1);
NorESM1_M_mrsosP1_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF1(:,i);
    new_temp2=NorESM1_M_new_mrsosP1(:,i);
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
    NorESM1_M_mrsosP1_new(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF2);
NorESM1_M_mrsosP2_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF2(:,i);
    new_temp2=NorESM1_M_new_mrsosP2(:,i);
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
    NorESM1_M_mrsosP2_new(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF3);
NorESM1_M_mrsosP3_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF3(:,i);
    new_temp2=NorESM1_M_new_mrsosP3(:,i);
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
    NorESM1_M_mrsosP3_new(:,count)=idx1;
    count=count+1;
end
% Differences of SSM_kw between models and observations over different time scales
NorESM1_M_SMPdif1=NorESM1_M_mrsosP1_new-SMP1;
NorESM1_M_SMPdif2=NorESM1_M_mrsosP2_new-SMP2;
NorESM1_M_SMPdif3=NorESM1_M_mrsosP3_new-SMP3;

% % ET_kw of NorESM1_M (hfls)
% Get ET_kw over different time scales
N=20440;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00002);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
[m1,n1]=size(NorESM1_M_hflsA);
[m2,n2]=size(NorESM1_M_mrsosF1);
NorESM1_M_hflsp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=NorESM1_M_hflsA(:,i);
    temp2=isnan(NorESM1_M_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        NorESM1_M_hflsp1(1,count)=nan;
    else
        A=NorESM1_M_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        NorESM1_M_hflsp1(1,count)=p1(1);
    end
    count=count+1;
end
hflsp11=reshape(NorESM1_M_hflsp1,[n2,m2]);
p11=rot90(hflsp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
NorESM1_M_hflsP1=[P2 P1];
NorESM1_M_hflsp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=NorESM1_M_hflsA(:,i);
    temp2=isnan(NorESM1_M_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        NorESM1_M_hflsp2(1,count)=nan;
    else
        A=NorESM1_M_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        NorESM1_M_hflsp2(1,count)=p2(1);
    end
    count=count+1;
end
hflsp22=reshape(NorESM1_M_hflsp2,[n2,m2]);
p22=rot90(hflsp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
NorESM1_M_hflsP2=[P2 P1];
NorESM1_M_hflsp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=NorESM1_M_hflsA(:,i);
    temp2=isnan(NorESM1_M_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        NorESM1_M_hflsp3(1,count)=nan;
    else
        A=NorESM1_M_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        NorESM1_M_hflsp3(1,count)=p3(1);
    end
    count=count+1;
end
hflsp33=reshape(NorESM1_M_hflsp3,[n2,m2]);
p33=rot90(hflsp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
NorESM1_M_hflsP3=[P2 P1];
% Change the spatial resolution of the model's ET_kw to the same as SMAP
[p,~]=size(NorESM1_M_lonData);
[q,~]=size(NorESM1_M_latData);
m=NorESM1_M_lonData(p);
NorESM1_M_lon=linspace(-m/2,m/2,p);
NorESM1_M_lon=NorESM1_M_lon.';
BNU_lon1=NorESM1_M_lon(1:p/2,1);
BNU_lon2=NorESM1_M_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    NorESM1_M_lat_D(count:count+(p-1),1)=NorESM1_M_latData(i,1);
    NorESM1_M_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
NorESM1_M_hflsp1=NorESM1_M_hflsp1.';
NorESM1_M_hflsp2=NorESM1_M_hflsp2.';
NorESM1_M_hflsp3=NorESM1_M_hflsp3.';
NorESM1_M_alpha1=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_hflsp1];
NorESM1_M_alpha2=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_hflsp2];
NorESM1_M_alpha3=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_hflsp3];
[m1,~]=size(NorESM1_M_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
NorESM1_M_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha1(i,1)-SMAP1(j,1))+abs(NorESM1_M_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P1(j,1)=NorESM1_M_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(NorESM1_M_new_P1,[964,406]);
NorESM1_M_new_hflsP1=new_P1.';
[m1,~]=size(NorESM1_M_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
NorESM1_M_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha2(i,1)-SMAP2(j,1))+abs(NorESM1_M_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P2(j,1)=NorESM1_M_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(NorESM1_M_new_P2,[964,406]);
NorESM1_M_new_hflsP2=new_mf2.';
[m1,~]=size(NorESM1_M_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
NorESM1_M_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha3(i,1)-SMAP3(j,1))+abs(NorESM1_M_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P3(j,1)=NorESM1_M_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(NorESM1_M_new_P3,[964,406]);
NorESM1_M_new_hflsP3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
NorESM1_M_hflsP1_new=NorESM1_M_new_hflsP1./temp;
NorESM1_M_hflsP2_new=NorESM1_M_new_hflsP2./temp;
NorESM1_M_hflsP3_new=NorESM1_M_new_hflsP3./temp;
% Differences of ET_kw between models and observations over different time scales
NorESM1_M_EPdif1=NorESM1_M_hflsP1_new-GLEAM_P1;
NorESM1_M_EPdif2=NorESM1_M_hflsP2_new-GLEAM_P2;
NorESM1_M_EPdif3=NorESM1_M_hflsP3_new-GLEAM_P3;

% % P_kw of NorESM1_M (pr)
% Get P_kw over different time scales
N=20440;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00002);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
[m1,n1]=size(NorESM1_M_prA);
[m2,n2]=size(NorESM1_M_mrsosF1);
NorESM1_M_prp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=NorESM1_M_prA(:,i);
    temp2=isnan(NorESM1_M_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        NorESM1_M_prp1(1,count)=nan;
    else
        A=NorESM1_M_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        NorESM1_M_prp1(1,count)=p1(1);
    end
    count=count+1;
end
prp11=reshape(NorESM1_M_prp1,[n2,m2]);
p11=rot90(prp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
NorESM1_M_prP1=[P2 P1];
NorESM1_M_prp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=NorESM1_M_prA(:,i);
    temp2=isnan(NorESM1_M_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        NorESM1_M_prp2(1,count)=nan;
    else
        A=NorESM1_M_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        NorESM1_M_prp2(1,count)=p2(1);
    end
    count=count+1;
end
prp22=reshape(NorESM1_M_prp2,[n2,m2]);
p22=rot90(prp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
NorESM1_M_prP2=[P2 P1];
NorESM1_M_prp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=NorESM1_M_prA(:,i);
    temp2=isnan(NorESM1_M_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        NorESM1_M_prp3(1,count)=nan;
    else
        A=NorESM1_M_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        NorESM1_M_prp3(1,count)=p3(1);
    end
    count=count+1;
end
prp33=reshape(NorESM1_M_prp3,[n2,m2]);
p33=rot90(prp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
NorESM1_M_prP3=[P2 P1];
[p,~]=size(NorESM1_M_lonData);
[q,~]=size(NorESM1_M_latData);
m=NorESM1_M_lonData(p);
NorESM1_M_lon=linspace(-m/2,m/2,p);
NorESM1_M_lon=NorESM1_M_lon.';
BNU_lon1=NorESM1_M_lon(1:p/2,1);
BNU_lon2=NorESM1_M_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    NorESM1_M_lat_D(count:count+(p-1),1)=NorESM1_M_latData(i,1);
    NorESM1_M_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
NorESM1_M_prp1=NorESM1_M_prp1.';
NorESM1_M_prp2=NorESM1_M_prp2.';
NorESM1_M_prp3=NorESM1_M_prp3.';
NorESM1_M_alpha1=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_prp1];
NorESM1_M_alpha2=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_prp2];
NorESM1_M_alpha3=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_prp3];
[m1,~]=size(NorESM1_M_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
NorESM1_M_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha1(i,1)-SMAP1(j,1))+abs(NorESM1_M_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P1(j,1)=NorESM1_M_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(NorESM1_M_new_P1,[964,406]);
NorESM1_M_new_prP1=new_P1.';
[m1,~]=size(NorESM1_M_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
NorESM1_M_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha2(i,1)-SMAP2(j,1))+abs(NorESM1_M_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P2(j,1)=NorESM1_M_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(NorESM1_M_new_P2,[964,406]);
NorESM1_M_new_prP2=new_mf2.';
[m1,~]=size(NorESM1_M_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
NorESM1_M_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_alpha3(i,1)-SMAP3(j,1))+abs(NorESM1_M_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_P3(j,1)=NorESM1_M_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(NorESM1_M_new_P3,[964,406]);
NorESM1_M_new_prP3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
NorESM1_M_prP1_new=NorESM1_M_new_prP1./temp;
NorESM1_M_prP2_new=NorESM1_M_new_prP2./temp;
NorESM1_M_prP3_new=NorESM1_M_new_prP3./temp;
% Differences of P_kw between models and observations over different time scales
NorESM1_M_PPdif1=NorESM1_M_prP1_new-ERA5_P1;
NorESM1_M_PPdif2=NorESM1_M_prP2_new-ERA5_P2;
NorESM1_M_PPdif3=NorESM1_M_prP3_new-ERA5_P3;