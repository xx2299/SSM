% % SSM_n of BNU_ESM (mrsos)
% Get SSM_n over different time scales
N=20440;
Fs=1;
A(20440,8192)=0;
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
Y=fft(dA);
BNU_ESM_Ayy=abs(Y).^2;
BNU_ESM_Ayy=BNU_ESM_Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00002);
b=find(abs(x-1/30)<=0.00002);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00002);
BNU_ESM_ybar1=sum(BNU_ESM_Ayy(b:a,:));
YBAR1=reshape(BNU_ESM_ybar1,[128,64]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
BNU_ESM_ybar2=sum(BNU_ESM_Ayy(c:b,:));
YBAR2=reshape(BNU_ESM_ybar2,[128,64]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
BNU_ESM_ybar3=sum(BNU_ESM_Ayy(d:c,:));
YBAR3=reshape(BNU_ESM_ybar3,[128,64]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
BNU_ESM_mrsosf1=BNU_ESM_ybar1./(BNU_ESM_ybar1+BNU_ESM_ybar2+BNU_ESM_ybar3);
BNU_ESM_mrsosf2=BNU_ESM_ybar2./(BNU_ESM_ybar1+BNU_ESM_ybar2+BNU_ESM_ybar3);
BNU_ESM_mrsosf3=BNU_ESM_ybar3./(BNU_ESM_ybar1+BNU_ESM_ybar2+BNU_ESM_ybar3);
BNU_ESM_mrsosF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
BNU_ESM_mrsosF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
BNU_ESM_mrsosF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(BNU_ESM_lonData);
[q,~]=size(BNU_ESM_latData);
m=BNU_ESM_lonData(p);
BNU_ESM_lon=linspace(-m/2,m/2,p);
BNU_ESM_lon=BNU_ESM_lon.';
BNU_lon1=BNU_ESM_lon(1:p/2,1);
BNU_lon2=BNU_ESM_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    BNU_ESM_lat_D(count:count+(p-1),1)=BNU_ESM_latData(i,1);
    BNU_ESM_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
BNU_ESM_mrsosf1=BNU_ESM_mrsosf1.';
BNU_ESM_mrsosf2=BNU_ESM_mrsosf2.';
BNU_ESM_mrsosf3=BNU_ESM_mrsosf3.';
BNU_ESM_alpha1=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_mrsosf1];
BNU_ESM_alpha2=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_mrsosf2];
BNU_ESM_alpha3=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_mrsosf3];
[m1,~]=size(BNU_ESM_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
BNU_ESM_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha1(i,1)-SMAP1(j,1))+abs(BNU_ESM_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P1(j,1)=BNU_ESM_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(BNU_ESM_new_P1,[964,406]);
BNU_ESM_new_mrsosf1=new_P1.';
[m1,~]=size(BNU_ESM_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
BNU_ESM_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha2(i,1)-SMAP2(j,1))+abs(BNU_ESM_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P2(j,1)=BNU_ESM_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(BNU_ESM_new_P2,[964,406]);
BNU_ESM_new_mrsosf2=new_mf2.';
[m1,~]=size(BNU_ESM_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BNU_ESM_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha3(i,1)-SMAP3(j,1))+abs(BNU_ESM_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P3(j,1)=BNU_ESM_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(BNU_ESM_new_P3,[964,406]);
BNU_ESM_new_mrsosf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
BNU_ESM_mrsosf1_new=BNU_ESM_new_mrsosf1./temp;
BNU_ESM_mrsosf2_new=BNU_ESM_new_mrsosf2./temp;
BNU_ESM_mrsosf3_new=BNU_ESM_new_mrsosf3./temp;

% % ET_n of BNU_ESM (hfls)
% Get ET_n over different time scales
N=20440;
Fs=1;
A(20440,8192)=0;
count=1;
for q=1:64
    for p=1:128
        C=hflsData(p,q,:);
        A(:,count)=C(:);
        count=count+1;
    end
end
[row,col] = size(A);
dA(row,col)=0;
for i=1:col
    dA(:,i)=detrend(A(:,i));
end
Y=fft(dA);
Ayy=abs(Y).^2;
Ayy=Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00002);
b=find(abs(x-1/30)<=0.00002);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00002);
BNU_ESM_ybar1=sum(Ayy(b:a,:));
YBAR1=reshape(BNU_ESM_ybar1,[128,64]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
BNU_ESM_ybar2=sum(Ayy(c:b,:));
YBAR2=reshape(BNU_ESM_ybar2,[128,64]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
BNU_ESM_ybar3=sum(Ayy(d:c,:));
YBAR3=reshape(BNU_ESM_ybar3,[128,64]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
BNU_ESM_hflsf1=BNU_ESM_ybar1./(BNU_ESM_ybar1+BNU_ESM_ybar2+BNU_ESM_ybar3);
BNU_ESM_hflsf2=BNU_ESM_ybar2./(BNU_ESM_ybar1+BNU_ESM_ybar2+BNU_ESM_ybar3);
BNU_ESM_hflsf3=BNU_ESM_ybar3./(BNU_ESM_ybar1+BNU_ESM_ybar2+BNU_ESM_ybar3);
BNU_ESM_hflsF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
BNU_ESM_hflsF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
BNU_ESM_hflsF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's ET_n to the same as SMAP
[p,~]=size(BNU_ESM_lonData);
[q,~]=size(BNU_ESM_latData);
m=BNU_ESM_lonData(p);
BNU_ESM_lon=linspace(-m/2,m/2,p);
BNU_ESM_lon=BNU_ESM_lon.';
BNU_lon1=BNU_ESM_lon(1:p/2,1);
BNU_lon2=BNU_ESM_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    BNU_ESM_lat_D(count:count+(p-1),1)=BNU_ESM_latData(i,1);
    BNU_ESM_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
BNU_ESM_hflsf1=BNU_ESM_hflsf1.';
BNU_ESM_hflsf2=BNU_ESM_hflsf2.';
BNU_ESM_hflsf3=BNU_ESM_hflsf3.';
BNU_ESM_alpha1=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_hflsf1];
BNU_ESM_alpha2=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_hflsf2];
BNU_ESM_alpha3=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_hflsf3];
[m1,~]=size(BNU_ESM_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
BNU_ESM_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha1(i,1)-SMAP1(j,1))+abs(BNU_ESM_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P1(j,1)=BNU_ESM_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(BNU_ESM_new_P1,[964,406]);
BNU_ESM_new_hflsf1=new_P1.';
[m1,~]=size(BNU_ESM_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
BNU_ESM_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha2(i,1)-SMAP2(j,1))+abs(BNU_ESM_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P2(j,1)=BNU_ESM_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(BNU_ESM_new_P2,[964,406]);
BNU_ESM_new_hflsf2=new_mf2.';
[m1,~]=size(BNU_ESM_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BNU_ESM_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha3(i,1)-SMAP3(j,1))+abs(BNU_ESM_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P3(j,1)=BNU_ESM_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(BNU_ESM_new_P3,[964,406]);
BNU_ESM_new_hflsf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
BNU_ESM_hflsf1_new=BNU_ESM_new_hflsf1./temp;
BNU_ESM_hflsf2_new=BNU_ESM_new_hflsf2./temp;
BNU_ESM_hflsf3_new=BNU_ESM_new_hflsf3./temp;

% % P_n of BCC-CSM1.1 (pr)
% Get P_n over different time scales
N=20440;
Fs=1;
A(20440,8192)=0;
count=1;
for q=1:64
    for p=1:128
        C=prData(p,q,:);
        A(:,count)=C(:);
        count=count+1;
    end
end
[row,col] = size(A);
dA(row,col)=0;
for i=1:col
    dA(:,i)=detrend(A(:,i));
end
A(A==0)=nan;
A=single(A);
AT=A.';
Z=mapminmax(AT,0,1);
A=Z.';
Y=fft(A);
Ayy=abs(Y).^2;
Ayy=Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00002);
b=find(abs(x-1/30)<=0.00002);
c=find(abs(x-1/90)<=0.00002);
d=find(abs(x-1/365)<=0.00001);
BNU_ESM_ybar1=sum(Ayy(b:a,:));
YBAR1=reshape(BNU_ESM_ybar1,[128,64]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
BNU_ESM_ybar2=sum(Ayy(c:b,:));
YBAR2=reshape(BNU_ESM_ybar2,[128,64]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
BNU_ESM_ybar3=sum(Ayy(d:c,:));
YBAR3=reshape(BNU_ESM_ybar3,[128,64]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
BNU_ESM_prf1=BNU_ESM_ybar1./(BNU_ESM_ybar1+BNU_ESM_ybar2+BNU_ESM_ybar3);
BNU_ESM_prf2=BNU_ESM_ybar2./(BNU_ESM_ybar1+BNU_ESM_ybar2+BNU_ESM_ybar3);
BNU_ESM_prf3=BNU_ESM_ybar3./(BNU_ESM_ybar1+BNU_ESM_ybar2+BNU_ESM_ybar3);
BNU_ESM_prF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
BNU_ESM_prF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
BNU_ESM_prF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's P_n to the same as SMAP
[p,~]=size(BNU_ESM_lonData);
[q,~]=size(BNU_ESM_latData);
m=BNU_ESM_lonData(p);
BNU_ESM_lon=linspace(-m/2,m/2,p);
BNU_ESM_lon=BNU_ESM_lon.';
BNU_lon1=BNU_ESM_lon(1:p/2,1);
BNU_lon2=BNU_ESM_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    BNU_ESM_lat_D(count:count+(p-1),1)=BNU_ESM_latData(i,1);
    BNU_ESM_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
BNU_ESM_prf1=BNU_ESM_prf1.';
BNU_ESM_prf2=BNU_ESM_prf2.';
BNU_ESM_prf3=BNU_ESM_prf3.';
BNU_ESM_alpha1=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_prf1];
BNU_ESM_alpha2=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_prf2];
BNU_ESM_alpha3=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_prf3];
[m1,~]=size(BNU_ESM_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
BNU_ESM_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha1(i,1)-SMAP1(j,1))+abs(BNU_ESM_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P1(j,1)=BNU_ESM_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(BNU_ESM_new_P1,[964,406]);
BNU_ESM_new_prf1=new_P1.';
[m1,~]=size(BNU_ESM_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
BNU_ESM_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha2(i,1)-SMAP2(j,1))+abs(BNU_ESM_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P2(j,1)=BNU_ESM_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(BNU_ESM_new_P2,[964,406]);
BNU_ESM_new_prf2=new_mf2.';
[m1,~]=size(BNU_ESM_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BNU_ESM_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha3(i,1)-SMAP3(j,1))+abs(BNU_ESM_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P3(j,1)=BNU_ESM_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(BNU_ESM_new_P3,[964,406]);
BNU_ESM_new_prf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
BNU_ESM_prf1_new=BNU_ESM_new_prf1./temp;
BNU_ESM_prf2_new=BNU_ESM_new_prf2./temp;
BNU_ESM_prf3_new=BNU_ESM_new_prf3./temp;

% Differences of SSM_n between models and observations over different time scales
BNU_ESM_mdif1=BNU_ESM_mF1-SMF1;
BNU_ESM_mdif2=BNU_ESM_mF2-SMF2;
BNU_ESM_mdif3=BNU_ESM_mF3-SMF3;
[m,n]=size(SMF1);
BNU_ESM_mdif1_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF1(:,i);
    temp2=BNU_ESM_mdif1(:,i);
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
    BNU_ESM_mdif1_try1(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF2);
BNU_ESM_mdif2_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF2(:,i);
    temp2=BNU_ESM_mdif2(:,i);
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
    BNU_ESM_mdif2_try1(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF3);
BNU_ESM_mdif3_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF3(:,i);
    temp2=BNU_ESM_mdif3(:,i);
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
    BNU_ESM_mdif3_try1(:,count)=idx1;
    count=count+1;
end

% % SSM_kw of BNU_ESM (mrsos)
% Get SSM_kw over different time scales
N=20440;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00002);
b=find(abs(x-1/30)<=0.00002);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00002);
[m1,n1]=size(BNU_ESM_mrsosA);
[m2,n2]=size(BNU_ESM_mrsosF1);
BNU_ESM_mrsosp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=BNU_ESM_mrsosA(:,i);
    temp2=isnan(BNU_ESM_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BNU_ESM_mrsosp1(1,count)=nan;
    else
        A=BNU_ESM_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        BNU_ESM_mrsosp1(1,count)=p1(1);
    end
    count=count+1;
end
mrsosp11=reshape(BNU_ESM_mrsosp1,[n2,m2]);
p11=rot90(mrsosp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
BNU_ESM_mrsosP1=[P2 P1];
BNU_ESM_mrsosp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=BNU_ESM_mrsosA(:,i);
    temp2=isnan(BNU_ESM_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BNU_ESM_mrsosp2(1,count)=nan;
    else
        A=BNU_ESM_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        BNU_ESM_mrsosp2(1,count)=p2(1);
    end
    count=count+1;
end
mrsosp22=reshape(BNU_ESM_mrsosp2,[n2,m2]);
p22=rot90(mrsosp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
BNU_ESM_mrsosP2=[P2 P1];
BNU_ESM_mrsosp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=BNU_ESM_mrsosA(:,i);
    temp2=isnan(BNU_ESM_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BNU_ESM_mrsosp3(1,count)=nan;
    else
        A=BNU_ESM_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        BNU_ESM_mrsosp3(1,count)=p3(1);
    end
    count=count+1;
end
mrsosp33=reshape(BNU_ESM_mrsosp3,[n2,m2]);
p33=rot90(mrsosp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
BNU_ESM_mrsosP3=[P2 P1];
% Change the spatial resolution of the model's SSM_kw to the same as SMAP
[p,~]=size(BNU_ESM_lonData);
[q,~]=size(BNU_ESM_latData);
m=BNU_ESM_lonData(p);
BNU_ESM_lon=linspace(-m/2,m/2,p);
BNU_ESM_lon=BNU_ESM_lon.';
BNU_lon1=BNU_ESM_lon(1:p/2,1);
BNU_lon2=BNU_ESM_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    BNU_ESM_lat_D(count:count+(p-1),1)=BNU_ESM_latData(i,1);
    BNU_ESM_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
BNU_ESM_mrsosp1=BNU_ESM_mrsosp1.';
BNU_ESM_mrsosp2=BNU_ESM_mrsosp2.';
BNU_ESM_mrsosp3=BNU_ESM_mrsosp3.';
BNU_ESM_alpha1=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_mrsosp1];
BNU_ESM_alpha2=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_mrsosp2];
BNU_ESM_alpha3=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_mrsosp3];
[m1,~]=size(BNU_ESM_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
BNU_ESM_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha1(i,1)-SMAP1(j,1))+abs(BNU_ESM_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P1(j,1)=BNU_ESM_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(BNU_ESM_new_P1,[964,406]);
BNU_ESM_new_mrsosP1=new_P1.';
[m1,~]=size(BNU_ESM_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
BNU_ESM_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha2(i,1)-SMAP2(j,1))+abs(BNU_ESM_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P2(j,1)=BNU_ESM_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(BNU_ESM_new_P2,[964,406]);
BNU_ESM_new_mrsosP2=new_mf2.';
[m1,~]=size(BNU_ESM_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BCC_CSM_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha3(i,1)-SMAP3(j,1))+abs(BNU_ESM_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P3(j,1)=BNU_ESM_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(BCC_CSM_new_P3,[964,406]);
BNU_ESM_new_mrsosP3=new_mf3.';
[m,n]=size(SMF1);
BNU_ESM_mrsosP1_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF1(:,i);
    new_temp2=BNU_ESM_new_mrsosP1(:,i);
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
    BNU_ESM_mrsosP1_new(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF2);
BNU_ESM_mrsosP2_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF2(:,i);
    new_temp2=BNU_ESM_new_mrsosP2(:,i);
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
    BNU_ESM_mrsosP2_new(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF3);
BNU_ESM_mrsosP3_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF3(:,i);
    new_temp2=BNU_ESM_new_mrsosP3(:,i);
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
    BNU_ESM_mrsosP3_new(:,count)=idx1;
    count=count+1;
end
% Differences of SSM_kw between models and observations over different time scales
BNU_ESM_SMPdif1=BNU_ESM_mrsosP1_new-SMP1;
BNU_ESM_SMPdif2=BNU_ESM_mrsosP2_new-SMP2;
BNU_ESM_SMPdif3=BNU_ESM_mrsosP3_new-SMP3;

% % ET_kw of BNU_ESM (hfls)
% Get ET_kw over different time scales
N=20440;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00002);
b=find(abs(x-1/30)<=0.00002);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00002);
[m1,n1]=size(BNU_ESM_hflsA);
[m2,n2]=size(BNU_ESM_mrsosF1);
BNU_ESM_hflsp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=BNU_ESM_hflsA(:,i);
    temp2=isnan(BNU_ESM_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BNU_ESM_hflsp1(1,count)=nan;
    else
        A=BNU_ESM_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        BNU_ESM_hflsp1(1,count)=p1(1);
    end
    count=count+1;
end
hflsp11=reshape(BNU_ESM_hflsp1,[n2,m2]);
p11=rot90(hflsp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
BNU_ESM_hflsP1=[P2 P1];
BNU_ESM_hflsp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=BNU_ESM_hflsA(:,i);
    temp2=isnan(BNU_ESM_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BNU_ESM_hflsp2(1,count)=nan;
    else
        A=BNU_ESM_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        BNU_ESM_hflsp2(1,count)=p2(1);
    end
    count=count+1;
end
hflsp22=reshape(BNU_ESM_hflsp2,[n2,m2]);
p22=rot90(hflsp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
BNU_ESM_hflsP2=[P2 P1];
BNU_ESM_hflsp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=BNU_ESM_hflsA(:,i);
    temp2=isnan(BNU_ESM_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BNU_ESM_hflsp3(1,count)=nan;
    else
        A=BNU_ESM_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        BNU_ESM_hflsp3(1,count)=p3(1);
    end
    count=count+1;
end
hflsp33=reshape(BNU_ESM_hflsp3,[n2,m2]);
p33=rot90(hflsp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
BNU_ESM_hflsP3=[P2 P1];
% Change the spatial resolution of the model's ET_kw to the same as SMAP
[p,~]=size(BNU_ESM_lonData);
[q,~]=size(BNU_ESM_latData);
m=BNU_ESM_lonData(p);
BNU_ESM_lon=linspace(-m/2,m/2,p);
BNU_ESM_lon=BNU_ESM_lon.';
BNU_lon1=BNU_ESM_lon(1:p/2,1);
BNU_lon2=BNU_ESM_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    BNU_ESM_lat_D(count:count+(p-1),1)=BNU_ESM_latData(i,1);
    BNU_ESM_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
BNU_ESM_hflsp1=BNU_ESM_hflsp1.';
BNU_ESM_hflsp2=BNU_ESM_hflsp2.';
BNU_ESM_hflsp3=BNU_ESM_hflsp3.';
BNU_ESM_alpha1=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_hflsp1];
BNU_ESM_alpha2=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_hflsp2];
BNU_ESM_alpha3=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_hflsp3];
[m1,~]=size(BNU_ESM_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
BNU_ESM_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha1(i,1)-SMAP1(j,1))+abs(BNU_ESM_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P1(j,1)=BNU_ESM_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(BNU_ESM_new_P1,[964,406]);
BNU_ESM_new_hflsP1=new_P1.';
[m1,~]=size(BNU_ESM_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
BNU_ESM_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha2(i,1)-SMAP2(j,1))+abs(BNU_ESM_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P2(j,1)=BNU_ESM_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(BNU_ESM_new_P2,[964,406]);
BNU_ESM_new_hflsP2=new_mf2.';
[m1,~]=size(BNU_ESM_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BNU_ESM_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha3(i,1)-SMAP3(j,1))+abs(BNU_ESM_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P3(j,1)=BNU_ESM_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(BNU_ESM_new_P3,[964,406]);
BNU_ESM_new_hflsP3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
BNU_ESM_hflsP1_new=BNU_ESM_new_hflsP1./temp;
BNU_ESM_hflsP2_new=BNU_ESM_new_hflsP2./temp;
BNU_ESM_hflsP3_new=BNU_ESM_new_hflsP3./temp;
% Differences of ET_kw between models and observations over different time scales
BNU_ESM_EPdif1=BNU_ESM_hflsP1_new-GLEAM_P1;
BNU_ESM_EPdif2=BNU_ESM_hflsP2_new-GLEAM_P2;
BNU_ESM_EPdif3=BNU_ESM_hflsP3_new-GLEAM_P3;

% % P_kw of BNU_ESM (pr)
% Get P_kw over different time scales
N=20440;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00002);
b=find(abs(x-1/30)<=0.00002);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00002);
[m1,n1]=size(BNU_ESM_prA);
[m2,n2]=size(BNU_ESM_mrsosF1);
BNU_ESM_prp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=BNU_ESM_prA(:,i);
    temp2=isnan(BNU_ESM_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BNU_ESM_prp1(1,count)=nan;
    else
        A=BNU_ESM_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        BNU_ESM_prp1(1,count)=p1(1);
    end
    count=count+1;
end
prp11=reshape(BNU_ESM_prp1,[n2,m2]);
p11=rot90(prp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
BNU_ESM_prP1=[P2 P1];
BNU_ESM_prp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=BNU_ESM_prA(:,i);
    temp2=isnan(BNU_ESM_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BNU_ESM_prp2(1,count)=nan;
    else
        A=BNU_ESM_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        BNU_ESM_prp2(1,count)=p2(1);
    end
    count=count+1;
end
prp22=reshape(BNU_ESM_prp2,[n2,m2]);
p22=rot90(prp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
BNU_ESM_prP2=[P2 P1];
BNU_ESM_prp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=BNU_ESM_prA(:,i);
    temp2=isnan(BNU_ESM_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BNU_ESM_prp3(1,count)=nan;
    else
        A=BNU_ESM_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        BNU_ESM_prp3(1,count)=p3(1);
    end
    count=count+1;
end
prp33=reshape(BNU_ESM_prp3,[n2,m2]);
p33=rot90(prp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
BNU_ESM_prP3=[P2 P1];
% Change the spatial resolution of the model's P_kw to the same as SMAP
[p,~]=size(BNU_ESM_lonData);
[q,~]=size(BNU_ESM_latData);
m=BNU_ESM_lonData(p);
BNU_ESM_lon=linspace(-m/2,m/2,p);
BNU_ESM_lon=BNU_ESM_lon.';
BNU_lon1=BNU_ESM_lon(1:p/2,1);
BNU_lon2=BNU_ESM_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    BNU_ESM_lat_D(count:count+(p-1),1)=BNU_ESM_latData(i,1);
    BNU_ESM_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
BNU_ESM_prp1=BNU_ESM_prp1.';
BNU_ESM_prp2=BNU_ESM_prp2.';
BNU_ESM_prp3=BNU_ESM_prp3.';
BNU_ESM_alpha1=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_prp1];
BNU_ESM_alpha2=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_prp2];
BNU_ESM_alpha3=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_prp3];
[m1,~]=size(BNU_ESM_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
BNU_ESM_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha1(i,1)-SMAP1(j,1))+abs(BNU_ESM_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P1(j,1)=BNU_ESM_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(BNU_ESM_new_P1,[964,406]);
BNU_ESM_new_prP1=new_P1.';
[m1,~]=size(BNU_ESM_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
BNU_ESM_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha2(i,1)-SMAP2(j,1))+abs(BNU_ESM_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P2(j,1)=BNU_ESM_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(BNU_ESM_new_P2,[964,406]);
BNU_ESM_new_prP2=new_mf2.';
[m1,~]=size(BNU_ESM_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BNU_ESM_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_alpha3(i,1)-SMAP3(j,1))+abs(BNU_ESM_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_P3(j,1)=BNU_ESM_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(BNU_ESM_new_P3,[964,406]);
BNU_ESM_new_prP3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
BNU_ESM_prP1_new=BNU_ESM_new_prP1./temp;
BNU_ESM_prP2_new=BNU_ESM_new_prP2./temp;
BNU_ESM_prP3_new=BNU_ESM_new_prP3./temp;
% Differences of P_kw between models and observations over different time scales
BNU_ESM_PPdif1=BNU_ESM_prP1_new-ERA5_P1;
BNU_ESM_PPdif2=BNU_ESM_prP2_new-ERA5_P2;
BNU_ESM_PPdif3=BNU_ESM_prP3_new-ERA5_P3;