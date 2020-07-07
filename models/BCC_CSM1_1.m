% % SSM_n of BCC-CSM1.1 (mrsos)
% Get SSM_n over different time scales
N=22995;
Fs=1;
A(22995,8192)=0;
count=1;
for q=1:64
    for p=1:128
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
Y=fft(dA);
BCC_CSM_Ayy=abs(Y).^2;
BCC_CSM_Ayy=BCC_CSM_Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00002);
b=767;
c=256;
d=find(abs(x-1/365)<=0.00002);
BCC_CSM_ybar1=sum(BCC_CSM_Ayy(b:a,:));
YBAR1=reshape(BCC_CSM_ybar1,[128,64]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
BCC_CSM_ybar2=sum(BCC_CSM_Ayy(c:b,:));
YBAR2=reshape(BCC_CSM_ybar2,[128,64]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
BCC_CSM_ybar3=sum(BCC_CSM_Ayy(d:c,:));
YBAR3=reshape(BCC_CSM_ybar3,[128,64]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
BCC_CSM_mrsosf1=BCC_CSM_ybar1./(BCC_CSM_ybar1+BCC_CSM_ybar2+BCC_CSM_ybar3);
BCC_CSM_mrsosf2=BCC_CSM_ybar2./(BCC_CSM_ybar1+BCC_CSM_ybar2+BCC_CSM_ybar3);
BCC_CSM_mrsosf3=BCC_CSM_ybar3./(BCC_CSM_ybar1+BCC_CSM_ybar2+BCC_CSM_ybar3);
BCC_CSM_mrsosF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
BCC_CSM_mrsosF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
BCC_CSM_mrsosF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(BCC_CSM_lonData);
[q,~]=size(BCC_CSM_latData);
m=BCC_CSM_lonData(p);
BCC_CSM_lon=linspace(-m/2,m/2,p);
BCC_CSM_lon=BCC_CSM_lon.';
BNU_lon1=BCC_CSM_lon(1:p/2,1);
BNU_lon2=BCC_CSM_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    BCC_CSM_lat_D(count:count+(p-1),1)=BCC_CSM_latData(i,1);
    BCC_CSM_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
BCC_CSM_mrsosf1=BCC_CSM_mrsosf1.';
BCC_CSM_mrsosf2=BCC_CSM_mrsosf2.';
BCC_CSM_mrsosf3=BCC_CSM_mrsosf3.';
BCC_CSM_alpha1=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_mrsosf1];
BCC_CSM_alpha2=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_mrsosf2];
BCC_CSM_alpha3=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_mrsosf3];
[m1,~]=size(BCC_CSM_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
BCC_CSM_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha1(i,1)-SMAP1(j,1))+abs(BCC_CSM_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P1(j,1)=BCC_CSM_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(BCC_CSM_new_P1,[964,406]);
BCC_CSM_new_mrsosf1=new_P1.';
[m1,~]=size(BCC_CSM_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
BCC_CSM_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha2(i,1)-SMAP2(j,1))+abs(BCC_CSM_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P2(j,1)=BCC_CSM_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(BCC_CSM_new_P2,[964,406]);
BCC_CSM_new_mrsosf2=new_mf2.';
[m1,~]=size(BCC_CSM_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BCC_CSM_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha3(i,1)-SMAP3(j,1))+abs(BCC_CSM_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P3(j,1)=BCC_CSM_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(BCC_CSM_new_P3,[964,406]);
BCC_CSM_new_mrsosf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
BCC_CSM_mrsosf1_new=BCC_CSM_new_mrsosf1./temp;
BCC_CSM_mrsosf2_new=BCC_CSM_new_mrsosf2./temp;
BCC_CSM_mrsosf3_new=BCC_CSM_new_mrsosf3./temp;

% % ET_n of BCC-CSM1.1 (hfls)
% Get ET_n over different time scales
N=22994;
Fs=1;
A(22994,8192)=0;
count=1;
for q=1:64
    for p=1:128
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
Y=fft(dA);
Ayy=abs(Y).^2;
Ayy=Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00002);
b=find(abs(x-1/30)<=0.000022);
c=find(abs(x-1/90)<=0.000022);
d=find(abs(x-1/365)<=0.00002);
BCC_CSM_ybar1=sum(Ayy(b:a,:));
YBAR1=reshape(BCC_CSM_ybar1,[128,64]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
BCC_CSM_ybar2=sum(Ayy(c:b,:));
YBAR2=reshape(BCC_CSM_ybar2,[128,64]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
BCC_CSM_ybar3=sum(Ayy(d:c,:));
YBAR3=reshape(BCC_CSM_ybar3,[128,64]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
BCC_CSM_hflsf1=BCC_CSM_ybar1./(BCC_CSM_ybar1+BCC_CSM_ybar2+BCC_CSM_ybar3);
BCC_CSM_hflsf2=BCC_CSM_ybar2./(BCC_CSM_ybar1+BCC_CSM_ybar2+BCC_CSM_ybar3);
BCC_CSM_hflsf3=BCC_CSM_ybar3./(BCC_CSM_ybar1+BCC_CSM_ybar2+BCC_CSM_ybar3);
BCC_CSM_hflsF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
BCC_CSM_hflsF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
BCC_CSM_hflsF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's ET_n to the same as SMAP
[p,~]=size(BCC_CSM_lonData);
[q,~]=size(BCC_CSM_latData);
m=BCC_CSM_lonData(p);
BCC_CSM_lon=linspace(-m/2,m/2,p);
BCC_CSM_lon=BCC_CSM_lon.';
BNU_lon1=BCC_CSM_lon(1:p/2,1);
BNU_lon2=BCC_CSM_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    BCC_CSM_lat_D(count:count+(p-1),1)=BCC_CSM_latData(i,1);
    BCC_CSM_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
BCC_CSM_hflsf1=BCC_CSM_hflsf1.';
BCC_CSM_hflsf2=BCC_CSM_hflsf2.';
BCC_CSM_hflsf3=BCC_CSM_hflsf3.';
BCC_CSM_alpha1=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_hflsf1];
BCC_CSM_alpha2=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_hflsf2];
BCC_CSM_alpha3=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_hflsf3];
[m1,~]=size(BCC_CSM_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
BCC_CSM_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha1(i,1)-SMAP1(j,1))+abs(BCC_CSM_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P1(j,1)=BCC_CSM_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(BCC_CSM_new_P1,[964,406]);
BCC_CSM_new_hflsf1=new_P1.';
[m1,~]=size(BCC_CSM_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
BCC_CSM_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha2(i,1)-SMAP2(j,1))+abs(BCC_CSM_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P2(j,1)=BCC_CSM_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(BCC_CSM_new_P2,[964,406]);
BCC_CSM_new_hflsf2=new_mf2.';
[m1,~]=size(BCC_CSM_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BCC_CSM_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha3(i,1)-SMAP3(j,1))+abs(BCC_CSM_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P3(j,1)=BCC_CSM_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(BCC_CSM_new_P3,[964,406]);
BCC_CSM_new_hflsf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
BCC_CSM_hflsf1_new=BCC_CSM_new_hflsf1./temp;
BCC_CSM_hflsf2_new=BCC_CSM_new_hflsf2./temp;
BCC_CSM_hflsf3_new=BCC_CSM_new_hflsf3./temp;

% % P_n of BCC-CSM1.1 (pr)
% Get P_n over different time scales
N=59494;
Fs=1;
A(59494,8192)=0;
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
Y=fft(dA);
Ayy=abs(Y).^2;
Ayy=Ayy/(N/2);
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
BCC_CSM_ybar1=sum(Ayy(b:a,:));
YBAR1=reshape(BCC_CSM_ybar1,[128,64]);
Fybar1=YBAR1.';
FYBAR1=flipud(Fybar1);
BCC_CSM_ybar2=sum(Ayy(c:b,:));
YBAR2=reshape(BCC_CSM_ybar2,[128,64]);
Fybar2=YBAR2.';
FYBAR2=flipud(Fybar2);
BCC_CSM_ybar3=sum(Ayy(d:c,:));
YBAR3=reshape(BCC_CSM_ybar3,[128,64]);
Fybar3=YBAR3.';
FYBAR3=flipud(Fybar3);
BCC_CSM_prf1=BCC_CSM_ybar1./(BCC_CSM_ybar1+BCC_CSM_ybar2+BCC_CSM_ybar3);
BCC_CSM_prf2=BCC_CSM_ybar2./(BCC_CSM_ybar1+BCC_CSM_ybar2+BCC_CSM_ybar3);
BCC_CSM_prf3=BCC_CSM_ybar3./(BCC_CSM_ybar1+BCC_CSM_ybar2+BCC_CSM_ybar3);
BCC_CSM_prF1=FYBAR1./(FYBAR1+FYBAR2+FYBAR3);
BCC_CSM_prF2=FYBAR2./(FYBAR1+FYBAR2+FYBAR3);
BCC_CSM_prF3=FYBAR3./(FYBAR1+FYBAR2+FYBAR3);
% Change the spatial resolution of the model's P_n to the same as SMAP
[p,~]=size(BCC_CSM_lonData);
[q,~]=size(BCC_CSM_latData);
m=BCC_CSM_lonData(p);
BCC_CSM_lon=linspace(-m/2,m/2,p);
BCC_CSM_lon=BCC_CSM_lon.';
BNU_lon1=BCC_CSM_lon(1:p/2,1);
BNU_lon2=BCC_CSM_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    BCC_CSM_lat_D(count:count+(p-1),1)=BCC_CSM_latData(i,1);
    BCC_CSM_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
BCC_CSM_prf1=BCC_CSM_prf1.';
BCC_CSM_prf2=BCC_CSM_prf2.';
BCC_CSM_prf3=BCC_CSM_prf3.';
BCC_CSM_alpha1=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_prf1];
BCC_CSM_alpha2=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_prf2];
BCC_CSM_alpha3=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_prf3];
[m1,~]=size(BCC_CSM_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
BCC_CSM_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha1(i,1)-SMAP1(j,1))+abs(BCC_CSM_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P1(j,1)=BCC_CSM_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(BCC_CSM_new_P1,[964,406]);
BCC_CSM_new_prf1=new_P1.';
[m1,~]=size(BCC_CSM_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
BCC_CSM_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha2(i,1)-SMAP2(j,1))+abs(BCC_CSM_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P2(j,1)=BCC_CSM_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(BCC_CSM_new_P2,[964,406]);
BCC_CSM_new_prf2=new_mf2.';
[m1,~]=size(BCC_CSM_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BCC_CSM_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha3(i,1)-SMAP3(j,1))+abs(BCC_CSM_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P3(j,1)=BCC_CSM_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(BCC_CSM_new_P3,[964,406]);
BCC_CSM_new_prf3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
BCC_CSM_prf1_new=BCC_CSM_new_prf1./temp;
BCC_CSM_prf2_new=BCC_CSM_new_prf2./temp;
BCC_CSM_prf3_new=BCC_CSM_new_prf3./temp;

% Differences of SSM_n between models and observations over different time scales
BCC_CSM_mdif1=BCC_CSM_mF1-SMF1;
BCC_CSM_mdif2=BCC_CSM_mF2-SMF2;
BCC_CSM_mdif3=BCC_CSM_mF3-SMF3;
[m,n]=size(SMF1);
BCC_CSM_mdif1_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF1(:,i);
    temp2=BCC_CSM_mdif1(:,i);
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
    BCC_CSM_mdif1_try1(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF2);
BCC_CSM_mdif2_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF2(:,i);
    temp2=BCC_CSM_mdif2(:,i);
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
    BCC_CSM_mdif2_try1(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF3);
BCC_CSM_mdif3_try1(m,n)=0;
count=1;
for i=1:n
    temp1=SMF3(:,i);
    temp2=BCC_CSM_mdif3(:,i);
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
    BCC_CSM_mdif3_try1(:,count)=idx1;
    count=count+1;
end

% % SSM_kw of BCC-CSM1.1 (mrsos)
% Get SSM_kw over different time scales
N=22995;
Fs=1;
F=((1:N)-1)*Fs/N;
x=F(1:(N+1)/2);
a=find(abs(x-1/7)<=0.00002);
b=767;
c=256;
d=find(abs(x-1/365)<=0.00002);
[m1,n1]=size(BCC_CSM_mrsosA);
[m2,n2]=size(BCC_CSM_mrsosF1);
BCC_CSM_mrsosp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=BCC_CSM_mrsosA(:,i);
    temp2=isnan(BCC_CSM_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BCC_CSM_mrsosp1(1,count)=nan;
    else
        A=BCC_CSM_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m1+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        BCC_CSM_mrsosp1(1,count)=p1(1);
    end
    count=count+1;
end
mrsosp11=reshape(BCC_CSM_mrsosp1,[n2,m2]);
p11=rot90(mrsosp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
BCC_CSM_mrsosP1=[P2 P1];
BCC_CSM_mrsosp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=BCC_CSM_mrsosA(:,i);
    temp2=isnan(BCC_CSM_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BCC_CSM_mrsosp2(1,count)=nan;
    else
        A=BCC_CSM_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m1+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        BCC_CSM_mrsosp2(1,count)=p2(1);
    end
    count=count+1;
end
mrsosp22=reshape(BCC_CSM_mrsosp2,[n2,m2]);
p22=rot90(mrsosp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
BCC_CSM_mrsosP2=[P2 P1];
BCC_CSM_mrsosp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=BCC_CSM_mrsosA(:,i);
    temp2=isnan(BCC_CSM_mrsosA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BCC_CSM_mrsosp3(1,count)=nan;
    else
        A=BCC_CSM_mrsosA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m1+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        BCC_CSM_mrsosp3(1,count)=p3(1);
    end
    count=count+1;
end
mrsosp33=reshape(BCC_CSM_mrsosp3,[n2,m2]);
p33=rot90(mrsosp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
BCC_CSM_mrsosP3=[P2 P1];
% Change the spatial resolution of the model's SSM_kw to the same as SMAP
[p,~]=size(BCC_CSM_lonData);
[q,~]=size(BCC_CSM_latData);
m=BCC_CSM_lonData(p);
BCC_CSM_lon=linspace(-m/2,m/2,p);
BCC_CSM_lon=BCC_CSM_lon.';
BNU_lon1=BCC_CSM_lon(1:p/2,1);
BNU_lon2=BCC_CSM_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    BCC_CSM_lat_D(count:count+(p-1),1)=BCC_CSM_latData(i,1);
    BCC_CSM_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
BCC_CSM_mrsosp1=BCC_CSM_mrsosp1.';
BCC_CSM_mrsosp2=BCC_CSM_mrsosp2.';
BCC_CSM_mrsosp3=BCC_CSM_mrsosp3.';
BCC_CSM_alpha1=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_mrsosp1];
BCC_CSM_alpha2=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_mrsosp2];
BCC_CSM_alpha3=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_mrsosp3];
[m1,~]=size(BCC_CSM_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
BCC_CSM_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha1(i,1)-SMAP1(j,1))+abs(BCC_CSM_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P1(j,1)=BCC_CSM_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(BCC_CSM_new_P1,[964,406]);
BCC_CSM_new_mrsosP1=new_P1.';
[m1,~]=size(BCC_CSM_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
BCC_CSM_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha2(i,1)-SMAP2(j,1))+abs(BCC_CSM_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P2(j,1)=BCC_CSM_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(BCC_CSM_new_P2,[964,406]);
BCC_CSM_new_mrsosP2=new_mf2.';
[m1,~]=size(BCC_CSM_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BCC_CSM_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha3(i,1)-SMAP3(j,1))+abs(BCC_CSM_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P3(j,1)=BCC_CSM_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(BCC_CSM_new_P3,[964,406]);
BCC_CSM_new_mrsosP3=new_mf3.';
[m,n]=size(SMF1);
BCC_CSM_mrsosP1_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF1(:,i);
    new_temp2=BCC_CSM_new_mrsosP1(:,i);
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
    BCC_CSM_mrsosP1_new(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF2);
BCC_CSM_mrsosP2_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF2(:,i);
    new_temp2=BCC_CSM_new_mrsosP2(:,i);
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
    BCC_CSM_mrsosP2_new(:,count)=idx1;
    count=count+1;
end
[m,n]=size(SMF3);
BCC_CSM_mrsosP3_new(m,n)=0;
count=1;
for i=1:n
    new_temp1=SMF3(:,i);
    new_temp2=BCC_CSM_new_mrsosP3(:,i);
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
    BCC_CSM_mrsosP3_new(:,count)=idx1;
    count=count+1;
end
% Differences of SSM_kw between models and observations over different time scales
BCC_CSM_SMPdif1=BCC_CSM_mrsosP1_new-SMP1;
BCC_CSM_SMPdif2=BCC_CSM_mrsosP2_new-SMP2;
BCC_CSM_SMPdif3=BCC_CSM_mrsosP3_new-SMP3;

% % ET_kw of BCC-CSM1.1 (hfls)
% Get ET_kw over different time scales
N=22994;
Fs=1;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00002);
b=find(abs(x-1/30)<=0.000022);
c=find(abs(x-1/90)<=0.000022);
d=find(abs(x-1/365)<=0.00002);
[m1,n1]=size(BCC_CSM_hflsA);
[m2,n2]=size(BCC_CSM_mrsosF1);
BCC_CSM_hflsp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=BCC_CSM_hflsA(:,i);
    temp2=isnan(BCC_CSM_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BCC_CSM_hflsp1(1,count)=nan;
    else
        A=BCC_CSM_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        BCC_CSM_hflsp1(1,count)=p1(1);
    end
    count=count+1;
end
hflsp11=reshape(BCC_CSM_hflsp1,[n2,m2]);
p11=rot90(hflsp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
BCC_CSM_hflsP1=[P2 P1];
BCC_CSM_hflsp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=BCC_CSM_hflsA(:,i);
    temp2=isnan(BCC_CSM_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BCC_CSM_hflsp2(1,count)=nan;
    else
        A=BCC_CSM_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        BCC_CSM_hflsp2(1,count)=p2(1);
    end
    count=count+1;
end
hflsp22=reshape(BCC_CSM_hflsp2,[n2,m2]);
p22=rot90(hflsp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
BCC_CSM_hflsP2=[P2 P1];
BCC_CSM_hflsp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=BCC_CSM_hflsA(:,i);
    temp2=isnan(BCC_CSM_hflsA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BCC_CSM_hflsp3(1,count)=nan;
    else
        A=BCC_CSM_hflsA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        BCC_CSM_hflsp3(1,count)=p3(1);
    end
    count=count+1;
end
hflsp33=reshape(BCC_CSM_hflsp3,[n2,m2]);
p33=rot90(hflsp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
BCC_CSM_hflsP3=[P2 P1];
% Change the spatial resolution of the model's ET_kw to the same as SMAP
[p,~]=size(BCC_CSM_lonData);
[q,~]=size(BCC_CSM_latData);
m=BCC_CSM_lonData(p);
BCC_CSM_lon=linspace(-m/2,m/2,p);
BCC_CSM_lon=BCC_CSM_lon.';
BNU_lon1=BCC_CSM_lon(1:p/2,1);
BNU_lon2=BCC_CSM_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    BCC_CSM_lat_D(count:count+(p-1),1)=BCC_CSM_latData(i,1);
    BCC_CSM_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
BCC_CSM_hflsp1=BCC_CSM_hflsp1.';
BCC_CSM_hflsp2=BCC_CSM_hflsp2.';
BCC_CSM_hflsp3=BCC_CSM_hflsp3.';
BCC_CSM_alpha1=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_hflsp1];
BCC_CSM_alpha2=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_hflsp2];
BCC_CSM_alpha3=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_hflsp3];
[m1,~]=size(BCC_CSM_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
BCC_CSM_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha1(i,1)-SMAP1(j,1))+abs(BCC_CSM_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P1(j,1)=BCC_CSM_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(BCC_CSM_new_P1,[964,406]);
BCC_CSM_new_hflsP1=new_P1.';
[m1,~]=size(BCC_CSM_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
BCC_CSM_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha2(i,1)-SMAP2(j,1))+abs(BCC_CSM_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P2(j,1)=BCC_CSM_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(BCC_CSM_new_P2,[964,406]);
BCC_CSM_new_hflsP2=new_mf2.';
[m1,~]=size(BCC_CSM_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BCC_CSM_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha3(i,1)-SMAP3(j,1))+abs(BCC_CSM_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P3(j,1)=BCC_CSM_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(BCC_CSM_new_P3,[964,406]);
BCC_CSM_new_hflsP3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
BCC_CSM_hflsP1_new=BCC_CSM_new_hflsP1./temp;
BCC_CSM_hflsP2_new=BCC_CSM_new_hflsP2./temp;
BCC_CSM_hflsP3_new=BCC_CSM_new_hflsP3./temp;
% Differences of ET_kw between models and observations over different time scales
BCC_CSM_EPdif1=BCC_CSM_hflsP1_new-GLEAM_P1;
BCC_CSM_EPdif2=BCC_CSM_hflsP2_new-GLEAM_P2;
BCC_CSM_EPdif3=BCC_CSM_hflsP3_new-GLEAM_P3;

% % P_kw of BCC-CSM1.1 (pr)
% Get P_kw over different time scales
N=59494;
Fs=1;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.00001);
b=find(abs(x-1/30)<=0.00001);
c=find(abs(x-1/90)<=0.00001);
d=find(abs(x-1/365)<=0.00001);
[m1,n1]=size(BCC_CSM_prA);
[m2,n2]=size(BCC_CSM_mrsosF1);
BCC_CSM_prp1(1,n1)=0;
count=1;
for i=1:n1
    temp1=BCC_CSM_prA(:,i);
    temp2=isnan(BCC_CSM_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BCC_CSM_prp1(1,count)=nan;
    else
        A=BCC_CSM_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        BCC_CSM_prp1(1,count)=p1(1);
    end
    count=count+1;
end
prp11=reshape(BCC_CSM_prp1,[n2,m2]);
p11=rot90(prp11);
P1=p11(:,1:m2);
P2=p11(:,(m2+1):n2);
BCC_CSM_prP1=[P2 P1];
BCC_CSM_prp2(1,n1)=0;
count=1;
for i=1:n1
    temp1=BCC_CSM_prA(:,i);
    temp2=isnan(BCC_CSM_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BCC_CSM_prp2(1,count)=nan;
    else
        A=BCC_CSM_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        BCC_CSM_prp2(1,count)=p2(1);
    end
    count=count+1;
end
prp22=reshape(BCC_CSM_prp2,[n2,m2]);
p22=rot90(prp22);
P1=p22(:,1:m2);
P2=p22(:,(m2+1):n2);
BCC_CSM_prP2=[P2 P1];
BCC_CSM_prp3(1,n1)=0;
count=1;
for i=1:n1
    temp1=BCC_CSM_prA(:,i);
    temp2=isnan(BCC_CSM_prA(:,i));
    if  max(temp1)-min(temp1)==0 || sum(temp2)==N
        BCC_CSM_prp3(1,count)=nan;
    else
        A=BCC_CSM_prA(:,i);
        var=sum((A(:,1)-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:(m1/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        BCC_CSM_prp3(1,count)=p3(1);
    end
    count=count+1;
end
prp33=reshape(BCC_CSM_prp3,[n2,m2]);
p33=rot90(prp33);
P1=p33(:,1:m2);
P2=p33(:,(m2+1):n2);
BCC_CSM_prP3=[P2 P1];
% Change the spatial resolution of the model's P_kw to the same as SMAP
[p,~]=size(BCC_CSM_lonData);
[q,~]=size(BCC_CSM_latData);
m=BCC_CSM_lonData(p);
BCC_CSM_lon=linspace(-m/2,m/2,p);
BCC_CSM_lon=BCC_CSM_lon.';
BNU_lon1=BCC_CSM_lon(1:p/2,1);
BNU_lon2=BCC_CSM_lon((p/2+1):p,1);
BNU_lon=[BNU_lon2;BNU_lon1];
count=1;
for i=1:q
    BCC_CSM_lat_D(count:count+(p-1),1)=BCC_CSM_latData(i,1);
    BCC_CSM_lon_D(count:count+(p-1),1)=BNU_lon;
    count=count+p;
end
BCC_CSM_prp1=BCC_CSM_prp1.';
BCC_CSM_prp2=BCC_CSM_prp2.';
BCC_CSM_prp3=BCC_CSM_prp3.';
BCC_CSM_alpha1=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_prp1];
BCC_CSM_alpha2=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_prp2];
BCC_CSM_alpha3=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_prp3];
[m1,~]=size(BCC_CSM_alpha1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
BCC_CSM_new_P1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha1(i,1)-SMAP1(j,1))+abs(BCC_CSM_alpha1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P1(j,1)=BCC_CSM_alpha1(idx,3);
    disp(j)
end
new_P1=reshape(BCC_CSM_new_P1,[964,406]);
BCC_CSM_new_prP1=new_P1.';
[m1,~]=size(BCC_CSM_alpha2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
BCC_CSM_new_P2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha2(i,1)-SMAP2(j,1))+abs(BCC_CSM_alpha2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P2(j,1)=BCC_CSM_alpha2(idx,3);
    disp(j)
end
new_mf2=reshape(BCC_CSM_new_P2,[964,406]);
BCC_CSM_new_prP2=new_mf2.';
[m1,~]=size(BCC_CSM_alpha3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BCC_CSM_new_P3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_alpha3(i,1)-SMAP3(j,1))+abs(BCC_CSM_alpha3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_P3(j,1)=BCC_CSM_alpha3(idx,3);
    disp(j)
end
new_mf3=reshape(BCC_CSM_new_P3,[964,406]);
BCC_CSM_new_prP3=new_mf3.';
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
BCC_CSM_prP1_new=BCC_CSM_new_prP1./temp;
BCC_CSM_prP2_new=BCC_CSM_new_prP2./temp;
BCC_CSM_prP3_new=BCC_CSM_new_prP3./temp;
% Differences of P_kw between models and observations over different time scales
BCC_CSM_PPdif1=BCC_CSM_prP1_new-ERA5_P1;
BCC_CSM_PPdif2=BCC_CSM_prP2_new-ERA5_P2;
BCC_CSM_PPdif3=BCC_CSM_prP3_new-ERA5_P3;