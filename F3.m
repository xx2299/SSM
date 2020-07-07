% % Noise color of SSM (SMAP)
N=1389;
Fs=1;
F=((1:N)-1)*Fs/N;
x=F(1:(N+1)/2);
a=find(abs(x-1/7)<=0.0004);
b=find(abs(x-1/30)<=0.0005);
c=find(abs(x-1/90)<=0.0004);
d=find(abs(x-1/365)<=0.0002);
e=2;
[m,n]=size(SM_Z);
SMp1(1,n)=0;
count=1;
for i=1:n
    temp=isnan(SM_Z(:,i));
    if sum(temp)>=1388
        SMp1(1,count)=nan;
    else
        A=SM_Z(:,i);
        var=sum((A-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        SMp1(1,count)=p1(1);
    end
    disp(i);
    count=count+1;
end
SMp1=single(SMp1);
SMp11=reshape(SMp1,[964,406]);
P1=imrotate(SMp11,-90);
SMP1=fliplr(P1);
[m,n]=size(SM_Z);
SMp2(1,n)=0;
count=1;
for i=1:n
    temp=isnan(SM_Z(:,i));
    if sum(temp)>=1388
        SMp2(1,count)=nan;
    else
        A=SM_Z(:,i);
        var=sum((A-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        SMp2(1,count)=p2(1);
    end
    disp(i);
    count=count+1;
end
SMp2=single(SMp2);
SMp22=reshape(SMp2,[964,406]);
P2=imrotate(SMp22,-90);
SMP2=fliplr(P2);
[m,n]=size(SM_Z);
SMp3(1,n)=0;
count=1;
for i=1:n
    temp=isnan(SM_Z(:,i));
    if sum(temp)>=1388
        SMp3(1,count)=nan;
    else
        A=SM_Z(:,i);
        var=sum((A-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        SMp3(1,count)=p3(1);
    end
    disp(i);
    count=count+1;
end
SMp3=single(SMp3);
SMp33=reshape(SMp3,[964,406]);
P3=imrotate(SMp33,-90);
SMP3=fliplr(P3);

% % Noise color of P (ERA5)
% Get P_kw over different time scales
[~,n1]=size(A);
N=8752;
Fs=8;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.0005);
b=find(abs(x-1/30)<=0.00045);
c=find(abs(x-1/90)<=0.0002);
d=find(abs(x-1/365)<=0.0002);
A1=A(:,1:(n1+1)/2);
[m,n2]=size(A1);
temp=isnan(A1);
temp=single(temp);
temp_sum=sum(temp);
Z(m,n2)=0;
count=1;
for i=1:n2
    if temp_sum(1,i)==0
        Z(:,count)=A1(:,i);
    else
        data=fillmissing(A1(:,i),'movmedian',10);
        data=double(data);
        Z(:,count)=fillgaps(data,80,40);
    end
    count=count+1;
    disp(i);
end
Z=single(Z);
[m,n]=size(Z);
Pp1_1(1,n)=0;
Pp1_2(1,n)=0;
Pp1_3(1,n)=0;
count=1;
for i=1:n
    A=Z(:,i);
    var=sum((A-mean(A)).^2)/length(A);
    Y=fft(A);
    Y=Y(1:(m/2),:);
    Ayy=abs(Y).^2;
    Ayy=Ayy/(N/2)/var;
    X=[ones(length(Ayy),1),x'];
    p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
    Pp1_1(1,count)=p1(1);
    p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
    Pp1_2(1,count)=p2(1);
    p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
    Pp1_3(1,count)=p3(1);
    disp(i);
    count=count+1;
end
Pp1_1=single(Pp1_1);
Pp1_2=single(Pp1_2);
Pp1_3=single(Pp1_3);
N=8776;
Fs=8;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.0005);
b=find(abs(x-1/30)<=0.00045);
c=find(abs(x-1/90)<=0.0002);
d=find(abs(x-1/365)<=0.0002);
A2=A(:,(n1/2+1):n1);
[m,n2]=size(A2);
temp=isnan(A2);
temp=single(temp);
temp_sum=sum(temp);
Z2(m,n2)=0;
count=1;
for i=1:n2
    if temp_sum(1,i)==0
        Z2(:,count)=A2(:,i);
    else
        data=fillmissing(A2(:,i),'movmedian',10);
        data=double(data);
        Z2(:,count)=fillgaps(data,80,40);
    end
    count=count+1;
    disp(i);
end
Z2=single(Z2);
[m,n]=size(Z2);
Pp2_1(1,n)=0;
Pp2_2(1,n)=0;
Pp2_3(1,n)=0;
count=1;
for i=1:n
    A=Z2(:,i);
    var=sum((A-mean(A)).^2)/length(A);
    Y=fft(A);
    Y=Y(1:(m/2),:);
    Ayy=abs(Y).^2;
    Ayy=Ayy/(N/2)/var;
    X=[ones(length(Ayy),1),x'];
    p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
    Pp2_1(1,count)=p1(1);
    p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
    Pp2_2(1,count)=p2(1);
    p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
    Pp2_3(1,count)=p3(1);
    disp(i);
    count=count+1;
end
Pp2_1=single(Pp2_1);
Pp2_2=single(Pp2_2);
Pp2_3=single(Pp2_3);
Pp_1=[Pp1_1,Pp2_1];
Pp_2=[Pp1_2,Pp2_2];
Pp_3=[Pp1_3,Pp2_3];
% Change the spatial resolution of P_kw to the same as SSM_kw
temp=isnan(A150401);
temp=single(temp);
temp_sum=sum(temp);
idx=find(temp_sum==120);
Pp1_tot(1,6483600)=0;
Pp1_tot(1,idx)=nan;
Pp1_tot=single(Pp1_tot);
nan_temp=isnan(Pp1_tot);
nan_temp=single(nan_temp);
idx1=find(nan_temp==0);
Pp1_tot(1,idx1)=Pp_1;
Pp2_tot(1,6483600)=0;
Pp2_tot(1,idx)=nan;
Pp2_tot=single(Pp2_tot);
nan_temp=isnan(Pp2_tot);
nan_temp=single(nan_temp);
idx2=find(nan_temp==0);
Pp2_tot(1,idx2)=Pp_2;
Pp3_tot(1,6483600)=0;
Pp3_tot(1,idx)=nan;
Pp3_tot=single(Pp3_tot);
nan_temp=isnan(Pp3_tot);
nan_temp=single(nan_temp);
idx3=find(nan_temp==0);
Pp3_tot(1,idx3)=Pp_3;
Pp4_tot(1,6483600)=0;
Pp4_tot(1,idx)=nan;
Pp4_tot=single(Pp4_tot);
nan_temp=isnan(Pp4_tot);
nan_temp=single(nan_temp);
idx4=find(nan_temp==0);
Pp4_tot(1,idx4)=Pp_4;
tp1_1=reshape(Pp1_tot,[3600 1801]);
tp1_2=flipud(tp1_1);
tp1=rot90(tp1_2,-1);
tp1(1502:1801,:)=nan;
a1=tp1(:,1:1800);
a2=tp1(:,1801:3600);
Pp_a3=[a2,a1];
tp2_1=reshape(Pp2_tot,[3600 1801]);
tp2_2=flipud(tp2_1);
tp2=rot90(tp2_2,-1);
tp2(1502:1801,:)=nan;
b1=tp2(:,1:1800);
b2=tp2(:,1801:3600);
Pp_b3=[b2,b1];
tp3_1=reshape(Pp3_tot,[3600 1801]);
tp3_2=flipud(tp3_1);
tp3=rot90(tp3_2,-1);
tp3(1502:1801,:)=nan;
c1=tp3(:,1:1800);
c2=tp3(:,1801:3600);
Pp_c3=[c2,c1];
tp4_1=reshape(Pp4_tot,[3600 1801]);
tp4_2=flipud(tp4_1);
tp4=rot90(tp4_2,-1);
tp4(1502:1801,:)=nan;
d1=tp4(:,1:1800);
d2=tp4(:,1801:3600);
Pp_d3=[d2,d1];
subplot(2,2,1)
map1=imagesc(Pp_a3);
set(map1,'alphadata',~isnan(Pp_a3))
colorbar
subplot(2,2,2)
map2=imagesc(Pp_b3);
set(map2,'alphadata',~isnan(Pp_b3))
colorbar
subplot(2,2,3)
map3=imagesc(Pp_c3);
set(map3,'alphadata',~isnan(Pp_c3))
colorbar
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
ERA5_tp1=Pp_a3.';
ERA5_PP1=reshape(ERA5_tp1,[6483600,1]);
ERA5_tp2=Pp_b3.';
ERA5_PP2=reshape(ERA5_tp2,[6483600,1]);
ERA5_tp3=Pp_c3.';
ERA5_PP3=reshape(ERA5_tp3,[6483600,1]);
ERAP1=[ERA5_lat,ERA5_lon,ERA5_PP1];
ERAP2=[ERA5_lat,ERA5_lon,ERA5_PP2];
ERAP3=[ERA5_lat,ERA5_lon,ERA5_PP3];
[m1,~]=size(ERAP1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
ERA5_new_f1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(ERAP1(i,1)-SMAP1(j,1))+abs(ERAP1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    ERA5_new_f1(j,1)=ERAP1(idx,3);
    disp(j)
end
new_f1=reshape(ERA5_new_f1,[964,406]);
ERA5_P1=new_f1.';
[m1,~]=size(ERAP2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
ERA5_new_f2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(ERAP2(i,1)-SMAP2(j,1))+abs(ERAP2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    ERA5_new_f2(j,1)=ERAP2(idx,3);
    disp(j)
end
new_f2=reshape(ERA5_new_f2,[964,406]);
ERA5_P2=new_f2.';
[m1,~]=size(ERAP3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
ERA5_new_f3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(ERAP3(i,1)-SMAP3(j,1))+abs(ERAP3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    ERA5_new_f3(j,1)=ERAP3(idx,3);
    disp(j)
end
new_f3=reshape(ERA5_new_f3,[964,406]);
ERA5_P3=new_f3.';

% % Noise color of ET (GLEAM)
% Get ET_kw over different time scales
EA11=EA1(:,1:518400);
EA12=EA1(:,518401:1036800);
EA21=EA2(:,1:518400);
EA22=EA2(:,518401:1036800);
EA1121=[EA11;EA21];
EA1222=[EA12;EA22];
N=1097;
Fs=1;
F=((1:N)-1)*Fs/N;
x=F(1:(N+1)/2);
a=find(abs(x-1/7)<=0.0003);
b=find(abs(x-1/30)<=0.0004);
c=find(abs(x-1/90)<=0.0002);
d=find(abs(x-1/365)<=0.00001);
[m,n]=size(EA1121);
Ep11(1,n)=0;
count=1;
for i=1:n
    temp=isnan(EA1121(:,i));
    if sum(temp)>=1096
        Ep11(1,count)=nan;
    else
        A=EA1121(:,i);
        var=sum((A-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        Ep11(1,count)=p1(1);
    end
    disp(i);
    count=count+1;
end
Ep12(1,n)=0;
count=1;
for i=1:n
    temp=isnan(EA1222(:,i));
    if sum(temp)>=1096
        Ep12(1,count)=nan;
    else
        A=EA1222(:,i);
        var=sum((A-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
        Ep12(1,count)=p1(1);
    end
    disp(i);
    count=count+1;
end
Ep11=single(Ep11);
Ep12=single(Ep12);
Ep1112=[Ep11,Ep12];
Ep1=reshape(Ep1112,[1440,720]);
E1=imrotate(Ep1,-90);
EP1=fliplr(E1);
[m,n]=size(EA1121);
Ep21(1,n)=0;
count=1;
for i=1:n
    temp=isnan(EA1121(:,i));
    if sum(temp)>=1096
        Ep21(1,count)=nan;
    else
        A=EA1121(:,i);
        var=sum((A-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        Ep21(1,count)=p2(1);
    end
    disp(i);
    count=count+1;
end
Ep22(1,n)=0;
count=1;
for i=1:n
    temp=isnan(EA1222(:,i));
    if sum(temp)>=1096
        Ep22(1,count)=nan;
    else
        A=EA1222(:,i);
        var=sum((A-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
        Ep22(1,count)=p1(1);
    end
    disp(i);
    count=count+1;
end
Ep21=single(Ep21);
Ep22=single(Ep22);
Ep2122=[Ep21,Ep22];
Ep2=reshape(Ep2122,[1440,720]);
E2=imrotate(Ep2,-90);
EP2=fliplr(E2);
[m,n]=size(EA1121);
Ep31(1,n)=0;
count=1;
for i=1:n
    temp=isnan(EA1121(:,i));
    if sum(temp)>=1096
        Ep31(1,count)=nan;
    else
        A=EA1121(:,i);
        var=sum((A-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        Ep31(1,count)=p3(1);
    end
    disp(i);
    count=count+1;
end
Ep32(1,n)=0;
count=1;
for i=1:n
    temp=isnan(EA1222(:,i));
    if sum(temp)>=1096
        Ep32(1,count)=nan;
    else
        A=EA1222(:,i);
        var=sum((A-mean(A)).^2)/length(A);
        Y=fft(A);
        Y=Y(1:((m+1)/2),:);
        Ayy=abs(Y).^2;
        Ayy=Ayy/(N/2)/var;
        X=[ones(length(Ayy),1),x'];
        p1=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
        Ep32(1,count)=p1(1);
    end
    disp(i);
    count=count+1;
end
Ep31=single(Ep31);
Ep32=single(Ep32);
Ep3132=[Ep31,Ep32];
Ep3=reshape(Ep3132,[1440,720]);
E3=imrotate(Ep3,-90);
EP3=fliplr(E3);
% Change the spatial resolution of ET_kw to the same as SSM_kw
EP1(601:720,:)=nan;
EP2(601:720,:)=nan;
EP3(601:720,:)=nan;
count=1;
for i=1:720
    GLEAM_lat(count:count+1439,1)=latGLEAM(i,1);
    GLEAM_lon(count:count+1439,1)=lonGLEAM;
    count=count+1440;
end
GLEAM_lat=single(GLEAM_lat);
GLEAM_lon=single(GLEAM_lon);
[m,n]=size(EP1);
EP1_t=EP1.';
EP2_t=EP2.';
EP3_t=EP3.';
EP1_reshape=reshape(EP1_t,[1,m*n]);
EP1_reshape=EP1_reshape.';
EP2_reshape=reshape(EP2_t,[1,m*n]);
EP2_reshape=EP2_reshape.';
EP3_reshape=reshape(EP3_t,[1,m*n]);
EP3_reshape=EP3_reshape.';
GLEAMP1=[GLEAM_lat,GLEAM_lon,EP1_reshape];
GLEAMP2=[GLEAM_lat,GLEAM_lon,EP2_reshape];
GLEAMP3=[GLEAM_lat,GLEAM_lon,EP3_reshape];
[m1,~]=size(GLEAMP1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
GLEAMP_new_f1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GLEAMP1(i,1)-SMAP1(j,1))+abs(GLEAMP1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GLEAMP_new_f1(j,1)=GLEAMP1(idx,3);
    disp(j)
end
new_f1=reshape(GLEAMP_new_f1,[964,406]);
GLEAM_P1=new_f1.';
[m1,~]=size(GLEAMP2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
GLEAMP_new_f2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GLEAMP2(i,1)-SMAP2(j,1))+abs(GLEAMP2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GLEAMP_new_f2(j,1)=GLEAMP2(idx,3);
    disp(j)
end
new_f2=reshape(GLEAMP_new_f2,[964,406]);
GLEAM_P2=new_f2.';
[m1,~]=size(GLEAMP3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
GLEAMP_new_f3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GLEAMP3(i,1)-SMAP3(j,1))+abs(GLEAMP3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GLEAMP_new_f3(j,1)=GLEAMP3(idx,3);
    disp(j)
end
new_f3=reshape(GLEAMP_new_f3,[964,406]);
GLEAM_P3=new_f3.';

% % Noise color of EP(ERA5)
% Get EP_kw over different time scales
N=8768;
Fs=8;
F=((1:N)-1)*Fs/N;
x=F(1:N/2);
a=find(abs(x-1/7)<=0.0005);
b=find(abs(x-1/30)<=0.00045);
c=find(abs(x-1/90)<=0.0002);
d=find(abs(x-1/365)<=0.0002);
[m,n]=size(A);
temp=isnan(A);
temp=single(temp);
temp_sum=sum(temp);
Z(m,n)=0;
count=1;
for i=1:n
    if temp_sum(1,i)==0
        Z(:,count)=A(:,i);
    else
        data=fillmissing(A(:,i),'movmedian',10);
        data=double(data);
        Z(:,count)=fillgaps(data,80,40);
    end
    count=count+1;
    disp(i);
end
Z=single(Z);
PEp_1(1,n)=0;
PEp_2(1,n)=0;
PEp_3(1,n)=0;
count=1;
for i=1:n
    A=A1(:,i);
    var=sum((A-mean(A)).^2)/length(A);
    Y=fft(A);
    Y=Y(1:(m/2),:);
    Ayy=abs(Y).^2;
    Ayy=Ayy/(N/2)/var;
    X=[ones(length(Ayy),1),x'];
    p1=polyfit(log(X(b:a,2)),log(Ayy(b:a)),1);
    PEp_1(1,count)=p1(1);
    p2=polyfit(log(X(c:b,2)),log(Ayy(c:b)),1);
    PEp_2(1,count)=p2(1);
    p3=polyfit(log(X(d:c,2)),log(Ayy(d:c)),1);
    PEp_3(1,count)=p3(1);
    disp(i);
    count=count+1;
end
PEp_1=single(PEp_1);
PEp_2=single(PEp_2);
PEp_3=single(PEp_3);
PEp_1(find(PEp_1==PEp_1(1,1)))=nan;
tPEp1_1=reshape(PEp_1,[1440 721]);
tPEp1_2=flipud(tPEp1_1);
tPEp1=rot90(tPEp1_2,-1);
a1=tPEp1(:,1:720);
a2=tPEp1(:,721:1440);
PEp_a3=[a2,a1];
PEp_2(find(PEp_2==PEp_2(1,1)))=nan;
tPEp2_1=reshape(PEp_2,[1440 721]);
tPEp2_2=flipud(tPEp2_1);
tPEp2=rot90(tPEp2_2,-1);
b1=tPEp2(:,1:720);
b2=tPEp2(:,721:1440);
PEp_b3=[b2,b1];
PEp_3(find(PEp_3==PEp_3(1,1)))=nan;
tPEp3_1=reshape(PEp_3,[1440 721]);
tPEp3_2=flipud(tPEp3_1);
tPEp3=rot90(tPEp3_2,-1);
c1=tPEp3(:,1:720);
c2=tPEp3(:,721:1440);
PEp_c3=[c2,c1];
% Change the spatial resolution of EP_kw to the same as SSM_kw
lon=linspace(-179.875,179.875,1440);
lon=lon.';
count=1;
for i=1:721
    ERA5_lat(count:count+1439,1)=latPE(i,1);
    ERA5_lon(count:count+1439,1)=lon;
    count=count+1440;
end
ERA5_lat=single(ERA5_lat);
ERA5_lon=single(ERA5_lon);
ERA5_tp1=PEp_a3.';
ERA5_PE1=reshape(ERA5_tp1,[1038240,1]);
ERA5_tp2=PEp_b3.';
ERA5_PE2=reshape(ERA5_tp2,[1038240,1]);
ERA5_tp3=PEp_c3.';
ERA5_PE3=reshape(ERA5_tp3,[1038240,1]);
ERAP1=[ERA5_lat,ERA5_lon,ERA5_PE1];
ERAP2=[ERA5_lat,ERA5_lon,ERA5_PE2];
ERAP3=[ERA5_lat,ERA5_lon,ERA5_PE3];
[m1,~]=size(ERAP1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
ERA5_new_f1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(ERAP1(i,1)-SMAP1(j,1))+abs(ERAP1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    ERA5_new_f1(j,1)=ERAP1(idx,3);
    disp(j)
end
new_f1=reshape(ERA5_new_f1,[964,406]);
ERA5_EpP1=new_f1.';
[m1,~]=size(ERAP2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
ERA5_new_f2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(ERAP2(i,1)-SMAP2(j,1))+abs(ERAP2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    ERA5_new_f2(j,1)=ERAP2(idx,3);
    disp(j)
end
new_f2=reshape(ERA5_new_f2,[964,406]);
ERA5_EpP2=new_f2.';
[m1,~]=size(ERAP3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
ERA5_new_f3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(ERAP3(i,1)-SMAP3(j,1))+abs(ERAP3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    ERA5_new_f3(j,1)=ERAP3(idx,3);
    disp(j)
end
new_f3=reshape(ERA5_new_f3,[964,406]);
ERA5_EpP3=new_f3.';

% % Noise color of SSM,ET,P,and PE over different time scales (Figure 3)
temp=~isnan(SMP1);
temp=single(temp);
temp(temp==0)=nan;
GLEAM_P1_new=GLEAM_P1./temp;
GLEAM_P2_new=GLEAM_P2./temp;
GLEAM_P3_new=GLEAM_P3./temp;
ERA5_P1_new=ERA5_P1./temp;
ERA5_P2_new=ERA5_P2./temp;
ERA5_P3_new=ERA5_P3./temp;
ERA5_EpP1_new=ERA5_EpP1./temp;
ERA5_EpP2_new=ERA5_EpP2./temp;
ERA5_EpP3_new=ERA5_EpP3./temp;
SM_mask11=abs(SMP1-2)<0.5;
SM_mask12=abs(SMP1-1)<0.5;
SM_mask13=abs(SMP1-0)<0.5;
SM_mask14=abs(SMP1-(-1))<0.5;
SM_mask15=abs(SMP1-(-2))<0.5;
SM_mask16=SMP1<=-2.5;
SM_mask21=abs(SMP2-2)<0.5;
SM_mask22=abs(SMP2-1)<0.5;
SM_mask23=abs(SMP2-0)<0.5;
SM_mask24=abs(SMP2-(-1))<0.5;
SM_mask25=abs(SMP2-(-2))<0.5;
SM_mask26=SMP2<=-2.5;
SM_mask31=abs(SMP3-2)<0.5;
SM_mask32=abs(SMP3-1)<0.5;
SM_mask33=abs(SMP3-0)<0.5;
SM_mask34=abs(SMP3-(-1))<0.5;
SM_mask35=abs(SMP3-(-2))<0.5;
SM_mask36=SMP3<=-2.5;
detalgx=lonSMAP;
detalgy=latSMAP;
[LON,LAT]=meshgrid(detalgx,detalgy);
subplot(2,3,1)
map1=pcolor(detalgx,detalgy,SMP1);
set(map1,'alphadata',~isnan(SMP1))
shading interp
colorbar
hold on
SM_h11=stipple(LON,LAT,SM_mask11,'density',500,'color',[0.93333 0.5098 0.93333],'marker','.','markersize',4);
SM_h12=stipple(LON,LAT,SM_mask12,'density',500,'color','b','marker','.','markersize',4);
SM_h13=stipple(LON,LAT,SM_mask13,'density',500,'color',[1 0.98039 0.94118],'marker','.','markersize',4);
SM_h14=stipple(LON,LAT,SM_mask14,'density',500,'color',[1 0.75294 0.79608],'marker','.','markersize',4);
SM_h15=stipple(LON,LAT,SM_mask15,'density',500,'color','r','marker','.','markersize',4);
SM_h16=stipple(LON,LAT,SM_mask16,'density',500,'color','k','marker','.','markersize',4);
hold off
set(gca,'Ydir','normal')
set(gca,'YTick',[-90,-60,-30,0,30,60,90])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(a)')
ylabel('SSM')
title('1/30<f<1/7 day^-1')
subplot(2,3,2)
map2=pcolor(detalgx,detalgy,SMP2);
set(map2,'alphadata',~isnan(SMP2))
shading interp
colorbar
hold on
SM_h21=stipple(LON,LAT,SM_mask21,'density',500,'color',[0.93333 0.5098 0.93333],'marker','.','markersize',4);
SM_h22=stipple(LON,LAT,SM_mask22,'density',500,'color','b','marker','.','markersize',4);
SM_h23=stipple(LON,LAT,SM_mask23,'density',500,'color',[1 0.98039 0.94118],'marker','.','markersize',4);
SM_h24=stipple(LON,LAT,SM_mask24,'density',500,'color',[1 0.75294 0.79608],'marker','.','markersize',4);
SM_h25=stipple(LON,LAT,SM_mask25,'density',500,'color','r','marker','.','markersize',4);
SM_h26=stipple(LON,LAT,SM_mask26,'density',500,'color','k','marker','.','markersize',4);
hold off
set(gca,'Ydir','normal')
set(gca,'YTick',[-90,-60,-30,0,30,60,90])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(b)')
title('1/90<f<1/30 day^-1')
subplot(2,3,3)
map3=pcolor(detalgx,detalgy,SMP3);
set(map3,'alphadata',~isnan(SMP3))
shading interp
colorbar
hold on
SM_h31=stipple(LON,LAT,SM_mask31,'density',500,'color',[0.93333 0.5098 0.93333],'marker','.','markersize',4);
SM_h32=stipple(LON,LAT,SM_mask32,'density',500,'color','b','marker','.','markersize',4);
SM_h33=stipple(LON,LAT,SM_mask33,'density',500,'color',[1 0.98039 0.94118],'marker','.','markersize',4);
SM_h34=stipple(LON,LAT,SM_mask34,'density',500,'color',[1 0.75294 0.79608],'marker','.','markersize',4);
SM_h35=stipple(LON,LAT,SM_mask35,'density',500,'color','r','marker','.','markersize',4);
SM_h36=stipple(LON,LAT,SM_mask36,'density',500,'color','k','marker','.','markersize',4);
hold off
set(gca,'Ydir','normal')
set(gca,'YTick',[-90,-60,-30,0,30,60,90])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(c)')
title('1/365<f<1/90 day^-1')
P_mask11=abs(ERA5_P1_new-2)<0.5;
P_mask12=abs(ERA5_P1_new-1)<0.5;
P_mask13=abs(ERA5_P1_new-0)<0.5;
P_mask14=abs(ERA5_P1_new-(-1))<0.5;
P_mask15=abs(ERA5_P1_new-(-2))<0.5;
P_mask16=ERA5_P1_new<=-2.5;
P_mask21=abs(ERA5_P2_new-2)<0.5;
P_mask22=abs(ERA5_P2_new-1)<0.5;
P_mask23=abs(ERA5_P2_new-0)<0.5;
P_mask24=abs(ERA5_P2_new-(-1))<0.5;
P_mask25=abs(ERA5_P2_new-(-2))<0.5;
P_mask26=ERA5_P2_new<=-2.5;
P_mask31=abs(ERA5_P3_new-2)<0.5;
P_mask32=abs(ERA5_P3_new-1)<0.5;
P_mask33=abs(ERA5_P3_new-0)<0.5;
P_mask34=abs(ERA5_P3_new-(-1))<0.5;
P_mask35=abs(ERA5_P3_new-(-2))<0.5;
P_mask36=ERA5_P3_new<=-2.5;
subplot(2,3,4)
map1=pcolor(detalgx,detalgy,ERA5_P1_new);
set(map1,'alphadata',~isnan(ERA5_P1_new))
shading interp
colorbar
hold on
P_h11=stipple(LON,LAT,P_mask11,'density',500,'color',[0.93333 0.5098 0.93333],'marker','.','markersize',4);
P_h12=stipple(LON,LAT,P_mask12,'density',500,'color','b','marker','.','markersize',4);
P_h13=stipple(LON,LAT,P_mask13,'density',500,'color',[1 0.98039 0.94118],'marker','.','markersize',4);
P_h14=stipple(LON,LAT,P_mask14,'density',500,'color',[1 0.75294 0.79608],'marker','.','markersize',4);
P_h15=stipple(LON,LAT,P_mask15,'density',500,'color','r','marker','.','markersize',4);
P_h16=stipple(LON,LAT,P_mask16,'density',500,'color','k','marker','.','markersize',4);
hold off
set(gca,'Ydir','normal')
set(gca,'YTick',[-90,-60,-30,0,30,60,90])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(d)')
ylabel('P')
subplot(2,3,5)
map2=pcolor(detalgx,detalgy,ERA5_P2_new);
set(map2,'alphadata',~isnan(ERA5_P2_new))
shading interp
colorbar
hold on
P_h21=stipple(LON,LAT,P_mask21,'density',500,'color',[0.93333 0.5098 0.93333],'marker','.','markersize',4);
P_h22=stipple(LON,LAT,P_mask22,'density',500,'color','b','marker','.','markersize',4);
P_h23=stipple(LON,LAT,P_mask23,'density',500,'color',[1 0.98039 0.94118],'marker','.','markersize',4);
P_h24=stipple(LON,LAT,P_mask24,'density',500,'color',[1 0.75294 0.79608],'marker','.','markersize',4);
P_h25=stipple(LON,LAT,P_mask25,'density',500,'color','r','marker','.','markersize',4);
P_h26=stipple(LON,LAT,P_mask26,'density',500,'color','k','marker','.','markersize',4);
hold off
set(gca,'Ydir','normal')
set(gca,'YTick',[-90,-60,-30,0,30,60,90])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(e)')
subplot(2,3,6)
map3=pcolor(detalgx,detalgy,ERA5_P3_new);
set(map3,'alphadata',~isnan(ERA5_P3_new))
shading interp
colorbar
hold on
P_h31=stipple(LON,LAT,P_mask31,'density',500,'color',[0.93333 0.5098 0.93333],'marker','.','markersize',4);
P_h32=stipple(LON,LAT,P_mask32,'density',500,'color','b','marker','.','markersize',4);
P_h33=stipple(LON,LAT,P_mask33,'density',500,'color',[1 0.98039 0.94118],'marker','.','markersize',4);
P_h34=stipple(LON,LAT,P_mask34,'density',500,'color',[1 0.75294 0.79608],'marker','.','markersize',4);
P_h35=stipple(LON,LAT,P_mask35,'density',500,'color','r','marker','.','markersize',4);
P_h36=stipple(LON,LAT,P_mask36,'density',500,'color','k','marker','.','markersize',4);
hold off
set(gca,'Ydir','normal')
set(gca,'YTick',[-90,-60,-30,0,30,60,90])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(f)')
E_mask11=abs(GLEAM_P1_new-2)<0.5;
E_mask12=abs(GLEAM_P1_new-1)<0.5;
E_mask13=abs(GLEAM_P1_new-0)<0.5;
E_mask14=abs(GLEAM_P1_new-(-1))<0.5;
E_mask15=abs(GLEAM_P1_new-(-2))<0.5;
E_mask16=GLEAM_P1_new<=-2.5;
E_mask21=abs(GLEAM_P2_new-2)<0.5;
E_mask22=abs(GLEAM_P2_new-1)<0.5;
E_mask23=abs(GLEAM_P2_new-0)<0.5;
E_mask24=abs(GLEAM_P2_new-(-1))<0.5;
E_mask25=abs(GLEAM_P2_new-(-2))<0.5;
E_mask26=GLEAM_P2_new<=-2.5;
E_mask31=abs(GLEAM_P3_new-2)<0.5;
E_mask32=abs(GLEAM_P3_new-1)<0.5;
E_mask33=abs(GLEAM_P3_new-0)<0.5;
E_mask34=abs(GLEAM_P3_new-(-1))<0.5;
E_mask35=abs(GLEAM_P3_new-(-2))<0.5;
E_mask36=GLEAM_P3_new<=-2.5;
subplot(2,3,1)
map1=pcolor(detalgx,detalgy,GLEAM_P1_new);
set(map1,'alphadata',~isnan(GLEAM_P1_new))
shading interp
colorbar
hold on
E_h11=stipple(LON,LAT,E_mask11,'density',500,'color',[0.93333 0.5098 0.93333],'marker','.','markersize',4);
E_h12=stipple(LON,LAT,E_mask12,'density',500,'color','b','marker','.','markersize',4);
E_h13=stipple(LON,LAT,E_mask13,'density',500,'color',[1 0.98039 0.94118],'marker','.','markersize',4);
E_h14=stipple(LON,LAT,E_mask14,'density',500,'color',[1 0.75294 0.79608],'marker','.','markersize',4);
E_h15=stipple(LON,LAT,E_mask15,'density',500,'color','r','marker','.','markersize',4);
E_h16=stipple(LON,LAT,E_mask16,'density',500,'color','k','marker','.','markersize',4);
hold off
set(gca,'Ydir','normal')
set(gca,'YTick',[-90,-60,-30,0,30,60,90])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(g)')
ylabel('ET')
subplot(2,3,2)
map2=pcolor(detalgx,detalgy,GLEAM_P2_new);
set(map2,'alphadata',~isnan(GLEAM_P2_new))
shading interp
colorbar
hold on
E_h21=stipple(LON,LAT,E_mask21,'density',500,'color',[0.93333 0.5098 0.93333],'marker','.','markersize',4);
E_h22=stipple(LON,LAT,E_mask22,'density',500,'color','b','marker','.','markersize',4);
E_h23=stipple(LON,LAT,E_mask23,'density',500,'color',[1 0.98039 0.94118],'marker','.','markersize',4);
E_h24=stipple(LON,LAT,E_mask24,'density',500,'color',[1 0.75294 0.79608],'marker','.','markersize',4);
E_h25=stipple(LON,LAT,E_mask25,'density',500,'color','r','marker','.','markersize',4);
E_h26=stipple(LON,LAT,E_mask26,'density',500,'color','k','marker','.','markersize',4);
hold off
set(gca,'Ydir','normal')
set(gca,'YTick',[-90,-60,-30,0,30,60,90])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(h)')
subplot(2,3,3)
map3=pcolor(detalgx,detalgy,GLEAM_P3_new);
set(map3,'alphadata',~isnan(GLEAM_P3_new))
shading interp
colorbar
hold on
E_h31=stipple(LON,LAT,E_mask31,'density',500,'color',[0.93333 0.5098 0.93333],'marker','.','markersize',4);
E_h32=stipple(LON,LAT,E_mask32,'density',500,'color','b','marker','.','markersize',4);
E_h33=stipple(LON,LAT,E_mask33,'density',500,'color',[1 0.98039 0.94118],'marker','.','markersize',4);
E_h34=stipple(LON,LAT,E_mask34,'density',500,'color',[1 0.75294 0.79608],'marker','.','markersize',4);
E_h35=stipple(LON,LAT,E_mask35,'density',500,'color','r','marker','.','markersize',4);
E_h36=stipple(LON,LAT,E_mask36,'density',500,'color','k','marker','.','markersize',4);
hold off
set(gca,'Ydir','normal')
set(gca,'YTick',[-90,-60,-30,0,30,60,90])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(i)')
PE_mask11=abs(ERA5_EpP1_new-2)<0.5;
PE_mask12=abs(ERA5_EpP1_new-1)<0.5;
PE_mask13=abs(ERA5_EpP1_new-0)<0.5;
PE_mask14=abs(ERA5_EpP1_new-(-1))<0.5;
PE_mask15=abs(ERA5_EpP1_new-(-2))<0.5;
PE_mask16=ERA5_EpP1_new<=-2.5;
PE_mask21=abs(ERA5_EpP2_new-2)<0.5;
PE_mask22=abs(ERA5_EpP2_new-1)<0.5;
PE_mask23=abs(ERA5_EpP2_new-0)<0.5;
PE_mask24=abs(ERA5_EpP2_new-(-1))<0.5;
PE_mask25=abs(ERA5_EpP2_new-(-2))<0.5;
PE_mask26=ERA5_EpP2_new<=-2.5;
PE_mask31=abs(ERA5_EpP3_new-2)<0.5;
PE_mask32=abs(ERA5_EpP3_new-1)<0.5;
PE_mask33=abs(ERA5_EpP3_new-0)<0.5;
PE_mask34=abs(ERA5_EpP3_new-(-1))<0.5;
PE_mask35=abs(ERA5_EpP3_new-(-2))<0.5;
PE_mask36=ERA5_EpP3_new<=-2.5;
subplot(2,3,4)
map1=pcolor(detalgx,detalgy,ERA5_EpP1_new);
set(map1,'alphadata',~isnan(ERA5_EpP1_new))
shading interp
colorbar
hold on
PE_h11=stipple(LON,LAT,PE_mask11,'density',500,'color',[0.93333 0.5098 0.93333],'marker','.','markersize',4);
PE_h12=stipple(LON,LAT,PE_mask12,'density',500,'color','b','marker','.','markersize',4);
PE_h13=stipple(LON,LAT,PE_mask13,'density',500,'color',[1 0.98039 0.94118],'marker','.','markersize',4);
PE_h14=stipple(LON,LAT,PE_mask14,'density',500,'color',[1 0.75294 0.79608],'marker','.','markersize',4);
PE_h15=stipple(LON,LAT,PE_mask15,'density',500,'color','r','marker','.','markersize',4);
PE_h16=stipple(LON,LAT,PE_mask16,'density',500,'color','k','marker','.','markersize',4);
hold off
set(gca,'Ydir','normal')
set(gca,'YTick',[-90,-60,-30,0,30,60,90])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(j)')
ylabel('E_p')
subplot(2,3,5)
map2=pcolor(detalgx,detalgy,ERA5_EpP2_new);
set(map2,'alphadata',~isnan(ERA5_EpP2_new))
shading interp
colorbar
hold on
PE_h21=stipple(LON,LAT,PE_mask21,'density',500,'color',[0.93333 0.5098 0.93333],'marker','.','markersize',4);
PE_h22=stipple(LON,LAT,PE_mask22,'density',500,'color','b','marker','.','markersize',4);
PE_h23=stipple(LON,LAT,PE_mask23,'density',500,'color',[1 0.98039 0.94118],'marker','.','markersize',4);
PE_h24=stipple(LON,LAT,PE_mask24,'density',500,'color',[1 0.75294 0.79608],'marker','.','markersize',4);
PE_h25=stipple(LON,LAT,PE_mask25,'density',500,'color','r','marker','.','markersize',4);
PE_h26=stipple(LON,LAT,PE_mask26,'density',500,'color','k','marker','.','markersize',4);
hold off
set(gca,'Ydir','normal')
set(gca,'YTick',[-90,-60,-30,0,30,60,90])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(k)')
subplot(2,3,6)
map3=pcolor(detalgx,detalgy,ERA5_EpP3_new);
set(map3,'alphadata',~isnan(ERA5_EpP3_new))
shading interp
colorbar
hold on
PE_h31=stipple(LON,LAT,PE_mask31,'density',500,'color',[0.93333 0.5098 0.93333],'marker','.','markersize',4);
PE_h32=stipple(LON,LAT,PE_mask32,'density',500,'color','b','marker','.','markersize',4);
PE_h33=stipple(LON,LAT,PE_mask33,'density',500,'color',[1 0.98039 0.94118],'marker','.','markersize',4);
PE_h34=stipple(LON,LAT,PE_mask34,'density',500,'color',[1 0.75294 0.79608],'marker','.','markersize',4);
PE_h35=stipple(LON,LAT,PE_mask35,'density',500,'color','r','marker','.','markersize',4);
PE_h36=stipple(LON,LAT,PE_mask36,'density',500,'color','k','marker','.','markersize',4);
hold off
set(gca,'Ydir','normal')
set(gca,'YTick',[-90,-60,-30,0,30,60,90])
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(l)')
axes('position',[0.81,0.14,0.15,0.76])
axis off
h=legend([PE_h31,PE_h32,PE_h33,PE_h34,PE_h35,PE_h36],'Violet (-2)','Blue (-1)','White (0)','Pink (1)','Red (2)','Black (>2)','Location','southeast','FontWeight','bold','FontSize',6);
title(h,'Color of Noise (\beta)')