% % Latitudinal average of SSM_n for models and SMAP observation over different time scales
% Latitudinal average of SSM_n for SMAP observation
SMF1(find(mean_SM<0.1))=nan;
SMF2(find(mean_SM<0.1))=nan;
SMF3(find(mean_SM<0.1))=nan;
[m,~]=size(SMF1);
SMAP_Avg1(m,1)=0;
count=1;
for i=1:m
    SMAP_Avg1(count,1)=nanmean(SMF1(i,:));
    count=count+1;
end
SMAP_Avg2(m,1)=0;
count=1;
for i=1:m
    SMAP_Avg2(count,1)=nanmean(SMF2(i,:));
    count=count+1;
end
SMAP_Avg3(m,1)=0;
count=1;
for i=1:m
    SMAP_Avg3(count,1)=nanmean(SMF3(i,:));
    count=count+1;
end
% % Latitudinal average of SSM_n for models
% BCC_CSM
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(BCC_CSM_lonData);
[q,~]=size(BCC_CSM_latData);
m=BCC_CSM_lonData(p);
BCC_CSM_lon=linspace(-m/2,m/2,p);
BCC_CSM_lon=BCC_CSM_lon.';
Can_lon1=BCC_CSM_lon(1:p/2,1);
Can_lon2=BCC_CSM_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    BCC_CSM_lat_D(count:count+(p-1),1)=BCC_CSM_latData(i,1);
    BCC_CSM_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
BCC_CSM_mrsosf1=BCC_CSM_mrsos_f1.';
BCC_CSM_mrsosf2=BCC_CSM_mrsos_f2.';
BCC_CSM_mrsosf3=BCC_CSM_mrsos_f3.';
BCC_CSM_m1=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_mrsosf1];
BCC_CSM_m2=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_mrsosf2];
BCC_CSM_m3=[BCC_CSM_lat_D,BCC_CSM_lon_D,BCC_CSM_mrsosf3];
[m1,~]=size(BCC_CSM_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
BCC_CSM_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_m1(i,1)-SMAP1(j,1))+abs(BCC_CSM_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_mf1(j,1)=BCC_CSM_m1(idx,3);
    disp(j)
end
new_mf1=reshape(BCC_CSM_new_mf1,[964,406]);
BCC_CSM_mF1=new_mf1.';
[m1,~]=size(BCC_CSM_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
BCC_CSM_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_m2(i,1)-SMAP2(j,1))+abs(BCC_CSM_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_mf2(j,1)=BCC_CSM_m2(idx,3);
    disp(j)
end
new_mf2=reshape(BCC_CSM_new_mf2,[964,406]);
BCC_CSM_mF2=new_mf2.';
[m1,~]=size(BCC_CSM_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BCC_CSM_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BCC_CSM_m3(i,1)-SMAP3(j,1))+abs(BCC_CSM_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BCC_CSM_new_mf3(j,1)=BCC_CSM_m3(idx,3);
    disp(j)
end
new_mf3=reshape(BCC_CSM_new_mf3,[964,406]);
BCC_CSM_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
BCC_CSM_mF1_new=BCC_CSM_mF1./temp;
BCC_CSM_mF2_new=BCC_CSM_mF2./temp;
BCC_CSM_mF3_new=BCC_CSM_mF3./temp;
BCC_CSM_mF1_new(find(mean_SM<0.1))=nan;
BCC_CSM_mF2_new(find(mean_SM<0.1))=nan;
BCC_CSM_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for BCC_CSM
BCC_CSM_AVG1_new(406,1)=0;
count=1;
for i=1:406
    BCC_CSM_AVG1_new(count,1)=nanmean(BCC_CSM_mF1_new(i,:));
    count=count+1;
end
BCC_CSM_AVG2_new(406,1)=0;
count=1;
for i=1:406
    BCC_CSM_AVG2_new(count,1)=nanmean(BCC_CSM_mF2_new(i,:));
    count=count+1;
end
BCC_CSM_AVG3_new(406,1)=0;
count=1;
for i=1:406
    BCC_CSM_AVG3_new(count,1)=nanmean(BCC_CSM_mF3_new(i,:));
    count=count+1;
end
% BNU_ESM
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(BNU_ESM_lonData);
[q,~]=size(BNU_ESM_latData);
m=BNU_ESM_lonData(p);
BNU_ESM_lon=linspace(-m/2,m/2,p);
BNU_ESM_lon=BNU_ESM_lon.';
Can_lon1=BNU_ESM_lon(1:p/2,1);
Can_lon2=BNU_ESM_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    BNU_ESM_lat_D(count:count+(p-1),1)=BNU_ESM_latData(i,1);
    BNU_ESM_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
BNU_ESM_mrsosf1=BNU_ESM_mrsos_f1.';
BNU_ESM_mrsosf2=BNU_ESM_mrsos_f2.';
BNU_ESM_mrsosf3=BNU_ESM_mrsos_f3.';
BNU_ESM_m1=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_mrsosf1];
BNU_ESM_m2=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_mrsosf2];
BNU_ESM_m3=[BNU_ESM_lat_D,BNU_ESM_lon_D,BNU_ESM_mrsosf3];
[m1,~]=size(BNU_ESM_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
BNU_ESM_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_m1(i,1)-SMAP1(j,1))+abs(BNU_ESM_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_mf1(j,1)=BNU_ESM_m1(idx,3);
    disp(j)
end
new_mf1=reshape(BNU_ESM_new_mf1,[964,406]);
BNU_ESM_mF1=new_mf1.';
[m1,~]=size(BNU_ESM_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
BNU_ESM_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_m2(i,1)-SMAP2(j,1))+abs(BNU_ESM_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_mf2(j,1)=BNU_ESM_m2(idx,3);
    disp(j)
end
new_mf2=reshape(BNU_ESM_new_mf2,[964,406]);
BNU_ESM_mF2=new_mf2.';
[m1,~]=size(BNU_ESM_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
BNU_ESM_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(BNU_ESM_m3(i,1)-SMAP3(j,1))+abs(BNU_ESM_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    BNU_ESM_new_mf3(j,1)=BNU_ESM_m3(idx,3);
    disp(j)
end
new_mf3=reshape(BNU_ESM_new_mf3,[964,406]);
BNU_ESM_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
BNU_ESM_mF1_new=BNU_ESM_mF1./temp;
BNU_ESM_mF2_new=BNU_ESM_mF2./temp;
BNU_ESM_mF3_new=BNU_ESM_mF3./temp;
BNU_ESM_mF1_new(find(mean_SM<0.1))=nan;
BNU_ESM_mF2_new(find(mean_SM<0.1))=nan;
BNU_ESM_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for BNU_ESM
BNU_ESM_AVG1_new(406,1)=0;
count=1;
for i=1:406
    BNU_ESM_AVG1_new(count,1)=nanmean(BNU_ESM_mF1_new(i,:));
    count=count+1;
end
BNU_ESM_AVG2_new(406,1)=0;
count=1;
for i=1:406
    BNU_ESM_AVG2_new(count,1)=nanmean(BNU_ESM_mF2_new(i,:));
    count=count+1;
end
BNU_ESM_AVG3_new(406,1)=0;
count=1;
for i=1:406
    BNU_ESM_AVG3_new(count,1)=nanmean(BNU_ESM_mF3_new(i,:));
    count=count+1;
end
% CanESM2
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(CanESM2_lonData);
[q,~]=size(CanESM2_latData);
m=CanESM2_lonData(p);
CanESM2_lon=linspace(-m/2,m/2,p);
CanESM2_lon=CanESM2_lon.';
Can_lon1=CanESM2_lon(1:p/2,1);
Can_lon2=CanESM2_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    CanESM2_lat_D(count:count+(p-1),1)=CanESM2_latData(i,1);
    CanESM2_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
CanESM2_mrsosf1=CanESM2_mrsos_f1.';
CanESM2_mrsosf2=CanESM2_mrsos_f2.';
CanESM2_mrsosf3=CanESM2_mrsos_f3.';
CanESM2_m1=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_mrsosf1];
CanESM2_m2=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_mrsosf2];
CanESM2_m3=[CanESM2_lat_D,CanESM2_lon_D,CanESM2_mrsosf3];
[m1,~]=size(CanESM2_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
CanESM2_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_m1(i,1)-SMAP1(j,1))+abs(CanESM2_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_mf1(j,1)=CanESM2_m1(idx,3);
    disp(j)
end
new_mf1=reshape(CanESM2_new_mf1,[964,406]);
CanESM2_mF1=new_mf1.';
[m1,~]=size(CanESM2_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
CanESM2_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_m2(i,1)-SMAP2(j,1))+abs(CanESM2_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_mf2(j,1)=CanESM2_m2(idx,3);
    disp(j)
end
new_mf2=reshape(CanESM2_new_mf2,[964,406]);
CanESM2_mF2=new_mf2.';
[m1,~]=size(CanESM2_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
CanESM2_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CanESM2_m3(i,1)-SMAP3(j,1))+abs(CanESM2_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CanESM2_new_mf3(j,1)=CanESM2_m3(idx,3);
    disp(j)
end
new_mf3=reshape(CanESM2_new_mf3,[964,406]);
CanESM2_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
CanESM2_mF1_new=CanESM2_mF1./temp;
CanESM2_mF2_new=CanESM2_mF2./temp;
CanESM2_mF3_new=CanESM2_mF3./temp;
CanESM2_mF1_new(find(mean_SM<0.1))=nan;
CanESM2_mF2_new(find(mean_SM<0.1))=nan;
CanESM2_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for CanESM2
CanESM2_AVG1_new(406,1)=0;
count=1;
for i=1:406
    CanESM2_AVG1_new(count,1)=nanmean(CanESM2_mF1_new(i,:));
    count=count+1;
end
CanESM2_AVG2_new(406,1)=0;
count=1;
for i=1:406
    CanESM2_AVG2_new(count,1)=nanmean(CanESM2_mF2_new(i,:));
    count=count+1;
end
CanESM2_AVG3_new(406,1)=0;
count=1;
for i=1:406
    CanESM2_AVG3_new(count,1)=nanmean(CanESM2_mF3_new(i,:));
    count=count+1;
end
% CNRM_CM5
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(CNRM_CM5_lonData);
[q,~]=size(CNRM_CM5_latData);
m=CNRM_CM5_lonData(p);
CNRM_CM5_lon=linspace(-m/2,m/2,p);
CNRM_CM5_lon=CNRM_CM5_lon.';
Can_lon1=CNRM_CM5_lon(1:p/2,1);
Can_lon2=CNRM_CM5_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    CNRM_CM5_lat_D(count:count+(p-1),1)=CNRM_CM5_latData(i,1);
    CNRM_CM5_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
CNRM_CM5_mrsosf1=CNRM_CM5_mrsos_f1.';
CNRM_CM5_mrsosf2=CNRM_CM5_mrsos_f2.';
CNRM_CM5_mrsosf3=CNRM_CM5_mrsos_f3.';
CNRM_CM5_m1=[CNRM_CM5_lat_D,CNRM_CM5_lon_D,CNRM_CM5_mrsosf1];
CNRM_CM5_m2=[CNRM_CM5_lat_D,CNRM_CM5_lon_D,CNRM_CM5_mrsosf2];
CNRM_CM5_m3=[CNRM_CM5_lat_D,CNRM_CM5_lon_D,CNRM_CM5_mrsosf3];
[m1,~]=size(CNRM_CM5_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
CNRM_CM5_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CNRM_CM5_m1(i,1)-SMAP1(j,1))+abs(CNRM_CM5_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CNRM_CM5_new_mf1(j,1)=CNRM_CM5_m1(idx,3);
    disp(j)
end
new_mf1=reshape(CNRM_CM5_new_mf1,[964,406]);
CNRM_CM5_mF1=new_mf1.';
[m1,~]=size(CNRM_CM5_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
CNRM_CM5_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CNRM_CM5_m2(i,1)-SMAP2(j,1))+abs(CNRM_CM5_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CNRM_CM5_new_mf2(j,1)=CNRM_CM5_m2(idx,3);
    disp(j)
end
new_mf2=reshape(CNRM_CM5_new_mf2,[964,406]);
CNRM_CM5_mF2=new_mf2.';
[m1,~]=size(CNRM_CM5_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
CNRM_CM5_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CNRM_CM5_m3(i,1)-SMAP3(j,1))+abs(CNRM_CM5_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CNRM_CM5_new_mf3(j,1)=CNRM_CM5_m3(idx,3);
    disp(j)
end
new_mf3=reshape(CNRM_CM5_new_mf3,[964,406]);
CNRM_CM5_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
CNRM_CM5_mF1_new=CNRM_CM5_mF1./temp;
CNRM_CM5_mF2_new=CNRM_CM5_mF2./temp;
CNRM_CM5_mF3_new=CNRM_CM5_mF3./temp;
CNRM_CM5_mF1_new(find(mean_SM<0.1))=nan;
CNRM_CM5_mF2_new(find(mean_SM<0.1))=nan;
CNRM_CM5_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for CNRM_CM5
CNRM_CM5_AVG1_new(406,1)=0;
count=1;
for i=1:406
    CNRM_CM5_AVG1_new(count,1)=nanmean(CNRM_CM5_mF1_new(i,:));
    count=count+1;
end
CNRM_CM5_AVG2_new(406,1)=0;
count=1;
for i=1:406
    CNRM_CM5_AVG2_new(count,1)=nanmean(CNRM_CM5_mF2_new(i,:));
    count=count+1;
end
CNRM_CM5_AVG3_new(406,1)=0;
count=1;
for i=1:406
    CNRM_CM5_AVG3_new(count,1)=nanmean(CNRM_CM5_mF3_new(i,:));
    count=count+1;
end
% CSIRO_Mk3.6
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(CSIRO_Mk_lonData);
[q,~]=size(CSIRO_Mk_latData);
m=CSIRO_Mk_lonData(p);
CSIRO_Mk_lon=linspace(-m/2,m/2,p);
CSIRO_Mk_lon=CSIRO_Mk_lon.';
Can_lon1=CSIRO_Mk_lon(1:p/2,1);
Can_lon2=CSIRO_Mk_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    CSIRO_Mk_lat_D(count:count+(p-1),1)=CSIRO_Mk_latData(i,1);
    CSIRO_Mk_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
CSIRO_Mk_mrsosf1=CSIRO_Mk_mrsos_f1.';
CSIRO_Mk_mrsosf2=CSIRO_Mk_mrsos_f2.';
CSIRO_Mk_mrsosf3=CSIRO_Mk_mrsos_f3.';
CSIRO_Mk_m1=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_mrsosf1];
CSIRO_Mk_m2=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_mrsosf2];
CSIRO_Mk_m3=[CSIRO_Mk_lat_D,CSIRO_Mk_lon_D,CSIRO_Mk_mrsosf3];
[m1,~]=size(CSIRO_Mk_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
CSIRO_Mk_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_m1(i,1)-SMAP1(j,1))+abs(CSIRO_Mk_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_mf1(j,1)=CSIRO_Mk_m1(idx,3);
    disp(j)
end
new_mf1=reshape(CSIRO_Mk_new_mf1,[964,406]);
CSIRO_Mk_mF1=new_mf1.';
[m1,~]=size(CSIRO_Mk_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
CSIRO_Mk_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_m2(i,1)-SMAP2(j,1))+abs(CSIRO_Mk_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_mf2(j,1)=CSIRO_Mk_m2(idx,3);
    disp(j)
end
new_mf2=reshape(CSIRO_Mk_new_mf2,[964,406]);
CSIRO_Mk_mF2=new_mf2.';
[m1,~]=size(CSIRO_Mk_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
CSIRO_Mk_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(CSIRO_Mk_m3(i,1)-SMAP3(j,1))+abs(CSIRO_Mk_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    CSIRO_Mk_new_mf3(j,1)=CSIRO_Mk_m3(idx,3);
    disp(j)
end
new_mf3=reshape(CSIRO_Mk_new_mf3,[964,406]);
CSIRO_Mk_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
CSIRO_Mk_mF1_new=CSIRO_Mk_mF1./temp;
CSIRO_Mk_mF2_new=CSIRO_Mk_mF2./temp;
CSIRO_Mk_mF3_new=CSIRO_Mk_mF3./temp;
CSIRO_Mk_mF1_new(find(mean_SM<0.1))=nan;
CSIRO_Mk_mF2_new(find(mean_SM<0.1))=nan;
CSIRO_Mk_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for CSIRO_Mk3.6
CSIRO_Mk_AVG1_new(406,1)=0;
count=1;
for i=1:406
    CSIRO_Mk_AVG1_new(count,1)=nanmean(CSIRO_Mk_mF1_new(i,:));
    count=count+1;
end
CSIRO_Mk_AVG2_new(406,1)=0;
count=1;
for i=1:406
    CSIRO_Mk_AVG2_new(count,1)=nanmean(CSIRO_Mk_mF2_new(i,:));
    count=count+1;
end
CSIRO_Mk_AVG3_new(406,1)=0;
count=1;
for i=1:406
    CSIRO_Mk_AVG3_new(count,1)=nanmean(CSIRO_Mk_mF3_new(i,:));
    count=count+1;
end
% GFDL_CM3
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(GFDL_CM3_lonData);
[q,~]=size(GFDL_CM3_latData);
m=GFDL_CM3_lonData(p);
GFDL_CM3_lon=linspace(-m/2,m/2,p);
GFDL_CM3_lon=GFDL_CM3_lon.';
Can_lon1=GFDL_CM3_lon(1:p/2,1);
Can_lon2=GFDL_CM3_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    GFDL_CM3_lat_D(count:count+(p-1),1)=GFDL_CM3_latData(i,1);
    GFDL_CM3_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
GFDL_CM3_mrsosf1=GFDL_CM3_mrsos_f1.';
GFDL_CM3_mrsosf2=GFDL_CM3_mrsos_f2.';
GFDL_CM3_mrsosf3=GFDL_CM3_mrsos_f3.';
GFDL_CM3_m1=[GFDL_CM3_lat_D,GFDL_CM3_lon_D,GFDL_CM3_mrsosf1];
GFDL_CM3_m2=[GFDL_CM3_lat_D,GFDL_CM3_lon_D,GFDL_CM3_mrsosf2];
GFDL_CM3_m3=[GFDL_CM3_lat_D,GFDL_CM3_lon_D,GFDL_CM3_mrsosf3];
[m1,~]=size(GFDL_CM3_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
GFDL_CM3_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_CM3_m1(i,1)-SMAP1(j,1))+abs(GFDL_CM3_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_CM3_new_mf1(j,1)=GFDL_CM3_m1(idx,3);
    disp(j)
end
new_mf1=reshape(GFDL_CM3_new_mf1,[964,406]);
GFDL_CM3_mF1=new_mf1.';
[m1,~]=size(GFDL_CM3_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
GFDL_CM3_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_CM3_m2(i,1)-SMAP2(j,1))+abs(GFDL_CM3_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_CM3_new_mf2(j,1)=GFDL_CM3_m2(idx,3);
    disp(j)
end
new_mf2=reshape(GFDL_CM3_new_mf2,[964,406]);
GFDL_CM3_mF2=new_mf2.';
[m1,~]=size(GFDL_CM3_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
GFDL_CM3_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_CM3_m3(i,1)-SMAP3(j,1))+abs(GFDL_CM3_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_CM3_new_mf3(j,1)=GFDL_CM3_m3(idx,3);
    disp(j)
end
new_mf3=reshape(GFDL_CM3_new_mf3,[964,406]);
GFDL_CM3_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
GFDL_CM3_mF1_new=GFDL_CM3_mF1./temp;
GFDL_CM3_mF2_new=GFDL_CM3_mF2./temp;
GFDL_CM3_mF3_new=GFDL_CM3_mF3./temp;
GFDL_CM3_mF1_new(find(mean_SM<0.1))=nan;
GFDL_CM3_mF2_new(find(mean_SM<0.1))=nan;
GFDL_CM3_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for GFDL_CM3
GFDL_CM3_AVG1_new(406,1)=0;
count=1;
for i=1:406
    GFDL_CM3_AVG1_new(count,1)=nanmean(GFDL_CM3_mF1_new(i,:));
    count=count+1;
end
GFDL_CM3_AVG2_new(406,1)=0;
count=1;
for i=1:406
    GFDL_CM3_AVG2_new(count,1)=nanmean(GFDL_CM3_mF2_new(i,:));
    count=count+1;
end
GFDL_CM3_AVG3_new(406,1)=0;
count=1;
for i=1:406
    GFDL_CM3_AVG3_new(count,1)=nanmean(GFDL_CM3_mF3_new(i,:));
    count=count+1;
end
% GFDL_ESM2G
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(GFDL_ESM2G_lonData);
[q,~]=size(GFDL_ESM2G_latData);
m=GFDL_ESM2G_lonData(p);
GFDL_ESM2G_lon=linspace(-m/2,m/2,p);
GFDL_ESM2G_lon=GFDL_ESM2G_lon.';
Can_lon1=GFDL_ESM2G_lon(1:p/2,1);
Can_lon2=GFDL_ESM2G_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    GFDL_ESM2G_lat_D(count:count+(p-1),1)=GFDL_ESM2G_latData(i,1);
    GFDL_ESM2G_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
GFDL_ESM2G_mrsosf1=GFDL_ESM2G_mrsos_f1.';
GFDL_ESM2G_mrsosf2=GFDL_ESM2G_mrsos_f2.';
GFDL_ESM2G_mrsosf3=GFDL_ESM2G_mrsos_f3.';
GFDL_ESM2G_m1=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_mrsosf1];
GFDL_ESM2G_m2=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_mrsosf2];
GFDL_ESM2G_m3=[GFDL_ESM2G_lat_D,GFDL_ESM2G_lon_D,GFDL_ESM2G_mrsosf3];
[m1,~]=size(GFDL_ESM2G_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
GFDL_ESM2G_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_m1(i,1)-SMAP1(j,1))+abs(GFDL_ESM2G_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_mf1(j,1)=GFDL_ESM2G_m1(idx,3);
    disp(j)
end
new_mf1=reshape(GFDL_ESM2G_new_mf1,[964,406]);
GFDL_ESM2G_mF1=new_mf1.';
[m1,~]=size(GFDL_ESM2G_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
GFDL_ESM2G_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_m2(i,1)-SMAP2(j,1))+abs(GFDL_ESM2G_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_mf2(j,1)=GFDL_ESM2G_m2(idx,3);
    disp(j)
end
new_mf2=reshape(GFDL_ESM2G_new_mf2,[964,406]);
GFDL_ESM2G_mF2=new_mf2.';
[m1,~]=size(GFDL_ESM2G_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
GFDL_ESM2G_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2G_m3(i,1)-SMAP3(j,1))+abs(GFDL_ESM2G_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2G_new_mf3(j,1)=GFDL_ESM2G_m3(idx,3);
    disp(j)
end
new_mf3=reshape(GFDL_ESM2G_new_mf3,[964,406]);
GFDL_ESM2G_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
GFDL_ESM2G_mF1_new=GFDL_ESM2G_mF1./temp;
GFDL_ESM2G_mF2_new=GFDL_ESM2G_mF2./temp;
GFDL_ESM2G_mF3_new=GFDL_ESM2G_mF3./temp;
GFDL_ESM2G_mF1_new(find(mean_SM<0.1))=nan;
GFDL_ESM2G_mF2_new(find(mean_SM<0.1))=nan;
GFDL_ESM2G_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for GFDL_ESM2G
GFDL_ESM2G_AVG1_new(406,1)=0;
count=1;
for i=1:406
    GFDL_ESM2G_AVG1_new(count,1)=nanmean(GFDL_ESM2G_mF1_new(i,:));
    count=count+1;
end
GFDL_ESM2G_AVG2_new(406,1)=0;
count=1;
for i=1:406
    GFDL_ESM2G_AVG2_new(count,1)=nanmean(GFDL_ESM2G_mF2_new(i,:));
    count=count+1;
end
GFDL_ESM2G_AVG3_new(406,1)=0;
count=1;
for i=1:406
    GFDL_ESM2G_AVG3_new(count,1)=nanmean(GFDL_ESM2G_mF3_new(i,:));
    count=count+1;
end
% GFDL_ESM2M
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(GFDL_ESM2M_lonData);
[q,~]=size(GFDL_ESM2M_latData);
m=GFDL_ESM2M_lonData(p);
GFDL_ESM2M_lon=linspace(-m/2,m/2,p);
GFDL_ESM2M_lon=GFDL_ESM2M_lon.';
Can_lon1=GFDL_ESM2M_lon(1:p/2,1);
Can_lon2=GFDL_ESM2M_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    GFDL_ESM2M_lat_D(count:count+(p-1),1)=GFDL_ESM2M_latData(i,1);
    GFDL_ESM2M_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
GFDL_ESM2M_mrsosf1=GFDL_ESM2M_mrsos_f1.';
GFDL_ESM2M_mrsosf2=GFDL_ESM2M_mrsos_f2.';
GFDL_ESM2M_mrsosf3=GFDL_ESM2M_mrsos_f3.';
GFDL_ESM2M_m1=[GFDL_ESM2M_lat_D,GFDL_ESM2M_lon_D,GFDL_ESM2M_mrsosf1];
GFDL_ESM2M_m2=[GFDL_ESM2M_lat_D,GFDL_ESM2M_lon_D,GFDL_ESM2M_mrsosf2];
GFDL_ESM2M_m3=[GFDL_ESM2M_lat_D,GFDL_ESM2M_lon_D,GFDL_ESM2M_mrsosf3];
[m1,~]=size(GFDL_ESM2M_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
GFDL_ESM2M_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2M_m1(i,1)-SMAP1(j,1))+abs(GFDL_ESM2M_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2M_new_mf1(j,1)=GFDL_ESM2M_m1(idx,3);
    disp(j)
end
new_mf1=reshape(GFDL_ESM2M_new_mf1,[964,406]);
GFDL_ESM2M_mF1=new_mf1.';
[m1,~]=size(GFDL_ESM2M_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
GFDL_ESM2M_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2M_m2(i,1)-SMAP2(j,1))+abs(GFDL_ESM2M_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2M_new_mf2(j,1)=GFDL_ESM2M_m2(idx,3);
    disp(j)
end
new_mf2=reshape(GFDL_ESM2M_new_mf2,[964,406]);
GFDL_ESM2M_mF2=new_mf2.';
[m1,~]=size(GFDL_ESM2M_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
GFDL_ESM2M_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(GFDL_ESM2M_m3(i,1)-SMAP3(j,1))+abs(GFDL_ESM2M_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    GFDL_ESM2M_new_mf3(j,1)=GFDL_ESM2M_m3(idx,3);
    disp(j)
end
new_mf3=reshape(GFDL_ESM2M_new_mf3,[964,406]);
GFDL_ESM2M_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
GFDL_ESM2M_mF1_new=GFDL_ESM2M_mF1./temp;
GFDL_ESM2M_mF2_new=GFDL_ESM2M_mF2./temp;
GFDL_ESM2M_mF3_new=GFDL_ESM2M_mF3./temp;
GFDL_ESM2M_mF1_new(find(mean_SM<0.1))=nan;
GFDL_ESM2M_mF2_new(find(mean_SM<0.1))=nan;
GFDL_ESM2M_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for GFDL_ESM2M
GFDL_ESM2M_AVG1_new(406,1)=0;
count=1;
for i=1:406
    GFDL_ESM2M_AVG1_new(count,1)=nanmean(GFDL_ESM2M_mF1_new(i,:));
    count=count+1;
end
GFDL_ESM2M_AVG2_new(406,1)=0;
count=1;
for i=1:406
    GFDL_ESM2M_AVG2_new(count,1)=nanmean(GFDL_ESM2M_mF2_new(i,:));
    count=count+1;
end
GFDL_ESM2M_AVG3_new(406,1)=0;
count=1;
for i=1:406
    GFDL_ESM2M_AVG3_new(count,1)=nanmean(GFDL_ESM2M_mF3_new(i,:));
    count=count+1;
end
% HadGEM2_CC
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(HadGEM2_CC_lonData);
[q,~]=size(HadGEM2_CC_latData);
m=HadGEM2_CC_lonData(p);
HadGEM2_CC_lon=linspace(-m/2,m/2,p);
HadGEM2_CC_lon=HadGEM2_CC_lon.';
Can_lon1=HadGEM2_CC_lon(1:p/2,1);
Can_lon2=HadGEM2_CC_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    HadGEM2_CC_lat_D(count:count+(p-1),1)=HadGEM2_CC_latData(i,1);
    HadGEM2_CC_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
HadGEM2_CC_mrsosf1=HadGEM2_CC_mrsos_f1.';
HadGEM2_CC_mrsosf2=HadGEM2_CC_mrsos_f2.';
HadGEM2_CC_mrsosf3=HadGEM2_CC_mrsos_f3.';
HadGEM2_CC_m1=[HadGEM2_CC_lat_D,HadGEM2_CC_lon_D,HadGEM2_CC_mrsosf1];
HadGEM2_CC_m2=[HadGEM2_CC_lat_D,HadGEM2_CC_lon_D,HadGEM2_CC_mrsosf2];
HadGEM2_CC_m3=[HadGEM2_CC_lat_D,HadGEM2_CC_lon_D,HadGEM2_CC_mrsosf3];
[m1,~]=size(HadGEM2_CC_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
HadGEM2_CC_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(HadGEM2_CC_m1(i,1)-SMAP1(j,1))+abs(HadGEM2_CC_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    HadGEM2_CC_new_mf1(j,1)=HadGEM2_CC_m1(idx,3);
    disp(j)
end
new_mf1=reshape(HadGEM2_CC_new_mf1,[964,406]);
HadGEM2_CC_mF1=new_mf1.';
[m1,~]=size(HadGEM2_CC_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
HadGEM2_CC_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(HadGEM2_CC_m2(i,1)-SMAP2(j,1))+abs(HadGEM2_CC_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    HadGEM2_CC_new_mf2(j,1)=HadGEM2_CC_m2(idx,3);
    disp(j)
end
new_mf2=reshape(HadGEM2_CC_new_mf2,[964,406]);
HadGEM2_CC_mF2=new_mf2.';
[m1,~]=size(HadGEM2_CC_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
HadGEM2_CC_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(HadGEM2_CC_m3(i,1)-SMAP3(j,1))+abs(HadGEM2_CC_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    HadGEM2_CC_new_mf3(j,1)=HadGEM2_CC_m3(idx,3);
    disp(j)
end
new_mf3=reshape(HadGEM2_CC_new_mf3,[964,406]);
HadGEM2_CC_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
HadGEM2_CC_mF1_new=HadGEM2_CC_mF1./temp;
HadGEM2_CC_mF2_new=HadGEM2_CC_mF2./temp;
HadGEM2_CC_mF3_new=HadGEM2_CC_mF3./temp;
HadGEM2_CC_mF1_new(find(mean_SM<0.1))=nan;
HadGEM2_CC_mF2_new(find(mean_SM<0.1))=nan;
HadGEM2_CC_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for HadGEM2_CC
HadGEM2_CC_AVG1_new(406,1)=0;
count=1;
for i=1:406
    HadGEM2_CC_AVG1_new(count,1)=nanmean(HadGEM2_CC_mF1_new(i,:));
    count=count+1;
end
HadGEM2_CC_AVG2_new(406,1)=0;
count=1;
for i=1:406
    HadGEM2_CC_AVG2_new(count,1)=nanmean(HadGEM2_CC_mF2_new(i,:));
    count=count+1;
end
HadGEM2_CC_AVG3_new(406,1)=0;
count=1;
for i=1:406
    HadGEM2_CC_AVG3_new(count,1)=nanmean(HadGEM2_CC_mF3_new(i,:));
    count=count+1;
end
% HadGEM2_ES
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(HadGEM2_ES_lonData);
[q,~]=size(HadGEM2_ES_latData);
m=HadGEM2_ES_lonData(p);
HadGEM2_ES_lon=linspace(-m/2,m/2,p);
HadGEM2_ES_lon=HadGEM2_ES_lon.';
Can_lon1=HadGEM2_ES_lon(1:p/2,1);
Can_lon2=HadGEM2_ES_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    HadGEM2_ES_lat_D(count:count+(p-1),1)=HadGEM2_ES_latData(i,1);
    HadGEM2_ES_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
HadGEM2_ES_mrsosf1=HadGEM2_ES_mrsos_f1.';
HadGEM2_ES_mrsosf2=HadGEM2_ES_mrsos_f2.';
HadGEM2_ES_mrsosf3=HadGEM2_ES_mrsos_f3.';
HadGEM2_ES_m1=[HadGEM2_ES_lat_D,HadGEM2_ES_lon_D,HadGEM2_ES_mrsosf1];
HadGEM2_ES_m2=[HadGEM2_ES_lat_D,HadGEM2_ES_lon_D,HadGEM2_ES_mrsosf2];
HadGEM2_ES_m3=[HadGEM2_ES_lat_D,HadGEM2_ES_lon_D,HadGEM2_ES_mrsosf3];
[m1,~]=size(HadGEM2_ES_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
HadGEM2_ES_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(HadGEM2_ES_m1(i,1)-SMAP1(j,1))+abs(HadGEM2_ES_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    HadGEM2_ES_new_mf1(j,1)=HadGEM2_ES_m1(idx,3);
    disp(j)
end
new_mf1=reshape(HadGEM2_ES_new_mf1,[964,406]);
HadGEM2_ES_mF1=new_mf1.';
[m1,~]=size(HadGEM2_ES_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
HadGEM2_ES_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(HadGEM2_ES_m2(i,1)-SMAP2(j,1))+abs(HadGEM2_ES_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    HadGEM2_ES_new_mf2(j,1)=HadGEM2_ES_m2(idx,3);
    disp(j)
end
new_mf2=reshape(HadGEM2_ES_new_mf2,[964,406]);
HadGEM2_ES_mF2=new_mf2.';
[m1,~]=size(HadGEM2_ES_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
HadGEM2_ES_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(HadGEM2_ES_m3(i,1)-SMAP3(j,1))+abs(HadGEM2_ES_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    HadGEM2_ES_new_mf3(j,1)=HadGEM2_ES_m3(idx,3);
    disp(j)
end
new_mf3=reshape(HadGEM2_ES_new_mf3,[964,406]);
HadGEM2_ES_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
HadGEM2_ES_mF1_new=HadGEM2_ES_mF1./temp;
HadGEM2_ES_mF2_new=HadGEM2_ES_mF2./temp;
HadGEM2_ES_mF3_new=HadGEM2_ES_mF3./temp;
HadGEM2_ES_mF1_new(find(mean_SM<0.1))=nan;
HadGEM2_ES_mF2_new(find(mean_SM<0.1))=nan;
HadGEM2_ES_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for HadGEM2_ES
HadGEM2_ES_AVG1_new(406,1)=0;
count=1;
for i=1:406
    HadGEM2_ES_AVG1_new(count,1)=nanmean(HadGEM2_ES_mF1_new(i,:));
    count=count+1;
end
HadGEM2_ES_AVG2_new(406,1)=0;
count=1;
for i=1:406
    HadGEM2_ES_AVG2_new(count,1)=nanmean(HadGEM2_ES_mF2_new(i,:));
    count=count+1;
end
HadGEM2_ES_AVG3_new(406,1)=0;
count=1;
for i=1:406
    HadGEM2_ES_AVG3_new(count,1)=nanmean(HadGEM2_ES_mF3_new(i,:));
    count=count+1;
end
% Institute_for_Numerical_Mathematics
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(inmcm4_lonData);
[q,~]=size(inmcm4_latData);
m=inmcm4_lonData(p);
inmcm4_lon=linspace(-m/2,m/2,p);
inmcm4_lon=inmcm4_lon.';
Can_lon1=inmcm4_lon(1:p/2,1);
Can_lon2=inmcm4_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    inmcm4_lat_D(count:count+(p-1),1)=inmcm4_latData(i,1);
    inmcm4_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
inmcm4_mrsosf1=inmcm4_mrsos_f1.';
inmcm4_mrsosf2=inmcm4_mrsos_f2.';
inmcm4_mrsosf3=inmcm4_mrsos_f3.';
inmcm4_m1=[inmcm4_lat_D,inmcm4_lon_D,inmcm4_mrsosf1];
inmcm4_m2=[inmcm4_lat_D,inmcm4_lon_D,inmcm4_mrsosf2];
inmcm4_m3=[inmcm4_lat_D,inmcm4_lon_D,inmcm4_mrsosf3];
[m1,~]=size(inmcm4_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
inmcm4_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(inmcm4_m1(i,1)-SMAP1(j,1))+abs(inmcm4_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    inmcm4_new_mf1(j,1)=inmcm4_m1(idx,3);
    disp(j)
end
new_mf1=reshape(inmcm4_new_mf1,[964,406]);
inmcm4_mF1=new_mf1.';
[m1,~]=size(inmcm4_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
inmcm4_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(inmcm4_m2(i,1)-SMAP2(j,1))+abs(inmcm4_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    inmcm4_new_mf2(j,1)=inmcm4_m2(idx,3);
    disp(j)
end
new_mf2=reshape(inmcm4_new_mf2,[964,406]);
inmcm4_mF2=new_mf2.';
[m1,~]=size(inmcm4_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
inmcm4_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(inmcm4_m3(i,1)-SMAP3(j,1))+abs(inmcm4_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    inmcm4_new_mf3(j,1)=inmcm4_m3(idx,3);
    disp(j)
end
new_mf3=reshape(inmcm4_new_mf3,[964,406]);
inmcm4_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
inmcm4_mF1_new=inmcm4_mF1./temp;
inmcm4_mF2_new=inmcm4_mF2./temp;
inmcm4_mF3_new=inmcm4_mF3./temp;
inmcm4_mF1_new(find(mean_SM<0.1))=nan;
inmcm4_mF2_new(find(mean_SM<0.1))=nan;
inmcm4_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for inmcm4
inmcm4_AVG1_new(406,1)=0;
count=1;
for i=1:406
    inmcm4_AVG1_new(count,1)=nanmean(inmcm4_mF1_new(i,:));
    count=count+1;
end
inmcm4_AVG2_new(406,1)=0;
count=1;
for i=1:406
    inmcm4_AVG2_new(count,1)=nanmean(inmcm4_mF2_new(i,:));
    count=count+1;
end
inmcm4_AVG3_new(406,1)=0;
count=1;
for i=1:406
    inmcm4_AVG3_new(count,1)=nanmean(inmcm4_mF3_new(i,:));
    count=count+1;
end
% MIROC5
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(MIROC5_lonData);
[q,~]=size(MIROC5_latData);
m=MIROC5_lonData(p);
MIROC5_lon=linspace(-m/2,m/2,p);
MIROC5_lon=MIROC5_lon.';
Can_lon1=MIROC5_lon(1:p/2,1);
Can_lon2=MIROC5_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    MIROC5_lat_D(count:count+(p-1),1)=MIROC5_latData(i,1);
    MIROC5_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
MIROC5_mrsosf1=MIROC5_mrsos_f1.';
MIROC5_mrsosf2=MIROC5_mrsos_f2.';
MIROC5_mrsosf3=MIROC5_mrsos_f3.';
MIROC5_m1=[MIROC5_lat_D,MIROC5_lon_D,MIROC5_mrsosf1];
MIROC5_m2=[MIROC5_lat_D,MIROC5_lon_D,MIROC5_mrsosf2];
MIROC5_m3=[MIROC5_lat_D,MIROC5_lon_D,MIROC5_mrsosf3];
[m1,~]=size(MIROC5_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
MIROC5_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MIROC5_m1(i,1)-SMAP1(j,1))+abs(MIROC5_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MIROC5_new_mf1(j,1)=MIROC5_m1(idx,3);
    disp(j)
end
new_mf1=reshape(MIROC5_new_mf1,[964,406]);
MIROC5_mF1=new_mf1.';
[m1,~]=size(MIROC5_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
MIROC5_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MIROC5_m2(i,1)-SMAP2(j,1))+abs(MIROC5_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MIROC5_new_mf2(j,1)=MIROC5_m2(idx,3);
    disp(j)
end
new_mf2=reshape(MIROC5_new_mf2,[964,406]);
MIROC5_mF2=new_mf2.';
[m1,~]=size(MIROC5_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
MIROC5_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MIROC5_m3(i,1)-SMAP3(j,1))+abs(MIROC5_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MIROC5_new_mf3(j,1)=MIROC5_m3(idx,3);
    disp(j)
end
new_mf3=reshape(MIROC5_new_mf3,[964,406]);
MIROC5_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
MIROC5_mF1_new=MIROC5_mF1./temp;
MIROC5_mF2_new=MIROC5_mF2./temp;
MIROC5_mF3_new=MIROC5_mF3./temp;
MIROC5_mF1_new(find(mean_SM<0.1))=nan;
MIROC5_mF2_new(find(mean_SM<0.1))=nan;
MIROC5_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for MIROC5
MIROC5_AVG1_new(406,1)=0;
count=1;
for i=1:406
    MIROC5_AVG1_new(count,1)=nanmean(MIROC5_mF1_new(i,:));
    count=count+1;
end
MIROC5_AVG2_new(406,1)=0;
count=1;
for i=1:406
    MIROC5_AVG2_new(count,1)=nanmean(MIROC5_mF2_new(i,:));
    count=count+1;
end
MIROC5_AVG3_new(406,1)=0;
count=1;
for i=1:406
    MIROC5_AVG3_new(count,1)=nanmean(MIROC5_mF3_new(i,:));
    count=count+1;
end
% MIROC_ESM
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(MIROC_ESM_lonData);
[q,~]=size(MIROC_ESM_latData);
m=MIROC_ESM_lonData(p);
MIROC_ESM_lon=linspace(-m/2,m/2,p);
MIROC_ESM_lon=MIROC_ESM_lon.';
Can_lon1=MIROC_ESM_lon(1:p/2,1);
Can_lon2=MIROC_ESM_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    MIROC_ESM_lat_D(count:count+(p-1),1)=MIROC_ESM_latData(i,1);
    MIROC_ESM_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
MIROC_ESM_mrsosf1=MIROC_ESM_mrsos_f1.';
MIROC_ESM_mrsosf2=MIROC_ESM_mrsos_f2.';
MIROC_ESM_mrsosf3=MIROC_ESM_mrsos_f3.';
MIROC_ESM_m1=[MIROC_ESM_lat_D,MIROC_ESM_lon_D,MIROC_ESM_mrsosf1];
MIROC_ESM_m2=[MIROC_ESM_lat_D,MIROC_ESM_lon_D,MIROC_ESM_mrsosf2];
MIROC_ESM_m3=[MIROC_ESM_lat_D,MIROC_ESM_lon_D,MIROC_ESM_mrsosf3];
[m1,~]=size(MIROC_ESM_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
MIROC_ESM_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MIROC_ESM_m1(i,1)-SMAP1(j,1))+abs(MIROC_ESM_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MIROC_ESM_new_mf1(j,1)=MIROC_ESM_m1(idx,3);
    disp(j)
end
new_mf1=reshape(MIROC_ESM_new_mf1,[964,406]);
MIROC_ESM_mF1=new_mf1.';
[m1,~]=size(MIROC_ESM_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
MIROC_ESM_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MIROC_ESM_m2(i,1)-SMAP2(j,1))+abs(MIROC_ESM_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MIROC_ESM_new_mf2(j,1)=MIROC_ESM_m2(idx,3);
    disp(j)
end
new_mf2=reshape(MIROC_ESM_new_mf2,[964,406]);
MIROC_ESM_mF2=new_mf2.';
[m1,~]=size(MIROC_ESM_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
MIROC_ESM_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MIROC_ESM_m3(i,1)-SMAP3(j,1))+abs(MIROC_ESM_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MIROC_ESM_new_mf3(j,1)=MIROC_ESM_m3(idx,3);
    disp(j)
end
new_mf3=reshape(MIROC_ESM_new_mf3,[964,406]);
MIROC_ESM_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
MIROC_ESM_mF1_new=MIROC_ESM_mF1./temp;
MIROC_ESM_mF2_new=MIROC_ESM_mF2./temp;
MIROC_ESM_mF3_new=MIROC_ESM_mF3./temp;
MIROC_ESM_mF1_new(find(mean_SM<0.1))=nan;
MIROC_ESM_mF2_new(find(mean_SM<0.1))=nan;
MIROC_ESM_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for MIROC_ESM
MIROC_ESM_AVG1_new(406,1)=0;
count=1;
for i=1:406
    MIROC_ESM_AVG1_new(count,1)=nanmean(MIROC_ESM_mF1_new(i,:));
    count=count+1;
end
MIROC_ESM_AVG2_new(406,1)=0;
count=1;
for i=1:406
    MIROC_ESM_AVG2_new(count,1)=nanmean(MIROC_ESM_mF2_new(i,:));
    count=count+1;
end
MIROC_ESM_AVG3_new(406,1)=0;
count=1;
for i=1:406
    MIROC_ESM_AVG3_new(count,1)=nanmean(MIROC_ESM_mF3_new(i,:));
    count=count+1;
end
% MIROC_ESM_CHEM
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(MIROC_ESM_CHEM_lonData);
[q,~]=size(MIROC_ESM_CHEM_latData);
m=MIROC_ESM_CHEM_lonData(p);
MIROC_ESM_CHEM_lon=linspace(-m/2,m/2,p);
MIROC_ESM_CHEM_lon=MIROC_ESM_CHEM_lon.';
Can_lon1=MIROC_ESM_CHEM_lon(1:p/2,1);
Can_lon2=MIROC_ESM_CHEM_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    MIROC_ESM_CHEM_lat_D(count:count+(p-1),1)=MIROC_ESM_CHEM_latData(i,1);
    MIROC_ESM_CHEM_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
MIROC_ESM_CHEM_mrsosf1=MIROC_ESM_CHEM_mrsos_f1.';
MIROC_ESM_CHEM_mrsosf2=MIROC_ESM_CHEM_mrsos_f2.';
MIROC_ESM_CHEM_mrsosf3=MIROC_ESM_CHEM_mrsos_f3.';
MIROC_ESM_CHEM_m1=[MIROC_ESM_CHEM_lat_D,MIROC_ESM_CHEM_lon_D,MIROC_ESM_CHEM_mrsosf1];
MIROC_ESM_CHEM_m2=[MIROC_ESM_CHEM_lat_D,MIROC_ESM_CHEM_lon_D,MIROC_ESM_CHEM_mrsosf2];
MIROC_ESM_CHEM_m3=[MIROC_ESM_CHEM_lat_D,MIROC_ESM_CHEM_lon_D,MIROC_ESM_CHEM_mrsosf3];
[m1,~]=size(MIROC_ESM_CHEM_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
MIROC_ESM_CHEM_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MIROC_ESM_CHEM_m1(i,1)-SMAP1(j,1))+abs(MIROC_ESM_CHEM_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MIROC_ESM_CHEM_new_mf1(j,1)=MIROC_ESM_CHEM_m1(idx,3);
    disp(j)
end
new_mf1=reshape(MIROC_ESM_CHEM_new_mf1,[964,406]);
MIROC_ESM_CHEM_mF1=new_mf1.';
[m1,~]=size(MIROC_ESM_CHEM_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
MIROC_ESM_CHEM_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MIROC_ESM_CHEM_m2(i,1)-SMAP2(j,1))+abs(MIROC_ESM_CHEM_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MIROC_ESM_CHEM_new_mf2(j,1)=MIROC_ESM_CHEM_m2(idx,3);
    disp(j)
end
new_mf2=reshape(MIROC_ESM_CHEM_new_mf2,[964,406]);
MIROC_ESM_CHEM_mF2=new_mf2.';
[m1,~]=size(MIROC_ESM_CHEM_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
MIROC_ESM_CHEM_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MIROC_ESM_CHEM_m3(i,1)-SMAP3(j,1))+abs(MIROC_ESM_CHEM_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MIROC_ESM_CHEM_new_mf3(j,1)=MIROC_ESM_CHEM_m3(idx,3);
    disp(j)
end
new_mf3=reshape(MIROC_ESM_CHEM_new_mf3,[964,406]);
MIROC_ESM_CHEM_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
MIROC_ESM_CHEM_mF1_new=MIROC_ESM_CHEM_mF1./temp;
MIROC_ESM_CHEM_mF2_new=MIROC_ESM_CHEM_mF2./temp;
MIROC_ESM_CHEM_mF3_new=MIROC_ESM_CHEM_mF3./temp;
MIROC_ESM_CHEM_mF1_new(find(mean_SM<0.1))=nan;
MIROC_ESM_CHEM_mF2_new(find(mean_SM<0.1))=nan;
MIROC_ESM_CHEM_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for MIROC_ESM_CHEM
MIROC_ESM_CHEM_AVG1_new(406,1)=0;
count=1;
for i=1:406
    MIROC_ESM_CHEM_AVG1_new(count,1)=nanmean(MIROC_ESM_CHEM_mF1_new(i,:));
    count=count+1;
end
MIROC_ESM_CHEM_AVG2_new(406,1)=0;
count=1;
for i=1:406
    MIROC_ESM_CHEM_AVG2_new(count,1)=nanmean(MIROC_ESM_CHEM_mF2_new(i,:));
    count=count+1;
end
MIROC_ESM_CHEM_AVG3_new(406,1)=0;
count=1;
for i=1:406
    MIROC_ESM_CHEM_AVG3_new(count,1)=nanmean(MIROC_ESM_CHEM_mF3_new(i,:));
    count=count+1;
end
% MRI_CGCM3
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(MRI_CGCM3_lonData);
[q,~]=size(MRI_CGCM3_latData);
m=MRI_CGCM3_lonData(p);
MRI_CGCM3_lon=linspace(-m/2,m/2,p);
MRI_CGCM3_lon=MRI_CGCM3_lon.';
Can_lon1=MRI_CGCM3_lon(1:p/2,1);
Can_lon2=MRI_CGCM3_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    MRI_CGCM3_lat_D(count:count+(p-1),1)=MRI_CGCM3_latData(i,1);
    MRI_CGCM3_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
MRI_CGCM3_mrsosf1=MRI_CGCM3_mrsos_f1.';
MRI_CGCM3_mrsosf2=MRI_CGCM3_mrsos_f2.';
MRI_CGCM3_mrsosf3=MRI_CGCM3_mrsos_f3.';
MRI_CGCM3_m1=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_mrsosf1];
MRI_CGCM3_m2=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_mrsosf2];
MRI_CGCM3_m3=[MRI_CGCM3_lat_D,MRI_CGCM3_lon_D,MRI_CGCM3_mrsosf3];
[m1,~]=size(MRI_CGCM3_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
MRI_CGCM3_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_m1(i,1)-SMAP1(j,1))+abs(MRI_CGCM3_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_mf1(j,1)=MRI_CGCM3_m1(idx,3);
    disp(j)
end
new_mf1=reshape(MRI_CGCM3_new_mf1,[964,406]);
MRI_CGCM3_mF1=new_mf1.';
[m1,~]=size(MRI_CGCM3_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
MRI_CGCM3_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_m2(i,1)-SMAP2(j,1))+abs(MRI_CGCM3_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_mf2(j,1)=MRI_CGCM3_m2(idx,3);
    disp(j)
end
new_mf2=reshape(MRI_CGCM3_new_mf2,[964,406]);
MRI_CGCM3_mF2=new_mf2.';
[m1,~]=size(MRI_CGCM3_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
MRI_CGCM3_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_CGCM3_m3(i,1)-SMAP3(j,1))+abs(MRI_CGCM3_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_CGCM3_new_mf3(j,1)=MRI_CGCM3_m3(idx,3);
    disp(j)
end
new_mf3=reshape(MRI_CGCM3_new_mf3,[964,406]);
MRI_CGCM3_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
MRI_CGCM3_mF1_new=MRI_CGCM3_mF1./temp;
MRI_CGCM3_mF2_new=MRI_CGCM3_mF2./temp;
MRI_CGCM3_mF3_new=MRI_CGCM3_mF3./temp;
MRI_CGCM3_mF1_new(find(mean_SM<0.1))=nan;
MRI_CGCM3_mF2_new(find(mean_SM<0.1))=nan;
MRI_CGCM3_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for MRI_CGCM3
MRI_CGCM3_AVG1_new(406,1)=0;
count=1;
for i=1:406
    MRI_CGCM3_AVG1_new(count,1)=nanmean(MRI_CGCM3_mF1_new(i,:));
    count=count+1;
end
MRI_CGCM3_AVG2_new(406,1)=0;
count=1;
for i=1:406
    MRI_CGCM3_AVG2_new(count,1)=nanmean(MRI_CGCM3_mF2_new(i,:));
    count=count+1;
end
MRI_CGCM3_AVG3_new(406,1)=0;
count=1;
for i=1:406
    MRI_CGCM3_AVG3_new(count,1)=nanmean(MRI_CGCM3_mF3_new(i,:));
    count=count+1;
end
% MRI_ESM1
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(MRI_ESM1_lonData);
[q,~]=size(MRI_ESM1_latData);
m=MRI_ESM1_lonData(p);
MRI_ESM1_lon=linspace(-m/2,m/2,p);
MRI_ESM1_lon=MRI_ESM1_lon.';
Can_lon1=MRI_ESM1_lon(1:p/2,1);
Can_lon2=MRI_ESM1_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    MRI_ESM1_lat_D(count:count+(p-1),1)=MRI_ESM1_latData(i,1);
    MRI_ESM1_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
MRI_ESM1_mrsosf1=MRI_ESM1_mrsos_f1.';
MRI_ESM1_mrsosf2=MRI_ESM1_mrsos_f2.';
MRI_ESM1_mrsosf3=MRI_ESM1_mrsos_f3.';
MRI_ESM1_m1=[MRI_ESM1_lat_D,MRI_ESM1_lon_D,MRI_ESM1_mrsosf1];
MRI_ESM1_m2=[MRI_ESM1_lat_D,MRI_ESM1_lon_D,MRI_ESM1_mrsosf2];
MRI_ESM1_m3=[MRI_ESM1_lat_D,MRI_ESM1_lon_D,MRI_ESM1_mrsosf3];
[m1,~]=size(MRI_ESM1_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
MRI_ESM1_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_ESM1_m1(i,1)-SMAP1(j,1))+abs(MRI_ESM1_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_ESM1_new_mf1(j,1)=MRI_ESM1_m1(idx,3);
    disp(j)
end
new_mf1=reshape(MRI_ESM1_new_mf1,[964,406]);
MRI_ESM1_mF1=new_mf1.';
[m1,~]=size(MRI_ESM1_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
MRI_ESM1_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_ESM1_m2(i,1)-SMAP2(j,1))+abs(MRI_ESM1_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_ESM1_new_mf2(j,1)=MRI_ESM1_m2(idx,3);
    disp(j)
end
new_mf2=reshape(MRI_ESM1_new_mf2,[964,406]);
MRI_ESM1_mF2=new_mf2.';
[m1,~]=size(MRI_ESM1_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
MRI_ESM1_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(MRI_ESM1_m3(i,1)-SMAP3(j,1))+abs(MRI_ESM1_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    MRI_ESM1_new_mf3(j,1)=MRI_ESM1_m3(idx,3);
    disp(j)
end
new_mf3=reshape(MRI_ESM1_new_mf3,[964,406]);
MRI_ESM1_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
MRI_ESM1_mF1_new=MRI_ESM1_mF1./temp;
MRI_ESM1_mF2_new=MRI_ESM1_mF2./temp;
MRI_ESM1_mF3_new=MRI_ESM1_mF3./temp;
MRI_ESM1_mF1_new(find(mean_SM<0.1))=nan;
MRI_ESM1_mF2_new(find(mean_SM<0.1))=nan;
MRI_ESM1_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for MRI_ESM1
MRI_ESM1_AVG1_new(406,1)=0;
count=1;
for i=1:406
    MRI_ESM1_AVG1_new(count,1)=nanmean(MRI_ESM1_mF1_new(i,:));
    count=count+1;
end
MRI_ESM1_AVG2_new(406,1)=0;
count=1;
for i=1:406
    MRI_ESM1_AVG2_new(count,1)=nanmean(MRI_ESM1_mF2_new(i,:));
    count=count+1;
end
MRI_ESM1_AVG3_new(406,1)=0;
count=1;
for i=1:406
    MRI_ESM1_AVG3_new(count,1)=nanmean(MRI_ESM1_mF3_new(i,:));
    count=count+1;
end
% NorESM1_M
% Change the spatial resolution of the model's SSM_n to the same as SMAP
[p,~]=size(NorESM1_M_lonData);
[q,~]=size(NorESM1_M_latData);
m=NorESM1_M_lonData(p);
NorESM1_M_lon=linspace(-m/2,m/2,p);
NorESM1_M_lon=NorESM1_M_lon.';
Can_lon1=NorESM1_M_lon(1:p/2,1);
Can_lon2=NorESM1_M_lon((p/2+1):p,1);
Can_lon=[Can_lon2;Can_lon1];
count=1;
for i=1:q
    NorESM1_M_lat_D(count:count+(p-1),1)=NorESM1_M_latData(i,1);
    NorESM1_M_lon_D(count:count+(p-1),1)=Can_lon;
    count=count+p;
end
NorESM1_M_mrsosf1=NorESM1_M_mrsos_f1.';
NorESM1_M_mrsosf2=NorESM1_M_mrsos_f2.';
NorESM1_M_mrsosf3=NorESM1_M_mrsos_f3.';
NorESM1_M_m1=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_mrsosf1];
NorESM1_M_m2=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_mrsosf2];
NorESM1_M_m3=[NorESM1_M_lat_D,NorESM1_M_lon_D,NorESM1_M_mrsosf3];
[m1,~]=size(NorESM1_M_m1);
[m2,~]=size(SMAP1);
dif(m1,1)=0;
NorESM1_M_new_mf1(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_m1(i,1)-SMAP1(j,1))+abs(NorESM1_M_m1(i,2)-SMAP1(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_mf1(j,1)=NorESM1_M_m1(idx,3);
    disp(j)
end
new_mf1=reshape(NorESM1_M_new_mf1,[964,406]);
NorESM1_M_mF1=new_mf1.';
[m1,~]=size(NorESM1_M_m2);
[m2,~]=size(SMAP2);
dif(m1,1)=0;
NorESM1_M_new_mf2(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_m2(i,1)-SMAP2(j,1))+abs(NorESM1_M_m2(i,2)-SMAP2(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_mf2(j,1)=NorESM1_M_m2(idx,3);
    disp(j)
end
new_mf2=reshape(NorESM1_M_new_mf2,[964,406]);
NorESM1_M_mF2=new_mf2.';
[m1,~]=size(NorESM1_M_m3);
[m2,~]=size(SMAP3);
dif(m1,1)=0;
NorESM1_M_new_mf3(m2,1)=0;
for j=1:m2
    count=1;
    for i=1:m1
        dif(count,1)=abs(NorESM1_M_m3(i,1)-SMAP3(j,1))+abs(NorESM1_M_m3(i,2)-SMAP3(j,2));
        count=count+1;
    end
    [~,idx]=min(dif);
    NorESM1_M_new_mf3(j,1)=NorESM1_M_m3(idx,3);
    disp(j)
end
new_mf3=reshape(NorESM1_M_new_mf3,[964,406]);
NorESM1_M_mF3=new_mf3.';
% Remove N/A regions and regions with SSM less than 0.1
temp=~isnan(SMF1);
temp=single(temp);
temp(temp==0)=nan;
NorESM1_M_mF1_new=NorESM1_M_mF1./temp;
NorESM1_M_mF2_new=NorESM1_M_mF2./temp;
NorESM1_M_mF3_new=NorESM1_M_mF3./temp;
NorESM1_M_mF1_new(find(mean_SM<0.1))=nan;
NorESM1_M_mF2_new(find(mean_SM<0.1))=nan;
NorESM1_M_mF3_new(find(mean_SM<0.1))=nan;
% Latitudinal average of SSM_n for NorESM1_M
NorESM1_M_AVG1_new(406,1)=0;
count=1;
for i=1:406
    NorESM1_M_AVG1_new(count,1)=nanmean(NorESM1_M_mF1_new(i,:));
    count=count+1;
end
NorESM1_M_AVG2_new(406,1)=0;
count=1;
for i=1:406
    NorESM1_M_AVG2_new(count,1)=nanmean(NorESM1_M_mF2_new(i,:));
    count=count+1;
end
NorESM1_M_AVG3_new(406,1)=0;
count=1;
for i=1:406
    NorESM1_M_AVG3_new(count,1)=nanmean(NorESM1_M_mF3_new(i,:));
    count=count+1;
end
% Latitudinal average of SSM_n for models
AVG1_new=(BCC_CSM_AVG1_new+BNU_ESM_AVG1_new+CanESM2_AVG1_new+CNRM_CM5_AVG1_new+CSIRO_Mk_AVG1_new+GFDL_CM3_AVG1_new+GFDL_ESM2G_AVG1_new+GFDL_ESM2M_AVG1_new+HadGEM2_CC_AVG1_new+HadGEM2_ES_AVG1_new+inmcm4_AVG1_new+MIROC5_AVG1_new+MIROC_ESM_AVG1_new+MIROC_ESM_CHEM_AVG1_new+MRI_CGCM3_AVG1_new+MRI_ESM1_AVG1_new+NorESM1_M_AVG1_new)/17;
AVG2_new=(BCC_CSM_AVG2_new+BNU_ESM_AVG2_new+CanESM2_AVG2_new+CNRM_CM5_AVG2_new+CSIRO_Mk_AVG2_new+GFDL_CM3_AVG2_new+GFDL_ESM2G_AVG2_new+GFDL_ESM2M_AVG2_new+HadGEM2_CC_AVG2_new+HadGEM2_ES_AVG2_new+inmcm4_AVG2_new+MIROC5_AVG2_new+MIROC_ESM_AVG2_new+MIROC_ESM_CHEM_AVG2_new+MRI_CGCM3_AVG2_new+MRI_ESM1_AVG2_new+NorESM1_M_AVG2_new)/17;
AVG3_new=(BCC_CSM_AVG3_new+BNU_ESM_AVG3_new+CanESM2_AVG3_new+CNRM_CM5_AVG3_new+CSIRO_Mk_AVG3_new+GFDL_CM3_AVG3_new+GFDL_ESM2G_AVG3_new+GFDL_ESM2M_AVG3_new+HadGEM2_CC_AVG3_new+HadGEM2_ES_AVG3_new+inmcm4_AVG3_new+MIROC5_AVG3_new+MIROC_ESM_AVG3_new+MIROC_ESM_CHEM_AVG3_new+MRI_CGCM3_AVG3_new+MRI_ESM1_AVG3_new+NorESM1_M_AVG3_new)/17;
% Latitudinal standard deviation of SSM_n across all models
count=1;
std_AVG1_new(406,1)=0;
for i=1:406
    std_f1(1,:)=BCC_CSM_AVG1_new(i,:);
    std_f1(2,:)=BNU_ESM_AVG1_new(i,:);
    std_f1(3,:)=CanESM2_AVG1_new(i,:);
    std_f1(4,:)=CNRM_CM5_AVG1_new(i,:);
    std_f1(5,:)=CSIRO_Mk_AVG1_new(i,:);
    std_f1(6,:)=GFDL_CM3_AVG1_new(i,:);
    std_f1(7,:)=GFDL_ESM2G_AVG1_new(i,:);
    std_f1(8,:)=GFDL_ESM2M_AVG1_new(i,:);
    std_f1(9,:)=HadGEM2_CC_AVG1_new(i,:);
    std_f1(10,:)=HadGEM2_ES_AVG1_new(i,:);
    std_f1(11,:)=inmcm4_AVG1_new(i,:);
    std_f1(12,:)=MIROC5_AVG1_new(i,:);
    std_f1(13,:)=MIROC_ESM_AVG1_new(i,:);
    std_f1(14,:)=MIROC_ESM_CHEM_AVG1_new(i,:);
    std_f1(15,:)=MRI_CGCM3_AVG1_new(i,:);
    std_f1(16,:)=MRI_ESM1_AVG1_new(i,:);
    std_f1(17,:)=NorESM1_M_AVG1_new(i,:);
    std_AVG1=std(std_f1);
    std_AVG1_new(count,:)=std_AVG1;
    count=count+1;
end
count=1;
std_AVG2_new(406,1)=0;
for i=1:406
    std_f2(1,:)=BCC_CSM_AVG2_new(i,:);
    std_f2(2,:)=BNU_ESM_AVG2_new(i,:);
    std_f2(3,:)=CanESM2_AVG2_new(i,:);
    std_f2(4,:)=CNRM_CM5_AVG2_new(i,:);
    std_f2(5,:)=CSIRO_Mk_AVG2_new(i,:);
    std_f2(6,:)=GFDL_CM3_AVG2_new(i,:);
    std_f2(7,:)=GFDL_ESM2G_AVG2_new(i,:);
    std_f2(8,:)=GFDL_ESM2M_AVG2_new(i,:);
    std_f2(9,:)=HadGEM2_CC_AVG2_new(i,:);
    std_f2(10,:)=HadGEM2_ES_AVG2_new(i,:);
    std_f2(11,:)=inmcm4_AVG2_new(i,:);
    std_f2(12,:)=MIROC5_AVG2_new(i,:);
    std_f2(13,:)=MIROC_ESM_AVG2_new(i,:);
    std_f2(14,:)=MIROC_ESM_CHEM_AVG2_new(i,:);
    std_f2(15,:)=MRI_CGCM3_AVG2_new(i,:);
    std_f2(16,:)=MRI_ESM1_AVG2_new(i,:);
    std_f2(17,:)=NorESM1_M_AVG2_new(i,:);
    std_AVG2=std(std_f2);
    std_AVG2_new(count,:)=std_AVG2;
    count=count+1;
end
count=1;
std_AVG3_new(406,1)=0;
for i=1:406
    std_f3(1,:)=BCC_CSM_AVG3_new(i,:);
    std_f3(2,:)=BNU_ESM_AVG3_new(i,:);
    std_f3(3,:)=CanESM2_AVG3_new(i,:);
    std_f3(4,:)=CNRM_CM5_AVG3_new(i,:);
    std_f3(5,:)=CSIRO_Mk_AVG3_new(i,:);
    std_f3(6,:)=GFDL_CM3_AVG3_new(i,:);
    std_f3(7,:)=GFDL_ESM2G_AVG3_new(i,:);
    std_f3(8,:)=GFDL_ESM2M_AVG3_new(i,:);
    std_f3(9,:)=HadGEM2_CC_AVG3_new(i,:);
    std_f3(10,:)=HadGEM2_ES_AVG3_new(i,:);
    std_f3(11,:)=inmcm4_AVG3_new(i,:);
    std_f3(12,:)=MIROC5_AVG3_new(i,:);
    std_f3(13,:)=MIROC_ESM_AVG3_new(i,:);
    std_f3(14,:)=MIROC_ESM_CHEM_AVG3_new(i,:);
    std_f3(15,:)=MRI_CGCM3_AVG3_new(i,:);
    std_f3(16,:)=MRI_ESM1_AVG3_new(i,:);
    std_f3(17,:)=NorESM1_M_AVG3_new(i,:);
    std_AVG3=std(std_f3);
    std_AVG3_new(count,:)=std_AVG3;
    count=count+1;
end
% Mean of renormalized global mean SSM for each latitude
mean_SM(find(mean_SM<0.1))=nan;
Flattened_mean_SM=mean_SM(:)';
MappedF_mean_SM=mapminmax(Flattened_mean_SM,0,1);
Mapped_mean_SM=reshape(MappedF_mean_SM,size(mean_SM));
[m,~]=size(Mapped_mean_SM);
Mapped_SMAP_Avg(m,1)=0;
count=1;
for i=1:m
    Mapped_SMAP_Avg(count,1)=nanmean(Mapped_mean_SM(i,:));
    count=count+1;
end
% Plot the diagram (Figure 6)
subplot(4,1,1)
plot(latSMAP,SMAP_Avg1,'Color',[0 0 0],'LineWidth',1)
hold on
shadedErrorBar(latSMAP,AVG1_new,std_AVG1_new,'lineprops',{'-r','MarkerFaceColor','r'});
hold off
grid on
xlim([-90,90])
ylim([0,1])
set(gca,'XTick',[-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])
set(gca,'YTick',[0,0.25,0.5,0.75,1])
xlabel('(a)')
title('1/30<f<1/7 day^{-1}')
legend('Observation (SMAP)','Model (CMIP5)')
subplot(4,1,2)
plot(latSMAP,SMAP_Avg2,'Color',[0 0 0],'LineWidth',1)
hold on
shadedErrorBar(latSMAP,AVG2_new,std_AVG2_new,'lineprops',{'-r','MarkerFaceColor','r'});
hold off
grid on
xlim([-90,90])
ylim([0,1])
set(gca,'XTick',[-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])
set(gca,'YTick',[0,0.25,0.5,0.75,1])
xlabel('(b)')
title('1/90<f<1/30 day^{-1}')
subplot(4,1,3)
plot(latSMAP,SMAP_Avg3,'Color',[0 0 0],'LineWidth',1)
hold on
shadedErrorBar(latSMAP,AVG3_new,std_AVG3_new,'lineprops',{'-r','MarkerFaceColor','r'});
hold off
grid on
xlim([-90,90])
ylim([0,1])
set(gca,'XTick',[-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])
set(gca,'YTick',[0,0.25,0.5,0.75,1])
xlabel('(c)')
title('1/365<f<1/90 day^{-1}')
h=text(-0.05,0.5,'Latitude-average SSM_n');
set(h,'fontsize',11,'rotation',90,'HorizontalAlignment','center')
subplot(4,1,4)
plot(latSMAP,Mapped_SMAP_Avg,'Color',[0 0 0],'LineWidth',1)
grid on
xlim([-90,90])
ylim([0,1])
set(gca,'XTick',[-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])
set(gca,'YTick',[0,0.25,0.5,0.75,1])
xlabel('(d)')
ylabel('Normalized mean SSM')
hAxis=axes('visible','off');
h=text(-0.05,0.5,'Latitude');
set(h,'fontsize',11,'VerticalAlignment','middle')