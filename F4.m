% % Model deficiencies of SSM_n over different time scales
% Averaged biases for SSM_n across all models
avg_mdif1_try1=(BCC_CSM_mdif1_try1+BNU_ESM_mdif1_try1+CanESM2_mdif1_try1+CNRM_CM5_mdif1_try1+CSIRO_Mk_mdif1_try1+GFDL_CM3_mdif1_try1+GFDL_ESM2G_mdif1_try1+GFDL_ESM2M_mdif1_try1+HadGEM2_CC_mdif1_try1+HadGEM2_ES_mdif1_try1+inmcm4_mdif1_try1+MIROC5_mdif1_try1+MIROC_ESM_mdif1_try1+MIROC_ESM_CHEM_mdif1_try1+MRI_CGCM3_mdif1_try1+MRI_ESM1_mdif1_try1+NorESM1_M_mdif1_try1)/17;
avg_mdif2_try1=(BCC_CSM_mdif2_try1+BNU_ESM_mdif2_try1+CanESM2_mdif2_try1+CNRM_CM5_mdif2_try1+CSIRO_Mk_mdif2_try1+GFDL_CM3_mdif2_try1+GFDL_ESM2G_mdif2_try1+GFDL_ESM2M_mdif2_try1+HadGEM2_CC_mdif2_try1+HadGEM2_ES_mdif2_try1+inmcm4_mdif2_try1+MIROC5_mdif2_try1+MIROC_ESM_mdif2_try1+MIROC_ESM_CHEM_mdif2_try1+MRI_CGCM3_mdif2_try1+MRI_ESM1_mdif2_try1+NorESM1_M_mdif2_try1)/17;
avg_mdif3_try1=(BCC_CSM_mdif3_try1+BNU_ESM_mdif3_try1+CanESM2_mdif3_try1+CNRM_CM5_mdif3_try1+CSIRO_Mk_mdif3_try1+GFDL_CM3_mdif3_try1+GFDL_ESM2G_mdif3_try1+GFDL_ESM2M_mdif3_try1+HadGEM2_CC_mdif3_try1+HadGEM2_ES_mdif3_try1+inmcm4_mdif3_try1+MIROC5_mdif3_try1+MIROC_ESM_mdif3_try1+MIROC_ESM_CHEM_mdif3_try1+MRI_CGCM3_mdif3_try1+MRI_ESM1_mdif3_try1+NorESM1_M_mdif3_try1)/17;
% Prepare for significance test (1st frequency range)
BCC_CSM_mdif1_s=BCC_CSM_mdif1_try1-avg_mdif1_try1;
BNU_ESM_mdif1_s=BNU_ESM_mdif1_try1-avg_mdif1_try1;
CanESM2_mdif1_s=CanESM2_mdif1_try1-avg_mdif1_try1;
CNRM_CM5_mdif1_s=CNRM_CM5_mdif1_try1-avg_mdif1_try1;
CSIRO_Mk_mdif1_s=CSIRO_Mk_mdif1_try1-avg_mdif1_try1;
GFDL_CM3_mdif1_s=GFDL_CM3_mdif1_try1-avg_mdif1_try1;
GFDL_ESM2G_mdif1_s=GFDL_ESM2G_mdif1_try1-avg_mdif1_try1;
GFDL_ESM2M_mdif1_s=GFDL_ESM2M_mdif1_try1-avg_mdif1_try1;
HadGEM2_CC_mdif1_s=HadGEM2_CC_mdif1_try1-avg_mdif1_try1;
HadGEM2_ES_mdif1_s=HadGEM2_ES_mdif1_try1-avg_mdif1_try1;
inmcm4_mdif1_s=inmcm4_mdif1_try1-avg_mdif1_try1;
MIROC5_mdif1_s=MIROC5_mdif1_try1-avg_mdif1_try1;
MIROC_ESM_mdif1_s=MIROC_ESM_mdif1_try1-avg_mdif1_try1;
MIROC_ESM_CHEM_mdif1_s=MIROC_ESM_CHEM_mdif1_try1-avg_mdif1_try1;
MRI_CGCM3_mdif1_s=MRI_CGCM3_mdif1_try1-avg_mdif1_try1;
MRI_ESM1_mdif1_s=MRI_ESM1_mdif1_try1-avg_mdif1_try1;
NorESM1_M_mdif1_s=NorESM1_M_mdif1_try1-avg_mdif1_try1;
BCC_CSM_mdif1_try1(BCC_CSM_mdif1_try1<0)=0;
BCC_CSM_mdif1_try1(BCC_CSM_mdif1_try1>0)=1;
BNU_ESM_mdif1_try1(BNU_ESM_mdif1_try1<0)=0;
BNU_ESM_mdif1_try1(BNU_ESM_mdif1_try1>0)=1;
CanESM2_mdif1_try1(CanESM2_mdif1_try1<0)=0;
CanESM2_mdif1_try1(CanESM2_mdif1_try1>0)=1;
CNRM_CM5_mdif1_try1(CNRM_CM5_mdif1_try1<0)=0;
CNRM_CM5_mdif1_try1(CNRM_CM5_mdif1_try1>0)=1;
CSIRO_Mk_mdif1_try1(CSIRO_Mk_mdif1_try1<0)=0;
CSIRO_Mk_mdif1_try1(CSIRO_Mk_mdif1_try1>0)=1;
GFDL_CM3_mdif1_try1(GFDL_CM3_mdif1_try1<0)=0;
GFDL_CM3_mdif1_try1(GFDL_CM3_mdif1_try1>0)=1;
GFDL_ESM2G_mdif1_try1(GFDL_ESM2G_mdif1_try1<0)=0;
GFDL_ESM2G_mdif1_try1(GFDL_ESM2G_mdif1_try1>0)=1;
GFDL_ESM2M_mdif1_try1(GFDL_ESM2M_mdif1_try1<0)=0;
GFDL_ESM2M_mdif1_try1(GFDL_ESM2M_mdif1_try1>0)=1;
HadGEM2_CC_mdif1_try1(HadGEM2_CC_mdif1_try1<0)=0;
HadGEM2_CC_mdif1_try1(HadGEM2_CC_mdif1_try1>0)=1;
HadGEM2_ES_mdif1_try1(HadGEM2_ES_mdif1_try1<0)=0;
HadGEM2_ES_mdif1_try1(HadGEM2_ES_mdif1_try1>0)=1;
inmcm4_mdif1_try1(inmcm4_mdif1_try1<0)=0;
inmcm4_mdif1_try1(inmcm4_mdif1_try1>0)=1;
MIROC5_mdif1_try1(MIROC5_mdif1_try1<0)=0;
MIROC5_mdif1_try1(MIROC5_mdif1_try1>0)=1;
MIROC_ESM_mdif1_try1(MIROC_ESM_mdif1_try1<0)=0;
MIROC_ESM_mdif1_try1(MIROC_ESM_mdif1_try1>0)=1;
MIROC_ESM_CHEM_mdif1_try1(MIROC_ESM_CHEM_mdif1_try1<0)=0;
MIROC_ESM_CHEM_mdif1_try1(MIROC_ESM_CHEM_mdif1_try1>0)=1;
MRI_CGCM3_mdif1_try1(MRI_CGCM3_mdif1_try1<0)=0;
MRI_CGCM3_mdif1_try1(MRI_CGCM3_mdif1_try1>0)=1;
MRI_ESM1_mdif1_try1(MRI_ESM1_mdif1_try1<0)=0;
MRI_ESM1_mdif1_try1(MRI_ESM1_mdif1_try1>0)=1;
NorESM1_M_mdif1_try1(NorESM1_M_mdif1_try1<0)=0;
NorESM1_M_mdif1_try1(NorESM1_M_mdif1_try1>0)=1;
[m,n]=size(BCC_CSM_mdif1_try1);
BCC_CSM_mdif1_r=reshape(BCC_CSM_mdif1_try1,[m*n,1]);
BNU_ESM_mdif1_r=reshape(BNU_ESM_mdif1_try1,[m*n,1]);
CanESM2_mdif1_r=reshape(CanESM2_mdif1_try1,[m*n,1]);
CNRM_CM5_mdif1_r=reshape(CNRM_CM5_mdif1_try1,[m*n,1]);
CSIRO_Mk_mdif1_r=reshape(CSIRO_Mk_mdif1_try1,[m*n,1]);
GFDL_CM3_mdif1_r=reshape(GFDL_CM3_mdif1_try1,[m*n,1]);
GFDL_ESM2G_mdif1_r=reshape(GFDL_ESM2G_mdif1_try1,[m*n,1]);
GFDL_ESM2M_mdif1_r=reshape(GFDL_ESM2M_mdif1_try1,[m*n,1]);
HadGEM2_CC_mdif1_r=reshape(HadGEM2_CC_mdif1_try1,[m*n,1]);
HadGEM2_ES_mdif1_r=reshape(HadGEM2_ES_mdif1_try1,[m*n,1]);
inmcm4_mdif1_r=reshape(inmcm4_mdif1_try1,[m*n,1]);
MIROC5_mdif1_r=reshape(MIROC5_mdif1_try1,[m*n,1]);
MIROC_ESM_mdif1_r=reshape(MIROC_ESM_mdif1_try1,[m*n,1]);
MIROC_ESM_CHEM_mdif1_r=reshape(MIROC_ESM_CHEM_mdif1_try1,[m*n,1]);
MRI_CGCM3_mdif1_r=reshape(MRI_CGCM3_mdif1_try1,[m*n,1]);
MRI_ESM1_mdif1_r=reshape(MRI_ESM1_mdif1_try1,[m*n,1]);
NorESM1_M_mdif1_r=reshape(NorESM1_M_mdif1_try1,[m*n,1]);
mdif1_r=[BCC_CSM_mdif1_r,BNU_ESM_mdif1_r,CanESM2_mdif1_r,CNRM_CM5_mdif1_r,CSIRO_Mk_mdif1_r,GFDL_CM3_mdif1_r,GFDL_ESM2G_mdif1_r,GFDL_ESM2M_mdif1_r,HadGEM2_CC_mdif1_r,HadGEM2_ES_mdif1_r,inmcm4_mdif1_r,MIROC5_mdif1_r,MIROC_ESM_mdif1_r,MIROC_ESM_CHEM_mdif1_r,MRI_CGCM3_mdif1_r,MRI_ESM1_mdif1_r,NorESM1_M_mdif1_r];
% The number of models with the same sign as averaged biaes in each pixel (1st frequency range)
mdif1_r_nan=isnan(mdif1_r);
mdif1_r_nan=single(mdif1_r_nan);
mdif1_r_nan_sum=sum(mdif1_r_nan,2);
mdif1_r_nan_sum_c=[mdif1_r,mdif1_r_nan_sum];
[row,~]=size(mdif1_r);
mdif1_count_0(row,1)=0;
mdif1_count_1(row,1)=0;
count=1;
for i=1:row
    if mdif1_r_nan_sum_c(i,18)==17
        num_0=nan;
        num_1=nan;
    else
        num_0=sum(mdif1_r(i,:)==0,2);
        num_1=sum(mdif1_r(i,:)~=0,2);
    end
    mdif1_count_0(count,:)=num_0;
    mdif1_count_1(count,:)=num_1;
    count=count+1;
    disp(i)
end
% Ratio of the number of models with the same sign as averaged biaes in each pixel to the number of total models (1st frequency range)
mdif1_count_0_real=mdif1_count_0/17;
mdif1_count_1_real=mdif1_count_1/17;
mdif1_positive_pct=reshape(mdif1_count_0_real,[m,n]);
mdif1_negative_pct=reshape(mdif1_count_1_real,[m,n]);
% Prepare for significance test (2nd frequency range)
BCC_CSM_mdif2_try1(BCC_CSM_mdif2_try1<0)=0;
BCC_CSM_mdif2_try1(BCC_CSM_mdif2_try1>0)=1;
BNU_ESM_mdif2_try1(BNU_ESM_mdif2_try1<0)=0;
BNU_ESM_mdif2_try1(BNU_ESM_mdif2_try1>0)=1;
CanESM2_mdif2_try1(CanESM2_mdif2_try1<0)=0;
CanESM2_mdif2_try1(CanESM2_mdif2_try1>0)=1;
CNRM_CM5_mdif2_try1(CNRM_CM5_mdif2_try1<0)=0;
CNRM_CM5_mdif2_try1(CNRM_CM5_mdif2_try1>0)=1;
CSIRO_Mk_mdif2_try1(CSIRO_Mk_mdif2_try1<0)=0;
CSIRO_Mk_mdif2_try1(CSIRO_Mk_mdif2_try1>0)=1;
GFDL_CM3_mdif2_try1(GFDL_CM3_mdif2_try1<0)=0;
GFDL_CM3_mdif2_try1(GFDL_CM3_mdif2_try1>0)=1;
GFDL_ESM2G_mdif2_try1(GFDL_ESM2G_mdif2_try1<0)=0;
GFDL_ESM2G_mdif2_try1(GFDL_ESM2G_mdif2_try1>0)=1;
GFDL_ESM2M_mdif2_try1(GFDL_ESM2M_mdif2_try1<0)=0;
GFDL_ESM2M_mdif2_try1(GFDL_ESM2M_mdif2_try1>0)=1;
HadGEM2_CC_mdif2_try1(HadGEM2_CC_mdif2_try1<0)=0;
HadGEM2_CC_mdif2_try1(HadGEM2_CC_mdif2_try1>0)=1;
HadGEM2_ES_mdif2_try1(HadGEM2_ES_mdif2_try1<0)=0;
HadGEM2_ES_mdif2_try1(HadGEM2_ES_mdif2_try1>0)=1;
inmcm4_mdif2_try1(inmcm4_mdif2_try1<0)=0;
inmcm4_mdif2_try1(inmcm4_mdif2_try1>0)=1;
MIROC5_mdif2_try1(MIROC5_mdif2_try1<0)=0;
MIROC5_mdif2_try1(MIROC5_mdif2_try1>0)=1;
MIROC_ESM_mdif2_try1(MIROC_ESM_mdif2_try1<0)=0;
MIROC_ESM_mdif2_try1(MIROC_ESM_mdif2_try1>0)=1;
MIROC_ESM_CHEM_mdif2_try1(MIROC_ESM_CHEM_mdif2_try1<0)=0;
MIROC_ESM_CHEM_mdif2_try1(MIROC_ESM_CHEM_mdif2_try1>0)=1;
MRI_CGCM3_mdif2_try1(MRI_CGCM3_mdif2_try1<0)=0;
MRI_CGCM3_mdif2_try1(MRI_CGCM3_mdif2_try1>0)=1;
MRI_ESM1_mdif2_try1(MRI_ESM1_mdif2_try1<0)=0;
MRI_ESM1_mdif2_try1(MRI_ESM1_mdif2_try1>0)=1;
NorESM1_M_mdif2_try1(NorESM1_M_mdif2_try1<0)=0;
NorESM1_M_mdif2_try1(NorESM1_M_mdif2_try1>0)=1;
[m,n]=size(BCC_CSM_mdif2_try1);
BCC_CSM_mdif2_r=reshape(BCC_CSM_mdif2_try1,[m*n,1]);
BNU_ESM_mdif2_r=reshape(BNU_ESM_mdif2_try1,[m*n,1]);
CanESM2_mdif2_r=reshape(CanESM2_mdif2_try1,[m*n,1]);
CNRM_CM5_mdif2_r=reshape(CNRM_CM5_mdif2_try1,[m*n,1]);
CSIRO_Mk_mdif2_r=reshape(CSIRO_Mk_mdif2_try1,[m*n,1]);
GFDL_CM3_mdif2_r=reshape(GFDL_CM3_mdif2_try1,[m*n,1]);
GFDL_ESM2G_mdif2_r=reshape(GFDL_ESM2G_mdif2_try1,[m*n,1]);
GFDL_ESM2M_mdif2_r=reshape(GFDL_ESM2M_mdif2_try1,[m*n,1]);
HadGEM2_CC_mdif2_r=reshape(HadGEM2_CC_mdif2_try1,[m*n,1]);
HadGEM2_ES_mdif2_r=reshape(HadGEM2_ES_mdif2_try1,[m*n,1]);
inmcm4_mdif2_r=reshape(inmcm4_mdif2_try1,[m*n,1]);
MIROC5_mdif2_r=reshape(MIROC5_mdif2_try1,[m*n,1]);
MIROC_ESM_mdif2_r=reshape(MIROC_ESM_mdif2_try1,[m*n,1]);
MIROC_ESM_CHEM_mdif2_r=reshape(MIROC_ESM_CHEM_mdif2_try1,[m*n,1]);
MRI_CGCM3_mdif2_r=reshape(MRI_CGCM3_mdif2_try1,[m*n,1]);
MRI_ESM1_mdif2_r=reshape(MRI_ESM1_mdif2_try1,[m*n,1]);
NorESM1_M_mdif2_r=reshape(NorESM1_M_mdif2_try1,[m*n,1]);
mdif2_r=[BCC_CSM_mdif2_r,BNU_ESM_mdif2_r,CanESM2_mdif2_r,CNRM_CM5_mdif2_r,CSIRO_Mk_mdif2_r,GFDL_CM3_mdif2_r,GFDL_ESM2G_mdif2_r,GFDL_ESM2M_mdif2_r,HadGEM2_CC_mdif2_r,HadGEM2_ES_mdif2_r,inmcm4_mdif2_r,MIROC5_mdif2_r,MIROC_ESM_mdif2_r,MIROC_ESM_CHEM_mdif2_r,MRI_CGCM3_mdif2_r,MRI_ESM1_mdif2_r,NorESM1_M_mdif2_r];
% The number of models with the same sign as averaged biaes in each pixel (2nd frequency range)
mdif2_r_nan=isnan(mdif2_r);
mdif2_r_nan=single(mdif2_r_nan);
mdif2_r_nan_sum=sum(mdif2_r_nan,2);
mdif2_r_nan_sum_c=[mdif2_r,mdif2_r_nan_sum];
[row,~]=size(mdif2_r);
mdif2_count_0(row,1)=0;
mdif2_count_1(row,1)=0;
count=1;
for i=1:row
    if mdif2_r_nan_sum_c(i,18)==17
        num_0=nan;
        num_1=nan;
    else
        num_0=sum(mdif2_r(i,:)==0,2);
        num_1=sum(mdif2_r(i,:)~=0,2);
    end
    mdif2_count_0(count,:)=num_0;
    mdif2_count_1(count,:)=num_1;
    count=count+1;
    disp(i)
end
% Ratio of the number of models with the same sign as averaged biaes in each pixel to the number of total models (2nd frequency range)
mdif2_count_0_real=mdif2_count_0/17;
mdif2_count_1_real=mdif2_count_1/17;
mdif2_positive_pct=reshape(mdif2_count_0_real,[m,n]);
mdif2_negative_pct=reshape(mdif2_count_1_real,[m,n]);
% Prepare for significance test (3rd frequency range)
BCC_CSM_mdif3_try1(BCC_CSM_mdif3_try1<0)=0;
BCC_CSM_mdif3_try1(BCC_CSM_mdif3_try1>0)=1;
BNU_ESM_mdif3_try1(BNU_ESM_mdif3_try1<0)=0;
BNU_ESM_mdif3_try1(BNU_ESM_mdif3_try1>0)=1;
CanESM2_mdif3_try1(CanESM2_mdif3_try1<0)=0;
CanESM2_mdif3_try1(CanESM2_mdif3_try1>0)=1;
CNRM_CM5_mdif3_try1(CNRM_CM5_mdif3_try1<0)=0;
CNRM_CM5_mdif3_try1(CNRM_CM5_mdif3_try1>0)=1;
CSIRO_Mk_mdif3_try1(CSIRO_Mk_mdif3_try1<0)=0;
CSIRO_Mk_mdif3_try1(CSIRO_Mk_mdif3_try1>0)=1;
GFDL_CM3_mdif3_try1(GFDL_CM3_mdif3_try1<0)=0;
GFDL_CM3_mdif3_try1(GFDL_CM3_mdif3_try1>0)=1;
GFDL_ESM2G_mdif3_try1(GFDL_ESM2G_mdif3_try1<0)=0;
GFDL_ESM2G_mdif3_try1(GFDL_ESM2G_mdif3_try1>0)=1;
GFDL_ESM2M_mdif3_try1(GFDL_ESM2M_mdif3_try1<0)=0;
GFDL_ESM2M_mdif3_try1(GFDL_ESM2M_mdif3_try1>0)=1;
HadGEM2_CC_mdif3_try1(HadGEM2_CC_mdif3_try1<0)=0;
HadGEM2_CC_mdif3_try1(HadGEM2_CC_mdif3_try1>0)=1;
HadGEM2_ES_mdif3_try1(HadGEM2_ES_mdif3_try1<0)=0;
HadGEM2_ES_mdif3_try1(HadGEM2_ES_mdif3_try1>0)=1;
inmcm4_mdif3_try1(inmcm4_mdif3_try1<0)=0;
inmcm4_mdif3_try1(inmcm4_mdif3_try1>0)=1;
MIROC5_mdif3_try1(MIROC5_mdif3_try1<0)=0;
MIROC5_mdif3_try1(MIROC5_mdif3_try1>0)=1;
MIROC_ESM_mdif3_try1(MIROC_ESM_mdif3_try1<0)=0;
MIROC_ESM_mdif3_try1(MIROC_ESM_mdif3_try1>0)=1;
MIROC_ESM_CHEM_mdif3_try1(MIROC_ESM_CHEM_mdif3_try1<0)=0;
MIROC_ESM_CHEM_mdif3_try1(MIROC_ESM_CHEM_mdif3_try1>0)=1;
MRI_CGCM3_mdif3_try1(MRI_CGCM3_mdif3_try1<0)=0;
MRI_CGCM3_mdif3_try1(MRI_CGCM3_mdif3_try1>0)=1;
MRI_ESM1_mdif3_try1(MRI_ESM1_mdif3_try1<0)=0;
MRI_ESM1_mdif3_try1(MRI_ESM1_mdif3_try1>0)=1;
NorESM1_M_mdif3_try1(NorESM1_M_mdif3_try1<0)=0;
NorESM1_M_mdif3_try1(NorESM1_M_mdif3_try1>0)=1;
[m,n]=size(BCC_CSM_mdif3_try1);
BCC_CSM_mdif3_r=reshape(BCC_CSM_mdif3_try1,[m*n,1]);
BNU_ESM_mdif3_r=reshape(BNU_ESM_mdif3_try1,[m*n,1]);
CanESM2_mdif3_r=reshape(CanESM2_mdif3_try1,[m*n,1]);
CNRM_CM5_mdif3_r=reshape(CNRM_CM5_mdif3_try1,[m*n,1]);
CSIRO_Mk_mdif3_r=reshape(CSIRO_Mk_mdif3_try1,[m*n,1]);
GFDL_CM3_mdif3_r=reshape(GFDL_CM3_mdif3_try1,[m*n,1]);
GFDL_ESM2G_mdif3_r=reshape(GFDL_ESM2G_mdif3_try1,[m*n,1]);
GFDL_ESM2M_mdif3_r=reshape(GFDL_ESM2M_mdif3_try1,[m*n,1]);
HadGEM2_CC_mdif3_r=reshape(HadGEM2_CC_mdif3_try1,[m*n,1]);
HadGEM2_ES_mdif3_r=reshape(HadGEM2_ES_mdif3_try1,[m*n,1]);
inmcm4_mdif3_r=reshape(inmcm4_mdif3_try1,[m*n,1]);
MIROC5_mdif3_r=reshape(MIROC5_mdif3_try1,[m*n,1]);
MIROC_ESM_mdif3_r=reshape(MIROC_ESM_mdif3_try1,[m*n,1]);
MIROC_ESM_CHEM_mdif3_r=reshape(MIROC_ESM_CHEM_mdif3_try1,[m*n,1]);
MRI_CGCM3_mdif3_r=reshape(MRI_CGCM3_mdif3_try1,[m*n,1]);
MRI_ESM1_mdif3_r=reshape(MRI_ESM1_mdif3_try1,[m*n,1]);
NorESM1_M_mdif3_r=reshape(NorESM1_M_mdif3_try1,[m*n,1]);
mdif3_r=[BCC_CSM_mdif3_r,BNU_ESM_mdif3_r,CanESM2_mdif3_r,CNRM_CM5_mdif3_r,CSIRO_Mk_mdif3_r,GFDL_CM3_mdif3_r,GFDL_ESM2G_mdif3_r,GFDL_ESM2M_mdif3_r,HadGEM2_CC_mdif3_r,HadGEM2_ES_mdif3_r,inmcm4_mdif3_r,MIROC5_mdif3_r,MIROC_ESM_mdif3_r,MIROC_ESM_CHEM_mdif3_r,MRI_CGCM3_mdif3_r,MRI_ESM1_mdif3_r,NorESM1_M_mdif3_r];
% The number of models with the same sign as averaged biaes in each pixel (3rd frequency range)
mdif3_r_nan=isnan(mdif3_r);
mdif3_r_nan=single(mdif3_r_nan);
mdif3_r_nan_sum=sum(mdif3_r_nan,2);
mdif3_r_nan_sum_c=[mdif3_r,mdif3_r_nan_sum];
[row,~]=size(mdif3_r);
mdif3_count_0(row,1)=0;
mdif3_count_1(row,1)=0;
count=1;
for i=1:row
    if mdif3_r_nan_sum_c(i,18)==17
        num_0=nan;
        num_1=nan;
    else
        num_0=sum(mdif3_r(i,:)==0,2);
        num_1=sum(mdif3_r(i,:)~=0,2);
    end
    mdif3_count_0(count,:)=num_0;
    mdif3_count_1(count,:)=num_1;
    count=count+1;
    disp(i)
end
% Ratio of the number of models with the same sign as averaged biaes in each pixel to the number of total models (3rd frequency range)
mdif3_count_0_real=mdif3_count_0/17;
mdif3_count_1_real=mdif3_count_1/17;
mdif3_positive_pct=reshape(mdif3_count_0_real,[m,n]);
mdif3_negative_pct=reshape(mdif3_count_1_real,[m,n]);
% Standard deviation of SSM_n across all models
count=1;
avg_mrsosf1_new(406,964)=0;
std_mrsosf1_new(406,964)=0;
for i=1:406
    std_f1(1,:)=BCC_CSM_mrsosf1_new(i,:);
    std_f1(2,:)=BNU_ESM_mrsosf1_new(i,:);
    std_f1(3,:)=CanESM2_mrsosf1_new(i,:);
    std_f1(4,:)=CNRM_CM5_mrsosf1_new(i,:);
    std_f1(5,:)=CSIRO_Mk_mrsosf1_new(i,:);
    std_f1(6,:)=GFDL_CM3_mrsosf1_new(i,:);
    std_f1(7,:)=GFDL_ESM2G_mrsosf1_new(i,:);
    std_f1(8,:)=GFDL_ESM2M_mrsosf1_new(i,:);
    std_f1(9,:)=HadGEM2_CC_mrsosf1_new(i,:);
    std_f1(10,:)=HadGEM2_ES_mrsosf1_new(i,:);
    std_f1(11,:)=inmcm4_mrsosf1_new(i,:);
    std_f1(12,:)=MIROC5_mrsosf1_new(i,:);
    std_f1(13,:)=MIROC_ESM_mrsosf1_new(i,:);
    std_f1(14,:)=MIROC_ESM_CHEM_mrsosf1_new(i,:);
    std_f1(15,:)=MRI_CGCM3_mrsosf1_new(i,:);
    std_f1(16,:)=MRI_ESM1_mrsosf1_new(i,:);
    std_f1(17,:)=NorESM1_M_mrsosf1_new(i,:);
    avg_mdif_f1=mean(std_f1);
    avg_mrsosf1_new(count,:)=avg_mdif_f1;
    std_mdif_f1=std(std_f1);
    std_mrsosf1_new(count,:)=std_mdif_f1;
    count=count+1;
end
count=1;
avg_mrsosf2_new(406,964)=0;
std_mrsosf2_new(406,964)=0;
for i=1:406
    std_f2(1,:)=BCC_CSM_mrsosf2_new(i,:);
    std_f2(2,:)=BNU_ESM_mrsosf2_new(i,:);
    std_f2(3,:)=CanESM2_mrsosf2_new(i,:);
    std_f2(4,:)=CNRM_CM5_mrsosf2_new(i,:);
    std_f2(5,:)=CSIRO_Mk_mrsosf2_new(i,:);
    std_f2(6,:)=GFDL_CM3_mrsosf2_new(i,:);
    std_f2(7,:)=GFDL_ESM2G_mrsosf2_new(i,:);
    std_f2(8,:)=GFDL_ESM2M_mrsosf2_new(i,:);
    std_f2(9,:)=HadGEM2_CC_mrsosf2_new(i,:);
    std_f2(10,:)=HadGEM2_ES_mrsosf2_new(i,:);
    std_f2(11,:)=inmcm4_mrsosf2_new(i,:);
    std_f2(12,:)=MIROC5_mrsosf2_new(i,:);
    std_f2(13,:)=MIROC_ESM_mrsosf2_new(i,:);
    std_f2(14,:)=MIROC_ESM_CHEM_mrsosf2_new(i,:);
    std_f2(15,:)=MRI_CGCM3_mrsosf2_new(i,:);
    std_f2(16,:)=MRI_ESM1_mrsosf2_new(i,:);
    std_f2(17,:)=NorESM1_M_mrsosf2_new(i,:);
    avg_mdif_f2=mean(std_f2);
    avg_mrsosf2_new(count,:)=avg_mdif_f2;
    std_mdif_f2=std(std_f2);
    std_mrsosf2_new(count,:)=std_mdif_f2;
    count=count+1;
end
count=1;
avg_mrsosf3_new(406,964)=0;
std_mrsosf3_new(406,964)=0;
for i=1:406
    std_f3(1,:)=BCC_CSM_mrsosf3_new(i,:);
    std_f3(2,:)=BNU_ESM_mrsosf3_new(i,:);
    std_f3(3,:)=CanESM2_mrsosf3_new(i,:);
    std_f3(4,:)=CNRM_CM5_mrsosf3_new(i,:);
    std_f3(5,:)=CSIRO_Mk_mrsosf3_new(i,:);
    std_f3(6,:)=GFDL_CM3_mrsosf3_new(i,:);
    std_f3(7,:)=GFDL_ESM2G_mrsosf3_new(i,:);
    std_f3(8,:)=GFDL_ESM2M_mrsosf3_new(i,:);
    std_f3(9,:)=HadGEM2_CC_mrsosf3_new(i,:);
    std_f3(10,:)=HadGEM2_ES_mrsosf3_new(i,:);
    std_f3(11,:)=inmcm4_mrsosf3_new(i,:);
    std_f3(12,:)=MIROC5_mrsosf3_new(i,:);
    std_f3(13,:)=MIROC_ESM_mrsosf3_new(i,:);
    std_f3(14,:)=MIROC_ESM_CHEM_mrsosf3_new(i,:);
    std_f3(15,:)=MRI_CGCM3_mrsosf3_new(i,:);
    std_f3(16,:)=MRI_ESM1_mrsosf3_new(i,:);
    std_f3(17,:)=NorESM1_M_mrsosf3_new(i,:);
    avg_mdif_f3=mean(std_f3);
    avg_mrsosf3_new(count,:)=avg_mdif_f3;
    std_mdif_f3=std(std_f3);
    std_mrsosf3_new(count,:)=std_mdif_f3;
    count=count+1;
end
% Coefficient of variation of SSM_n across all models
cv_nsm1=std_mrsosf1_new./avg_mrsosf1_new;
cv_nsm2=std_mrsosf2_new./avg_mrsosf2_new;
cv_nsm3=std_mrsosf3_new./avg_mrsosf3_new;
nm_max1=max(max(cv_nsm1));
nm_max2=max(max(cv_nsm2));
nm_max3=max(max(cv_nsm3));
nm_maxn=[nm_max1,nm_max2,nm_max3];
nm_MAX=max(nm_maxn);
% Remove regions with SSM less than 0.1
avg_mdif1_try1(find(mean_SM<0.1))=-1;
avg_mdif2_try1(find(mean_SM<0.1))=-1;
avg_mdif3_try1(find(mean_SM<0.1))=-1;
cv_nsm1(find(mean_SM<0.1))=-0.1;
cv_nsm2(find(mean_SM<0.1))=-0.1;
cv_nsm3(find(mean_SM<0.1))=-0.1;
mdif1_positive_pct(find(mean_SM<0.1))=nan;
mdif1_negative_pct(find(mean_SM<0.1))=nan;
% Significance test by stippling for averaged biases for SSM_n
mask11=(mdif1_positive_pct==1) | (mdif1_negative_pct==1);
mask12=((mdif1_positive_pct>=14/17) & (mdif1_positive_pct<1)) | ((mdif1_negative_pct>=14/17) & (mdif1_negative_pct<1));
mdif2_positive_pct(find(mean_SM<0.1))=nan;
mdif2_negative_pct(find(mean_SM<0.1))=nan;
mask21=(mdif2_positive_pct==1) | (mdif2_negative_pct==1);
mask22=((mdif2_positive_pct>=14/17) & (mdif2_positive_pct<1)) | ((mdif2_negative_pct>=14/17) & (mdif2_negative_pct<1));
mdif3_positive_pct(find(mean_SM<0.1))=nan;
mdif3_negative_pct(find(mean_SM<0.1))=nan;
mask31=(mdif3_positive_pct==1) | (mdif3_negative_pct==1);
mask32=((mdif3_positive_pct>=14/17) & (mdif3_positive_pct<1)) | ((mdif3_negative_pct>=14/17) & (mdif3_negative_pct<1));
% Global map of averaged biased and coefficient of variation for SSM_n over different time scales (Figure 4)
x_axis=lonSMAP;
y_axis=latSMAP;
[LON,LAT]=meshgrid(x_axis,y_axis);
subplot(2,3,1)
map1=pcolor(x_axis,y_axis,avg_mdif1_try1);
set(map1,'alphadata',~isnan(avg_mdif1_try1))
colorbar
caxis([-1,1])
shading interp
hold on
h11=stipple(LON,LAT,mask11,'density',100,'color','k','marker','+','markersize',4);
h12=stipple(LON,LAT,mask12,'density',100,'color','k','marker','.','markersize',4);
hold off
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(a)')
ylabel('SSM_n diff AVG')
title('1/30<f<1/7 day^{-1}')
subplot(2,3,2)
map2=pcolor(x_axis,y_axis,avg_mdif2_try1);
set(map2,'alphadata',~isnan(avg_mdif2_try1))
colorbar
caxis([-1,1])
shading interp
hold on
h21=stipple(LON,LAT,mask21,'density',100,'color','k','marker','+','markersize',4);
h22=stipple(LON,LAT,mask22,'density',100,'color','k','marker','.','markersize',4);
hold off
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(b)')
title('1/90<f<1/30 day^{-1}')
subplot(2,3,3)
map3=pcolor(x_axis,y_axis,avg_mdif3_try1);
set(map3,'alphadata',~isnan(avg_mdif3_try1))
colorbar
caxis([-1,1])
shading interp
hold on
h31=stipple(LON,LAT,mask31,'density',100,'color','k','marker','+','markersize',4);
h32=stipple(LON,LAT,mask32,'density',100,'color','k','marker','.','markersize',4);
hold off
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(c)')
title('1/365<f<1/90 day^{-1}')
subplot(2,3,4)
map4=pcolor(x_axis,y_axis,cv_nsm1);
set(map4,'alphadata',~isnan(cv_nsm1))
colorbar
caxis([0,1.5])
shading interp
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(d)')
ylabel('SSM_n CV')
subplot(2,3,5)
map5=pcolor(x_axis,y_axis,cv_nsm2);
set(map5,'alphadata',~isnan(cv_nsm2))
colorbar
caxis([0,1.5])
shading interp
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(e)')
subplot(2,3,6)
map6=pcolor(x_axis,y_axis,cv_nsm3);
set(map6,'alphadata',~isnan(cv_nsm3))
colorbar
caxis([0,1.5])
shading interp
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(f)')
axes('position',[0.81,0.14,0.15,0.34])
axis off
colorbar('eastoutside')
caxis([-1,1])
hold on
axes('position',[0.41,0.14,0.15,0.34])
axis off
colorbar('eastoutside')
caxis([0,0.4])
hold off
h3=legend([h31,h32],'100%','>80%','Location','SouthEast','FontWeight','bold','FontSize',6);
title(h3,'Significance Test')