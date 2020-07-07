% % Average biases of spectral slope biases for SSM, ET, and P between models and observation over different time scales
% Averaged biases of SSM_kw across all models
avg_SMPdif1=(BCC_CSM_SMPdif1+BNU_ESM_SMPdif1+CanESM2_SMPdif1+CNRM_CM5_SMPdif1+CSIRO_Mk_SMPdif1+GFDL_CM3_SMPdif1+GFDL_ESM2G_SMPdif1+GFDL_ESM2M_SMPdif1+HadGEM2_CC_SMPdif1+HadGEM2_ES_SMPdif1+inmcm4_SMPdif1+MIROC5_SMPdif1+MIROC_ESM_SMPdif1+MIROC_ESM_CHEM_SMPdif1+MRI_CGCM3_SMPdif1+MRI_ESM1_SMPdif1+NorESM1_M_SMPdif1)/17;
avg_SMPdif2=(BCC_CSM_SMPdif2+BNU_ESM_SMPdif2+CanESM2_SMPdif2+CNRM_CM5_SMPdif2+CSIRO_Mk_SMPdif2+GFDL_CM3_SMPdif2+GFDL_ESM2G_SMPdif2+GFDL_ESM2M_SMPdif2+HadGEM2_CC_SMPdif2+HadGEM2_ES_SMPdif2+inmcm4_SMPdif2+MIROC5_SMPdif2+MIROC_ESM_SMPdif2+MIROC_ESM_CHEM_SMPdif2+MRI_CGCM3_SMPdif2+MRI_ESM1_SMPdif2+NorESM1_M_SMPdif2)/17;
avg_SMPdif3=(BCC_CSM_SMPdif3+BNU_ESM_SMPdif3+CanESM2_SMPdif3+CNRM_CM5_SMPdif3+CSIRO_Mk_SMPdif3+GFDL_CM3_SMPdif3+GFDL_ESM2G_SMPdif3+GFDL_ESM2M_SMPdif3+HadGEM2_CC_SMPdif3+HadGEM2_ES_SMPdif3+inmcm4_SMPdif3+MIROC5_SMPdif3+MIROC_ESM_SMPdif3+MIROC_ESM_CHEM_SMPdif3+MRI_CGCM3_SMPdif3+MRI_ESM1_SMPdif3+NorESM1_M_SMPdif3)/17;
% Prepare for significance test for averaged biases of SSM_kw (1st frequency range)
BCC_CSM_SMPdif1(BCC_CSM_SMPdif1<0)=0;
BCC_CSM_SMPdif1(BCC_CSM_SMPdif1>0)=1;
BNU_ESM_SMPdif1(BNU_ESM_SMPdif1<0)=0;
BNU_ESM_SMPdif1(BNU_ESM_SMPdif1>0)=1;
CanESM2_SMPdif1(CanESM2_SMPdif1<0)=0;
CanESM2_SMPdif1(CanESM2_SMPdif1>0)=1;
CNRM_CM5_SMPdif1(CNRM_CM5_SMPdif1<0)=0;
CNRM_CM5_SMPdif1(CNRM_CM5_SMPdif1>0)=1;
CSIRO_Mk_SMPdif1(CSIRO_Mk_SMPdif1<0)=0;
CSIRO_Mk_SMPdif1(CSIRO_Mk_SMPdif1>0)=1;
GFDL_CM3_SMPdif1(GFDL_CM3_SMPdif1<0)=0;
GFDL_CM3_SMPdif1(GFDL_CM3_SMPdif1>0)=1;
GFDL_ESM2G_SMPdif1(GFDL_ESM2G_SMPdif1<0)=0;
GFDL_ESM2G_SMPdif1(GFDL_ESM2G_SMPdif1>0)=1;
GFDL_ESM2M_SMPdif1(GFDL_ESM2M_SMPdif1<0)=0;
GFDL_ESM2M_SMPdif1(GFDL_ESM2M_SMPdif1>0)=1;
HadGEM2_CC_SMPdif1(HadGEM2_CC_SMPdif1<0)=0;
HadGEM2_CC_SMPdif1(HadGEM2_CC_SMPdif1>0)=1;
HadGEM2_ES_SMPdif1(HadGEM2_ES_SMPdif1<0)=0;
HadGEM2_ES_SMPdif1(HadGEM2_ES_SMPdif1>0)=1;
inmcm4_SMPdif1(inmcm4_SMPdif1<0)=0;
inmcm4_SMPdif1(inmcm4_SMPdif1>0)=1;
MIROC5_SMPdif1(MIROC5_SMPdif1<0)=0;
MIROC5_SMPdif1(MIROC5_SMPdif1>0)=1;
MIROC_ESM_SMPdif1(MIROC_ESM_SMPdif1<0)=0;
MIROC_ESM_SMPdif1(MIROC_ESM_SMPdif1>0)=1;
MIROC_ESM_CHEM_SMPdif1(MIROC_ESM_CHEM_SMPdif1<0)=0;
MIROC_ESM_CHEM_SMPdif1(MIROC_ESM_CHEM_SMPdif1>0)=1;
MRI_CGCM3_SMPdif1(MRI_CGCM3_SMPdif1<0)=0;
MRI_CGCM3_SMPdif1(MRI_CGCM3_SMPdif1>0)=1;
MRI_ESM1_SMPdif1(MRI_ESM1_SMPdif1<0)=0;
MRI_ESM1_SMPdif1(MRI_ESM1_SMPdif1>0)=1;
NorESM1_M_SMPdif1(NorESM1_M_SMPdif1<0)=0;
NorESM1_M_SMPdif1(NorESM1_M_SMPdif1>0)=1;
[m,n]=size(BCC_CSM_SMPdif1);
BCC_CSM_SMPdif1_r=reshape(BCC_CSM_SMPdif1,[m*n,1]);
BNU_ESM_SMPdif1_r=reshape(BNU_ESM_SMPdif1,[m*n,1]);
CanESM2_SMPdif1_r=reshape(CanESM2_SMPdif1,[m*n,1]);
CNRM_CM5_SMPdif1_r=reshape(CNRM_CM5_SMPdif1,[m*n,1]);
CSIRO_Mk_SMPdif1_r=reshape(CSIRO_Mk_SMPdif1,[m*n,1]);
GFDL_CM3_SMPdif1_r=reshape(GFDL_CM3_SMPdif1,[m*n,1]);
GFDL_ESM2G_SMPdif1_r=reshape(GFDL_ESM2G_SMPdif1,[m*n,1]);
GFDL_ESM2M_SMPdif1_r=reshape(GFDL_ESM2M_SMPdif1,[m*n,1]);
HadGEM2_CC_SMPdif1_r=reshape(HadGEM2_CC_SMPdif1,[m*n,1]);
HadGEM2_ES_SMPdif1_r=reshape(HadGEM2_ES_SMPdif1,[m*n,1]);
inmcm4_SMPdif1_r=reshape(inmcm4_SMPdif1,[m*n,1]);
MIROC5_SMPdif1_r=reshape(MIROC5_SMPdif1,[m*n,1]);
MIROC_ESM_SMPdif1_r=reshape(MIROC_ESM_SMPdif1,[m*n,1]);
MIROC_ESM_CHEM_SMPdif1_r=reshape(MIROC_ESM_CHEM_SMPdif1,[m*n,1]);
MRI_CGCM3_SMPdif1_r=reshape(MRI_CGCM3_SMPdif1,[m*n,1]);
MRI_ESM1_SMPdif1_r=reshape(MRI_ESM1_SMPdif1,[m*n,1]);
NorESM1_M_SMPdif1_r=reshape(NorESM1_M_SMPdif1,[m*n,1]);
SMPdif1_r=[BCC_CSM_SMPdif1_r,BNU_ESM_SMPdif1_r,CanESM2_SMPdif1_r,CNRM_CM5_SMPdif1_r,CSIRO_Mk_SMPdif1_r,GFDL_CM3_SMPdif1_r,GFDL_ESM2G_SMPdif1_r,GFDL_ESM2M_SMPdif1_r,HadGEM2_CC_SMPdif1_r,HadGEM2_ES_SMPdif1_r,inmcm4_SMPdif1_r,MIROC5_SMPdif1_r,MIROC_ESM_SMPdif1_r,MIROC_ESM_CHEM_SMPdif1_r,MRI_CGCM3_SMPdif1_r,MRI_ESM1_SMPdif1_r,NorESM1_M_SMPdif1_r];
% The number of models with the same sign as averaged biaes in each pixel (SSM_kw, 1st frequency range)
SMPdif1_r_nan=isnan(SMPdif1_r);
SMPdif1_r_nan=single(SMPdif1_r_nan);
SMPdif1_r_nan_sum=sum(SMPdif1_r_nan,2);
SMPdif1_r_nan_sum_c=[SMPdif1_r,SMPdif1_r_nan_sum];
[row,~]=size(SMPdif1_r);
SMPdif1_count_0(row,1)=0;
SMPdif1_count_1(row,1)=0;
count=1;
for i=1:row
    if SMPdif1_r_nan_sum_c(i,18)==17
        num_0=nan;
        num_1=nan;
    else
        num_0=sum(SMPdif1_r(i,:)==0,2);
        num_1=sum(SMPdif1_r(i,:)~=0,2);
    end
    SMPdif1_count_0(count,:)=num_0;
    SMPdif1_count_1(count,:)=num_1;
    count=count+1;
    disp(i)
end
% Ratio of the number of models with the same sign as averaged biaes in each pixel to the number of total models (SSM_kw, 1st frequency range)
SMPdif1_count_0_real=SMPdif1_count_0/17;
SMPdif1_count_1_real=SMPdif1_count_1/17;
SMPdif1_positive_pct=reshape(SMPdif1_count_0_real,[m,n]);
SMPdif1_negative_pct=reshape(SMPdif1_count_1_real,[m,n]);
% Prepare for significance test for averaged biases of SSM_kw (2nd frequency range)
BCC_CSM_SMPdif2(BCC_CSM_SMPdif2<0)=0;
BCC_CSM_SMPdif2(BCC_CSM_SMPdif2>0)=1;
BNU_ESM_SMPdif2(BNU_ESM_SMPdif2<0)=0;
BNU_ESM_SMPdif2(BNU_ESM_SMPdif2>0)=1;
CanESM2_SMPdif2(CanESM2_SMPdif2<0)=0;
CanESM2_SMPdif2(CanESM2_SMPdif2>0)=1;
CNRM_CM5_SMPdif2(CNRM_CM5_SMPdif2<0)=0;
CNRM_CM5_SMPdif2(CNRM_CM5_SMPdif2>0)=1;
CSIRO_Mk_SMPdif2(CSIRO_Mk_SMPdif2<0)=0;
CSIRO_Mk_SMPdif2(CSIRO_Mk_SMPdif2>0)=1;
GFDL_CM3_SMPdif2(GFDL_CM3_SMPdif2<0)=0;
GFDL_CM3_SMPdif2(GFDL_CM3_SMPdif2>0)=1;
GFDL_ESM2G_SMPdif2(GFDL_ESM2G_SMPdif2<0)=0;
GFDL_ESM2G_SMPdif2(GFDL_ESM2G_SMPdif2>0)=1;
GFDL_ESM2M_SMPdif2(GFDL_ESM2M_SMPdif2<0)=0;
GFDL_ESM2M_SMPdif2(GFDL_ESM2M_SMPdif2>0)=1;
HadGEM2_CC_SMPdif2(HadGEM2_CC_SMPdif2<0)=0;
HadGEM2_CC_SMPdif2(HadGEM2_CC_SMPdif2>0)=1;
HadGEM2_ES_SMPdif2(HadGEM2_ES_SMPdif2<0)=0;
HadGEM2_ES_SMPdif2(HadGEM2_ES_SMPdif2>0)=1;
inmcm4_SMPdif2(inmcm4_SMPdif2<0)=0;
inmcm4_SMPdif2(inmcm4_SMPdif2>0)=1;
MIROC5_SMPdif2(MIROC5_SMPdif2<0)=0;
MIROC5_SMPdif2(MIROC5_SMPdif2>0)=1;
MIROC_ESM_SMPdif2(MIROC_ESM_SMPdif2<0)=0;
MIROC_ESM_SMPdif2(MIROC_ESM_SMPdif2>0)=1;
MIROC_ESM_CHEM_SMPdif2(MIROC_ESM_CHEM_SMPdif2<0)=0;
MIROC_ESM_CHEM_SMPdif2(MIROC_ESM_CHEM_SMPdif2>0)=1;
MRI_CGCM3_SMPdif2(MRI_CGCM3_SMPdif2<0)=0;
MRI_CGCM3_SMPdif2(MRI_CGCM3_SMPdif2>0)=1;
MRI_ESM1_SMPdif2(MRI_ESM1_SMPdif2<0)=0;
MRI_ESM1_SMPdif2(MRI_ESM1_SMPdif2>0)=1;
NorESM1_M_SMPdif2(NorESM1_M_SMPdif2<0)=0;
NorESM1_M_SMPdif2(NorESM1_M_SMPdif2>0)=1;
[m,n]=size(BCC_CSM_SMPdif2);
BCC_CSM_SMPdif2_r=reshape(BCC_CSM_SMPdif2,[m*n,1]);
BNU_ESM_SMPdif2_r=reshape(BNU_ESM_SMPdif2,[m*n,1]);
CanESM2_SMPdif2_r=reshape(CanESM2_SMPdif2,[m*n,1]);
CNRM_CM5_SMPdif2_r=reshape(CNRM_CM5_SMPdif2,[m*n,1]);
CSIRO_Mk_SMPdif2_r=reshape(CSIRO_Mk_SMPdif2,[m*n,1]);
GFDL_CM3_SMPdif2_r=reshape(GFDL_CM3_SMPdif2,[m*n,1]);
GFDL_ESM2G_SMPdif2_r=reshape(GFDL_ESM2G_SMPdif2,[m*n,1]);
GFDL_ESM2M_SMPdif2_r=reshape(GFDL_ESM2M_SMPdif2,[m*n,1]);
HadGEM2_CC_SMPdif2_r=reshape(HadGEM2_CC_SMPdif2,[m*n,1]);
HadGEM2_ES_SMPdif2_r=reshape(HadGEM2_ES_SMPdif2,[m*n,1]);
inmcm4_SMPdif2_r=reshape(inmcm4_SMPdif2,[m*n,1]);
MIROC5_SMPdif2_r=reshape(MIROC5_SMPdif2,[m*n,1]);
MIROC_ESM_SMPdif2_r=reshape(MIROC_ESM_SMPdif2,[m*n,1]);
MIROC_ESM_CHEM_SMPdif2_r=reshape(MIROC_ESM_CHEM_SMPdif2,[m*n,1]);
MRI_CGCM3_SMPdif2_r=reshape(MRI_CGCM3_SMPdif2,[m*n,1]);
MRI_ESM1_SMPdif2_r=reshape(MRI_ESM1_SMPdif2,[m*n,1]);
NorESM1_M_SMPdif2_r=reshape(NorESM1_M_SMPdif2,[m*n,1]);
SMPdif2_r=[BCC_CSM_SMPdif2_r,BNU_ESM_SMPdif2_r,CanESM2_SMPdif2_r,CNRM_CM5_SMPdif2_r,CSIRO_Mk_SMPdif2_r,GFDL_CM3_SMPdif2_r,GFDL_ESM2G_SMPdif2_r,GFDL_ESM2M_SMPdif2_r,HadGEM2_CC_SMPdif2_r,HadGEM2_ES_SMPdif2_r,inmcm4_SMPdif2_r,MIROC5_SMPdif2_r,MIROC_ESM_SMPdif2_r,MIROC_ESM_CHEM_SMPdif2_r,MRI_CGCM3_SMPdif2_r,MRI_ESM1_SMPdif2_r,NorESM1_M_SMPdif2_r];
% The number of models with the same sign as averaged biaes in each pixel (SSM_kw, 2nd frequency range)
SMPdif2_r_nan=isnan(SMPdif2_r);
SMPdif2_r_nan=single(SMPdif2_r_nan);
SMPdif2_r_nan_sum=sum(SMPdif2_r_nan,2);
SMPdif2_r_nan_sum_c=[SMPdif2_r,SMPdif2_r_nan_sum];
[row,~]=size(SMPdif2_r);
SMPdif2_count_0(row,1)=0;
SMPdif2_count_1(row,1)=0;
count=1;
for i=1:row
    if SMPdif2_r_nan_sum_c(i,18)==17
        num_0=nan;
        num_1=nan;
    else
        num_0=sum(SMPdif2_r(i,:)==0,2);
        num_1=sum(SMPdif2_r(i,:)~=0,2);
    end
    SMPdif2_count_0(count,:)=num_0;
    SMPdif2_count_1(count,:)=num_1;
    count=count+1;
    disp(i)
end
% Ratio of the number of models with the same sign as averaged biaes in each pixel to the number of total models (SSM_kw, 2nd frequency range)
SMPdif2_count_0_real=SMPdif2_count_0/17;
SMPdif2_count_1_real=SMPdif2_count_1/17;
SMPdif2_positive_pct=reshape(SMPdif2_count_0_real,[m,n]);
SMPdif2_negative_pct=reshape(SMPdif2_count_1_real,[m,n]);
% Prepare for significance test for averaged biases of SSM_kw (3rd frequency range)
BCC_CSM_SMPdif3(BCC_CSM_SMPdif3<0)=0;
BCC_CSM_SMPdif3(BCC_CSM_SMPdif3>0)=1;
BNU_ESM_SMPdif3(BNU_ESM_SMPdif3<0)=0;
BNU_ESM_SMPdif3(BNU_ESM_SMPdif3>0)=1;
CanESM2_SMPdif3(CanESM2_SMPdif3<0)=0;
CanESM2_SMPdif3(CanESM2_SMPdif3>0)=1;
CNRM_CM5_SMPdif3(CNRM_CM5_SMPdif3<0)=0;
CNRM_CM5_SMPdif3(CNRM_CM5_SMPdif3>0)=1;
CSIRO_Mk_SMPdif3(CSIRO_Mk_SMPdif3<0)=0;
CSIRO_Mk_SMPdif3(CSIRO_Mk_SMPdif3>0)=1;
GFDL_CM3_SMPdif3(GFDL_CM3_SMPdif3<0)=0;
GFDL_CM3_SMPdif3(GFDL_CM3_SMPdif3>0)=1;
GFDL_ESM2G_SMPdif3(GFDL_ESM2G_SMPdif3<0)=0;
GFDL_ESM2G_SMPdif3(GFDL_ESM2G_SMPdif3>0)=1;
GFDL_ESM2M_SMPdif3(GFDL_ESM2M_SMPdif3<0)=0;
GFDL_ESM2M_SMPdif3(GFDL_ESM2M_SMPdif3>0)=1;
HadGEM2_CC_SMPdif3(HadGEM2_CC_SMPdif3<0)=0;
HadGEM2_CC_SMPdif3(HadGEM2_CC_SMPdif3>0)=1;
HadGEM2_ES_SMPdif3(HadGEM2_ES_SMPdif3<0)=0;
HadGEM2_ES_SMPdif3(HadGEM2_ES_SMPdif3>0)=1;
inmcm4_SMPdif3(inmcm4_SMPdif3<0)=0;
inmcm4_SMPdif3(inmcm4_SMPdif3>0)=1;
MIROC5_SMPdif3(MIROC5_SMPdif3<0)=0;
MIROC5_SMPdif3(MIROC5_SMPdif3>0)=1;
MIROC_ESM_SMPdif3(MIROC_ESM_SMPdif3<0)=0;
MIROC_ESM_SMPdif3(MIROC_ESM_SMPdif3>0)=1;
MIROC_ESM_CHEM_SMPdif3(MIROC_ESM_CHEM_SMPdif3<0)=0;
MIROC_ESM_CHEM_SMPdif3(MIROC_ESM_CHEM_SMPdif3>0)=1;
MRI_CGCM3_SMPdif3(MRI_CGCM3_SMPdif3<0)=0;
MRI_CGCM3_SMPdif3(MRI_CGCM3_SMPdif3>0)=1;
MRI_ESM1_SMPdif3(MRI_ESM1_SMPdif3<0)=0;
MRI_ESM1_SMPdif3(MRI_ESM1_SMPdif3>0)=1;
NorESM1_M_SMPdif3(NorESM1_M_SMPdif3<0)=0;
NorESM1_M_SMPdif3(NorESM1_M_SMPdif3>0)=1;
[m,n]=size(BCC_CSM_SMPdif3);
BCC_CSM_SMPdif3_r=reshape(BCC_CSM_SMPdif3,[m*n,1]);
BNU_ESM_SMPdif3_r=reshape(BNU_ESM_SMPdif3,[m*n,1]);
CanESM2_SMPdif3_r=reshape(CanESM2_SMPdif3,[m*n,1]);
CNRM_CM5_SMPdif3_r=reshape(CNRM_CM5_SMPdif3,[m*n,1]);
CSIRO_Mk_SMPdif3_r=reshape(CSIRO_Mk_SMPdif3,[m*n,1]);
GFDL_CM3_SMPdif3_r=reshape(GFDL_CM3_SMPdif3,[m*n,1]);
GFDL_ESM2G_SMPdif3_r=reshape(GFDL_ESM2G_SMPdif3,[m*n,1]);
GFDL_ESM2M_SMPdif3_r=reshape(GFDL_ESM2M_SMPdif3,[m*n,1]);
HadGEM2_CC_SMPdif3_r=reshape(HadGEM2_CC_SMPdif3,[m*n,1]);
HadGEM2_ES_SMPdif3_r=reshape(HadGEM2_ES_SMPdif3,[m*n,1]);
inmcm4_SMPdif3_r=reshape(inmcm4_SMPdif3,[m*n,1]);
MIROC5_SMPdif3_r=reshape(MIROC5_SMPdif3,[m*n,1]);
MIROC_ESM_SMPdif3_r=reshape(MIROC_ESM_SMPdif3,[m*n,1]);
MIROC_ESM_CHEM_SMPdif3_r=reshape(MIROC_ESM_CHEM_SMPdif3,[m*n,1]);
MRI_CGCM3_SMPdif3_r=reshape(MRI_CGCM3_SMPdif3,[m*n,1]);
MRI_ESM1_SMPdif3_r=reshape(MRI_ESM1_SMPdif3,[m*n,1]);
NorESM1_M_SMPdif3_r=reshape(NorESM1_M_SMPdif3,[m*n,1]);
SMPdif3_r=[BCC_CSM_SMPdif3_r,BNU_ESM_SMPdif3_r,CanESM2_SMPdif3_r,CNRM_CM5_SMPdif3_r,CSIRO_Mk_SMPdif3_r,GFDL_CM3_SMPdif3_r,GFDL_ESM2G_SMPdif3_r,GFDL_ESM2M_SMPdif3_r,HadGEM2_CC_SMPdif3_r,HadGEM2_ES_SMPdif3_r,inmcm4_SMPdif3_r,MIROC5_SMPdif3_r,MIROC_ESM_SMPdif3_r,MIROC_ESM_CHEM_SMPdif3_r,MRI_CGCM3_SMPdif3_r,MRI_ESM1_SMPdif3_r,NorESM1_M_SMPdif3_r];
% The number of models with the same sign as averaged biaes in each pixel (SSM_kw, 3rd frequency range)
SMPdif3_r_nan=isnan(SMPdif3_r);
SMPdif3_r_nan=single(SMPdif3_r_nan);
SMPdif3_r_nan_sum=sum(SMPdif3_r_nan,2);
SMPdif3_r_nan_sum_c=[SMPdif3_r,SMPdif3_r_nan_sum];
[row,~]=size(SMPdif3_r);
SMPdif3_count_0(row,1)=0;
SMPdif3_count_1(row,1)=0;
count=1;
for i=1:row
    if SMPdif3_r_nan_sum_c(i,18)==17
        num_0=nan;
        num_1=nan;
    else
        num_0=sum(SMPdif3_r(i,:)==0,2);
        num_1=sum(SMPdif3_r(i,:)~=0,2);
    end
    SMPdif3_count_0(count,:)=num_0;
    SMPdif3_count_1(count,:)=num_1;
    count=count+1;
    disp(i)
end
% Ratio of the number of models with the same sign as averaged biaes in each pixel to the number of total models (SSM_kw, 3rd frequency range)
SMPdif3_count_0_real=SMPdif3_count_0/17;
SMPdif3_count_1_real=SMPdif3_count_1/17;
SMPdif3_positive_pct=reshape(SMPdif3_count_0_real,[m,n]);
SMPdif3_negative_pct=reshape(SMPdif3_count_1_real,[m,n]);
% Averaged biases of ET_kw across all models
avg_EPdif1=(BCC_CSM_EPdif1+BNU_ESM_EPdif1+CanESM2_EPdif1+CNRM_CM5_EPdif1+CSIRO_Mk_EPdif1+GFDL_CM3_EPdif1+GFDL_ESM2G_EPdif1+GFDL_ESM2M_EPdif1+HadGEM2_CC_EPdif1+HadGEM2_ES_EPdif1+inmcm4_EPdif1+MIROC5_EPdif1+MIROC_ESM_EPdif1+MIROC_ESM_CHEM_EPdif1+MRI_CGCM3_EPdif1+MRI_ESM1_EPdif1+NorESM1_M_EPdif1)/17;
avg_EPdif2=(BCC_CSM_EPdif2+BNU_ESM_EPdif2+CanESM2_EPdif2+CNRM_CM5_EPdif2+CSIRO_Mk_EPdif2+GFDL_CM3_EPdif2+GFDL_ESM2G_EPdif2+GFDL_ESM2M_EPdif2+HadGEM2_CC_EPdif2+HadGEM2_ES_EPdif2+inmcm4_EPdif2+MIROC5_EPdif2+MIROC_ESM_EPdif2+MIROC_ESM_CHEM_EPdif2+MRI_CGCM3_EPdif2+MRI_ESM1_EPdif2+NorESM1_M_EPdif2)/17;
avg_EPdif3=(BCC_CSM_EPdif3+BNU_ESM_EPdif3+CanESM2_EPdif3+CNRM_CM5_EPdif3+CSIRO_Mk_EPdif3+GFDL_CM3_EPdif3+GFDL_ESM2G_EPdif3+GFDL_ESM2M_EPdif3+HadGEM2_CC_EPdif3+HadGEM2_ES_EPdif3+inmcm4_EPdif3+MIROC5_EPdif3+MIROC_ESM_EPdif3+MIROC_ESM_CHEM_EPdif3+MRI_CGCM3_EPdif3+MRI_ESM1_EPdif3+NorESM1_M_EPdif3)/17;
% Define the range of average biases of ET_kw
EPdif_max1=max(max(avg_EPdif1));
EPdif_max2=max(max(avg_EPdif2));
EPdif_max3=max(max(avg_EPdif3));
EPdif_maxn=[EPdif_max1,EPdif_max2,EPdif_max3];
EPdif_MAX=max(EPdif_maxn);
EPdif_min1=min(min(avg_EPdif1));
EPdif_min2=min(min(avg_EPdif2));
EPdif_min3=min(min(avg_EPdif3));
EPdif_minn=[EPdif_min1,EPdif_min2,EPdif_min3];
EPdif_MIN=min(EPdif_minn);
% Prepare for significance test for averaged biases of ET_kw (1st frequency range)
BCC_CSM_EPdif1(BCC_CSM_EPdif1<0)=0;
BCC_CSM_EPdif1(BCC_CSM_EPdif1>0)=1;
BNU_ESM_EPdif1(BNU_ESM_EPdif1<0)=0;
BNU_ESM_EPdif1(BNU_ESM_EPdif1>0)=1;
CanESM2_EPdif1(CanESM2_EPdif1<0)=0;
CanESM2_EPdif1(CanESM2_EPdif1>0)=1;
CNRM_CM5_EPdif1(CNRM_CM5_EPdif1<0)=0;
CNRM_CM5_EPdif1(CNRM_CM5_EPdif1>0)=1;
CSIRO_Mk_EPdif1(CSIRO_Mk_EPdif1<0)=0;
CSIRO_Mk_EPdif1(CSIRO_Mk_EPdif1>0)=1;
GFDL_CM3_EPdif1(GFDL_CM3_EPdif1<0)=0;
GFDL_CM3_EPdif1(GFDL_CM3_EPdif1>0)=1;
GFDL_ESM2G_EPdif1(GFDL_ESM2G_EPdif1<0)=0;
GFDL_ESM2G_EPdif1(GFDL_ESM2G_EPdif1>0)=1;
GFDL_ESM2M_EPdif1(GFDL_ESM2M_EPdif1<0)=0;
GFDL_ESM2M_EPdif1(GFDL_ESM2M_EPdif1>0)=1;
HadGEM2_CC_EPdif1(HadGEM2_CC_EPdif1<0)=0;
HadGEM2_CC_EPdif1(HadGEM2_CC_EPdif1>0)=1;
HadGEM2_ES_EPdif1(HadGEM2_ES_EPdif1<0)=0;
HadGEM2_ES_EPdif1(HadGEM2_ES_EPdif1>0)=1;
inmcm4_EPdif1(inmcm4_EPdif1<0)=0;
inmcm4_EPdif1(inmcm4_EPdif1>0)=1;
MIROC5_EPdif1(MIROC5_EPdif1<0)=0;
MIROC5_EPdif1(MIROC5_EPdif1>0)=1;
MIROC_ESM_EPdif1(MIROC_ESM_EPdif1<0)=0;
MIROC_ESM_EPdif1(MIROC_ESM_EPdif1>0)=1;
MIROC_ESM_CHEM_EPdif1(MIROC_ESM_CHEM_EPdif1<0)=0;
MIROC_ESM_CHEM_EPdif1(MIROC_ESM_CHEM_EPdif1>0)=1;
MRI_CGCM3_EPdif1(MRI_CGCM3_EPdif1<0)=0;
MRI_CGCM3_EPdif1(MRI_CGCM3_EPdif1>0)=1;
MRI_ESM1_EPdif1(MRI_ESM1_EPdif1<0)=0;
MRI_ESM1_EPdif1(MRI_ESM1_EPdif1>0)=1;
NorESM1_M_EPdif1(NorESM1_M_EPdif1<0)=0;
NorESM1_M_EPdif1(NorESM1_M_EPdif1>0)=1;
[m,n]=size(BCC_CSM_EPdif1);
BCC_CSM_EPdif1_r=reshape(BCC_CSM_EPdif1,[m*n,1]);
BNU_ESM_EPdif1_r=reshape(BNU_ESM_EPdif1,[m*n,1]);
CanESM2_EPdif1_r=reshape(CanESM2_EPdif1,[m*n,1]);
CNRM_CM5_EPdif1_r=reshape(CNRM_CM5_EPdif1,[m*n,1]);
CSIRO_Mk_EPdif1_r=reshape(CSIRO_Mk_EPdif1,[m*n,1]);
GFDL_CM3_EPdif1_r=reshape(GFDL_CM3_EPdif1,[m*n,1]);
GFDL_ESM2G_EPdif1_r=reshape(GFDL_ESM2G_EPdif1,[m*n,1]);
GFDL_ESM2M_EPdif1_r=reshape(GFDL_ESM2M_EPdif1,[m*n,1]);
HadGEM2_CC_EPdif1_r=reshape(HadGEM2_CC_EPdif1,[m*n,1]);
HadGEM2_ES_EPdif1_r=reshape(HadGEM2_ES_EPdif1,[m*n,1]);
inmcm4_EPdif1_r=reshape(inmcm4_EPdif1,[m*n,1]);
MIROC5_EPdif1_r=reshape(MIROC5_EPdif1,[m*n,1]);
MIROC_ESM_EPdif1_r=reshape(MIROC_ESM_EPdif1,[m*n,1]);
MIROC_ESM_CHEM_EPdif1_r=reshape(MIROC_ESM_CHEM_EPdif1,[m*n,1]);
MRI_CGCM3_EPdif1_r=reshape(MRI_CGCM3_EPdif1,[m*n,1]);
MRI_ESM1_EPdif1_r=reshape(MRI_ESM1_EPdif1,[m*n,1]);
NorESM1_M_EPdif1_r=reshape(NorESM1_M_EPdif1,[m*n,1]);
EPdif1_r=[BCC_CSM_EPdif1_r,BNU_ESM_EPdif1_r,CanESM2_EPdif1_r,CNRM_CM5_EPdif1_r,CSIRO_Mk_EPdif1_r,GFDL_CM3_EPdif1_r,GFDL_ESM2G_EPdif1_r,GFDL_ESM2M_EPdif1_r,HadGEM2_CC_EPdif1_r,HadGEM2_ES_EPdif1_r,inmcm4_EPdif1_r,MIROC5_EPdif1_r,MIROC_ESM_EPdif1_r,MIROC_ESM_CHEM_EPdif1_r,MRI_CGCM3_EPdif1_r,MRI_ESM1_EPdif1_r,NorESM1_M_EPdif1_r];
% The number of models with the same sign as averaged biaes in each pixel (ET_kw, 1st frequency range)
EPdif1_r_nan=isnan(EPdif1_r);
EPdif1_r_nan=single(EPdif1_r_nan);
EPdif1_r_nan_sum=sum(EPdif1_r_nan,2);
EPdif1_r_nan_sum_c=[EPdif1_r,EPdif1_r_nan_sum];
[row,~]=size(EPdif1_r);
EPdif1_count_0(row,1)=0;
EPdif1_count_1(row,1)=0;
count=1;
for i=1:row
    if EPdif1_r_nan_sum_c(i,18)==17
        num_0=nan;
        num_1=nan;
    else
        num_0=sum(EPdif1_r(i,:)==0,2);
        num_1=sum(EPdif1_r(i,:)~=0,2);
    end
    EPdif1_count_0(count,:)=num_0;
    EPdif1_count_1(count,:)=num_1;
    count=count+1;
    disp(i)
end
% Ratio of the number of models with the same sign as averaged biaes in each pixel to the number of total models (ET_kw, 1st frequency range)
EPdif1_count_0_real=EPdif1_count_0/17;
EPdif1_count_1_real=EPdif1_count_1/17;
EPdif1_positive_pct=reshape(EPdif1_count_0_real,[m,n]);
EPdif1_negative_pct=reshape(EPdif1_count_1_real,[m,n]);
% Prepare for significance test for averaged biases of ET_kw (2nd frequency range)
BCC_CSM_EPdif2(BCC_CSM_EPdif2<0)=0;
BCC_CSM_EPdif2(BCC_CSM_EPdif2>0)=1;
BNU_ESM_EPdif2(BNU_ESM_EPdif2<0)=0;
BNU_ESM_EPdif2(BNU_ESM_EPdif2>0)=1;
CanESM2_EPdif2(CanESM2_EPdif2<0)=0;
CanESM2_EPdif2(CanESM2_EPdif2>0)=1;
CNRM_CM5_EPdif2(CNRM_CM5_EPdif2<0)=0;
CNRM_CM5_EPdif2(CNRM_CM5_EPdif2>0)=1;
CSIRO_Mk_EPdif2(CSIRO_Mk_EPdif2<0)=0;
CSIRO_Mk_EPdif2(CSIRO_Mk_EPdif2>0)=1;
GFDL_CM3_EPdif2(GFDL_CM3_EPdif2<0)=0;
GFDL_CM3_EPdif2(GFDL_CM3_EPdif2>0)=1;
GFDL_ESM2G_EPdif2(GFDL_ESM2G_EPdif2<0)=0;
GFDL_ESM2G_EPdif2(GFDL_ESM2G_EPdif2>0)=1;
GFDL_ESM2M_EPdif2(GFDL_ESM2M_EPdif2<0)=0;
GFDL_ESM2M_EPdif2(GFDL_ESM2M_EPdif2>0)=1;
HadGEM2_CC_EPdif2(HadGEM2_CC_EPdif2<0)=0;
HadGEM2_CC_EPdif2(HadGEM2_CC_EPdif2>0)=1;
HadGEM2_ES_EPdif2(HadGEM2_ES_EPdif2<0)=0;
HadGEM2_ES_EPdif2(HadGEM2_ES_EPdif2>0)=1;
inmcm4_EPdif2(inmcm4_EPdif2<0)=0;
inmcm4_EPdif2(inmcm4_EPdif2>0)=1;
MIROC5_EPdif2(MIROC5_EPdif2<0)=0;
MIROC5_EPdif2(MIROC5_EPdif2>0)=1;
MIROC_ESM_EPdif2(MIROC_ESM_EPdif2<0)=0;
MIROC_ESM_EPdif2(MIROC_ESM_EPdif2>0)=1;
MIROC_ESM_CHEM_EPdif2(MIROC_ESM_CHEM_EPdif2<0)=0;
MIROC_ESM_CHEM_EPdif2(MIROC_ESM_CHEM_EPdif2>0)=1;
MRI_CGCM3_EPdif2(MRI_CGCM3_EPdif2<0)=0;
MRI_CGCM3_EPdif2(MRI_CGCM3_EPdif2>0)=1;
MRI_ESM1_EPdif2(MRI_ESM1_EPdif2<0)=0;
MRI_ESM1_EPdif2(MRI_ESM1_EPdif2>0)=1;
NorESM1_M_EPdif2(NorESM1_M_EPdif2<0)=0;
NorESM1_M_EPdif2(NorESM1_M_EPdif2>0)=1;
[m,n]=size(BCC_CSM_EPdif2);
BCC_CSM_EPdif2_r=reshape(BCC_CSM_EPdif2,[m*n,1]);
BNU_ESM_EPdif2_r=reshape(BNU_ESM_EPdif2,[m*n,1]);
CanESM2_EPdif2_r=reshape(CanESM2_EPdif2,[m*n,1]);
CNRM_CM5_EPdif2_r=reshape(CNRM_CM5_EPdif2,[m*n,1]);
CSIRO_Mk_EPdif2_r=reshape(CSIRO_Mk_EPdif2,[m*n,1]);
GFDL_CM3_EPdif2_r=reshape(GFDL_CM3_EPdif2,[m*n,1]);
GFDL_ESM2G_EPdif2_r=reshape(GFDL_ESM2G_EPdif2,[m*n,1]);
GFDL_ESM2M_EPdif2_r=reshape(GFDL_ESM2M_EPdif2,[m*n,1]);
HadGEM2_CC_EPdif2_r=reshape(HadGEM2_CC_EPdif2,[m*n,1]);
HadGEM2_ES_EPdif2_r=reshape(HadGEM2_ES_EPdif2,[m*n,1]);
inmcm4_EPdif2_r=reshape(inmcm4_EPdif2,[m*n,1]);
MIROC5_EPdif2_r=reshape(MIROC5_EPdif2,[m*n,1]);
MIROC_ESM_EPdif2_r=reshape(MIROC_ESM_EPdif2,[m*n,1]);
MIROC_ESM_CHEM_EPdif2_r=reshape(MIROC_ESM_CHEM_EPdif2,[m*n,1]);
MRI_CGCM3_EPdif2_r=reshape(MRI_CGCM3_EPdif2,[m*n,1]);
MRI_ESM1_EPdif2_r=reshape(MRI_ESM1_EPdif2,[m*n,1]);
NorESM1_M_EPdif2_r=reshape(NorESM1_M_EPdif2,[m*n,1]);
EPdif2_r=[BCC_CSM_EPdif2_r,BNU_ESM_EPdif2_r,CanESM2_EPdif2_r,CNRM_CM5_EPdif2_r,CSIRO_Mk_EPdif2_r,GFDL_CM3_EPdif2_r,GFDL_ESM2G_EPdif2_r,GFDL_ESM2M_EPdif2_r,HadGEM2_CC_EPdif2_r,HadGEM2_ES_EPdif2_r,inmcm4_EPdif2_r,MIROC5_EPdif2_r,MIROC_ESM_EPdif2_r,MIROC_ESM_CHEM_EPdif2_r,MRI_CGCM3_EPdif2_r,MRI_ESM1_EPdif2_r,NorESM1_M_EPdif2_r];
% The number of models with the same sign as averaged biaes in each pixel (ET_kw, 2nd frequency range)
EPdif2_r_nan=isnan(EPdif2_r);
EPdif2_r_nan=single(EPdif2_r_nan);
EPdif2_r_nan_sum=sum(EPdif2_r_nan,2);
EPdif2_r_nan_sum_c=[EPdif2_r,EPdif2_r_nan_sum];
[row,~]=size(EPdif2_r);
EPdif2_count_0(row,1)=0;
EPdif2_count_1(row,1)=0;
count=1;
for i=1:row
    if EPdif2_r_nan_sum_c(i,18)==17
        num_0=nan;
        num_1=nan;
    else
        num_0=sum(EPdif2_r(i,:)==0,2);
        num_1=sum(EPdif2_r(i,:)~=0,2);
    end
    EPdif2_count_0(count,:)=num_0;
    EPdif2_count_1(count,:)=num_1;
    count=count+1;
    disp(i)
end
% Ratio of the number of models with the same sign as averaged biaes in each pixel to the number of total models (SSM_kw, 2nd frequency range)
EPdif2_count_0_real=EPdif2_count_0/17;
EPdif2_count_1_real=EPdif2_count_1/17;
EPdif2_positive_pct=reshape(EPdif2_count_0_real,[m,n]);
EPdif2_negative_pct=reshape(EPdif2_count_1_real,[m,n]);
% Prepare for significance test for averaged biases of ET_kw (3rd frequency range)
BCC_CSM_EPdif3(BCC_CSM_EPdif3<0)=0;
BCC_CSM_EPdif3(BCC_CSM_EPdif3>0)=1;
BNU_ESM_EPdif3(BNU_ESM_EPdif3<0)=0;
BNU_ESM_EPdif3(BNU_ESM_EPdif3>0)=1;
CanESM2_EPdif3(CanESM2_EPdif3<0)=0;
CanESM2_EPdif3(CanESM2_EPdif3>0)=1;
CNRM_CM5_EPdif3(CNRM_CM5_EPdif3<0)=0;
CNRM_CM5_EPdif3(CNRM_CM5_EPdif3>0)=1;
CSIRO_Mk_EPdif3(CSIRO_Mk_EPdif3<0)=0;
CSIRO_Mk_EPdif3(CSIRO_Mk_EPdif3>0)=1;
GFDL_CM3_EPdif3(GFDL_CM3_EPdif3<0)=0;
GFDL_CM3_EPdif3(GFDL_CM3_EPdif3>0)=1;
GFDL_ESM2G_EPdif3(GFDL_ESM2G_EPdif3<0)=0;
GFDL_ESM2G_EPdif3(GFDL_ESM2G_EPdif3>0)=1;
GFDL_ESM2M_EPdif3(GFDL_ESM2M_EPdif3<0)=0;
GFDL_ESM2M_EPdif3(GFDL_ESM2M_EPdif3>0)=1;
HadGEM2_CC_EPdif3(HadGEM2_CC_EPdif3<0)=0;
HadGEM2_CC_EPdif3(HadGEM2_CC_EPdif3>0)=1;
HadGEM2_ES_EPdif3(HadGEM2_ES_EPdif3<0)=0;
HadGEM2_ES_EPdif3(HadGEM2_ES_EPdif3>0)=1;
inmcm4_EPdif3(inmcm4_EPdif3<0)=0;
inmcm4_EPdif3(inmcm4_EPdif3>0)=1;
MIROC5_EPdif3(MIROC5_EPdif3<0)=0;
MIROC5_EPdif3(MIROC5_EPdif3>0)=1;
MIROC_ESM_EPdif3(MIROC_ESM_EPdif3<0)=0;
MIROC_ESM_EPdif3(MIROC_ESM_EPdif3>0)=1;
MIROC_ESM_CHEM_EPdif3(MIROC_ESM_CHEM_EPdif3<0)=0;
MIROC_ESM_CHEM_EPdif3(MIROC_ESM_CHEM_EPdif3>0)=1;
MRI_CGCM3_EPdif3(MRI_CGCM3_EPdif3<0)=0;
MRI_CGCM3_EPdif3(MRI_CGCM3_EPdif3>0)=1;
MRI_ESM1_EPdif3(MRI_ESM1_EPdif3<0)=0;
MRI_ESM1_EPdif3(MRI_ESM1_EPdif3>0)=1;
NorESM1_M_EPdif3(NorESM1_M_EPdif3<0)=0;
NorESM1_M_EPdif3(NorESM1_M_EPdif3>0)=1;
[m,n]=size(BCC_CSM_EPdif3);
BCC_CSM_EPdif3_r=reshape(BCC_CSM_EPdif3,[m*n,1]);
BNU_ESM_EPdif3_r=reshape(BNU_ESM_EPdif3,[m*n,1]);
CanESM2_EPdif3_r=reshape(CanESM2_EPdif3,[m*n,1]);
CNRM_CM5_EPdif3_r=reshape(CNRM_CM5_EPdif3,[m*n,1]);
CSIRO_Mk_EPdif3_r=reshape(CSIRO_Mk_EPdif3,[m*n,1]);
GFDL_CM3_EPdif3_r=reshape(GFDL_CM3_EPdif3,[m*n,1]);
GFDL_ESM2G_EPdif3_r=reshape(GFDL_ESM2G_EPdif3,[m*n,1]);
GFDL_ESM2M_EPdif3_r=reshape(GFDL_ESM2M_EPdif3,[m*n,1]);
HadGEM2_CC_EPdif3_r=reshape(HadGEM2_CC_EPdif3,[m*n,1]);
HadGEM2_ES_EPdif3_r=reshape(HadGEM2_ES_EPdif3,[m*n,1]);
inmcm4_EPdif3_r=reshape(inmcm4_EPdif3,[m*n,1]);
MIROC5_EPdif3_r=reshape(MIROC5_EPdif3,[m*n,1]);
MIROC_ESM_EPdif3_r=reshape(MIROC_ESM_EPdif3,[m*n,1]);
MIROC_ESM_CHEM_EPdif3_r=reshape(MIROC_ESM_CHEM_EPdif3,[m*n,1]);
MRI_CGCM3_EPdif3_r=reshape(MRI_CGCM3_EPdif3,[m*n,1]);
MRI_ESM1_EPdif3_r=reshape(MRI_ESM1_EPdif3,[m*n,1]);
NorESM1_M_EPdif3_r=reshape(NorESM1_M_EPdif3,[m*n,1]);
EPdif3_r=[BCC_CSM_EPdif3_r,BNU_ESM_EPdif3_r,CanESM2_EPdif3_r,CNRM_CM5_EPdif3_r,CSIRO_Mk_EPdif3_r,GFDL_CM3_EPdif3_r,GFDL_ESM2G_EPdif3_r,GFDL_ESM2M_EPdif3_r,HadGEM2_CC_EPdif3_r,HadGEM2_ES_EPdif3_r,inmcm4_EPdif3_r,MIROC5_EPdif3_r,MIROC_ESM_EPdif3_r,MIROC_ESM_CHEM_EPdif3_r,MRI_CGCM3_EPdif3_r,MRI_ESM1_EPdif3_r,NorESM1_M_EPdif3_r];
% The number of models with the same sign as averaged biaes in each pixel (ET_kw, 3rd frequency range)
EPdif3_r_nan=isnan(EPdif3_r);
EPdif3_r_nan=single(EPdif3_r_nan);
EPdif3_r_nan_sum=sum(EPdif3_r_nan,2);
EPdif3_r_nan_sum_c=[EPdif3_r,EPdif3_r_nan_sum];
[row,~]=size(EPdif3_r);
EPdif3_count_0(row,1)=0;
EPdif3_count_1(row,1)=0;
count=1;
for i=1:row
    if EPdif3_r_nan_sum_c(i,18)==17
        num_0=nan;
        num_1=nan;
    else
        num_0=sum(EPdif3_r(i,:)==0,2);
        num_1=sum(EPdif3_r(i,:)~=0,2);
    end
    EPdif3_count_0(count,:)=num_0;
    EPdif3_count_1(count,:)=num_1;
    count=count+1;
    disp(i)
end
% Ratio of the number of models with the same sign as averaged biaes in each pixel to the number of total models (SSM_kw, 3rd frequency range)
EPdif3_count_0_real=EPdif3_count_0/17;
EPdif3_count_1_real=EPdif3_count_1/17;
EPdif3_positive_pct=reshape(EPdif3_count_0_real,[m,n]);
EPdif3_negative_pct=reshape(EPdif3_count_1_real,[m,n]);
% Averaged biases of P_kw across all models
avg_PPdif1=(BCC_CSM_PPdif1+BNU_ESM_PPdif1+CanESM2_PPdif1+CNRM_CM5_PPdif1+CSIRO_Mk_PPdif1+GFDL_CM3_PPdif1+GFDL_ESM2G_PPdif1+GFDL_ESM2M_PPdif1+HadGEM2_CC_PPdif1+HadGEM2_ES_PPdif1+inmcm4_PPdif1+MIROC5_PPdif1+MIROC_ESM_PPdif1+MIROC_ESM_CHEM_PPdif1+MRI_CGCM3_PPdif1+MRI_ESM1_PPdif1+NorESM1_M_PPdif1)/17;
avg_PPdif2=(BCC_CSM_PPdif2+BNU_ESM_PPdif2+CanESM2_PPdif2+CNRM_CM5_PPdif2+CSIRO_Mk_PPdif2+GFDL_CM3_PPdif2+GFDL_ESM2G_PPdif2+GFDL_ESM2M_PPdif2+HadGEM2_CC_PPdif2+HadGEM2_ES_PPdif2+inmcm4_PPdif2+MIROC5_PPdif2+MIROC_ESM_PPdif2+MIROC_ESM_CHEM_PPdif2+MRI_CGCM3_PPdif2+MRI_ESM1_PPdif2+NorESM1_M_PPdif2)/17;
avg_PPdif3=(BCC_CSM_PPdif3+BNU_ESM_PPdif3+CanESM2_PPdif3+CNRM_CM5_PPdif3+CSIRO_Mk_PPdif3+GFDL_CM3_PPdif3+GFDL_ESM2G_PPdif3+GFDL_ESM2M_PPdif3+HadGEM2_CC_PPdif3+HadGEM2_ES_PPdif3+inmcm4_PPdif3+MIROC5_PPdif3+MIROC_ESM_PPdif3+MIROC_ESM_CHEM_PPdif3+MRI_CGCM3_PPdif3+MRI_ESM1_PPdif3+NorESM1_M_PPdif3)/17;
% Define the range of average biases of P_kw
PPdif_max1=max(max(avg_PPdif1));
PPdif_max2=max(max(avg_PPdif2));
PPdif_max3=max(max(avg_PPdif3));
PPdif_maxn=[PPdif_max1,PPdif_max2,PPdif_max3];
PPdif_MAX=max(PPdif_maxn);
PPdif_min1=min(min(avg_PPdif1));
PPdif_min2=min(min(avg_PPdif2));
PPdif_min3=min(min(avg_PPdif3));
PPdif_minn=[PPdif_min1,PPdif_min2,PPdif_min3];
PPdif_MIN=min(PPdif_minn);
% Prepare for significance test for averaged biases of P_kw (1st frequency range)
BCC_CSM_PPdif1(BCC_CSM_PPdif1<0)=0;
BCC_CSM_PPdif1(BCC_CSM_PPdif1>0)=1;
BNU_ESM_PPdif1(BNU_ESM_PPdif1<0)=0;
BNU_ESM_PPdif1(BNU_ESM_PPdif1>0)=1;
CanESM2_PPdif1(CanESM2_PPdif1<0)=0;
CanESM2_PPdif1(CanESM2_PPdif1>0)=1;
CNRM_CM5_PPdif1(CNRM_CM5_PPdif1<0)=0;
CNRM_CM5_PPdif1(CNRM_CM5_PPdif1>0)=1;
CSIRO_Mk_PPdif1(CSIRO_Mk_PPdif1<0)=0;
CSIRO_Mk_PPdif1(CSIRO_Mk_PPdif1>0)=1;
GFDL_CM3_PPdif1(GFDL_CM3_PPdif1<0)=0;
GFDL_CM3_PPdif1(GFDL_CM3_PPdif1>0)=1;
GFDL_ESM2G_PPdif1(GFDL_ESM2G_PPdif1<0)=0;
GFDL_ESM2G_PPdif1(GFDL_ESM2G_PPdif1>0)=1;
GFDL_ESM2M_PPdif1(GFDL_ESM2M_PPdif1<0)=0;
GFDL_ESM2M_PPdif1(GFDL_ESM2M_PPdif1>0)=1;
HadGEM2_CC_PPdif1(HadGEM2_CC_PPdif1<0)=0;
HadGEM2_CC_PPdif1(HadGEM2_CC_PPdif1>0)=1;
HadGEM2_ES_PPdif1(HadGEM2_ES_PPdif1<0)=0;
HadGEM2_ES_PPdif1(HadGEM2_ES_PPdif1>0)=1;
inmcm4_PPdif1(inmcm4_PPdif1<0)=0;
inmcm4_PPdif1(inmcm4_PPdif1>0)=1;
MIROC5_PPdif1(MIROC5_PPdif1<0)=0;
MIROC5_PPdif1(MIROC5_PPdif1>0)=1;
MIROC_ESM_PPdif1(MIROC_ESM_PPdif1<0)=0;
MIROC_ESM_PPdif1(MIROC_ESM_PPdif1>0)=1;
MIROC_ESM_CHEM_PPdif1(MIROC_ESM_CHEM_PPdif1<0)=0;
MIROC_ESM_CHEM_PPdif1(MIROC_ESM_CHEM_PPdif1>0)=1;
MRI_CGCM3_PPdif1(MRI_CGCM3_PPdif1<0)=0;
MRI_CGCM3_PPdif1(MRI_CGCM3_PPdif1>0)=1;
MRI_ESM1_PPdif1(MRI_ESM1_PPdif1<0)=0;
MRI_ESM1_PPdif1(MRI_ESM1_PPdif1>0)=1;
NorESM1_M_PPdif1(NorESM1_M_PPdif1<0)=0;
NorESM1_M_PPdif1(NorESM1_M_PPdif1>0)=1;
[m,n]=size(BCC_CSM_PPdif1);
BCC_CSM_PPdif1_r=reshape(BCC_CSM_PPdif1,[m*n,1]);
BNU_ESM_PPdif1_r=reshape(BNU_ESM_PPdif1,[m*n,1]);
CanESM2_PPdif1_r=reshape(CanESM2_PPdif1,[m*n,1]);
CNRM_CM5_PPdif1_r=reshape(CNRM_CM5_PPdif1,[m*n,1]);
CSIRO_Mk_PPdif1_r=reshape(CSIRO_Mk_PPdif1,[m*n,1]);
GFDL_CM3_PPdif1_r=reshape(GFDL_CM3_PPdif1,[m*n,1]);
GFDL_ESM2G_PPdif1_r=reshape(GFDL_ESM2G_PPdif1,[m*n,1]);
GFDL_ESM2M_PPdif1_r=reshape(GFDL_ESM2M_PPdif1,[m*n,1]);
HadGEM2_CC_PPdif1_r=reshape(HadGEM2_CC_PPdif1,[m*n,1]);
HadGEM2_ES_PPdif1_r=reshape(HadGEM2_ES_PPdif1,[m*n,1]);
inmcm4_PPdif1_r=reshape(inmcm4_PPdif1,[m*n,1]);
MIROC5_PPdif1_r=reshape(MIROC5_PPdif1,[m*n,1]);
MIROC_ESM_PPdif1_r=reshape(MIROC_ESM_PPdif1,[m*n,1]);
MIROC_ESM_CHEM_PPdif1_r=reshape(MIROC_ESM_CHEM_PPdif1,[m*n,1]);
MRI_CGCM3_PPdif1_r=reshape(MRI_CGCM3_PPdif1,[m*n,1]);
MRI_ESM1_PPdif1_r=reshape(MRI_ESM1_PPdif1,[m*n,1]);
NorESM1_M_PPdif1_r=reshape(NorESM1_M_PPdif1,[m*n,1]);
PPdif1_r=[BCC_CSM_PPdif1_r,BNU_ESM_PPdif1_r,CanESM2_PPdif1_r,CNRM_CM5_PPdif1_r,CSIRO_Mk_PPdif1_r,GFDL_CM3_PPdif1_r,GFDL_ESM2G_PPdif1_r,GFDL_ESM2M_PPdif1_r,HadGEM2_CC_PPdif1_r,HadGEM2_ES_PPdif1_r,inmcm4_PPdif1_r,MIROC5_PPdif1_r,MIROC_ESM_PPdif1_r,MIROC_ESM_CHEM_PPdif1_r,MRI_CGCM3_PPdif1_r,MRI_ESM1_PPdif1_r,NorESM1_M_PPdif1_r];
% The number of models with the same sign as averaged biaes in each pixel (P_kw, 1st frequency range)
PPdif1_r_nan=isnan(PPdif1_r);
PPdif1_r_nan=single(PPdif1_r_nan);
PPdif1_r_nan_sum=sum(PPdif1_r_nan,2);
PPdif1_r_nan_sum_c=[PPdif1_r,PPdif1_r_nan_sum];
[row,~]=size(PPdif1_r);
PPdif1_count_0(row,1)=0;
PPdif1_count_1(row,1)=0;
count=1;
for i=1:row
    if PPdif1_r_nan_sum_c(i,18)==17
        num_0=nan;
        num_1=nan;
    else
        num_0=sum(PPdif1_r(i,:)==0,2);
        num_1=sum(PPdif1_r(i,:)~=0,2);
    end
    PPdif1_count_0(count,:)=num_0;
    PPdif1_count_1(count,:)=num_1;
    count=count+1;
    disp(i)
end
% Ratio of the number of models with the same sign as averaged biaes in each pixel to the number of total models (P_kw, 1st frequency range)
PPdif1_count_0_real=PPdif1_count_0/17;
PPdif1_count_1_real=PPdif1_count_1/17;
PPdif1_positive_pct=reshape(PPdif1_count_0_real,[m,n]);
PPdif1_negative_pct=reshape(PPdif1_count_1_real,[m,n]);
% Prepare for significance test for averaged biases of P_kw (2nd frequency range)
BCC_CSM_PPdif2(BCC_CSM_PPdif2<0)=0;
BCC_CSM_PPdif2(BCC_CSM_PPdif2>0)=1;
BNU_ESM_PPdif2(BNU_ESM_PPdif2<0)=0;
BNU_ESM_PPdif2(BNU_ESM_PPdif2>0)=1;
CanESM2_PPdif2(CanESM2_PPdif2<0)=0;
CanESM2_PPdif2(CanESM2_PPdif2>0)=1;
CNRM_CM5_PPdif2(CNRM_CM5_PPdif2<0)=0;
CNRM_CM5_PPdif2(CNRM_CM5_PPdif2>0)=1;
CSIRO_Mk_PPdif2(CSIRO_Mk_PPdif2<0)=0;
CSIRO_Mk_PPdif2(CSIRO_Mk_PPdif2>0)=1;
GFDL_CM3_PPdif2(GFDL_CM3_PPdif2<0)=0;
GFDL_CM3_PPdif2(GFDL_CM3_PPdif2>0)=1;
GFDL_ESM2G_PPdif2(GFDL_ESM2G_PPdif2<0)=0;
GFDL_ESM2G_PPdif2(GFDL_ESM2G_PPdif2>0)=1;
GFDL_ESM2M_PPdif2(GFDL_ESM2M_PPdif2<0)=0;
GFDL_ESM2M_PPdif2(GFDL_ESM2M_PPdif2>0)=1;
HadGEM2_CC_PPdif2(HadGEM2_CC_PPdif2<0)=0;
HadGEM2_CC_PPdif2(HadGEM2_CC_PPdif2>0)=1;
HadGEM2_ES_PPdif2(HadGEM2_ES_PPdif2<0)=0;
HadGEM2_ES_PPdif2(HadGEM2_ES_PPdif2>0)=1;
inmcm4_PPdif2(inmcm4_PPdif2<0)=0;
inmcm4_PPdif2(inmcm4_PPdif2>0)=1;
MIROC5_PPdif2(MIROC5_PPdif2<0)=0;
MIROC5_PPdif2(MIROC5_PPdif2>0)=1;
MIROC_ESM_PPdif2(MIROC_ESM_PPdif2<0)=0;
MIROC_ESM_PPdif2(MIROC_ESM_PPdif2>0)=1;
MIROC_ESM_CHEM_PPdif2(MIROC_ESM_CHEM_PPdif2<0)=0;
MIROC_ESM_CHEM_PPdif2(MIROC_ESM_CHEM_PPdif2>0)=1;
MRI_CGCM3_PPdif2(MRI_CGCM3_PPdif2<0)=0;
MRI_CGCM3_PPdif2(MRI_CGCM3_PPdif2>0)=1;
MRI_ESM1_PPdif2(MRI_ESM1_PPdif2<0)=0;
MRI_ESM1_PPdif2(MRI_ESM1_PPdif2>0)=1;
NorESM1_M_PPdif2(NorESM1_M_PPdif2<0)=0;
NorESM1_M_PPdif2(NorESM1_M_PPdif2>0)=1;
[m,n]=size(BCC_CSM_PPdif2);
BCC_CSM_PPdif2_r=reshape(BCC_CSM_PPdif2,[m*n,1]);
BNU_ESM_PPdif2_r=reshape(BNU_ESM_PPdif2,[m*n,1]);
CanESM2_PPdif2_r=reshape(CanESM2_PPdif2,[m*n,1]);
CNRM_CM5_PPdif2_r=reshape(CNRM_CM5_PPdif2,[m*n,1]);
CSIRO_Mk_PPdif2_r=reshape(CSIRO_Mk_PPdif2,[m*n,1]);
GFDL_CM3_PPdif2_r=reshape(GFDL_CM3_PPdif2,[m*n,1]);
GFDL_ESM2G_PPdif2_r=reshape(GFDL_ESM2G_PPdif2,[m*n,1]);
GFDL_ESM2M_PPdif2_r=reshape(GFDL_ESM2M_PPdif2,[m*n,1]);
HadGEM2_CC_PPdif2_r=reshape(HadGEM2_CC_PPdif2,[m*n,1]);
HadGEM2_ES_PPdif2_r=reshape(HadGEM2_ES_PPdif2,[m*n,1]);
inmcm4_PPdif2_r=reshape(inmcm4_PPdif2,[m*n,1]);
MIROC5_PPdif2_r=reshape(MIROC5_PPdif2,[m*n,1]);
MIROC_ESM_PPdif2_r=reshape(MIROC_ESM_PPdif2,[m*n,1]);
MIROC_ESM_CHEM_PPdif2_r=reshape(MIROC_ESM_CHEM_PPdif2,[m*n,1]);
MRI_CGCM3_PPdif2_r=reshape(MRI_CGCM3_PPdif2,[m*n,1]);
MRI_ESM1_PPdif2_r=reshape(MRI_ESM1_PPdif2,[m*n,1]);
NorESM1_M_PPdif2_r=reshape(NorESM1_M_PPdif2,[m*n,1]);
PPdif2_r=[BCC_CSM_PPdif2_r,BNU_ESM_PPdif2_r,CanESM2_PPdif2_r,CNRM_CM5_PPdif2_r,CSIRO_Mk_PPdif2_r,GFDL_CM3_PPdif2_r,GFDL_ESM2G_PPdif2_r,GFDL_ESM2M_PPdif2_r,HadGEM2_CC_PPdif2_r,HadGEM2_ES_PPdif2_r,inmcm4_PPdif2_r,MIROC5_PPdif2_r,MIROC_ESM_PPdif2_r,MIROC_ESM_CHEM_PPdif2_r,MRI_CGCM3_PPdif2_r,MRI_ESM1_PPdif2_r,NorESM1_M_PPdif2_r];
% The number of models with the same sign as averaged biaes in each pixel (P_kw, 2nd frequency range)
PPdif2_r_nan=isnan(PPdif2_r);
PPdif2_r_nan=single(PPdif2_r_nan);
PPdif2_r_nan_sum=sum(PPdif2_r_nan,2);
PPdif2_r_nan_sum_c=[PPdif2_r,PPdif2_r_nan_sum];
[row,~]=size(PPdif2_r);
PPdif2_count_0(row,1)=0;
PPdif2_count_1(row,1)=0;
count=1;
for i=1:row
    if PPdif2_r_nan_sum_c(i,18)==17
        num_0=nan;
        num_1=nan;
    else
        num_0=sum(PPdif2_r(i,:)==0,2);
        num_1=sum(PPdif2_r(i,:)~=0,2);
    end
    PPdif2_count_0(count,:)=num_0;
    PPdif2_count_1(count,:)=num_1;
    count=count+1;
    disp(i)
end
% Ratio of the number of models with the same sign as averaged biaes in each pixel to the number of total models (P_kw, 2nd frequency range)
PPdif2_count_0_real=PPdif2_count_0/17;
PPdif2_count_1_real=PPdif2_count_1/17;
PPdif2_positive_pct=reshape(PPdif2_count_0_real,[m,n]);
PPdif2_negative_pct=reshape(PPdif2_count_1_real,[m,n]);
% Prepare for significance test for averaged biases of P_kw (3rd frequency range)
BCC_CSM_PPdif3(BCC_CSM_PPdif3<0)=0;
BCC_CSM_PPdif3(BCC_CSM_PPdif3>0)=1;
BNU_ESM_PPdif3(BNU_ESM_PPdif3<0)=0;
BNU_ESM_PPdif3(BNU_ESM_PPdif3>0)=1;
CanESM2_PPdif3(CanESM2_PPdif3<0)=0;
CanESM2_PPdif3(CanESM2_PPdif3>0)=1;
CNRM_CM5_PPdif3(CNRM_CM5_PPdif3<0)=0;
CNRM_CM5_PPdif3(CNRM_CM5_PPdif3>0)=1;
CSIRO_Mk_PPdif3(CSIRO_Mk_PPdif3<0)=0;
CSIRO_Mk_PPdif3(CSIRO_Mk_PPdif3>0)=1;
GFDL_CM3_PPdif3(GFDL_CM3_PPdif3<0)=0;
GFDL_CM3_PPdif3(GFDL_CM3_PPdif3>0)=1;
GFDL_ESM2G_PPdif3(GFDL_ESM2G_PPdif3<0)=0;
GFDL_ESM2G_PPdif3(GFDL_ESM2G_PPdif3>0)=1;
GFDL_ESM2M_PPdif3(GFDL_ESM2M_PPdif3<0)=0;
GFDL_ESM2M_PPdif3(GFDL_ESM2M_PPdif3>0)=1;
HadGEM2_CC_PPdif3(HadGEM2_CC_PPdif3<0)=0;
HadGEM2_CC_PPdif3(HadGEM2_CC_PPdif3>0)=1;
HadGEM2_ES_PPdif3(HadGEM2_ES_PPdif3<0)=0;
HadGEM2_ES_PPdif3(HadGEM2_ES_PPdif3>0)=1;
inmcm4_PPdif3(inmcm4_PPdif3<0)=0;
inmcm4_PPdif3(inmcm4_PPdif3>0)=1;
MIROC5_PPdif3(MIROC5_PPdif3<0)=0;
MIROC5_PPdif3(MIROC5_PPdif3>0)=1;
MIROC_ESM_PPdif3(MIROC_ESM_PPdif3<0)=0;
MIROC_ESM_PPdif3(MIROC_ESM_PPdif3>0)=1;
MIROC_ESM_CHEM_PPdif3(MIROC_ESM_CHEM_PPdif3<0)=0;
MIROC_ESM_CHEM_PPdif3(MIROC_ESM_CHEM_PPdif3>0)=1;
MRI_CGCM3_PPdif3(MRI_CGCM3_PPdif3<0)=0;
MRI_CGCM3_PPdif3(MRI_CGCM3_PPdif3>0)=1;
MRI_ESM1_PPdif3(MRI_ESM1_PPdif3<0)=0;
MRI_ESM1_PPdif3(MRI_ESM1_PPdif3>0)=1;
NorESM1_M_PPdif3(NorESM1_M_PPdif3<0)=0;
NorESM1_M_PPdif3(NorESM1_M_PPdif3>0)=1;
[m,n]=size(BCC_CSM_PPdif3);
BCC_CSM_PPdif3_r=reshape(BCC_CSM_PPdif3,[m*n,1]);
BNU_ESM_PPdif3_r=reshape(BNU_ESM_PPdif3,[m*n,1]);
CanESM2_PPdif3_r=reshape(CanESM2_PPdif3,[m*n,1]);
CNRM_CM5_PPdif3_r=reshape(CNRM_CM5_PPdif3,[m*n,1]);
CSIRO_Mk_PPdif3_r=reshape(CSIRO_Mk_PPdif3,[m*n,1]);
GFDL_CM3_PPdif3_r=reshape(GFDL_CM3_PPdif3,[m*n,1]);
GFDL_ESM2G_PPdif3_r=reshape(GFDL_ESM2G_PPdif3,[m*n,1]);
GFDL_ESM2M_PPdif3_r=reshape(GFDL_ESM2M_PPdif3,[m*n,1]);
HadGEM2_CC_PPdif3_r=reshape(HadGEM2_CC_PPdif3,[m*n,1]);
HadGEM2_ES_PPdif3_r=reshape(HadGEM2_ES_PPdif3,[m*n,1]);
inmcm4_PPdif3_r=reshape(inmcm4_PPdif3,[m*n,1]);
MIROC5_PPdif3_r=reshape(MIROC5_PPdif3,[m*n,1]);
MIROC_ESM_PPdif3_r=reshape(MIROC_ESM_PPdif3,[m*n,1]);
MIROC_ESM_CHEM_PPdif3_r=reshape(MIROC_ESM_CHEM_PPdif3,[m*n,1]);
MRI_CGCM3_PPdif3_r=reshape(MRI_CGCM3_PPdif3,[m*n,1]);
MRI_ESM1_PPdif3_r=reshape(MRI_ESM1_PPdif3,[m*n,1]);
NorESM1_M_PPdif3_r=reshape(NorESM1_M_PPdif3,[m*n,1]);
PPdif3_r=[BCC_CSM_PPdif3_r,BNU_ESM_PPdif3_r,CanESM2_PPdif3_r,CNRM_CM5_PPdif3_r,CSIRO_Mk_PPdif3_r,GFDL_CM3_PPdif3_r,GFDL_ESM2G_PPdif3_r,GFDL_ESM2M_PPdif3_r,HadGEM2_CC_PPdif3_r,HadGEM2_ES_PPdif3_r,inmcm4_PPdif3_r,MIROC5_PPdif3_r,MIROC_ESM_PPdif3_r,MIROC_ESM_CHEM_PPdif3_r,MRI_CGCM3_PPdif3_r,MRI_ESM1_PPdif3_r,NorESM1_M_PPdif3_r];
% The number of models with the same sign as averaged biaes in each pixel (P_kw, 3rd frequency range)
PPdif3_r_nan=isnan(PPdif3_r);
PPdif3_r_nan=single(PPdif3_r_nan);
PPdif3_r_nan_sum=sum(PPdif3_r_nan,2);
PPdif3_r_nan_sum_c=[PPdif3_r,PPdif3_r_nan_sum];
[row,~]=size(PPdif3_r);
PPdif3_count_0(row,1)=0;
PPdif3_count_1(row,1)=0;
count=1;
for i=1:row
    if PPdif3_r_nan_sum_c(i,18)==17
        num_0=nan;
        num_1=nan;
    else
        num_0=sum(PPdif3_r(i,:)==0,2);
        num_1=sum(PPdif3_r(i,:)~=0,2);
    end
    PPdif3_count_0(count,:)=num_0;
    PPdif3_count_1(count,:)=num_1;
    count=count+1;
    disp(i)
end
% Ratio of the number of models with the same sign as averaged biaes in each pixel to the number of total models (P_kw, 3rd frequency range)
PPdif3_count_0_real=PPdif3_count_0/17;
PPdif3_count_1_real=PPdif3_count_1/17;
PPdif3_positive_pct=reshape(PPdif3_count_0_real,[m,n]);
PPdif3_negative_pct=reshape(PPdif3_count_1_real,[m,n]);
% Define the range of average biases of SSM_kw
SMPdif_max1=max(max(avg_SMPdif1));
SMPdif_max2=max(max(avg_SMPdif2));
SMPdif_max3=max(max(avg_SMPdif3));
SMPdif_maxn=[SMPdif_max1,SMPdif_max2,SMPdif_max3];
SMPdif_MAX=max(SMPdif_maxn);
SMPdif_min1=min(min(avg_SMPdif1));
SMPdif_min2=min(min(avg_SMPdif2));
SMPdif_min3=min(min(avg_SMPdif3));
SMPdif_minn=[SMPdif_min1,SMPdif_min2,SMPdif_min3];
SMPdif_MIN=min(SMPdif_minn);
% Change ET_kw and P_kw have the same architecture as SSM_kw
temp=~isnan(avg_SMPdif1);
temp=single(temp);
temp(temp==0)=nan;
avg_EPdif1_new=avg_EPdif1./temp;
avg_EPdif2_new=avg_EPdif2./temp;
avg_EPdif3_new=avg_EPdif3./temp;
avg_PPdif1_new=avg_PPdif1./temp;
avg_PPdif2_new=avg_PPdif2./temp;
avg_PPdif3_new=avg_PPdif3./temp;
% Remove regions with SSM less than 0.1 for SSM_kw
avg_SMPdif1(find(mean_SM<0.1))=-9;
avg_SMPdif2(find(mean_SM<0.1))=-9;
avg_SMPdif3(find(mean_SM<0.1))=-9;
SMPdif1_positive_pct(find(mean_SM<0.1))=nan;
SMPdif1_negative_pct(find(mean_SM<0.1))=nan;
SMPdif2_positive_pct(find(mean_SM<0.1))=nan;
SMPdif2_negative_pct(find(mean_SM<0.1))=nan;
SMPdif3_positive_pct(find(mean_SM<0.1))=nan;
SMPdif3_negative_pct(find(mean_SM<0.1))=nan;
% Remove N/A regions for SSM_kw
avg_SMP_m1(avg_SMP_m1>0)=1;
avg_SMP_m1(avg_SMP_m1<0)=-1;
avg_SMP_m1(find(isnan(avg_SMP_m1)==1))=0;
avg_SMP_m1=single(avg_SMP_m1);
avg_SMP_m2(avg_SMP_m2>0)=1;
avg_SMP_m2(avg_SMP_m2<0)=-1;
avg_SMP_m2(find(isnan(avg_SMP_m2)==1))=0;
avg_SMP_m2=single(avg_SMP_m2);
avg_SMP_m3(avg_SMP_m3>0)=1;
avg_SMP_m3(avg_SMP_m3<0)=-1;
avg_SMP_m3(find(isnan(avg_SMP_m3)==1))=0;
avg_SMP_m3=single(avg_SMP_m3);
SMP1(SMP1>0)=1;
SMP1(SMP1<0)=-1;
SMP1(find(isnan(SMP1)==1))=0;
SMP2(SMP2>0)=1;
SMP2(SMP2<0)=-1;
SMP2(find(isnan(SMP2)==1))=0;
SMP3(SMP3>0)=1;
SMP3(SMP3<0)=-1;
SMP3(find(isnan(SMP3)==1))=0;
avg_SMPdif1(find(avg_SMP_m1~=SMP1))=9;
avg_SMPdif2(find(avg_SMP_m2~=SMP2))=9;
avg_SMPdif3(find(avg_SMP_m1~=SMP3))=9;
% Remove N/A regions for ET_kw
avg_EP_m1(avg_EP_m1>0)=1;
avg_EP_m1(avg_EP_m1<0)=-1;
avg_EP_m1(find(isnan(avg_EP_m1)==1))=0;
avg_EP_m1=single(avg_EP_m1);
avg_EP_m2(avg_EP_m2>0)=1;
avg_EP_m2(avg_EP_m2<0)=-1;
avg_EP_m2(find(isnan(avg_EP_m2)==1))=0;
avg_EP_m2=single(avg_EP_m2);
avg_EP_m3(avg_EP_m3>0)=1;
avg_EP_m3(avg_EP_m3<0)=-1;
avg_EP_m3(find(isnan(avg_EP_m3)==1))=0;
avg_EP_m3=single(avg_EP_m3);
GLEAM_P1_new(GLEAM_P1_new>0)=1;
GLEAM_P1_new(GLEAM_P1_new<0)=-1;
GLEAM_P1_new(find(isnan(GLEAM_P1_new)==1))=0;
GLEAM_P2_new(GLEAM_P2_new>0)=1;
GLEAM_P2_new(GLEAM_P2_new<0)=-1;
GLEAM_P2_new(find(isnan(GLEAM_P2_new)==1))=0;
GLEAM_P3_new(GLEAM_P3_new>0)=1;
GLEAM_P3_new(GLEAM_P3_new<0)=-1;
GLEAM_P3_new(find(isnan(GLEAM_P3_new)==1))=0;
% Remove N/A regions for P_kw
avg_EPdif1_new(find(avg_EP_m1~=GLEAM_P1_new))=9;
avg_EPdif2_new(find(avg_EP_m2~=GLEAM_P2_new))=9;
avg_EPdif3_new(find(avg_EP_m1~=GLEAM_P3_new))=9;
avg_PP_m1(avg_PP_m1>0)=1;
avg_PP_m1(avg_PP_m1<0)=-1;
avg_PP_m1(find(isnan(avg_PP_m1)==1))=0;
avg_PP_m1=single(avg_PP_m1);
avg_PP_m2(avg_PP_m2>0)=1;
avg_PP_m2(avg_PP_m2<0)=-1;
avg_PP_m2(find(isnan(avg_PP_m2)==1))=0;
avg_PP_m2=single(avg_PP_m2);
avg_PP_m3(avg_PP_m3>0)=1;
avg_PP_m3(avg_PP_m3<0)=-1;
avg_PP_m3(find(isnan(avg_PP_m3)==1))=0;
avg_PP_m3=single(avg_PP_m3);
ERA5_P1_new(ERA5_P1_new>0)=1;
ERA5_P1_new(ERA5_P1_new<0)=-1;
ERA5_P1_new(find(isnan(ERA5_P1_new)==1))=0;
ERA5_P2_new(ERA5_P2_new>0)=1;
ERA5_P2_new(ERA5_P2_new<0)=-1;
ERA5_P2_new(find(isnan(ERA5_P2_new)==1))=0;
ERA5_P3_new(ERA5_P3_new>0)=1;
ERA5_P3_new(ERA5_P3_new<0)=-1;
ERA5_P3_new(find(isnan(ERA5_P3_new)==1))=0;
avg_PPdif1_new(find(avg_PP_m1~=ERA5_P1_new))=9;
avg_PPdif2_new(find(avg_PP_m2~=ERA5_P2_new))=9;
avg_PPdif3_new(find(avg_PP_m1~=ERA5_P3_new))=9;
% Significance test by stippling for averaged biases for SSM_kw
SM_mask11=(SMPdif1_positive_pct==1) | (SMPdif1_negative_pct==1);
SM_mask12=((SMPdif1_positive_pct>=14/17) & (SMPdif1_positive_pct<1)) | ((SMPdif1_negative_pct>=14/17) & (SMPdif1_negative_pct<1));
SM_mask21=(SMPdif2_positive_pct==1) | (SMPdif2_negative_pct==1);
SM_mask22=((SMPdif2_positive_pct>=14/17) & (SMPdif2_positive_pct<1)) | ((SMPdif2_negative_pct>=14/17) & (SMPdif2_negative_pct<1));
SM_mask31=(SMPdif3_positive_pct==1) | (SMPdif3_negative_pct==1);
SM_mask32=((SMPdif3_positive_pct>=14/17) & (SMPdif3_positive_pct<1)) | ((SMPdif3_negative_pct>=14/17) & (SMPdif3_negative_pct<1));
% Significance test by stippling for averaged biases for ET_kw
E_mask11=(EPdif1_positive_pct==1) | (EPdif1_negative_pct==1);
E_mask12=((EPdif1_positive_pct>=14/17) & (EPdif1_positive_pct<1)) | ((EPdif1_negative_pct>=14/17) & (EPdif1_negative_pct<1));
E_mask21=(EPdif2_positive_pct==1) | (EPdif2_negative_pct==1);
E_mask22=((EPdif2_positive_pct>=14/17) & (EPdif2_positive_pct<1)) | ((EPdif2_negative_pct>=14/17) & (EPdif2_negative_pct<1));
E_mask31=(EPdif3_positive_pct==1) | (EPdif3_negative_pct==1);
E_mask32=((EPdif3_positive_pct>=14/17) & (EPdif3_positive_pct<1)) | ((EPdif3_negative_pct>=14/17) & (EPdif3_negative_pct<1));
% Significance test by stippling for averaged biases for P_kw
P_mask11=(PPdif1_positive_pct==1) | (PPdif1_negative_pct==1);
P_mask12=((PPdif1_positive_pct>=14/17) & (PPdif1_positive_pct<1)) | ((PPdif1_negative_pct>=14/17) & (PPdif1_negative_pct<1));
P_mask21=(PPdif2_positive_pct==1) | (PPdif2_negative_pct==1);
P_mask22=((PPdif2_positive_pct>=14/17) & (PPdif2_positive_pct<1)) | ((PPdif2_negative_pct>=14/17) & (PPdif2_negative_pct<1));
P_mask31=(PPdif3_positive_pct==1) | (PPdif3_negative_pct==1);
P_mask32=((PPdif3_positive_pct>=14/17) & (PPdif3_positive_pct<1)) | ((PPdif3_negative_pct>=14/17) & (PPdif3_negative_pct<1));
% Global map of averaged biased for SSM_kw, ET_kw, and P_kw over different time scales (Figure 9)
detalgx=lonSMAP;
detalgy=latSMAP;
[LON,LAT]=meshgrid(detalgx,detalgy);
subplot(3,3,1)
map1=pcolor(detalgx,detalgy,avg_SMPdif1);
set(map1,'alphadata',~isnan(avg_SMPdif1))
shading interp
colorbar
caxis([-8,8])
hold on
SM_h11=stipple(LON,LAT,SM_mask11,'density',80,'color','k','marker','+','markersize',2);
SM_h12=stipple(LON,LAT,SM_mask12,'density',80,'color','k','marker','.','markersize',2);
hold off
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(a)')
ylabel('SSM_k_\Omega diff')
title('1/30<f<1/7 day^{-1}')
subplot(3,3,2)
map2=pcolor(detalgx,detalgy,avg_SMPdif2);
set(map2,'alphadata',~isnan(avg_SMPdif2))
shading interp
colorbar
caxis([-8,8])
hold on
SM_h21=stipple(LON,LAT,SM_mask21,'density',80,'color','k','marker','+','markersize',2);
SM_h22=stipple(LON,LAT,SM_mask22,'density',80,'color','k','marker','.','markersize',2);
hold off
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(b)')
title('1/90<f<1/30 day^{-1}')
subplot(3,3,3)
map3=pcolor(detalgx,detalgy,avg_SMPdif3);
set(map3,'alphadata',~isnan(avg_SMPdif3))
shading interp
colorbar
caxis([-8,8])
hold on
SM_h31=stipple(LON,LAT,SM_mask31,'density',80,'color','k','marker','+','markersize',2);
SM_h32=stipple(LON,LAT,SM_mask32,'density',80,'color','k','marker','.','markersize',2);
hold off
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(c)')
title('1/365<f<1/90 day^{-1}')
subplot(3,3,4)
map1=pcolor(detalgx,detalgy,avg_EPdif1_new);
set(map1,'alphadata',~isnan(avg_EPdif1_new))
shading interp
colorbar
caxis([-8,8])
hold on
E_h11=stipple(LON,LAT,E_mask11,'density',80,'color','k','marker','+','markersize',2);
E_h12=stipple(LON,LAT,E_mask12,'density',80,'color','k','marker','.','markersize',2);
hold off
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(d)')
ylabel('ET_k_\Omega diff')
subplot(3,3,5)
map2=pcolor(detalgx,detalgy,avg_EPdif2_new);
set(map2,'alphadata',~isnan(avg_EPdif2_new))
shading interp
colorbar
caxis([-8,8])
hold on
E_h21=stipple(LON,LAT,E_mask21,'density',80,'color','k','marker','+','markersize',2);
E_h22=stipple(LON,LAT,E_mask22,'density',80,'color','k','marker','.','markersize',2);
hold off
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(e)')
subplot(3,3,6)
map3=pcolor(detalgx,detalgy,avg_EPdif3_new);
set(map3,'alphadata',~isnan(avg_EPdif3_new))
shading interp
colorbar
caxis([-8,8])
hold on
E_h31=stipple(LON,LAT,E_mask31,'density',80,'color','k','marker','+','markersize',2);
E_h32=stipple(LON,LAT,E_mask32,'density',80,'color','k','marker','.','markersize',2);
hold off
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(f)')
subplot(3,3,7)
map4=pcolor(detalgx,detalgy,avg_PPdif1_new);
set(map4,'alphadata',~isnan(avg_PPdif1_new))
shading interp
colorbar
caxis([-8,8])
hold on
P_h11=stipple(LON,LAT,P_mask11,'density',80,'color','k','marker','+','markersize',2);
P_h12=stipple(LON,LAT,P_mask12,'density',80,'color','k','marker','.','markersize',2);
hold off
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(g)')
ylabel('P_k_\Omega diff')
subplot(3,3,8)
map5=pcolor(detalgx,detalgy,avg_PPdif2_new);
set(map5,'alphadata',~isnan(avg_PPdif2_new))
shading interp
colorbar
caxis([-8,8])
hold on
P_h21=stipple(LON,LAT,P_mask21,'density',80,'color','k','marker','+','markersize',2);
P_h22=stipple(LON,LAT,P_mask22,'density',80,'color','k','marker','.','markersize',2);
hold off
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(h)')
subplot(3,3,9)
map6=pcolor(detalgx,detalgy,avg_PPdif3_new);
set(map6,'alphadata',~isnan(avg_PPdif3_new))
shading interp
colorbar
caxis([-8,8])
hold on
P_h31=stipple(LON,LAT,P_mask31,'density',80,'color','k','marker','+','markersize',2);
P_h32=stipple(LON,LAT,P_mask32,'density',80,'color','k','marker','.','markersize',2);
hold off
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
xlabel('(i)')
axes('position',[0.81,0.14,0.15,0.76])
axis off
colorbar('eastoutside')
caxis([-8,8])
P_h3=legend([P_h31,P_h32],'100%','>80%','Location','SouthEast','FontWeight','bold','FontSize',6);
title(P_h3,'Significance Test')