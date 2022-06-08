
m1=0.54; m2=0.24; m3=0.2; m4=0.2; m5=0.3; m6=0.1; m7=0.1; m8=0.5; m9=0.2; m10=0.1; m11=1; m12=1; m13=0.32; m14=0.4; m15=0.7; 
m16=0.5; m17=0.5; m18=3.4; m19=0.9; m20=0.6; m21=0.08; m22=0.1; m23=0.1; m24=0.1; m25=0.9; m26=0.5; m27=0.1; m28=28; 
m29=0.3; m30=1; m31=0.1; m32=0.2; m33=13; m34=0.6; m35=0.3; m36=0.3; m37=0.4; m38=0.3; m39=0.2;

n1=2.6; n2=0.35; n3=0.29; n4=0.04; n5=0.4; n6=20; n7=0.1; n8=0.5; n9=0.6; n10=0.3; n11=0.6; n12=9; n13=2; 
n14=0.1; n15=2; n16=0.1; n17=0.5; n18=0.5; n19=0.2;

p1=0.13; p2=0.27; p3=0.1; p4=0.5; p5=1; p6=0.2; p7=0.3; p8=0.6; p9=0.8; p10=0.54; p11=0.5; p12=10; p13=0.1; p14=0.14; p15=2; p16=0.62; 
p17=17; p18=4; p19=1; p20=0.1; p21=1; p22=0.5; p23=0.37; p24=11; p25=2; p26=0.3; p27=0.8; p28=2; p29=0.1; p30=0.9; p31=0.1; p32=0.1;  p33=0.2;

g1=0.1; g2=0.01; g3=0.6; g4=0.005; g5=0.2; g6=0.3; g7=1; g8=0.04; g9=0.3; g10=0.5; g11=0.7; g12=0.1; g13=1; g14=0.02; g15=0.4; g16=0.3; 
g17=0.6; g18=0.4; g19=0.4; g20=0.03; g21=0.4; g22=0.1; g23=0.4; g24=0.3; g25=0.4; g26=0.3; g27=0.2; g28=0.1; g29=1;

a=2; b=2; c=2; d=2; e=2; f=2; g=2; h=2; i=2; j=2;
q1=1; q2=1.56; q3=3; 

A0=1; n16=0; m24=0.5; m23=0.5; g23=0.6; g25=0.5; g4=0.006; m6=0.2;

period=24;
dusk=12;
day_numb=3;

Loffset=0;
Lamplitude=1;

L=[];

t1=[0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36];
cLm6=[3.79 1.36 0.4 0.055 0.009275 0.004335 0.022135 0.025 0.155 1.4 4.81 5.49 4.735 1.17 0.395 0.065 0.005155 0.002805 0.0016465]/4.67;
cLm12=[4.67 1.9 0.69 0.1 0.01 0.0032 0.0025 0.001285 0.00579 0.065 0.715 2.785 4.285 1.48 0.665 0.135 0.009405 0.00225 0.00154]/4.67;
cLm18=[3.07 2.02 0.835 0.23 0.045 0.006315 0.00675 0.009195 0.03 0.115 0.09 0.595 2.5 2.44 0.93 0.285 0.045 0.00931 0.00487]/4.67;


t2=[0 4.0190 8.2687 12 16.26350935 20.299 23.773];
Cl_m=[1 0.1221 0.004 0.001 0.004 0.1386 1.0098]; % wt LD LHY Farre05
%Cl_m_79=[0.633634374 0.2475 0.2574 0.118983192 0.066 0.0726 0.227915116 0.633634374];

Cl_m_Kieron=[1 0.406852248 0.147751606 0.021413276 0.002141328 0.000685225 0.000535332 0.000275161 0.001239829 0.01391863 0.153104925 0.596359743 0.86509636]; % wt LD LHY Kieron
t3=[0 2 4 6 8 10 12 14 16 18 20 22 24];

%Cl_m_Megan=[0.330026161 0.287523197 0.170062831 0.098909046 0.046149859
%0.004802833 0.002650503 0.00887496 0.007361477 0.004616688 0.030730877 0.195102266 0.260808779];
%t3m=[0.25 1 2 4 7 9.75 11.75 12.33 13.5 15 18 21 23.75];

%Cl_p1=[0.474074074 0.518518519 1.185555498 0.681481481 0.342026559 0.266666667 0.297774092 0.237037037 0.15108206 0.133333333 0.15108206 0.312518861 0.62292729];
%t4=[0 1.41 3.76 5.64 8.46 11.28 12.22 13.16 15.06934637 17.86 20.68 23.97 24.44];
%Cl_p2=[0.62292729 0.814814815 0.711111111 0.342026559 0.079780219 0.122166093 0.044444444 0.103703704 0.755555556];
%t5=[0.44 1.85 4.67 8.9 11.72 12.20220573 13.6 16.42 24.88];

Cl_p_cca1=[0.31 0.72 1 0.74 0.54 0.336 0.18];
t55=[0 0.24 1.12 2.04 4.24 6 22];

%Ct_m_S=[0.096806814 0.123255233 0.430110785 0.99993022 0.75740988 0.57862019 0.478412376 0.378418307];
%t6=[1 4 7 10 13 16 19 22];

%Ct_m_M1=[0.212306955 0.135260076 0.388658699 0.999897271 0.672876074 0.477690648 0.296202445];
%Ct_m_M2=[0.296202445 0.130123617 0.452008355 0.782453857 0.600965654 0.409204534 0.236277095];
%Ct_m_M1_cl=[0.311611821 1.121460124 0.674588227 0.6283601 0.505085094 0.392083005 0.359552101];
%Ct_m_M2_cl=[0.359552101 1.308084786 0.763620176 0.648905934 0.475978495 0.337294114 0.294490292];
%t7=[0 4 8 12 16 20 24];

Cg_m=[0.49 0.49132948 1 0.034682081 0 0]; % wt LD GI
t8=[1.1 5 9.2 13.20340865 17.30115603 21];

Cp5_m=[0.020408163 0.624489796 1 0.628624435 0.057142857 0]; %wt LD PRR5
t9=[5.0641 7.0556 9.0471 11.0386 13.0301 15.0785];

%Cp5=[0.0082 0.01478276 0.0123 0.0451 0.4305 1.004533469 0.774997613 0.3854 0.1148 0.0082 0.0123 0.01478276];
%t10=[1.083 3.080110388 5.074280737 7.182904705 9.006 11.1155846 13.11 15.10596789 17.1 19.152 21.147 23.08528148];

%Cg=[0.139784946 0.279569892 0.752995329 1 0.483870968 0.172043011];
%t11=[1.085 5.27 9.3 13.02 17.205 21.235];

Cp5_m_M=[0.060606061 0.078787879 1 0.878787879 0.272727273 0.115787716 0.145958722 0.090909091 0.067759636]; % wt LD PRR5
t12=[0 3.36 6.24 9.12 12.24 15.12 18.12158933 21.12136359 24.00479952];

%Ct=[0.091 0.0455 0.30085101 1.001 1.0101 0.655452729];
%t13=[1.16 5.075 9.135 13.05 17.11245745 21.025];

Ct_m_K=[0.361904762 0.371428571 0.685714286 1 0.514285714 0.49047619]; % wt LD TOC1
t14=[0.962 4.884 8.954 12.95 16.946 20.942];

Cz=[0.157894737 0.184210526 0.450453757 1 0.738719413 0.526315789]; % wt LD ZTL
t15=[1.068 4.99669811 9.256 13.172 17.088 21.182];

%C9_m=[0.06 1.000449899 0.16 0.09 0.05 0.23 0.1 0.08];
%t16=[-1.644 1.294 7.17 10.334 22.54199969 28.64 34.29273929 48.07799996];
C9_mm=[0.06 1.000449899 0.16 0.09 0.05]; %wt LD PRR9
t16m=[-1.644 1.294 7.17 10.334 22.54199969];

C9_m=[1 0.55 0.15]; % Matsushika 00
t16=[1 3.3 14.4];

Cg_m_L=[0.037974684 0.177666694 0.348101266 0.297468354 0.265822785 1 0.670886076 0.456091601 0.076996994 0.040028831 0.056962025 0.044303797 0.044303797];
t17=[0 0.735391052 1.3 2.08 4.42 7.28 10.14 12 13.82 15.38 17.98 21.1 23.96];

%C9_m_D=[0.05 0.23 0.1 0.08];
%t18=[-1.458 4.64 10.293 24.08];

%C7_m=[0.113402062 0.463917526 1.00021254 0.845612233 0.432989691 0.350515464 0.092783505 0.030927835 0.092783505 0.546780644 0.649811732 0.403249633 0.382000151 0.216494845 0.20721393 0.248280301];
%t19=[-1.656 1.704 4.168 7.304 9.992 13.8 16.49314875 19.4 22.32091878 27.46729394 31.27492798 34.184 37.32 40.23722248 42.92491681 47.84997353];
%C7_mm=[0.113402062 0.463917526 1.00021254 0.845612233 0.432989691 0.350515464 0.092783505 0.030927835 0.092783505]; %wt LD PRR7
%t19m=[-1.656 1.704 4.168 7.304 9.992 13.8 16.49314875 19.4 22.32091878];

C7_m=[0.18 0.32 0.52 0.71 1 0.72 0.39 0.21 0.14 0.18];
t19n=[0 2.18 3.27 4.91 6.54 7.63 11.45 13.63 16.08 24];

%C7_m_D=[0.092783505 0.546780644 0.649811732 0.403249633 0.382000151 0.216494845 0.20721393 0.248280301];
%t20=[-1.679081217 3.46729394 7.27492798 10.184 13.32 16.23722248 18.92491681 23.84997353];


%C7=[0.064 0.312409987 1 0.85614952 0.793009458 0.224];
%t21=[0 4.342 8.183 12.358 16.366 20.207];

Ct_m_N=[0.186440678 0.169491525 0.279661017 1.000574383 0.974576271 0.644569385 0.56779661 0.610404851]; % wt LD TOC1
t22=[0 2.99 6.095 9.085 12.305 15.41 18.515 21.85];

%Ct_m_N1=[0.185339259 0.172413793 0.24137931 1.000264201 0.82790543 0.505747126 0.563218391 0.482758621 0.195402299];
%t23=[0 3.77 6.96 9.57 13.34 16.24 19.14 22.62743468 25.81];

%Cl_m_579=[0.4025 0.598 0.575459816 0.586950807 0.632918044 0.289562601 0.449089356 0.46];
%t24=[0 3.08 6.44 9.52 12.32 15.4 20.44 23.8];

%Cl_m_57=[0.2972859 0.2265625 0.1953125 0.296875 0.3359375 0.188959166 0.2578125 0.3671875];
%t25=[0 2.99 6.095 9.085 12.305 15.41 18.515 21.85];

%Cc_m_79=[1.3205 0.693760225 0.883959558 0.665 0.427922014 0.1805 0.57];
%t26=[0 4 8 12 12 16 20];

%Cl_mLL=[1 0.7115 0.1170 0.125 0.1154 0.3942 0.4712 0.5192 0.3365 0.1058 0.1827 0.1741 0.4712 0.5292 0.3077 0.1741 0.1837 0.5 0.4135 0.4912 0.2421];
%t27=[0 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80];

%Ct_mLL=[0.24 0.259 0.318 1.23 0.606 0.44 0.20 0.17 0.39 0.587 0.596 0.337 0.23 0.289 0.50 0.933 0.577 0.394  0.251 0.25 0.529];
%t28=[0 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80];

%Cl_sk2=[0.003316018 0.151467606 0.585378948 0.055234325 0.80849913 0.691075701 0.000164909 3.275441672 1.452356518 1.386231399 0.969255298];
%Cl_sk2=0.2/3.275*Cl_sk2;
%t29=[9.5 10.5 11.5 12.5 13.5 14.5 15.5 16.5 17.5 18.5 19.5];

%t30=[0.96 3.36 5.28 7.68 13.44 13.47424209 15.36 17.28 20.64 24 24.48 25.44 27.36 31.68 36 37.03784011 37.92 40.8 49.44931951];
%Cl_pox=[0.469 0.642 0.901 0.975 0.704 0.556 0.568 0.272 0.111 0.222 0.185 0.654 0.975 0.988 1 0.951 0.605 0.309 0.704];
%Cl_mox=[0.835826533 0.924050633 1 0.987341772 0.721518987 0.759493671 0.544303797 0.569620253 0.620253165 0.746835443 0.759493671 0.671363575 0.772566811 0.646065831 0.936708861 0.671959972 0.87378455 0.860759494 0.632911392];

%t31=[-1.052 1.064 3.272 5.296978696 7.872 9.712779647 12.012];
%C9m=[0.118199288 1.000328677 0.358974359 0.397435897 0.294871795 0.115384615 0.08974359];

%t32=[-1.052 3.273108394 5.48 7.504 11.92070768 12 14.116 16.324 22.948];
%C9p=[0.116703094 0.620769612 0.670886076 1 0.455696203 0.342708516 0.166492993 0.113924051 0.116703094];

t33=[1.22 3.22 5.22 7.22 9.22 11.22 13.23 15.23 17.23 19.23 21.23 23.23 25.23 27.23 29.24 31.24 33.24 35.24 37.24 39.24 41.24 43.25 45.25 47.25 49.25 51.25 53.25 55.25 57.26 59.26 61.26 63.26 65.26 67.26 69.26 71.27 73.27 75.27 77.27 79.27 81.27 83.27 85.28 87.28 89.28 91.28 93.28 95.28 97.28 99.28 101.29 103.29 105.29 107.29 109.29 111.29 113.29 115.30 117.30 119.30 121.30  123.30 125.30 127.30 129.31 131.31 133.31 135.31 137.31 139.31 141.31 143.32 145.32 147.32 149.32]; 
CTlc=[1.45E+07 2.02E+07 2.70E+07 3.07E+07 3.15E+07 3.09E+07 2.90E+07 2.58E+07 2.24E+07 1.93E+07  1.83E+07 2.03E+07 2.27E+07 2.49E+07 2.53E+07 2.39E+07 2.10E+07 1.81E+07 1.65E+07 1.69E+07 1.82E+07 1.95E+07 2.03E+07 1.99E+07 1.85E+07 1.67E+07 1.53E+07 1.52E+07 1.59E+07 1.67E+07 1.74E+07 1.76E+07 1.68E+07 1.55E+07 1.43E+07 1.39E+07 1.43E+07 1.49E+07 1.51E+07 1.50E+07 1.48E+07 1.39E+07 1.32E+07 1.28E+07 1.28E+07 1.30E+07 1.32E+07 1.32E+07 1.28E+07 1.27E+07 1.23E+07 1.21E+07 1.21E+07 1.21E+07 1.22E+07 1.21E+07 1.20E+07 1.18E+07 1.19E+07 1.18E+07 1.16E+07 1.17E+07 1.15E+07 1.13E+07 1.11E+07 1.11E+07 1.10E+07 1.09E+07 1.11E+07 1.12E+07 1.11E+07 1.11E+07 1.11E+07 1.09E+07 1.09E+07]; % TOC1:LUC in lhy/cca1 in LL; raw data 
CTlcg=[4.54E+06 5.71E+06 8.42E+06 1.08E+07 1.15E+07 1.13E+07 1.00E+07 8.64E+06 7.39E+06 6.38E+06 5.79E+06 5.64E+06 5.91E+06 6.16E+06 6.19E+06 5.98E+06 5.64E+06 5.29E+06 4.98E+06 4.74E+06 4.51E+06 4.36E+06 4.19E+06 4.08E+06 3.94E+06 3.80E+06 3.64E+06 3.50E+06 3.45E+06 3.40E+06 3.36E+06 3.38E+06 3.30E+06 3.20E+06 3.10E+06 3.02E+06 2.95E+06 2.92E+06 2.91E+06 2.89E+06 2.87E+06 2.84E+06 2.81E+06 2.78E+06 2.75E+06 2.68E+06 2.63E+06 2.58E+06 2.52E+06 2.54E+06 2.49E+06 2.46E+06 2.42E+06 2.38E+06 2.36E+06 2.31E+06 2.30E+06 2.26E+06 2.25E+06 2.24E+06 2.22E+06 2.22E+06 2.18E+06 2.16E+06 2.15E+06 2.13E+06 2.09E+06 2.10E+06 2.07E+06 2.08E+06 2.05E+06 2.05E+06 2.04E+06 2.05E+06 2.07E+06]; % TOC1:LUC in lhy/cca1/gi in LL; raw data 

t34=[0 4 4 8 8 12 12 16 16 20 20 24 24]; % LUX Hazen 2005
CLUX=[0.13 0.11 0.2 0.95 0.25 0.22 0.69 0.38 1 0.1 0.66 0.2 0.43];

t34=[]; % LUX Hazen 2005
CLUX=[];
t34=[0 4 8 12 16 20 24]; % LUX Hazen 2005
CLUX=[0.13 0.11 0.95 0.22 0.38 0.1 0.2];
CLUX2=[0.17 0.13 0.71 1 0.33 0.25 0.3]; % LUX Helfer 2011
C92=[0 1 0.25 0.04 0.02 0 0.03]; % PRR9 Helfer 2011

t35=[0 2 4 6 8 10 12 14 16 18 20 22];  % ELF3 Dixon 2011
CE3=[0.51 0.24 0.27 0.49 0.83 1.39 1.13 1.03 1.84 1.56 1.25 1.14];
CE3=CE3/3.;
CP7=[0.06 0.85 1.83 2.85 3.19 1.95 0.59 0.09 0.02 0.01 0.02 0.03];
CP7=CP7/3.19;
CP9=[0.16 3.73 3.88 2.73 0.7 0.11 0.04 0.01 0.01 0.02 0.04 0.09];
CP9=CP9/3.88;

n16_ar=[0 0.03 0.05 0.1 0.2 0.3];
per_n16=[24.2 24.9 25.5 25.9 27 34];

t = [0 day_numb*period];

%Solving the ODEs and draw traectories
options=odeset('stat','off');

    y0=[0.9548 0.956 0.0768 0.0206 0.5005 0.0656 0.0251 0.1502 0.0699 0.0137 0.0873 0.2505 0.1419 0.1554 0.0663 0.1811 0.0849 0.2527 0.4923 0.0811 0.856 1.3143 0.2893 0.1485 0.2234 0.8445 0.4068 0.0995 0.6628 0.4027 0.2362 0.2843 0.1342 0.4764 0.0319];
    %y0=[0.76 0.96 0.08 0.02 0.45 0.06 0.03 0.15 0.07 0.02 0.16 0.25 0.12 0.19 0.06 0.09 0.05 0.25 0.43 0.08 0.73 1.31 0.29 0.14 0.22 0.84 0.41 0.1 0.58 0.44 0.21 0.3 0.12 0.48 0.03]; % toc1-ox n16=0.1    
%    y0=[0.83 0.96 0.08 0.02 0.53 0 0.03 0.22 0.1 0.01 0 0.25 0.18 0.12 0.07 0.48 0.26 0.24 0.7 0.09 1.06 1.31 0.26 0.14 0.16 0.84 0.41 0.1 0.94 0.27 0.36 0.21 0.13 0.48 0.03]; %toc1
%    y0=[0 0.96 0.08 0.08 0 0.04 0.01 0 0 0.02 0.05 0.25 0.13 0.36 0.06 0 0 0.36 0.25 0 0 1.31 0.58 0.27 0.89 0.84 0.41 0.14 0.26 1 0 0.41 0 0.48 0.03]; %lhy/cca1
%    y0=[0 0.96 0.1 0.05 0 0 0.01 0 0 0.02 0 0.26 0.21 0.24 0.05 0 0 0.26 0.33 0 0 1.31 0.58 0.27 0.54 0.84 0.41 0.1 0.43 1 0 0.41 0 0.48 0.03]; %lhy/cca1/toc1
%    y0=[0.49 0.96 0.08 0.04 0.34 0.14 0 0.23 0.21 0.06 0.16 0.25 0.04 0.49 0 0 0 0.61 0.91 0.1 0.74 1.31 0.3 0.08 0.04 0.84 0.41 0.3 1.2 0.33 0.33 0.22 0.13 0.48 0.03]; %prr7/prr9
    
    
[T, Y] = ode23s(@Model_ABA_art_fin,t,y0,options,q1,q2,q3,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21,g22,g23,g24,g25,g26,g27,g28,g29,a,b,c,d,e,f,g,h,i,j,dusk,period,Loffset,Lamplitude,A0);

%T=T-24*3;

GN=p28*Y(:,10)./(p29+m19+p17*Y(:,25));
EGN=(p18*Y(:,4)+p17*Y(:,25).*GN)./(m10*Y(:,26)+m9*Y(:,27)+p22);
E34=p25*Y(:,19).*Y(:,25)./(p26*Y(:,29)+p21+m10*Y(:,26)+m9*Y(:,27));
AR=A0*Y(:,21)./(Y(:,21)+g29);
AR=0.5*(A0+Y(:,21)+g29-sqrt((A0+Y(:,21)+g29).^2-4*A0*Y(:,21)));

%pip=g26^2./(g26^2.+Y(:,31).*Y(:,31));
%mim=(n19/m40+n20/m40.*L);
%st=mim'.*pip;

figure (11)
%plot(t13,Ct,'LineStyle','none','Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);
%hold on;
%plot(t6,Ct_m_S,'LineStyle','none','Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);
%hold on;
%title('TOC1 mRNA Somers04 - black; Nikamuchi-blue+red; LHY prot-magenta; PRR5 prot-green');
plot(t14,Ct_m_K,'LineStyle','none','Marker','square','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',4);
hold on;
plot(t22,Ct_m_N,'LineStyle','none','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',4);
hold on;
%plot(t28,Ct_mLL,'LineStyle','none','Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);
%hold on;
%plot(t23,Ct_m_N1,'LineStyle','none','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',4);
%hold on;
plot(T,Y(:,11),'k');
hold on;
plot(T,Y(:,2),'r:');
hold on;
plot(T,Y(:,6),'b');
hold on;
plot(T,Y(:,5),'m');
hold on;
plot(T,Y(:,13),'g');
hold on;

title({'TOC1 mRNA-blue; TOC1 prot-black (mod-dot; X-dash-dot)'; 'EC-green'});

figure(2)
plot(T,Y(:,23),'b');
hold on;
plot(T,Y(:,11),'b:');
hold on;
plot(T,Y(:,24),'k:');
hold on;
plot(T,Y(:,25),'k-.');
hold on;
plot(T,Y(:,4),'k');     
hold on;
plot(T,EGN,'k--');     
hold on;
plot(T,Y(:,24)+Y(:,25)+Y(:,4)+Y(:,13)+EGN+E34,'r');     
hold on;
%plot(T,Y(:,5),'r:');
%hold on;
plot(T,Y(:,2),'r:');
hold on;
plot(T,Y(:,26),'g');
hold on;
plot(T,Y(:,22),'g:');
hold on;
plot(T,Y(:,27),'m');
hold on;
plot(T,Y(:,13),'c');
hold on;
plot(t35,CE3,'LineStyle','none','Marker','diamond','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',4), 
hold on;
title({'LHY prot-magenta; COP1-green (cyt-dot); CUL4-mag; ELF3 mRNA-blue (TOC1dot)';'ELF3 prot-black (cyt-dot; nucl-dash-dot; ELF3-GI cyt-solid;  ELF3-GInucl-dash); tot-red; EC-cyan'});

figure(3)
plot(t15,Cz,'LineStyle','none','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',4);
hold on;
plot(T,Y(:,12),'k');     
hold on;
plot(T,Y(:,3),'k:');     
hold on;
plot(T,Y(:,2),'r:');     
hold on;
plot(T,Y(:,12)+Y(:,3),'r');     
hold on;
plot(T,Y(:,14),'b');     
hold on;
plot(T,Y(:,10),'b:');     
hold on;
plot(T,Y(:,13),'g');     
hold on;
plot(T,Y(:,4),'g:');     
hold on;
title({'ZTL prot black (straight-free; dot-with GI); sum-red;';'GI-blue (mRNA, prot); EC-green (dot-ELF3 with GI)'});

figure(4)
plot(t16m,C9_mm,'LineStyle','none','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',8);
hold on;
plot(t16,C9_m,'LineStyle','none','Marker','diamond','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',8);
hold on;
plot(t34,C92,'LineStyle','none','Marker','diamond','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);
hold on;
%plot(t35,CP9,'LineStyle','none','Marker','diamond','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',8);
%hold on;
%plot(t31,C9m,'LineStyle','none','Marker','diamond','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',8);
%hold on;
%plot(t32,C9p,'LineStyle','none','Marker','pentagram','MarkerEdgeColor','r','MarkerFaceColor','k','MarkerSize',8);
%hold on;

%plot(t18,C9_m_D,'LineStyle','none','Marker','diamond','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
%hold on;
plot(T,Y(:,15),'b');     
hold on;
plot(T,Y(:,2),'r:');     
hold on;
plot(T,Y(:,5),'m');     
hold on;
%plot(T,Y(:,25),'g');     
%hold on;
plot(T,Y(:,7),'k');
hold on;
plot(T,Y(:,13),'g');
hold on;

title({'PRR9 mRNA-blue, proteins-black ; magenta-LHY, ELF3 prot.-green'; 'EC-green, X pr-red'});
%title('PRR9 mRNA-blue, proteins-black ; magenta-LHY, ELF3 prot.-green, X pr-red');

figure(44)
%plot(t19m,C7_mm,'LineStyle','none','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',8);
%hold on;
plot(t19n,C7_m,'LineStyle','none','Marker','diamond','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',8);
hold on;
%plot(t35,CP7,'LineStyle','none','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',8);
%hold on;
%plot(t20,C7_m_D,'LineStyle','none','Marker','diamond','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
%hold on;
plot(T,Y(:,16),'b');     
hold on;
plot(T,Y(:,5)+Y(:,20),'m');     
hold on;
plot(T,Y(:,7),'g');     
hold on;
plot(T,Y(:,2),'r:');     
hold on;
%plot(t21,C7,'LineStyle','none','Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);
%hold on;
plot(T,Y(:,17),'k');     
hold on;
title('PRR7 mRNA-blue (Nakamichi03); prot-black; magenta-LHY+LHYm; green-PRR9');

figure(5)
%plot(t8,Cg_m,'LineStyle','none','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',4);
%hold on;
%plot(t17,Cg_m_L,'LineStyle','none','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',4);
%hold on;
plot(T,Y(:,14),'b');     
hold on;
plot(T,Y(:,2),'r:');     
hold on;
plot(T,Y(:,5),'m');     
hold on;
%plot(T,Y(:,20),'m:');     
%hold on;
plot(T,Y(:,13),'g');
hold on;
plot(T,Y(:,10),'k:');     
hold on;
%plot(T,Y(:,3),'k--');     
%hold on;
plot(T,Y(:,4),'r--');     
hold on;
%plot(T,GN,'k');     
%hold on;
plot(T,Y(:,10)+Y(:,3)+Y(:,4)+GN+EGN,'r');     
hold on;
plot(T,EGN,'r:');     
hold on;
%plot(t11,Cg,'LineStyle','none','Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);
%hold on;
title({'GI mRNA-blue; LHY prot-magenta; EC-green);GI prot black (cyt-dot; nucl-solid;';'ZTL-GI-dash);GI tot-red (ELF3-GI cyt-dash; ELF3-GInucl-dot)'});

figure(6)
plot(t9,Cp5_m,'LineStyle','none','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',4);
hold on;
plot(t12,Cp5_m_M,'LineStyle','none','Marker','diamond','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',4);
hold on;
plot(T,Y(:,8),'b');     
hold on;
plot(T,Y(:,9),'r');     
hold on;
plot(T,Y(:,5),'m:');     
hold on;
plot(T,Y(:,17),'g');     
hold on;
%plot(t10,Cp5,'LineStyle','none','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','k','MarkerSize',4);
%hold on;
title('NA mRNA-blue; NA prot-red; PRR7-green; LHY prot-magenta');


figure(7)
%plot(t30,Cl_pox,'LineStyle','none','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',4);
%hold on;
%plot(t30,Cl_mox,'LineStyle','none','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',4);
%hold on;

%plot(t4,Cl_p1,'LineStyle','none','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',4);
%hold on;
%plot(t5,Cl_p2,'LineStyle','none','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',4);
%hold on;
%plot(T,Y(:,5)+Y(:,20),'b');     
%hold on;
%plot(T,Y(:,2),'r:');
%hold on;
%title('LHY mRNA');
%axis([0 24 0 1.2]);
%plot(t2,Cl_m,'LineStyle','none','Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);
%hold on;
%plot(t2,Cl_m_79,'LineStyle','none','Marker','square','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',4);
%hold on;
%plot(t26,Cc_m_79,'LineStyle','none','Marker','square','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',4);
%hold on;
%plot(t24,Cl_m_579,'LineStyle','none','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',4);
%hold on;
%plot(t25,Cl_m_57,'LineStyle','none','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',4);
%hold on;
%plot(t27,Cl_mLL,'LineStyle','none','Marker','diamond','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',4);
%hold on;
plot(T,Y(:,8),'b');
hold on;
plot(T,Y(:,9),'b:');
hold on;
plot(T,Y(:,16),'g');
hold on;
plot(T,Y(:,17),'g:');     
hold on;
plot(T,Y(:,15),'m');
hold on;
plot(T,Y(:,7),'m:');
hold on;
%plot(T,Y(:,5),'r:');     
%hold on;
%plot(T,Y(:,20),'r--');
%hold on;
%plot(T,Y(:,5)+Y(:,20),'r');     
%hold on;
%plot(T,Y(:,18),'g');
%hold on;
%plot(T,Y(:,13),'c');
%hold on;
plot(T,Y(:,6),'r');
hold on;
plot(T,Y(:,11),'r:');
hold on;
%plot(T,Y(:,21),'k:');
%hold on;
%plot(T,L(:),'r');
%hold on;

%plot(T,n1*Y(:,18).^b./(Y(:,18).^b+g2^b),'g');
%hold on;
%plot(T,n0,'g--');
%hold on;

%plot(t3,Cl_m_Kieron,'LineStyle','none','Marker','diamond','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);
%hold on;
%plot(t3m,Cl_m_Megan./0.33003,'LineStyle','none','Marker','pentagram','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);
%hold on;
%plot(t1,cLm6,'LineStyle','none','Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8), 
%hold on, 
%plot(t1,cLm18,'LineStyle','none','Marker','diamond','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8), 
%hold on;
plot(T,Y(:,1),'k');
hold on;
plot(T,Y(:,13),'g');
hold on;

%title({'LHY mRNA black (squares for Farre05, diamond - for Kieron)';'dot: green-PRR7; magenta-PRR9; NI-black'});title({'LHY mRNA black (squares for Farre05, diamond - for Kieron)';'dot: green-PRR7; magenta-PRR9; PRR5-black'});
title('LHY mRNA-black, PRR9-magenta; PRR7-green; NI-blue; EC-cyan');

figure (9)
plot(T,Y(:,28),'b');
hold on;
%plot(T,Y(:,29),'k');
%hold on;
%plot(t34,CLUX,'LineStyle','none','Marker','diamond','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4), 
%hold on;
%plot(t34,CLUX2,'LineStyle','none','Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4), 
%hold on;
plot(T,Y(:,18),'k');
hold on;
%plot(T,Y(:,19),'k:');
%hold on;
%plot(T,E34,'k--');
%hold on;
plot(T,Y(:,13),'g');
hold on;
%plot(T,Y(:,13)+Y(:,29),'r');
%hold on;
%plot(T,Y(:,13)+Y(:,19)+E34,'m');
%hold on;
%plot(T,Y(:,2),'r:');
%hold on;
plot(T,Y(:,5),'m');
hold on;
%title({'LUX mRNA-blue; ELF4 mRNA-blue dot; LUX prot-black; ELF4-dot(with ELF3-dash)';' EC-green; LUXtot-red; ELF4tot-mag'});

figure (1)
plot(T,Y(:,21),'k:');
hold on;
plot(T,Y(:,5),'m:');
hold on;
plot(T,Y(:,30),'b');
hold on;
plot(T,Y(:,6),'m');
hold on;
plot(T,Y(:,11),'c');
hold on;
plot(T,Y(:,2),'r:');
hold on;
plot(T,Y(:,31),'r');
hold on;
plot(T,Y(:,33),'r:');
hold on;
plot(T,Y(:,32),'g');
hold on;
plot(T,AR,'k');
hold on;
title('AR black (ABAR mRNA-dot); PP2C-blue; S-red (mRNA-dot); LHY prot-magenda; TOC1 prot-cyan; green-stomata');




