% This file initializes all the parameter values and initial conditions of the model.
% Last modified: 10/28/2013
% Author: Sridevi Nagaraja

% Assignment of parameter values as well as parameter identifying number (as comments)

% Platelets
kdP = 0.69;%0.69;%P1 
kT_P =1.25e-8;%P2 

%Neutrophils
kN_in = 3e2;%3e2; %P3 
kN_deact = 0.1; %P4
k1_ingest = 2; %P5 
k2_ingest = 4.71E6;%P6 

%Macrophages
kM_in =4e2;%4e2; %P7		
k1_ingest_tilde = 0.1; %P8 
k2_ingest_tilde = 1.0E6;% P9 
kd_M = 8.3e-3; %P10
				
%TGF
kTGF_pro = 1.88e-6; %P11 
kTGF_anti = 1.6e-8; %P12
kTGF_F = 1.6e-7; %P98 
kd_TGF = 0.693;  %P13 
 
%PDGF				 
kPDGF_pro = 6.00E-8; %P14
kPDGF_anti = 0.1*kPDGF_pro;%P15	
kd_PDGF	= 1.73E-1	; %P16

%TNF	
kTNF_pro_N = 1.66E-08;%*0.05;%P44
kTNF_pro = 3.46E-07;%P17
kTNF_anti = 4.29E-08;%0.05; %P18
kd_TNF = 0.5331;%P19

%IL-1
kIL1_pro_N = 1.70E-08;%P45
kIL1_pro = 1.23E-06;%P20
kIL1_anti = 2.45E-07;%P21
kd_IL1 = 0.1732;%P22		
				
%IL-6	
kIL6_pro_N = 8.30E-10;%P46
kIL6_pro = 1.18E-6;%P23
kIL6_anti = 0.1*kIL6_pro;%P24
kIL6_F = 2e-8;%P102 
kIL6_EC = 8.33e-7; %P145
kd_IL6 = 0.462;%P25

%IL-10				
kIL10_pro= 7.60E-08;%P26
kIL10_anti=1.55E-07;%P27
kIL10_F = 1.55e-7;%P89
kd_IL10= 1.93E-1;%P28

%IL-8
kIL8_pro = 2.5e-6;%P29
kIL8_anti = 6.94e-7;%P30
kIL8_F = 2.5e-7;%P90
kIL8_EC = 6.25e-7; %P141
kd_IL8= 0.693;%P31

%IL-12
kIL12_pro= 8.33E-07;%P32
kd_IL12= 0.05775;%P33

%MIP-1alpha
kMIP1_pro= 3.61E-08;%P34
kMIP1_anti= 4.25E-08;%P35
kd_MIP1= 0.385;%P36

%MIP-2
kMIP2_pro = 1.84E-07;%P37
kMIP2_anti =1.92E-07;%P38
kd_MIP2=0.2772 ;%P39

%IP-10
kIP10_pro = 7.50E-07;%P40
kIP10_anti =1.56E-06;%P41
kd_IP10 = 0.1732;%P42


% Feedback loop parameters
%fdn_TNF_IL10 
a1 = 0.4666; %P47
b1 = -1.528; %P48
c1 = 0.5332; %P49

%fdn_IL6_IL10 
a2 = 0.3298; %P50
b2 = -1.189; %P51
c2 = 0.6695; %P52

%fdn_IL1_IL10 
a3 = 0.6334; %P53
b3 = -1.794; %P54
c3 = 0.3667; %P55

%fdn_TNF_tgf 
a5 = 0.6211; %P59
b5 = -0.8305; %P60
c5 = 0.4466; %P61

%fdn_IL1_tgf 
a6 = 0.69; %P62
b6 = -20.37; %P63
c6 = 0.31; %P64

%fdn_IL1_IL6 
a7 = 4.459; %P65
b7 = 0.1571; %P66

%fdn_TNF_IL6 
a8 = 4.488; %P67
b8 = 0.1541; %P68

%fdn_IL12_TNF 
a9 = 0.8671; %P69
b9 = -2.794; %P70
c9 = 0.1307; %P71

%fup_IL6_tgf 
a10 = 0.9821; %P72

%fup_IL10_tgf 
a11 = 274.5; %P73

%fup_coll_MCP1
a12 = -2e-5;%P107
a13 = 0.0123;%P108

%fup_ANG2_VEGF
a14 = 5.804; %P139
a15 = 1.135; %P140

% Fibroblasts
kF_in = 100;%  P56
k_prolif_F = 0.0385; %P57 
k_apop_F =0.0375;%P58 

% Myofibroblasts
k_prolif_myoF = 0.5*k_prolif_F;%P106
kapop_myoF = 0.02;%P74 

% fibronectin
kfibnec_pro = 1e-5;%P75
kfibnec_F = 9.6515e-5;%P76 
kd_fibnec = 0.0330;%P77

% b-FGF
kFGF_pro = 1.30e-7;%P78
k_FGF_anti = 1.87e-7; %P43
kd_FGF = 0.0912;%P79

% MMP9
kMMP9_pro = 1.3e-3;%P80
kMMP9_F = 9.03e-6;%P81
kd_MMP9 = 0.0365;%P82

% MMP1
kMMP1_pro = 2e-8;%2e-3;%P83
kd_MMP1 = 0.3456;%P84

% MMP2
kMMP2_pro = 2e-7;%P85
kd_MMP2 = 0.3456; % P86

% TIMP-1
kTIMP_pro = 2.5e-6;%P103
kTIMP_F = 1.98E-05;%P87
kd_TIMP = 0.6300;%P88

% Cell density and crowding effect parameters
M_max = 0.002*1e-3;%P91
F_max = 0.0025*1e-3; %P92
Coll_max = 0.0004*1e-6;%P93 
myof_max = 2*F_max;%P94 

%Tropocollagen kinetics
kprod_F = 833;% P95 
kdeg_MMP = 2.18e-5;%P96 
kdeg_F = 6.25e-7;% P97 % 

k_chemotaxis_tgfpdgf = 0.1;%P99

%Collagen fiber kinetics
k_polym_low = 0.06;%P100
k_polym_high = 1.2;%P101

%MCP1
kMCP1_pro = 4.16e-8;%;%P104 
kMCP1_res = 2.081;%2.081;%P105

%VEGF
kVEGF_pro = 3e-9;%P109 2.4e-9
kVEGF_anti = 3.0e-7;%P110 *Change for impaired angiogenesis reduce 3-fold*
kVEGF_F = 8.75e-8;%*3.5;%P111  
kVEGF_EC = 5.2e-10; %P146 
kVEGF_K = 1.46e-7; %P158
kdeg_VEGF = 0.693;%P112 

%ANG-1
kANG1_EC = 1.25e-5;%P113
kANG1_F = 6.94e-10;%P114
kdeg_ANG1 = 0.0042;%P115

%ANG-2
kANG2_EC = 2.4e-5;%P116
kANG2_F = 5.5e-10;%P117
kdeg_ANG2 = 0.1667;%P118 

%TSP-1
kTSP1_P = 0;%P119 %used as HBOT parameter
kTSP1_pro = 6e-5;%P120
kTSP1_F = 6.54e-4;%P121
kTSP1_EC = 5.56e-5;%P122 
kdeg_TSP1 = 1.18;%1.18;%P123

%Endostatin
kendo_EC = 4.6e-9;%P124
kdeg_endo = 7e-3;%P125  1.25e-4**** CHANGED**


%Endothelial cells on capillary tips
kin_EC = 300;%P142
ECmax = 1e4; % cells/mL P136 

ktip_growth = 0.0345;%h-1 %P126 
kprolif_EC = 0.03; %P137 0.03
kapop_EC = 0.09;  %h-1%P138 change for impaired angiogenesis, multiple by 1.2


%capillary sprouts
kmax_cap = 1e5;% %P129 cell/ml 1e10
kremod = 0.0542;%/1.15; %h-1 P130 

ktip_VEGF = 0.0021;% h-1 P143 
ktipECdiff = 4.1e-2; %P144  4.1667e-05 *CHANGED**

kanasttips = 3.7e-8; %mL cell-1 h-1%P127  
kanasttipsprout = 3.7e-10;% mL cell-1 h-1%P128  

kapop_ecsprout = 0.5;%P147 %of EC apoptosis in vessels 


%Oxygen
koxy_bv = 2.8;% mmHg cell-1 h-1 P131   2.8e-4
kdeg_O = 1.3;%h-1 P132 1e-3
kdegO_M = 6.8e-6;%ml h-1 cell-1P133  
kdegO_F =4.5e-6;%ml h-1 cell-1 P134  
O2max = 5400; %ngmL-1%M P135 %NOT USED

%PEDF
kPEDF_mpro = 1.3e-7;%P148
kPEDF_EC = 2.5e-6;%P149
kd_PEDF = 1.386;%P150

%Keratinocytes
kin_K = 75; %P151
kprolif_K = 486.48;%P152;
kmax_K =1e4;%P153;

%KGF
kKGF_F = 3e-8;%P154
kdeg_KGF = 0.299;% P155;

%CXCL1
kCXCL1_Mpro = 6e-9;%*0.05; % P156
kCXCL1_F = 2.16e-7;%*0.05; %P159 **** changed from 2.16e-6
kdeg_CXCL1 = 0.2772;%P157


% Parameter vector 
Param  = [kdP kT_P kN_in kN_deact k1_ingest k2_ingest kM_in k1_ingest_tilde k2_ingest_tilde kd_M kTGF_pro kTGF_anti kd_TGF kPDGF_pro kPDGF_anti kd_PDGF kTNF_pro kTNF_anti kd_TNF kIL1_pro kIL1_anti kd_IL1 kIL6_pro kIL6_anti ...
    kd_IL6 kIL10_pro kIL10_anti kd_IL10 kIL8_pro kIL8_anti kd_IL8 kIL12_pro kd_IL12 kMIP1_pro kMIP1_anti kd_MIP1 kMIP2_pro kMIP2_anti kd_MIP2 kIP10_pro kIP10_anti kd_IP10 k_FGF_anti kTNF_pro_N kIL1_pro_N kIL6_pro_N a1 b1 c1 a2 b2 c2 a3 b3 c3 ...
    kF_in k_prolif_F k_apop_F a5 b5 c5 a6 b6 c6 a7 b7 a8 b8 a9 b9 c9 a10 a11 kapop_myoF kfibnec_pro kfibnec_F kd_fibnec kFGF_pro kd_FGF kMMP9_pro kMMP9_F kd_MMP9 kMMP1_pro kd_MMP1 kMMP2_pro kd_MMP2 kTIMP_F kd_TIMP kIL10_F kIL8_F...
   M_max F_max Coll_max myof_max kprod_F kdeg_MMP kdeg_F kTGF_F k_chemotaxis_tgfpdgf k_polym_low k_polym_high kIL6_F kTIMP_pro kMCP1_pro kMCP1_res k_prolif_myoF a12 a13 kVEGF_pro kVEGF_anti kVEGF_F kdeg_VEGF kANG1_EC kANG1_F kdeg_ANG1...
   kANG2_EC kANG2_F kdeg_ANG2 kTSP1_P kTSP1_pro kTSP1_F kTSP1_EC kdeg_TSP1 kendo_EC kdeg_endo ktip_growth kanasttips kanasttipsprout kmax_cap kremod koxy_bv kdeg_O kdegO_M kdegO_F O2max ECmax kprolif_EC kapop_EC a14 a15 kIL8_EC kin_EC ktip_VEGF ktipECdiff kIL6_EC kVEGF_EC kapop_ecsprout ...
   kPEDF_mpro kPEDF_EC kd_PEDF kin_K kprolif_K kmax_K kKGF_F kdeg_KGF kCXCL1_Mpro kdeg_CXCL1 kVEGF_K kCXCL1_F];


% Initial conditions
P_init =2e8;% Platelets/mL; % INITIAL PLATELET CONCENTRATION;*** CHANGE MODEL INPUT HERE *****
Nact_init =0;% cells mL-1
Napop_init = 0;% cells mL-1
Mpro_init = 0;%cells mL-1
Manti_init = 0;% cells mL-1
tgf_init =0;%ng mL-1 
pdgf_init =0;%ng mL-1
tnf_init = 0;%ng mL-1
IL1_init = 0;%ng mL-10
IL6_init= 0;%ng mL-1
IL10_init = 0;%ng mL-1
IL8_init = 0;%ng mL-1
IL12_init = 0;%ng mL-1
MIP1_init = 0;%ng mL-1
MIP2_init = 0;%ng mL-1
IP10_init = 0;%ng mL-1
F_init = 0;
myoF_init = 0;
fibnec_init = 0;%ng mL-1
fgf2_init = 0;%ng mL-1
mmp9_init = 0;%ng mL-1
mmp1_init = 0;
mmp2_init = 0;
timp_init = 0;%ng mL-1
coll_init = 0;%ng mL-1
mmp9coll_init = 0;
mmptimp_init = 0;
col1_fib_init = 0;
mcp1_init = 0;
intermed_init = 0;
VEGF_init = 0;
ANG1_init = 0;
ANG2_fib_init = 0;
TSP1_init = 0;
endo_init = 0;
captip_init = 0;
capsprout_init = 0;
O_init = 0;
EC_init = 0;
PEDF_init = 0;
K_init = 0;
KGF_init = 0;
CXCL1_init = 0;

% Initial values vector
yinit = [Nact_init Napop_init Mpro_init Manti_init tgf_init pdgf_init tnf_init IL1_init IL6_init IL10_init P_init IL8_init IL12_init MIP1_init MIP2_init IP10_init F_init myoF_init fibnec_init fgf2_init mmp9_init timp_init coll_init mmp1_init mmp2_init col1_fib_init mcp1_init intermed_init...
  VEGF_init EC_init ANG1_init ANG2_fib_init TSP1_init endo_init capsprout_init O_init PEDF_init K_init KGF_init CXCL1_init];
