% Last modified:10/28/2013
% Author:Sridevi Nagaraja
% Main function of the cytokine inhibiton model. This routine initializes the parameter values and initial conditions (Parameters.m), defines the time span of the simulation, 
% calls the in-built MATLAB delay differential equation solver DDE23, saves the predicted values of all model variables and plots the results.

clear all
close all
clc

Parameters_inhibition; % Initializing all parameter values and initial conditions

% Defining the simulation time points 
tinit = 0; % Initial time
tday =42; % Maximum time (days)
tstep  = 1;% Time step (one hour)
tmax2 = tday*24;
delaytime = 12;
dosetime = [24];
delaytime = 12;
delaytime2 = 48;

% Param1  = [kdP kT_P kN_in kN_deact k1_ingest k2_ingest kM_in k1_ingest_tilde k2_ingest_tilde kd_M kTGF_pro kTGF_anti kd_TGF kPDGF_pro kPDGF_anti kd_PDGF kTNF_pro kTNF_anti kd_TNF kIL1_pro kIL1_anti kd_IL1 kIL6_pro kIL6_anti ...
%     kd_IL6 kIL10_pro kIL10_anti kd_IL10 kIL8_pro kIL8_anti kd_IL8 kIL12_pro kd_IL12 kMIP1_pro kMIP1_anti kd_MIP1 kMIP2_pro kMIP2_anti kd_MIP2 kIP10_pro kIP10_anti kd_IP10 k_FGF_anti kTNF_pro_N kIL1_pro_N kIL6_pro_N a1 b1 c1 a2 b2 c2 a3 b3 c3 ...
%     kF_in k_prolif_F k_apop_F a5 b5 c5 a6 b6 c6 a7 b7 a8 b8 a9 b9 c9 a10 a11 kapop_myoF kfibnec_pro kfibnec_F kd_fibnec kFGF_pro kd_FGF kMMP9_pro kMMP9_F kd_MMP9 kMMP1_pro kd_MMP1 kMMP2_pro kd_MMP2 kTIMP_F kd_TIMP kIL10_F kIL8_F...
%    M_max F_max Coll_max myof_max kprod_F kdeg_MMP kdeg_F kTGF_F k_chemotaxis_tgfpdgf k_polym_low k_polym_high kIL6_F kTIMP_pro kMCP1_pro kMCP1_res k_prolif_myoF a12 a13 kVEGF_pro kVEGF_anti kVEGF_F kdeg_VEGF kANG1_EC kANG1_F kdeg_ANG1...
%    kANG2_EC kANG2_F kdeg_ANG2 kTSP1_P kTSP1_pro kTSP1_F kTSP1_EC kdeg_TSP1 kendo_EC kdeg_endo ktip_growth kanasttips kanasttipsprout kmax_cap kremod koxy_bv kdeg_O kdegO_M kdegO_F O2max ECmax kprolif_EC kapop_EC a14 a15 kIL8_EC kin_EC ktip_VEGF ktipECdiff kIL6_EC kVEGF_EC kapop_ecsprout ...
%    kPEDF_mpro kPEDF_EC kd_PEDF kin_K kprolif_K kmax_K kKGF_F kdeg_KGF kCXCL1_Mpro kdeg_CXCL1 kVEGF_K kCXCL1_F kontgf ki_tgf k_add_FGF k_add_ANG2];

for k = 1:length(dosetime)
            k
tmax1 = dosetime(k); % h
tmax2 = tday*24;
tspan1 = (tinit:tstep:tmax1);
tspan2 = (tmax1:tstep:tmax2);

yinit1 = [Nact_init Napop_init Mpro_init Manti_init tgf_init pdgf_init tnf_init IL1_init IL6_init IL10_init P_init IL8_init IL12_init MIP1_init MIP2_init IP10_init F_init myoF_init fibnec_init fgf2_init mmp9_init timp_init coll_init mmp1_init mmp2_init col1_fib_init mcp1_init intermed_init...
  VEGF_init EC_init ANG1_init ANG2_fib_init TSP1_init endo_init capsprout_init O_init PEDF_init K_init KGF_init CXCL1_init Itgf_init I_tgf_init];

% Calling ode
sol1 = dde23(@woundhealingequations2_inhibition,[delaytime delaytime2],yinit,tspan1,[],Param);  
sol = deval(sol1,tspan1);
Ntot = sol(1,:)+sol(2,:); % Calculating total neutrophils 
Mtot = sol(3,:)+sol(4,:); % Calculating total macrophages 
Colltot = sol(23,:)+sol(26,:)+sol(28,:); % Calculating total collagen
bloodvesselden = sol(30,:)+sol(35,:);
Y = [tspan1; Ntot; Mtot; Colltot; sol(1,:); sol(2,:); sol(3,:); sol(4,:); sol(5,:); sol(6,:); sol(7,:); sol(8,:); sol(9,:); sol(10,:);sol(11,:); sol(12,:); sol(13,:); sol(14,:); sol(15,:); sol(16,:);sol(17,:);sol(18,:);sol(19,:);sol(20,:);sol(21,:);sol(22,:);sol(23,:);sol(24,:);sol(25,:);sol(26,:);sol(27,:);sol(28,:);sol(29,:);sol(30,:);sol(31,:);sol(32,:);sol(33,:);sol(34,:);sol(35,:);sol(36,:); sol(37,:);sol(38,:);sol(39,:);sol(40,:); sol(41,:);sol(42,:);bloodvesselden];
g1 = transpose(Y);% Saving the outputs
yinit2 = g1(end,5:46);
yinit2(32)=50; %Uncomment to introduce TGF-beta inhibitor
yinit2(41) =0; %Uncomment to introduce TNF-alpha inhibitor
clear sol1 sol Y

% Calling ode
sol1 = dde23(@woundhealingequations2_inhibition,[delaytime delaytime2],yinit2,tspan2,[],Param);
sol = deval(sol1,tspan2);
Ntot = sol(1,:)+sol(2,:);
Mtot = sol(3,:)+sol(4,:); 
Colltot = sol(23,:)+sol(26,:)+sol(28,:); % Calculating total collagen
bloodvesselden = sol(30,:)+sol(35,:);
Y = [tspan2; Ntot; Mtot; Colltot; sol(1,:); sol(2,:); sol(3,:); sol(4,:); sol(5,:); sol(6,:); sol(7,:); sol(8,:); sol(9,:); sol(10,:);sol(11,:); sol(12,:); sol(13,:); sol(14,:); sol(15,:); sol(16,:);sol(17,:);sol(18,:);sol(19,:);sol(20,:);sol(21,:);sol(22,:);sol(23,:);sol(24,:);sol(25,:);sol(26,:);sol(27,:);sol(28,:);sol(29,:);sol(30,:);sol(31,:);sol(32,:);sol(33,:);sol(34,:);sol(35,:);sol(36,:); sol(37,:);sol(38,:);sol(39,:);sol(40,:); sol(41,:);sol(42,:);bloodvesselden];
g2 = transpose(Y);
% yinit3 = g2(end,5:46);
% yinit3(20) = 20; 
F_FGF{k} = [g1(1:dosetime(k),:);g2];%Uncomment when TGF-beta inhibitor is introduced
%F_ANG2{k} = [g1(1:dosetime(k),:);g2]; %Uncomment when TNF-alpha inhibitor is introduced
%F_CXCL8{k} = [g1(1:dosetime(k),:);g2];%Uncomment when CXCL8 inhibitor is introduced
clear g1 g2 yinit2
end


