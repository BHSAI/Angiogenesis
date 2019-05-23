function dYout = woundhealingequations2_inhibition(t,yin,Z,Param)

% Created on 05/29/2012
% Author: Sridevi Nagaraja
% Last updated on 01/04/2016
% Function with ODE description for the 16 variables in inflammation model
% This file describes the 27 ordinary differential equation and one delay differential equation used in the model to describe the 16 model variables. 
% The file also contains descriptions for the 11 dimensionless feedback loops of cytokine regulation as well as chemotaxis functions for TGF-beta,TNf-alpha, PDGF,CXCL8 and MIP-1alpha

% Assigning initial values; Comments show the identifying number assigned to each output variable
%Inflammation components
Nact  = yin(1);%5
Napop = yin(2);%6
Mpro  = yin(3);%7
Manti = yin(4);%8
tgf   = yin(5);%9
pdgf  = yin(6);%10
tnf   = yin(7);%11
IL1   = yin(8);%12
IL6   = yin(9);%13
IL10  = yin(10);%14
P     = yin(11);%15
IL8   = yin(12);%16
IL12  = yin(13);%17
MIP1  = yin(14);%18
MIP2  = yin(15);%19
IP10  = yin(16);%20
% proliferation components
F     = yin(17);%21
myoF  = yin(18);%22
fibnec= yin(19);%23
fgf2  = yin(20);%24
mmp9  = yin(21);%25
timp  = yin(22);%26
coll  = yin(23);%27
mmp1  = yin(24);%28
mmp2  = yin(25);%29
col1_fib = yin(26);%30
MCP1 = yin(27);%31
intermed = yin(28);%32
%Angiogenesis components
VEGF = yin(29);%33
EC = yin(30);%34
ANG1  = yin(31);%35
ANG2 = yin(32);%36
TSP1 = yin(33);%37
endo = yin(34);%38
%captip = yin(35);
capsprout = yin(35);%39
O = yin(36);%40
PEDF = yin(37);%41
K = yin(38);%42
KGF = yin(39);%43
CXCL1 = yin(40);%44
Itgf = yin(41);%45
I_tgf = yin(42);%46


% Outputs identifying number 1 and 2 are used to identify the total neutrophils and total macrophage concentrations calculated in the main
% file by adding the concentration of their respective individual phenotypes

%Assigning delay to TGF-beta
tgfdelay = Z(5,1);
tgfdelayf = Z(5,2);

%Mechanical stress corresponding to collagen fiber concentration and mechanical strain due to mechanical stress
Mstress = 3e-9*col1_fib;% 
Mstrain = 0.203*Mstress;%

%Intializing the output matrix
 dYout = zeros(36,1);

% Chemotaxis functions

chemotaxis_N = tgfcurve_N_Brandes(tgf)+pdgfcurve_N_Deuel(pdgf)+CXCL1curve_neutrophils(CXCL1);%+IL8_curve_N(IL8);%

chemotaxis_M = tgfcurve_M_Wahl_shifted(tgfdelay)+pdgfcurve_M_Deuel(pdgf)+MIP1_curve_M(MIP1)+TNF_curve_M(tnf)+ MCPcurve_M(MCP1);

chemotaxis_F = fgfcurve_F(fgf2)+fibneccurve_F(fibnec)+tnfcurve_F(tnf)+Param(99)*(pdgfcurve_F(pdgf)+tgfcurve_F(tgf));

chemotaxis_EC = VEGFcurve_EC(VEGF)+FGFcurve_EC(fgf2)+ANG2curve_EC(ANG2);

chemotaxis_K = VEGFcurve_K(VEGF);

% Negative feedback functions for cytokine (dimensionless)
%IL-10 mediated downregulation
fdn_TNF_IL10 = Param(47)*exp(Param(48)*IL10)+Param(49); % supression of TNF by IL-10

fdn_IL6_IL10 = Param(50)*exp(Param(51)*IL10)+Param(52);  %supression of IL-6 by IL-10

fdn_IL1_IL10 = Param(53)*exp(Param(54)*IL10)+Param(55);%supression of IL-1 by IL-10


%TGF mediated downregulation
fdn_TNF_tgf = Param(59)*exp(Param(60)*tgf)+Param(61); %supression of TNF by TGF

fdn_IL1_tgf = Param(62)*exp(Param(63)*tgf)+Param(64); %supression of IL-1 by TGF


%IL-6 mediated downregulation
fdn_IL1_IL6 = Param(65)/(Param(65)+IL6^Param(66)); %supression of IL-1 by IL-6

fdn_TNF_IL6 = Param(67)/(Param(67)+IL6^Param(68)); %supression of TNF by IL-6


% TNF mediated downregulation
fdn_IL12_TNF = Param(69)*exp(Param(70)*IL12)+Param(71);% supression of IL-12 by TNF


% Positive Feedback function for  cytokines (dimensionless)
%TGF-beta mediated upregulation
fup_IL6_tgf = 1+(Param(72)*tgf/(1+tgf));

fup_IL10_tgf = 1+(Param(73)*tgf/(1+Param(73)*tgf));

fup_ANG2_VEGF = 1+((Param(139)*VEGF)/(1+ VEGF^Param(140)));

%MCP mediated collagen upregulation by fibroblasts
if MCP1<400
    fup_coll_MCP1 = Param(107)*(MCP1^2)+Param(108)*(MCP1);
else
    fup_coll_MCP1 = 0;
end

% Regulation of EC proliferation by VEGF and endostatin
%f_up_VEGF_ECprolif = Param()/(1+((Param()/VEGF)*(1+(endo/Param()))));


amin= P>1e-10; % To avoid TGF-beta increase when platelet concentration value becomes very low
fib1 = fibnec>0; % To implement the requirement of Fibronectin for TGF-beta mediated fibroblast conversion to myofibroblast
o1 = O>0;
% 
% % %TGF degradation manipulation
%  tintervention_TGF = 24;
%  tinterventionend_TGF = 1008;
%  if t>tintervention_TGF && t<tinterventionend_TGF
%      Param(112) = Param(112)/2;
%  end
% % % % % %Macrophage efflux rate
%  tintervention_Mac =12;
%  tinterventionend_Mac =120;
%  if t>tintervention_Mac && t<tinterventionend_Mac
%      Param(10)= Param(10)/3;
%  end
% % % %  %Fibroblast apoptosis rate
%  tintervention_apop =12;
%  tinterventionend_apop = 120;
%  if t>tintervention_apop && t<tinterventionend_apop
%      Param(58) = Param(58)/3;
%  end
% % %  %Fibroblast migration rate
%  tintervention_Fib =96;
%  tinterventionend_Fib =360;
%  if t>tintervention_Fib && t<tinterventionend_Fib
%      Param(56) = Param(56)*3;
%  end
% %  %EC migration rate
%  tintervention_EC = 24;
%  tinterventionend_EC = 1008;
%  if t>tintervention_EC && t<tinterventionend_EC
%      Param(142) = Param(142)*2;
%  end
% Calculation of differentials
dPdt = -Param(1)*P;

dNactdt = Param(3)*amin*chemotaxis_N-Param(4)*Nact; %cells mL-1 hr-1

dNapopdt = Param(4)*Nact-((Param(5)*Napop)/(Param(6)+Napop))*Mpro;%cells mL-1 hr-1

dMprodt =Param(7)*amin*chemotaxis_M-((Param(8)*Napop)/(Param(9)+Napop))*Mpro-Param(10)*(1-(4*(Param(8)*Napop)/(Param(9)+Napop)))*Mpro;%cells mL-1 hr-1

dMantidt = ((Param(8)*Napop)/(Param(9)+Napop))*Mpro -Param(10)*Manti;%cells mL-1 hr-1

dFdt = Param(56)*chemotaxis_F+mechstrain_Fibprolif(Mstrain)*tgfcurve_F1(tgf)*Param(57)*F*(1-(Param(91)*(Mpro+Manti))-(Param(92)*F)-(Param(93)*col1_fib)-(Param(94)*myoF))-fib1*tgfcurve_myoF(tgf)*F-Param(58)*F;%tgfcurve_F1(tgf)*

dmyoFdt = fib1*tgfcurve_myoF(tgf)*F +Param(106)*myoF*(1-(Param(91)*(Mpro+Manti))-(Param(92)*F)-(Param(93)*col1_fib)-(Param(94)*myoF))-Param(74)*myoF;%

dcolldt = Param(95)*(1+fup_coll_MCP1)*tgfcurve_coll(tgf)*mechstrain_coll(Mstrain)*crowdingcurve_coll(col1_fib)*(fup_collfib_O(O)*F+2*fup_coll_O(O)*myoF)-Param(100)*coll+Param(96)*(mmp1+mmp2)*col1_fib;%

dintermeddt = Param(100)*coll-Param(101)*intermed;

dcoll1_fibdt = Param(101)*intermed-Param(96)*(mmp1+mmp2)*col1_fib-Param(97)*(F+2*myoF)*col1_fib;

dtgfdt = Param(2)*P+Param(11)*Mpro+Param(12)*Manti+Param(98)*mechstrain_tgfbeta(Mstrain)*(F+2*myoF)-Param(13)*tgf-Param(160)*Itgf*tgf+Param(161)*I_tgf; %ng mL-1 hr-1

dpdgfdt =Param(14)*Mpro+Param(15)*Manti-Param(16)*pdgf;%ng mL-1 hr-1

dtnfdt = Param(44)*Nact+Param(17)*fdn_TNF_IL10*fdn_TNF_tgf*fdn_TNF_IL6*Mpro+Param(18)*Manti- Param(19)*tnf;%ng mL-1 hr-1%%% influence of IL-10

dIL1dt = Param(45)*Nact+Param(20)*fdn_IL1_IL10*fdn_IL1_tgf*fdn_IL1_IL6*Mpro+Param(21)*Manti-Param(22)*IL1;%ng mL-1 hr-1

dIL6dt = Param(46)*Nact+Param(23)*fup_IL6_tgf*fdn_IL6_IL10*Mpro+Param(24)*Manti+Param(102)*mechstrain_IL6(Mstrain)*F+Param(145)*EC-Param(25)*IL6;%ng mL-1 hr-1

dIL10dt =Param(26)*fup_IL10_tgf*Mpro+Param(27)*Manti+Param(89)*F-Param(28)*IL10;%ng mL-1 hr-1

dIL8dt = Param(29)*Mpro+Param(30)*Manti+Param(90)* mechstrain_IL8(Mstrain)*F+fup_CXCL8_ANG1(ANG1)*Param(141)*EC-Param(31)*IL8;%ng mL-1 hr-1

dIL12dt =Param(32)*fdn_IL12_TNF*Mpro-Param(33)*IL12;%ng mL-1 hr-1

dMIP1dt =Param(34)*Mpro+Param(35)*Manti-Param(36)*MIP1;%ng mL-1 hr-1

dMIP2dt =Param(37)*Mpro+Param(38)*Manti-Param(39)*MIP2;%ng mL-1 hr-1

dIP10dt =Param(40)*Mpro+Param(41)*Manti-Param(42)*IP10;%ng mL-1 hr-1

dfibnecdt = Param(75)*Mpro+Param(76)*F-Param(77)*fibnec;%ng mL-1 hr-1

dFGF2dt =Param(78)*Mpro+Param(43)*Manti-Param(79)*fgf2+Param(162);%ng mL-1 hr-1

dMMP9dt = Param(80)*Mpro+Param(81)*F-Param(82)*mmp9;%ng mL-1 hr-1

dMMP1dt = Param(83)*Mpro-Param(84)*mmp1;% ng mL-1 hr-1 +1.4e-10*(F+2*myoF)*(coll>0)

dMMP2dt = Param(85)*Mpro-Param(86)*mmp2;% ng mL-1 hr-1

dTIMPdt =Param(103)*Mpro+Param(87)*F-Param(88)*timp; %ng mL-1 hr-1

dMCP1dt = Param(104)*Mpro-Param(105)*MCP1;%ng mL-1 hr-1

dVEGFdt = Param(109)*fup_VEGF_O2(O)*Mpro+Param(110)*fup_VEGF_O2(O)*Manti+fup_VEGF_O_fib(O)*Param(111)*F+Param(146)*EC+Param(158)*K-Param(112)*VEGF;%ng mL-1 hr-1-Param(109)*VEGF*EC+Param(109)*P+
 
dANG1dt = Param(113)*EC+Param(114)*F-Param(115)*ANG1;%ng mL-1 hr-1

dANG2dt = Param(116)*fdn_ANG2_IL6(IL6)*fdn_ANG2_TNF(tnf)*EC+Param(117)*F-Param(118)*ANG2+Param(163);%ng mL-1 hr-1*fup_ANG2_VEGF*fdn_ANG2_IL6(IL6)*fdn_ANG2_TNF(tnf)
 
dTSP1dt = Param(120)*Mpro+Param(121)*F+Param(122)*EC-Param(123)*TSP1;%ng mL-1 hr-1Param(119)*P+

dendodt = Param(124)*EC-Param(125)*endo;%ng mL-1 hr-1

dECdt = heaviside(Param(136)-EC)*(Param(142)*ECchemocurve_TSP1(TSP1)*chemotaxis_EC+Param(137)*fup_ECprolif_PEDF(PEDF)*ECprolifcurve_endo(endo)*fup_ECprolif_VEGF(VEGF)*EC)-Param(138)*fup_EC_PEDF(PEDF)*fup_ECapop_VEGF(VEGF)*EC;%  +Param(126)*capsprout

dcapsproutdt = heaviside(Param(129)-capsprout)*(((Param(127)*EC)+(Param(128)*capsprout))*EC+((Param(144)+Param(143)*VEGF)*EC))-Param(147)*fup_ECapop_VEGF(VEGF)*Param(138)*capsprout+Param(130)*capsprout*(1-(capsprout/Param(129)));% 

dOdt = Param(131)*capsprout*(1+Param(119))-((Param(134)*F)+(Param(133)*(Mpro+Manti)))*O-Param(132)*O;%ng mL-1 h-1

dPEDFdt = Param(148)*fup_PEDF_O(O)*F+Param(149)*fdn_PEDF_VEGF(VEGF)*K-Param(150)*PEDF;

dKdt = Param(151)*chemotaxis_K + fup_Kprolif_TGF(tgf)*fup_Kprolif_IL1(IL1)*Param(152)*(1-K/Param(153));

dKGFdt  = Param(154)*F-Param(155)*KGF;

dCXCL1dt =  Param(156)*Mpro+Param(159)*F-Param(157)*CXCL1;

dItgf = -Param(160)*Itgf*tgf+Param(161)*I_tgf;

dI_tgf = Param(160)*Itgf*tgf-Param(161)*I_tgf;



%Assiging output values in the differential matrix dYout
dYout(1)  = dNactdt;
dYout(2)  = dNapopdt;
dYout(3)  = dMprodt;
dYout(4)  = dMantidt;
dYout(5)  = dtgfdt;
dYout(6)  = dpdgfdt;
dYout(7)  = dtnfdt;
dYout(8)  = dIL1dt;
dYout(9)  = dIL6dt;
dYout(10) = dIL10dt;
dYout(11) = dPdt;
dYout(12) = dIL8dt;
dYout(13) = dIL12dt;
dYout(14) = dMIP1dt;
dYout(15) = dMIP2dt;
dYout(16) = dIP10dt;
dYout(17) = dFdt;
dYout(18) = dmyoFdt;
dYout(19) = dfibnecdt;
dYout(20) = dFGF2dt;
dYout(21) = dMMP9dt;
dYout(22) = dTIMPdt;
dYout(23) = dcolldt;
dYout(24) = dMMP1dt;
dYout(25) = dMMP2dt;
dYout(26) = dcoll1_fibdt;
dYout(27) = dMCP1dt;
dYout(28) = dintermeddt;
dYout(29) = dVEGFdt;
dYout(30) = dECdt;
dYout(31) = dANG1dt;
dYout(32) = dANG2dt;
dYout(33) = dTSP1dt;
dYout(34) = dendodt;
%dYout(35) = dcaptipdt;
dYout(35) = dcapsproutdt;
dYout(36) = dOdt;
dYout(37) = dPEDFdt;
dYout(38) = dKdt;
dYout(39) = dKGFdt;
dYout(40) = dCXCL1dt;
dYout(41) = dItgf;
dYout(42) = dI_tgf;

% Code to exit if simulation takes too long
MAXTIME = 60;
MINSTEP = 1e-12; %Minimum step

persistent tprev elapsedtime 
if isempty(tprev)
    tprev = -inf;
end
if isempty(elapsedtime)
    elapsedtime = tic;
end
timestep = t - tprev;
tprev = t;
if (t>0.01)&&(timestep>= 0) && (timestep < MINSTEP)
    error('Stopped. Time step is too small')
elseif (toc(elapsedtime) > MAXTIME)
    error('Stopped. Taking too long.')
end
clear elapsedtime tprev
%--------- Mechanical stress functions------------------------------------------- 
%Bellego 2009 allergy
% function x24 = mechstrain_MCP1(Mstrain)
% x24 = 1;
% if Mstrain<30
%     x24 = 0.0481*Mstrain+1;
% end

function x22 = mechstrain_IL8(Mstrain)
x22 = 1;
if Mstrain<30
    x22 = 0.0481*Mstrain+1;
end

function x21 = mechstrain_IL6(Mstrain)
x21 = 1;
if Mstrain<30
    x21 = 0.0429*Mstrain+1;
end

%Yang 2004 Jbiomech engg
function x20 = mechstrain_tgfbeta(Mstrain)
x20 = 1;
if Mstrain<10
    x20 = 0.08*Mstrain+0.9767;
end

function x19 = mechstrain_coll(Mstrain)
x19 = 1;
if Mstrain<10
    x19 = 0.0313*Mstrain+0.995;
end

function x18 = mechstrain_Fibprolif(Mstrain)
x18 = 1;
if Mstrain<10
    x18 = 0.0099*Mstrain+0.9928;
end
%---------------------------------------------------------------------------
%Function to calculate the effect of PEDF on EC proliferation
% Input: PEDF
%Output: Multiplier for EC proliferation rate

function x45= CXCL1curve_neutrophils(CXCL1)%Yao 2013 Vascular Biology
f45p1 = -4.2303;%-0.423;
f45p2 = 632.37;%63.237; 

x45= 0;
if CXCL1<100
    x45 = f45p1*CXCL1^2+f45p2*CXCL1;
end

%---------------------------------------------------------------------------
%Function to calculate the effect of PEDF on VEGF production 
% Input: PEDF
%Output: Multiplier for EC proliferation rate

function x44= fdn_PEDF_VEGF(VEGF)%Wang 2012 PLOS one
f44p1 = -0.0173;
f44p2 = 1; 

x44= 1;
if VEGF<40
    x44 = f44p1*VEGF+f44p2;
end
%----------------------------------------------------------------------------
%Function to calculate the effect of PEDF on EC proliferation
% Input: PEDF
%Output: Multiplier for EC proliferation rate

function x43= fup_ECprolif_PEDF(PEDF)%Wang 2012 PLOS one
f43p1 = -3e-5;
f43p2 = 1; 

x43= 1;
if PEDF<13600
    x43 = f43p1*PEDF+f43p2;
end
%---------------------------------------------------------------------------
%Function to calculate the effect of PEDF on Ec apoptosis
% Input: PEDF
%Output: Multiplier for EC apoptosis rate

function x42= fup_EC_PEDF(PEDF)%Wang 2012 PLOS one
f42p1 = -4e-7;
f42p2 = 0.0018;
f42p3 =1;

x42= 1;
if PEDF<3400
    x42 = f42p1*PEDF^2+f42p2*PEDF+f42p3;
end

%---------------------------------------------------------------------------
%Function to calculate the effect of hypoxia on PEDF prodcution by
%keratinocytes
% Input: oxygen
%Output: Multiplier for PEDF production rate by keratinocytes

function x41= fup_PEDF_O(O)%Wang 2012 PLOS one
f41p1 = -2e-8;
f41p2 = 0.0003;
f41p3 =0.1444;

x41 = 1;
if O<9198
    x41 = f41p1*O^2+f41p2*O+f41p3;
end

%---------------------------------------------------------------------------
%Function to calculate the effect of TGF-beta level on keratinocyte
% proliferation
% Input: IL-1beta
%Output: Multiplier for myofibroblast proliferation rate

function x40= fup_Kprolif_TGF(tgf)%Wang 2012 PLOS one
f40p1 = 0.0722;
f40p2 = 0.8333;

x40 = 1;
if tgf<40
    x40 = f40p1*tgf+f40p2;
end
%---------------------------------------------------------------------------
%Function to calculate the effect of IL-1beta level on keratonicyte
% proliferation
% Input: IL-1beta
%Output: Multiplier for myofibroblast proliferation rate

function x39 = fup_Kprolif_IL1(IL1)%Wang 2012 PLOS one
f39p1 = 0.0678;
f39p2 = 0.8556;

x39 = 1;
if IL1<40
    x39 = f39p1*IL1+f39p2;
end

%---------------------------------------------------------------------------
% %Function to calculate chemotaxis due to VEGF 
% %Input: VEGF concentration
% %Output: migrating keratinocytes

function x38 = VEGFcurve_K(VEGF)%Yongman 2006 mol med

 f38p1 = -0.1321;
 f38p2 =18.271;
 
 x38  =1;
 if VEGF<100
     x38 = f38p1*VEGF^2+f38p2*VEGF;
 end
%-----------------------------------------------------------------------------

% %Function to calculate chemotaxis due to VEGF 
% %Input: VEGF concentration
% %Output: migrating keratinocytes

% function x37 = fup_fibprolif_O(O)
% f37p1 = -0.0022;
% f37p2 = 20.662;
% 
% x37  =1;
% if O>876 && O<9000
%     x37 = f37p1*O+f37p2;
% end
%-------------------------------------------------------------------------
% %Function to calculate the effect of oxygen level on collagen production by
% %fibroblasts
% % Input: Oxygen concentration
% %Output: Multiplier for myofibroblast proliferation rate

function x36 = fup_collfib_O(O)%Siddiqui 1996 WRR
f36p1 = -6e-5;
f36p2 = 1.5539;

x36 = 1;
if O>876 && O<9000
    x36 = f36p1*O+f36p2;
end

% ---------------------------------------------------------------------------
% %Function to calculate the effect of oxygen level on VEGF production by
% %fibroblasts
% %Input: current oxygen concentration
% %Output: Multiplier for VEGF production rate

function x35  = fup_VEGF_O_fib(O)%van vilmmeren 2010 J apply Physiol
f35p1 = 8e-7;
f35p2 = -0.0035;
f35p3 = 4.4965;
x35 = 1;
if O<3066
    x35 = f35p1*O^2+f35p2*O+f35p3;
end

% ----------------------------------------------------------------------------
% %Function to calculate effect of oxygen on collagen synthesis from
% %myofibroblasts
% %Input: current oxygen concentration
% %Output: Multiplier for collagen synthesis rate

function x34 = fup_coll_O(O)%van vilmmeren 2010 J apply Physiol
% f34p1 = -0.0001;
% f34p2 = 1.811;

f34p1 = 1e-4;
f34p2 = 0.0688;
%x34 = 1;
% if O<9198
%     x34 =f34p1*O+f34p2;
%if O<10000
    x34 =f34p1*O+f34p2;
%end
%------------------------------------------------------------------------------ 
% %function to calculate effect of VEGF on EC proliferation 
% %Input: current VEGF concentration
% %Output: Multiplier for EC proliferation term

function x33 = fup_ECprolif_VEGF(VEGF)%Yoo 2005 Invest Opthalmol vis sci
f33p1 = -3e-5;
f33p2 = 0.0147;

x33 = 1;
if VEGF<500
    x33 = f33p1*VEGF^2+f33p2*VEGF+1;
end

% %---------------------------------------------------------------------------
% %function to calculate effect of oxygen level on monocyte production of VEGF 
% %Input: current oxygen concentration
% %Output: Multiplier for VEGF prodcution rate
% 
function x32 = fup_VEGF_O2(O)%Xiong 1998_ amJ pathol
% f32p1 = 6e-6;
% f32p2 = 0.0109;
% f32p3 = 6.0239;

f32p1 = -0.0014;
f32p2 = 1.3528;


x32 = 1;
%if O<1172
if O<250
%x32 = f32p1*O^2-f32p2*O+f32p3;
x32 = f32p1*O+f32p2;
end

%---------------------------------------------------------------------------
%function to calculate VEGF effect on EC apoptosis rate
%Input: current vEGF concentration
%Output: Multiplier for reduction in EC apoptosis rate

function x31 = fup_ECapop_VEGF(VEGF)%Gerber 1998 JBC
f31p1 = 0.002;
f31p2 = 0.0633;

x31 = 1;
if VEGF<30
x31 = f31p1*VEGF^2-f31p2*VEGF+1;
end

%---------------Chemotaxis functions--------------------------------------------------
%Function to calculate the positive feedback of ANG-1 on EC production of
%CXCL8

function x30 = fup_CXCL8_ANG1(ANG1)
f30p1 = 0.0062;

x30 = 1;
if ANG1<300
x30 =  f30p1*ANG1+1;
end

%------------------------------------------------------------------------------------
%Function to calculate the negative feedback of IL-6 alpha on ANG2
%production by EC
% Input: current IL-6 concentration
% Output: Reduction in EC production of ANG 2

function x29 = fdn_ANG2_IL6(IL6)%Orfans 2007 critical care med
f29p1 = 0.000817;
f29p2 = -8.163;
f29p3 = 0.1719;
f29p4 = -0.2825;

x29 = 1;
if IL6<150
    x29 = f29p1*exp(f29p2*IL6)+f29p3*exp(f29p4*IL6);
end
% -----------------------------------------------------------------------------
% %Function to calculate the negative feedback of TNF alpha on ANG2
% %production by EC
% % Input: current TNF-alpha concentration
% % Output: Reduction in EC production of ANG 2
% 
function x28 = fdn_ANG2_TNF(tnf)%Orfans 2007 critical care med
f28p1 = 0.000222;
f28p2 = -9.743;
f28p3 = 0.1529;
f28p4 = -0.1252;

x28 = 1;
if tnf<150
    x28 = f28p1*exp(f28p2*tnf)+f28p3*exp(f28p4*tnf);
end
%-----------------------------------------------------------------------------     
%Function to calculate the effect of TSP-1 on EC migration
% Input: current TSP-1 concentration
% Output: Reduction in EC chemotaxis
function x27 = ECchemocurve_TSP1(TSP1)%Short JCB 2005
f27p1 = -0.0003;
f27p2 = 0.8692;

x27 = 1;
if TSP1<1000
    x27 = f27p1*TSP1+f27p2;
end
% ----------------------------------------------------------------------------    
% %Function to calculate the downregulation  EC proliferation by endostatin 
% % Input: current endostatin concentration
% % Output: fraction of reduction in EC proliferation rate
% 
function x26 = ECprolifcurve_endo(endo)%Oreilly 1997 Cell

f26p1 = 2e-6;
f26p2 = 0.0023;
f26p3 = 1.0073;

x26=1;
if endo<1000
x26 = f26p1*endo^2-f26p2*endo+f26p3;
end
%-----------------------------------------------------------------------------     
% Function to calculate the effect of  VEGF on EC migration
% Input: current  VEGF concentration
% Output: Number of migrated endothelial cells
function x25 = VEGFcurve_EC(VEGF)%pascal 1999 JBC

f25p1 = -0.6704;%2.391e5;
f25p2 = 32.811; %0.1081;
% f25p3 = -7142;
% f25p4 = -3.554;

x25=0;
if VEGF<38
%x25 = f25p1*exp(f25p2*VEGF) + f25p3*exp(f25p4*VEGF);
x25 = f25p1*VEGF^2 + f25p2*VEGF;
end
% -----------------------------------------------------------------------------
% Function to calculate the effect of  FGF on EC migration
% Input: current  ANG concentration
% Output: Number of migrated endothelial cells
function x24 = ANG2curve_EC(ANG2)%Murdoch JI 2007

f24p1 = -0.06;%60.12
f24p2 = 7.03;%-0.0002908

% f24p1 = 60.12;
% f24p2 = -0.0002908;
% f24p3 = -33.16;
% f24p4 = -0.4761;

x24 =0;
if ANG2<100
%x24 = f24p1*exp(f24p2*ANG2) + f24p3*exp(f24p4*ANG2);
x24 = f24p1*ANG2^2 + f24p2*ANG2;
end
%-------------------------------------------------------------------------
function x23 = FGFcurve_EC(fgf2)%Grant et al 1992 Invest opthalmol and vis sci
% Function to calculate the effect of  FGF on EC migration
% Input: current  FGF concentration
% Output: Number of migrated endothelial cells

f23p1 =-0.002;%-5.893
f23p2 = 0.8119;%32.57
%f23p3 = 48.11;

x23 =0;
if fgf2<200
%x23 = f23p1*(fgf2^2)+f23p2*(fgf2)+f23p3;
x23 = f23p1*fgf2^2+f23p2*(fgf2);
end
%-------------------------------------------------------------------------------------
function x17 = tgfcurve_F1(tgf) %Jin et al
% Function to calculate the upregulating effect of TGF-beta on fibroblast proliferation
% Input: current TGF-beta concentration
% Output: Multiplication factor for fibroblast proliferation rate
   
  f17p1 =     -1.313;
  f17p2 =     9.42; 

x17=1;
if tgf<10
x17 = f17p1*(tgf^2)+f17p2*(tgf);
end

%--------------------------------------------------------------------------------------
function x16 = tgfcurve_coll(tgf)
% Function to calculate the upregulating effect of TGF-beta on collagen production
% Input: current TGF-beta concentration
% Output: Multiplication factor for collagen production rate by fibroblasts and myofibroblasts

 f16p1 =      0.0105;  %Grotendorst et al
 f16p2 =     -0.1671;  
 f16p3 =      0.6417 ; 
 f16p4 =      0.2722;  
 
x16 = 1;
if tgf<=10
x16 = f16p1*(tgf^3)+f16p2*(tgf^2)+f16p3*(tgf)+f16p4;
end

%-------------------------------------------------------------------------------------
function x15 = crowdingcurve_coll(coll)
% Function to calculate the effect of cellular crowding on collagen production
% Input: current tropocollagen concentration
% Output: Multiplication factor for collagen production rate by fibroblasts and myofibroblasts

coll = coll*1e-3; % convert to ug/mL

f15p1 = -4.33e-10;
f15p2 = 9e-7;
f15p3 = 0.00055;
f15p4 = 0.13;

x15=1;
if coll<=1250
x15 = f15p1*(coll^3)+f15p2*(coll^2)- f15p3*coll + f15p4;
end
 
%------------------------------------------------------------------------------------ 
function x14 = tnfcurve_F(tnf)
% Function to calculate  chemotactic activity of fibroblasts due to TNF-alpha
% Input: current TNF-alpha concentration
% Output: Number of migrated fibroblasts

tnf =tnf*1e3;

       f14p1 =        2.52;  
       f14p2 =  -0.0006002;  
       f14p3 =      -2.453;  
       f14p4 =    -0.00446;  

       x14 = 0;
       if tnf<4000
           x14 = (f14p1*exp(f14p2*tnf) + f14p3*exp(f14p4*tnf));     
       end   
       
%--------------------------------------------------------------------------------------       
function x13 = fibneccurve_F(fibnec)
% Function to calculate  chemotactic activity of fibroblasts due to fibronectin
% Input: current fibronectin concentration
% Output: Number of migrated fibroblasts

       f13p1 =       15.05;  
       f13p2 =  -1.391e-05;  
       f13p3 =      -15.65;  
       f13p4 =   -0.001837;  
       
       x13 =0;
       if fibnec>=25
       x13 = f13p1*exp(f13p2*fibnec) + f13p3*exp(f13p4*fibnec);
       end
       
%--------------------------------------------------------------------------------------          
function x12 = fgfcurve_F(fgf2)
% Function to calculate  chemotactic activity of fibroblasts due to FGF
% Input: current FGF concentration
% Output: Number of migrated fibroblasts

       f12p1 = -0.007759;	
       f12p2 =  0.8591;
       
 x12 = 0;
    if fgf2<=110
       x12 = f12p1*fgf2^2+f12p2*fgf2;	 
    end      
  
%--------------------------------------------------------------------------------------
function x11 = tgfcurve_myoF(tgf)
% Function to calculate differentiation of fibroblasts into myofibroblasts due to TGF
% Input: current TGF concentration
% Output: % of differentiated fibroblasts

       f11p1 =    0.003405  ;
       f11p2 =    0.005605 ;
       f11p3 =   -0.001993  ;
       f11p4 =      -4.582  ;
       
x11 = 0;
if tgf<=25  
x11 = (f11p1*exp(f11p2*tgf) + f11p3*exp(f11p4*tgf));
end

%---------------------------------------------------------------------------------------
function x10 = pdgfcurve_F(pdgf)
% Function to calculate  chemotactic activity of fibroblasts due to PDGF
% Input: current PDGF concentration
% Output: Number of migrated fibroblasts

      f10p1 =        5.56;
      f10p2 =   -0.001388;
      f10p3 =      -3.167;
      f10p4 =     -0.5014;
       
x10 = 0;
if pdgf<=300
    
x10 = (f10p1*exp(f10p2*pdgf) + f10p3*exp(f10p4*pdgf));

end

%-------------------------------------------------------------------------------------
function x9 = tgfcurve_F(tgf)
% Function to calculate chemotactic activity of fibroblasts due to TGF
% Input: current TGF concentration
% Output: Number of migrated fibroblasts

tgf = tgf*1e3;
       f9p1 =       110.9;  
       f9p2 =    -0.01954;  
       f9p3 =      -119.5;  
       f9p4 =     -0.1376;  
       
 x9=0;
  if tgf>=0.7 && tgf<=200
 x9 = (f9p1*exp(f9p2*tgf) + f9p3*exp(f9p4*tgf));
 end
 
%------------------------------------------------------------------------------------
function x8 = MCPcurve_M(MCP1)

% Function to calculate  chemotactic activity of monocytes due to MCP-1
% Input: Current MCP-1 concentration
% Output: Number of migrated macrophages
% Two curves are used to fit the data(linear and quadratic). 
% The curve shifts form quadratic to linear at MCP-1 concentration of 1 pg ml-1

% constants for linear section of TGF-beta curve
f8p1_linear =-0.058678;
f8p2_linear = 126.1;

f8p1_quad = -8.152;
f8p2_quad = 110;

mcp_transition = 12;%ng mL-1

x8=0;
if MCP1<mcp_transition
    x8 = f8p1_quad*MCP1^2+f8p2_quad*MCP1;%quadratic equation 
elseif MCP1>=mcp_transition 
    x8 = f8p1_linear*MCP1+f8p2_linear;% linear   
end

%-------------------------------------------------------------------------------------
function x7 = tgfcurve_N_Brandes(tgf)

% Function to calculate  chemotactic activity of neutrophils due to TGF-beta
% Input: Current TGF-beta concentration
% Output: Number of migrated neutrophils 
% Curve fitting of data from Brandes et al. 
% Two curves are used to fit the data from Brandes et al. (linear and
% quadratic). The curve shifts form quadratic to linear to TGF-betaconcentration of 0.958 pg ml-1

%Convert TGF-beta concentration in pg ml-1
tgf = tgf*1e3;%pg ml-1

% constants for quadratic section of TGF-beta curve
f7p1 =-537.3;
f7p2 = 698.05;

% constants for linear section of TGF-beta curve
f7p1_linear =-12.778;
f7p2_linear = 187.78;

tgf_transition = 0.958;%pg mL-1

x7=0;
if tgf<tgf_transition
    x7 = (f7p1*tgf^2+f7p2*tgf); % quadratic curve
elseif tgf>=tgf_transition && tgf<=10
x7 = (f7p1_linear*tgf+f7p2_linear);% linear   
end

%-------------------------------------------------------------------------%
function x6 = tgfcurve_M_Wahl_shifted(tgf)

% Function to calculate chemotactic activity of macrophages produced by TGF-beta
% Input: Current TGF-beta concentration
% Output: Number of migrated macrophages
% Curve fitting of data from Wahl et al. 
% Two curves are used to fit the data from Wahl et al. (linear and quadratic). 
% The curve shifts form quadratic to linear to TGF-beta concentration of 1 pg ml-1

%Convert TGF-beta concentration in pg ml-1
tgf = tgf*1e3;%pg mL-1

%Constants for quadratic section of TGF-beta curve
f6p1 =-240.29;%-240.29;
f6p2 = 293.93;%298.93;

%Constants for linear section of TGF-beta curve
f6p1_linear =-0.5926;
f6p2_linear = 60.593;

tgf_transition = 1;%pg mL-1

x6=0;
if tgf<tgf_transition
    x6 = (f6p1*tgf^2+f6p2*tgf); % quadratic curve
elseif tgf>=tgf_transition && tgf<=10
    x6 = (f6p1_linear*tgf+f6p2_linear);% linear   
end

%-------------------------------------------------------------------------%
function x5 = pdgfcurve_N_Deuel(pdgf)

% Function to calculate chemotactic activity of neutrophils due to PDGF
% Curve fitting on data from Deuel et al. 1991(p1x^2+p2x)
% Input: Current PDGF concentration
% Output: Number of migrated neutrophils 

% Quadratic coefficients for PDGF curve
f5p1 = -1.0204;%-6.1227; 
f5p2 = 5.9525;%35.71;

c5 = [f5p1 f5p2 0];
r5 = roots(c5);

x5 = 0;
if pdgf <=r5(2)
x5 = f5p1*pdgf^2+f5p2*pdgf; % cells quadratic curve
end

%-------------------------------------------------------------------------%
function x4 = pdgfcurve_M_Deuel(pdgf)

% Function to calculate chemotactic activity of macrophages produced by PDGF
% Curve fitting on data from Deuel 1982 J clinical invest(p1x^2+p2x)
% Input: Current PDGF concentration
% Output: Number of migrated macrophages

% Quadratic coefficients for PDGF curve
f4p1 =-0.3538;
f4p2 = 11.978;

c4 = [f4p1 f4p2 0];
r4 = roots(c4);

x4=0;
if pdgf <=r4(2)
x4 = f4p1*pdgf^2+f4p2*pdgf; % cells quadratic curve
end
%-------------------------------------------------------------------------%
function x3 = IL8_curve_N(IL8)

% Function to calculate  chemotactic activity of neutrophils due to CXCL8
% Curve fitting of data from Wang et al. 1993_JI (t1x^2+t2x)
% Input: Current CXCL8 concentration
% Output: Number of migrated neutrophils

% Quadratic coefficients for IL-8 curve
f3p1 =-0.1045;%
f3p2 = 13.678;%

c3 = [f3p1 f3p2 0];
r3 = roots(c3);
 
x3 = 0;
if IL8<=r3(2)
x3 = (f3p1*IL8^2+f3p2*IL8); % quadratic curve, zero intercept
end

%-------------------------------------------------------------------------%
function x2 = TNF_curve_M(tnf)

% Function to calculate  chemotactic activity of macrophages due to TNF-alpha
% Curve fitting of data from Pai et al.1996_J Am soc Neph (t1x^2+t2x)
% Input: Current TNF-alpha concentration
% Output: Number of migrated macrophages

% Quadratic coefficients for TNF-alpha curve
f2p1 = -0.3164;
f2p2 = 10.708;

c2 = [f2p1 f2p2 0];
r2 = roots(c2);

x2 = 0;
if tnf<=r2(2)
x2 = (f2p1*tnf^2+f2p2*tnf); % quadratic curve, zero intercept
end

%-------------------------------------------------------------------------%
function x1 = MIP1_curve_M(MIP1)

% Function to calculate  chemotactic activity of macrophages due to MIP-1alpha
% Curve fitting of data from Wang et al. 1993_JI (t1x^3+t2x^2+t3x)
% Input: Current MIP-1alpha concentration
% Output: Number of migrated macrophages

% Cubic coefficients 
f1p1 =0.0006;
f1p2 = -0.1481;
f1p3= 11.36;

x1 = (f1p1*MIP1^3+f1p2*MIP1^2+f1p3*MIP1); % cubic curve, zero intercept
