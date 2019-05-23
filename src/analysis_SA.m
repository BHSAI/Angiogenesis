% Thsi routine calculates fold changes in output variables with respect to simulation
% with default parameter set. Please load the default output variable "Y"
% from "outputvars.mat"

for i = 1:size(Y_main_global,2)
    F = Y_main_global{i};
    F1 = F./Y;
    foldchanges{i} = F1;
end

for i=1:40
[AT(i) PH(i) RI(i) RP(i)] = timescalculation(Y(i,:),tspan);
end

% Calculate the fold changes in each index with respect to its default
% value
AT_foldchanges = ATf./AT;
PH_foldchanges = PHf./PH;
RI_foldchanges = DTf./RI;
RP_foldchanges = RPf./RP;


%Outputs of interest
P = [29 30 35 36 38]; %Output numbers int he code
Pname = ['VEGF', 'Endothelial cells', 'Blood vessel density','Oxygen','Keratinocytes']; %Out put names


% Filter the fold change values in output of interest at final simulation time-point
for i=1:size(foldchanges,2)
    S = foldchanges{i};
    S2 = S(P,end);
    foldchanges_fiveoutputs(:,i) = S2;
end

PH_foldchanges_fiveoutputs = PH_foldchanges(:,P);
PH_foldchanges_fiveoutputs = PH_foldchanges_fiveoutputs';


%Separation of wound healing simulations based on peak height of model outputs
impairedcases = find(PH_foldchanges_fiveoutputs(1,:)<=0.2 & PH_foldchanges_fiveoutputs(2,:)<=0.2);%impaired angiogenesis
normalcases = find(PH_foldchanges_fiveoutputs(1,:)>=0.8&PH_foldchanges_fiveoutputs(1,:)<=1.5&PH_foldchanges_fiveoutputs(2,:)>=0.8&PH_foldchanges_fiveoutputs(2,:)<=1.5);% normal angiogenesis


%Separation of parameter sets into normal and impaired sets
Parametersets_impaired = Param1(impairedcases,:);
Parametersets_normal = Param1(normalcases,:);

%Separation of PH values into normal and impaired sets
PH_values_normal = PH_foldchanges(normalcases,P);
PH_values_impaired = PH_foldchanges(impairedcases,P);


%Second filter
% impaired angiogenesis with model-predicted mechanism disruption
impairedcases_parameter = find(Parametersets_impaired(:,110)<=2e-7 & Parametersets_impaired(:,138)>=0.099);
%normal angiogenesis cases with model-predicted mechanism disruption
normalcases_parameters = find(Parametersets_normal(:,110)>2e-7 & Parametersets_normal(:,138)<0.099);

%Separation of parameter sets into normal and impaired sets after second
%filter
Parametersets_impaired_secondfilter =Parametersets_impaired(impairedcases_parameter,:);
Parametersets_normal_secondfilter = Parametersets_normal(normalcases_parameters,:);

%Separation of PH values into normal and impaired sets after second filter
PH_values_normal_secondfilter = PH_values_normal(normalcases_parameters,:);
PH_values_impaired_secondfilter = PH_values_impaired(impairedcases_parameter,:);

%----------------------------------------
% [p1 p2] = find(errornumber(:,1)==0);
% k=1;
% for i=1:size(p1,1)
% G{k} = F_Paramseti{p1(i)};
% k=k+1;
% end
% for i=1:size(G,2)
%     F = G{i};
%     D = F{1};
%     VEGF(:,i) = D(:,33);
%     EC(:,i) = D(:,34);
% end
% Mean_VEGF = mean(VEGF');
% std_VEGF = std(VEGF');
% Mean_EC = mean(EC');
% std_EC = std(EC');
% 
% Mean_VEGF = Mean_VEGF';
% std_VEGF = std_VEGF';
% Mean_EC = Mean_EC';
% std_EC = std_EC';

% Mean and standard deviation in normal and impaired model outputs
Mean_normal = mean(PH_values_normal_secondfilter);
std_normal = std(PH_values_normal_secondfilter);
Mean_impaired = mean(PH_values_impaired_secondfilter);
std_impaired = std(PH_values_impaired_secondfilter);

figure
hold on
bar(1:5,Mean_normal)
errorbar(1:5,Mean_normal,std_normal,'.')
hold on
bar(1:5, Mean_impaired)
errorbar(1:5,Mean_impaired,std_impaired,'.')

% This part of the routine performs the knockout of one protein at a time by assiging the
% value zero to all the production and degradaiton rates of that protein.
% It then simulates the time courses of all the model output variables
% after the protein knockout

%TGF KO
P_TGF = [11 12 98 13];
Parametersets_impaired_secondfilter_TGFKO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_TGFKO(:,P_TGF) = 0;
[ATf_TGF PHf_TGF DTf_TGF RPf_TGF Y_output_TGF] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_TGFKO = PHf_TGF(:,P);
Mean_impaired_TGFKO = mean(PH_values_impaired_TGFKO);
std_impaired_TGFKO = std(PH_values_impaired_TGFKO);

%PDGF KO
P_PDGF = [14 15 16];
Parametersets_impaired_secondfilter_PDGFKO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_PDGFKO(:,P_PDGF) = 0;
[ATf_PDGF PHf_PDGF DTf_PDGF RPf_PDGF Y_output_PDGF] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_PDGFKO = PHf_PDGF(:,P);
Mean_impaired_PDGFKO = mean(PH_values_impaired_PDGFKO);
std_impaired_PDGFKO = std(PH_values_impaired_PDGFKO);

%TNF KO
P_TNF = [44 17 18 19];
Parametersets_impaired_secondfilter_TNFKO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_TNFKO(:,P_TNF) = 0;
[ATf_TNF PHf_TNF DTf_TNF RPf_TNF Y_output_TNF] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_TNFKO = PHf_TNF(:,P);
Mean_impaired_TNFKO = mean(PH_values_impaired_TNFKO);
std_impaired_TNFKO = std(PH_values_impaired_TNFKO);

%IL1 KO
P_IL1 = [45 20 21 22];
Parametersets_impaired_secondfilter_IL1KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_IL1KO(:,P_IL1) = 0;
[ATf_IL1 PHf_IL1 DTf_IL1 RPf_IL1 Y_output_IL1] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_IL1KO = PHf_IL1(:,P);
Mean_impaired_IL1KO = mean(PH_values_impaired_IL1KO);
std_impaired_IL1KO = std(PH_values_impaired_IL1KO);

%IL6 KO
P_IL6 = [46 23 24 102 145 25];
Parametersets_impaired_secondfilter_IL6KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_IL6KO(:,P_IL6) = 0;
[ATf_IL6 PHf_IL6 DTf_IL6 RPf_IL6 Y_output_IL6] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_IL6KO = PHf_IL6(:,P);
Mean_impaired_IL6KO = mean(PH_values_impaired_IL6KO);
std_impaired_IL6KO = std(PH_values_impaired_IL6KO);

%IL10 KO
P_IL10 = [26 27 89 28];
Parametersets_impaired_secondfilter_IL10KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_IL10KO(:,P_IL10) = 0;
[ATf_IL10 PHf_IL10 DTf_IL10 RPf_IL10 Y_output_IL10] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_IL10KO = PHf_IL10(:,P);
Mean_impaired_IL10KO = mean(PH_values_impaired_IL10KO);
std_impaired_IL10KO = std(PH_values_impaired_IL10KO);

%IL8 KO
P_IL8 = [29 30 90 141 31];
Parametersets_impaired_secondfilter_IL8KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_IL8KO(:,P_IL8) = 0;
[ATf_IL8 PHf_IL8 DTf_IL8 RPf_IL8 Y_output_IL8] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_IL8KO = PHf_IL8(:,P);
Mean_impaired_IL8KO = mean(PH_values_impaired_IL8KO);
std_impaired_IL8KO = std(PH_values_impaired_IL8KO);

%IL12 KO
P_IL12 = [32 33];
Parametersets_impaired_secondfilter_IL12KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_IL12KO(:,P_IL12) = 0;
[ATf_IL12 PHf_IL12 DTf_IL12 RPf_IL12 Y_output_IL12] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_IL12KO = PHf_IL12(:,P);
Mean_impaired_IL12KO = mean(PH_values_impaired_IL12KO);
std_impaired_IL12KO = std(PH_values_impaired_IL12KO);

%MIP1 KO
P_MIP1 = [34 35 36];
Parametersets_impaired_secondfilter_MIP1KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_MIP1KO(:,P_MIP1) = 0;
[ATf_MIP1 PHf_MIP1 DTf_MIP1 RPf_MIP1 Y_output_MIP1] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_MIP1KO = PHf_MIP1(:,P);
Mean_impaired_MIP1KO = mean(PH_values_impaired_MIP1KO);
std_impaired_MIP1KO = std(PH_values_impaired_MIP1KO);

%MIP2 KO
P_MIP2 = [37 38 39];
Parametersets_impaired_secondfilter_MIP2KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_MIP2KO(:,P_MIP2) = 0;
[ATf_MIP2 PHf_MIP2 DTf_MIP2 RPf_MIP2 Y_output_MIP2] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_MIP2KO = PHf_MIP2(:,P);
Mean_impaired_MIP2KO = mean(PH_values_impaired_MIP2KO);
std_impaired_MIP2KO = std(PH_values_impaired_MIP2KO);

%IP10 KO
P_IP10 = [40 41 42];
Parametersets_impaired_secondfilter_IP10KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_IP10KO(:,P_IP10) = 0;
[ATf_IP10 PHf_IP10 DTf_IP10 RPf_IP10 Y_output_IP10] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_IP10KO = PHf_IP10(:,P);
Mean_impaired_IP10KO = mean(PH_values_impaired_IP10KO);
std_impaired_IP10KO = std(PH_values_impaired_IP10KO);

% fibronectin KO
P_fibronectin = [75 76 77];
Parametersets_impaired_secondfilter_fibronectinKO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_fibronectinKO(:,P_fibronectin) = 0;
[ATf_fib PHf_fib DTf_fib RPf_fib Y_output_fib] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_fibronectinKO = PHf_fibronectin(:,P);
Mean_impaired_fibronectinKO = mean(PH_values_impaired_fibronectinKO);
std_impaired_fibronectinKO = std(PH_values_impaired_fibronectinKO);

% FGF KO
P_FGF = [78 43 79];
Parametersets_impaired_secondfilter_FGFKO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_FGFKO(:,P_FGF) = 0;
[ATf_FGF PHf_FGF DTf_FGF RPf_FGF Y_output_FGF] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_FGFKO = PHf_FGF(:,P);
Mean_impaired_FGFKO = mean(PH_values_impaired_FGFKO);
std_impaired_FGFKO = std(PH_values_impaired_FGFKO);

% MMP9 KO
P_MMP9 = [80 81 82];
Parametersets_impaired_secondfilter_MMP9KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_MMP9KO(:,P_MMP9) = 0;
[ATf_MMP9 PHf_MMP9 DTf_MMP9 RPf_MMP9 Y_output_MMP9] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_MMP9KO = PHf_MMP9(:,P);
Mean_impaired_MMP9KO = mean(PH_values_impaired_MMP9KO);
std_impaired_MMP9KO = std(PH_values_impaired_MMP9KO);

% MMP1 KO
P_MMP1 = [83 84];
Parametersets_impaired_secondfilter_MMP1KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_MMP1KO(:,P_MMP1) = 0;
[ATf_MMP1 PHf_MMP1 DTf_MMP1 RPf_MMP1 Y_output_MMP1] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_MMP1KO = PHf_MMP1(:,P);
Mean_impaired_MMP1KO = mean(PH_values_impaired_MMP1KO);
std_impaired_MMP1KO = std(PH_values_impaired_MMP1KO);

% MMP2 KO
P_MMP2 = [85 86];
Parametersets_impaired_secondfilter_MMP2KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_MMP2KO(:,P_MMP2) = 0;
[ATf_MMP2 PHf_MMP2 DTf_MMP2 RPf_MMP2 Y_output_MMP2] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delayetime2); 

PH_values_impaired_MMP2KO = PHf_MMP2(:,P);
Mean_impaired_MMP2KO = mean(PH_values_impaired_MMP2KO);
std_impaired_MMP2KO = std(PH_values_impaired_MMP2KO);

% TIMP1 KO
P_TIMP1 = [103 87 88];
Parametersets_impaired_secondfilter_TIMP1KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_TIMP1KO(:,P_TIMP1) = 0;
[ATf_TIMP1 PHf_TIMP1 DTf_TIMP1 RPf_TIMP1 Y_output_TIMP1] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_TIMP1KO = PHf_TIMP1(:,P);
Mean_impaired_TIMP1KO = mean(PH_values_impaired_TIMP1KO);
std_impaired_TIMP1KO = std(PH_values_impaired_TIMP1KO);

% MCP1 KO
P_MCP1 = [104 105];
Parametersets_impaired_secondfilter_MCP1KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_MCP1KO(:,P_MCP1) = 0;
[ATf_MCP1 PHf_MCP1 DTf_MCP1 RPf_MCP1 Y_output_MCP1] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_MCP1KO = PHf_MCP1(:,P);
Mean_impaired_MCP1KO = mean(PH_values_impaired_MCP1KO);
std_impaired_MCP1KO = std(PH_values_impaired_MCP1KO);
    
% VEGF KO
P_VEGF = [109 110 111 146 158 112];
Parametersets_impaired_secondfilter_VEGFKO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_VEGFKO(:,P_VEGF) = 0;
[ATf_VEGF PHf_VEGF DTf_VEGF RPf_VEGF Y_output_VEGF] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_VEGFKO = PHf_VEGF(:,P);
Mean_impaired_VEGFKO = mean(PH_values_impaired_VEGFKO);
std_impaired_VEGFKO = std(PH_values_impaired_VEGFKO);

% ANG1 KO
P_ANG1 = [113 114 115];
Parametersets_impaired_secondfilter_ANG1KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_ANG1KO(:,P_ANG1) = 0;
[ATf_ANG1 PHf_ANG1 DTf_ANG1 RPf_ANG1 Y_output_ANG1] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_ANG1KO = PHf_ANG1(:,P);
Mean_impaired_ANG1KO = mean(PH_values_impaired_ANG1KO);
std_impaired_ANG1KO = std(PH_values_impaired_ANG1KO);

% ANG2 KO
P_ANG2 = [116 117 118];
Parametersets_impaired_secondfilter_ANG2KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_ANG2KO(:,P_ANG2) = 0;
[ATf_ANG2 PHf_ANG2 DTf_ANG2 RPf_ANG2 Y_output_ANG2] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_ANG2KO = PHf_ANG2(:,P);
Mean_impaired_ANG2KO = mean(PH_values_impaired_ANG2KO);
std_impaired_ANG2KO = std(PH_values_impaired_ANG2KO);

% TSP1 KO
P_TSP1 = [119 120 121 122 123];
Parametersets_impaired_secondfilter_TSP1KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_TSP1KO(:,P_TSP1) = 0;
[ATf_TSP1 PHf_TSP1 DTf_TSP1 RPf_TSP1 Y_output_TSP1] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_TSP1KO = PHf_TSP1(:,P);
Mean_impaired_TSP1KO = mean(PH_values_impaired_TSP1KO);
std_impaired_TSP1KO = std(PH_values_impaired_TSP1KO);

% endostatin KO
P_endostatin = [124 125];
Parametersets_impaired_secondfilter_endostatinKO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_endostatinKO(:,P_endostatin) = 0;
[ATf_endo PHf_endo DTf_endo RPf_endo Y_output_endo] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_endostatinKO = PHf_endostatin(:,P);
Mean_impaired_endostatinKO = mean(PH_values_impaired_endostatinKO);
std_impaired_endostatinKO = std(PH_values_impaired_endostatinKO);

% Oxygen KO
P_O = [131 132 1333 134 135];
Parametersets_impaired_secondfilter_OKO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_OKO(:,P_O) = 0;
[ATf_O PHf_O DTf_O RPf_O Y_output_O] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_OKO = PHf_O(:,P);
Mean_impaired_OKO = mean(PH_values_impaired_OKO);
std_impaired_OKO = std(PH_values_impaired_OKO);

% PEDF KO
P_PEDF = [148 149 150];
Parametersets_impaired_secondfilter_PEDFKO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_PEDFKO(:,P_PEDF) = 0;
[ATf_PEDF PHf_PEDF DTf_PEDF RPf_PEDF Y_output_PEDF] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impaired_PEDFKO = PHf_PEDF(:,P);
Mean_impaired_PEDFKO = mean(PH_values_impaired_PEDFKO);
std_impaired_PEDFKO = std(PH_values_impaired_PEDFKO);

% CXCL1 KO
P_CXCL1 = [156 159 157];
Parametersets_impaired_secondfilter_CXCL1KO = Parametersets_impaired_secondfilter;
Parametersets_impaired_secondfilter_CXCL1KO(:,P_CXCL1) = 0;
[ATf_CXCL1 PHf_CXCL1 DTf_CXCL1 RPf_CXCL1 Y_output_CXCL1] = knockoutcalculation(Param1,tmax,tspan,yinit,delaytime,delaytime2); 

PH_values_impairedCXCL1KO = PHf_CXCL1(:,P);
Mean_impaired_CXCL1KO = mean(PH_values_impaired_CXCL1KO);
std_impaired_CXCL1KO = std(PH_values_impaired_CXCL1KO);
        
%--------------------------------------------------------------------
% list_proteins = {'TGF','PDGF','TNF','IL1', 'IL6','IL10','IL8','IL12','MIP1','MIP2','IP10','fibnec','fgf2','mmp9','timp','mmp1',
%     'mmp2','MCP1','VEGF','ANG1','ANG2','TSP1','endo','O','PEDF','KGF'};
% 
% i=1;
% load ('PEDFKO.mat');   
% AT_impaired{i} = ATf_KO;
% PH_impaired{i} = PHf_KO;
% DT_impaired{i} = DTf_KO;
% RP_impaired{i} = RPF_KO;
% i=i+1;
% clear ATf_KO PHf_KO DTf_KO RPF_KO
% 
% for i = 1:size(AT_impaired,2)
%     [ATmeanKO(i,:) ATstdKO(i,:) PHmeanKO(i,:) PHstdKO(i,:) DTmeanKO(i,:) DTstdKO(i,:) RPmeanKO(i,:) RPstdKO(i,:)] = PHimpairedKO(AT_impaired{i}, PH_impaired{i},DT_impaired{i},RP_impaired{i}, P) ;
% end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
