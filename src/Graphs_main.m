% This file contains algorithms to plot the variable kinetics for the acute inflammation model
% Author: Sridevi Nagaraja
% Last modified: August 27, 2012
%A = {'N_t_o_t','M_t_o_t','Platelets','N_a_c_t','N_a_p_o_p','M_p_r_o','M_a_n_t_i','TGF-\beta','PDGF','TNF-\alpha','IL-1\beta','IL-6','IL-10','CXCL8','IL12','MIP-1\alpha','MIP-2','IP-10'};
% A = {'N_t_o_t','M_t_o_t','N_a_c_t','N_a_p_o_p','M_p_r_o','M_a_n_t_i','PDGF','TNF-\alpha','IL-1\beta','IL-6','IL-10','CXCL8','IL12','MIP-1\alpha','MIP-2','IP-10'};
% B= {'P1','P10','P20','P30','P40','P50','P60','P70'};
% set(gca,'XTick',1:10:70, 'Ytick',1:16)
% format_ticks(gca,B,A)
% array 'g' contains all the outputs 
% cytokines are multiplied by 1000 to plot in pgml-1

% Extracting outputs from output array
tspan    = g(:,1); % time points
Ntot     = g(:,2); % Total neutrophils
Mtot     = g(:,3); % Total neutrophils
Colltot  = g(:,4); % Total collagen
Nact     = g(:,5); % Active neutrophils
Napop    = g(:,6); % Apoptotic neutrophils
Mpro     = g(:,7); % Pro-inflammatory macrophages
Manti    = g(:,8); % Anti-inflammatory macrophages
TGF      = g(:,9); % TGF-beta
PDGF     = g(:,10); % PDGF
TNF      = g(:,11); % TNF-alpha
IL1      = g(:,12); % IL-1beta
IL6      = g(:,13); % IL-6
IL10     = g(:,14); % IL-10
P        = g(:,15); % Platelets
IL8      = g(:,16); % IL-8
IL12     = g(:,17); % IL-12
MIP1     = g(:,18); % MIP-1alpha
MIP2     = g(:,19); % MIP-2
IP10     = g(:,20); % IP-10
F        = g(:,21);
myoF     = g(:,22);
fibnec   = g(:,23);%23
fgf2     = g(:,24);%24
mmp9     = g(:,25);%25
timp     = g(:,26);%26
Coll     = g(:,27);%29
mmp1     = g(:,28);%28
mmp2     = g(:,29);%29
col1_fib = g(:,30);%30
MCP1 = g(:,31);%31
intermed = g(:,32);%32
VEGF = g(:,33);%33
EC = g(:,34);%34
ANG1  = g(:,35);%35
ANG2 = g(:,36);%36
TSP1 = g(:,37);%37
endo = g(:,38);%38
%captip = g(:,39);%39
capsprout = g(:,39);%40
O = g(:,40);%41
PEDF = g(:,41);
K = g(:,42);
KGF = g(:,43);
CXCL1 = g(:,44);
totdensity = g(:,45);


figure(1)
h1 = plot(tspan./24, Ntot,'g'); %tspan is divided by 24 to plot time in days
hold on;
h2 = plot(tspan./24,Mtot,'b');
hold on
h3 = plot(tspan./24,F,'r');
hold on
h4 = plot(tspan./24,myoF,'k');
hold on
h5 = plot(tspan./24,EC,'c');
legend([h1 h2 h3 h4 h5],{'Neutrophils','Macrophages','Fibroblasts','Myofibroblasts','Endothelial cells'});
legend boxoff
box off
xlabel('Time  (d)')
ylabel('Cell concentration (cells mL^{-1})')
set(gca,'FontSize',20,'FontName','Arial','FontWeight','Bold','LineWidth',1.5 )
hold off

figure(2)
h1 = plot(tspan./24, Colltot,'g'); %tspan is divided by 24 to plot time in days
hold on;
h2 = plot(tspan./24,col1_fib,'b');
hold on
legend([h1 h2],{'total collagen','Collagen fiber'});
legend boxoff
box off
xlabel('Time  (d)')
ylabel('(cells mL^{-1})')
set(gca,'FontSize',20,'FontName','Arial','FontWeight','Bold','LineWidth',1.5 )
hold off

figure(3)
h1 = plot(tspan./24, EC,'g'); %tspan is divided by 24 to plot time in days
hold on;
h2 = plot(tspan./24,capsprout,'r');
hold on
h3 = plot(tspan./24,totdensity,'k');
legend([h1 h2 h3],{'ECtip','sprout','total density'});
legend boxoff
box off
xlabel('Time  (d)')
ylabel('(cells mL^{-1})')
set(gca,'FontSize',20,'FontName','Arial','FontWeight','Bold','LineWidth',1.5 )
hold off

figure(4)
h1 = plot(tspan./24, VEGF,'g'); %tspan is divided by 24 to plot time in days
hold on;
h2 = plot(tspan./24,ANG1,'b');
hold on
h3 = plot(tspan./24,ANG2,'r');
legend([h1 h2 h3],{'VEGF','ANG1','ANG2'});
legend boxoff
box off
xlabel('Time  (d)')
ylabel('(ng mL^{-1})')
set(gca,'FontSize',20,'FontName','Arial','FontWeight','Bold','LineWidth',1.5 )
hold off

figure(5)
h1 = plot(tspan./24, O,'g'); %tspan is divided by 24 to plot time in days
legend([h1],{'Oxygen'});
legend boxoff
box off
xlabel('Time  (d)')
ylabel('(nM)')
set(gca,'FontSize',20,'FontName','Arial','FontWeight','Bold','LineWidth',1.5 )
hold off

figure(6)
h1 = plot(tspan./24, PEDF,'g'); %tspan is divided by 24 to plot time in days
hold on;
h2 = plot(tspan./24,KGF,'r');
legend([h1 h2],{'PEDF','KGF'});
legend boxoff
box off
xlabel('Time  (d)')
ylabel('(ng mL^{-1})')
set(gca,'FontSize',20,'FontName','Arial','FontWeight','Bold','LineWidth',1.5 )
hold off

figure(7)
h1 = plot(tspan./24, K,'g'); %tspan is divided by 24 to plot time in days
legend([h1],{'Keratinocyte'});
legend boxoff
box off
xlabel('Time  (d)')
ylabel('(cells mL^{-1})')
set(gca,'FontSize',20,'FontName','Arial','FontWeight','Bold','LineWidth',1.5 )
hold off




