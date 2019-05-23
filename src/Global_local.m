% Simulation routine for extended sensitivity analysis (in the vicinities of 10,000 random parameter sets)
% Author: Sridevi Nagaraja
% Last modified: August 27, 2012
% This file contains code for the extended sensitivity analysis. 
% This routine initializes parameters values and initial conditions and generates 'iter' number of random parameter sets
% It performs local sensitivity analysis on the 'iter' number of parameter sets using function Param_var_local.m
clear all
close all
clc

%*** CHANGE HERE TO SPECIFY NUMBER OF GENERATED PARAMETER SETS *** 
iter = 60000; % NUMBER OF RANDOM PARAMETER SETS TO BE GENERATED 

Parameters;%Initializes parameters and initial conditions

rangefactor = 2.5; % Limit for uniform distribution variation for creation of random parameter sets

Param1 = Latinhypercube(Param,rangefactor,iter); % function to create 'iter' number of random parameters sets using Latin Hypercube Sampling technique

% Time span
%Initializing time parameters
tinit = 0; % Initial time
tday = 42; % Maximum time (days)
tstep  = 1;% Time step (one hour)
tstep_sens = 24;
tmax = tday*24; %in hours
tspan = (tinit:tstep:tmax);
tspan_sens = (tinit:tstep_sens:tmax);
delaytime = 12;
delaytime2 = 48;
tic;
tstart=tic;

AT_def=zeros(1,40);
PH_def= zeros(1,40);
DT_def = zeros(1,40);
RP_def = zeros(1,40);
G_sen_def=zeros(40,tday);
Y_main_def=zeros(tmax,40);
E = zeros(iter,1);

%Assigning parallel processors, comment if not using parallel processing
 parfor_progress(iter)
 parpool('local',12)

% Function to calculate the raw and logarithmic sensitivities for 'iter' number of parameter sets and extract the critical triggers of chronic inflammation 
parfor i = 1:iter % if you are not using parallel processing, comment this line and uncomment the line below
%for i = 1:iter
    i
    try
    Gsen= Param_var_local(tmax,tspan,tspan_sens,yinit,Param1(i,:),delaytime,delaytime2);
    [AT PH DT RP Y_main]= Curvecharacteristic(tmax,tspan,yinit,Param1(i,:),delaytime,delaytime2);
    
    ATf(i,:) = AT;
    PHf(i,:) = PH;
    DTf(i,:) = DT;
    RPf(i,:) = RP;
    G_global{i} = Gsen;
    Y_main_global{i} = Y_main;
    parfor_progress %****comment if not using parallel processing****
     catch
    E(i) = i;
    ATf(i,:) = AT_def;
    PHf(i,:) = PH_def;
    DTf(i,:) = DT_def;
    RPf(i,:) = RP_def;
    Y_main_global{i} = Y_main_def;
    G_global{i} = G_RP_def;
    parfor_progress %****comment if not using parallel processing****
    end
end
parfor_progress(0); %****comment if not using parallel processing****
delete(gcp('nocreate')); %****comment if not using parallel processing****

runtime = toc(tstart)/60;

save('outputvars','tspan','iter','Param1','ATf','PHf','DTf','RPf');
save('G_global','-v7.3');
save('Y_main_global','-v7.3');


