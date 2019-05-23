% Last modified:10/28/2013
% Author:Sridevi Nagaraja
% Main function of the acute inflammation model. This routine initializes the parameter values and initial conditions (Parameters.m), defines the time span of the simulation, 
% calls the in-built MATLAB delay differential equation solver DDE23, saves the predicted values of all model variables and plots the results.
% Results are plotted using Graphs_main

clear all
close all
clc

tic;
tstart=tic;

Parameters;% Initializing all parameter values and initial conditions
%Param = Param1(8,:);
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
 
% Calling ode
sol1 = dde23(@woundhealingequations2,[delaytime delaytime2],yinit,tspan,[],Param);  
sol = deval(sol1,tspan);
Ntot = sol(1,:)+sol(2,:); % Calculating total neutrophils 
Mtot = sol(3,:)+sol(4,:); % Calculating total macrophages 
Colltot = sol(23,:)+sol(26,:)+sol(28,:); % Calculating total collagen
bloodvesselden = sol(30,:)+sol(35,:);
Y = [tspan; Ntot; Mtot; Colltot; sol(1,:); sol(2,:); sol(3,:); sol(4,:); sol(5,:); sol(6,:); sol(7,:); sol(8,:); sol(9,:); sol(10,:);sol(11,:); sol(12,:); sol(13,:); sol(14,:); sol(15,:); sol(16,:);sol(17,:);sol(18,:);sol(19,:);sol(20,:);sol(21,:);sol(22,:);sol(23,:);sol(24,:);sol(25,:);sol(26,:);sol(27,:);sol(28,:);sol(29,:);sol(30,:);sol(31,:);sol(32,:);sol(33,:);sol(34,:);sol(35,:);sol(36,:); sol(37,:);sol(38,:);sol(39,:);sol(40,:);bloodvesselden];
Y1 = [sol(1,:); sol(2,:); sol(3,:); sol(4,:); sol(5,:); sol(6,:); sol(7,:); sol(8,:); sol(9,:); sol(10,:);sol(11,:); sol(12,:); sol(13,:); sol(14,:); sol(15,:); sol(16,:);sol(17,:);sol(18,:);sol(19,:);sol(20,:);sol(21,:);sol(22,:);sol(23,:);sol(24,:);sol(25,:);sol(26,:);sol(27,:);sol(28,:);sol(29,:);sol(30,:);sol(31,:);sol(32,:);sol(33,:);sol(34,:);sol(35,:);sol(36,:); sol(37,:);sol(38,:);sol(39,:);sol(40,:)];
g = transpose(Y);% Saving the outputs
gmin = min(g);
gmax = max(g);

%Calculate normalized values of model outputs
for  i=1:45
gnorm(:,i) = (g(:,i)-gmin(i))./(gmax(i)-gmin(i));
end
runtime = toc(tstart)/60;

%Plotting results
Graphs_main
save('outputvars','gnorm','g');

%---------------Local sensitivity analysis---------------------------------
% tspan_sens = 0:24:tmax; % increase time step to one day for sensitivity analysis
% [Gsen_local]= Param_var_local(tmax,tspan,tspan_sens,yinit,Param,delaytime,delaytime2);

