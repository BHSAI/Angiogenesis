%function [AT PH DT RP Ymain_default G_PHfinal G_ATfinal G_DTfinal G_RPfinal]= Curvecharacteristic(tmax,tspan, Yinit, Param,delaytime, delaytime2)
function [AT PH DT RP Ymain_default1]= Curvecharacteristic(tmax,tspan, Yinit, Param,delaytime, delaytime2)

% Function for calculation of inflammation indices and sensitivity values of inflammation indices by implementing finite difference formula
% Author: Sridevi Nagaraja

%Inputs: tmax: maximum time of simulation
%        tspan: simulation time steps
%        Yinit: initial value vector of model variables for each iteration
%        Param: array of 60,000 default parameter sets created by Latin Hypercube Sampling
%        delaytime: lag in DDE for pro-inflammatory macrophages

%Outputs: AT, PH, DT, RP:Activation time, Peak height, Resolution interval, and Resolution plateau for all model outputs using default parameter set
%        Ymain_default: Time courses of each model output variable fromt he
%        60,000 parameter sets

%Intializing output matrices
Ymain_add = cell([],[]);% Output matrix for paramter varied by adding 1% default value
Ymain_subs= cell([],[]);% Output matrix for paramter varied by substracting 1% default value
Ymain_default = cell([],[]);% Output matrix for default parameter set 

% Calling DDE with default parameters
sol = dde23(@woundhealingequations2,[delaytime delaytime2],Yinit,tspan,[],Param);
sol1 = deval(sol,tspan); % calculating  solution at particular time points
Ymain_default1 = [sol1(1,:); sol1(2,:); sol1(3,:); sol1(4,:); sol1(5,:); sol1(6,:); sol1(7,:); sol1(8,:); sol1(9,:); sol1(10,:);sol1(11,:); sol1(12,:); sol1(13,:); sol1(14,:); sol1(15,:); sol1(16,:);sol1(17,:);sol1(18,:);sol1(19,:);sol1(20,:);sol1(21,:);sol1(22,:);sol1(23,:);sol1(24,:);sol1(25,:);sol1(26,:);sol1(27,:);sol1(28,:);sol1(29,:);sol1(30,:);sol1(31,:);sol1(32,:);sol1(33,:);sol1(34,:);sol1(35,:);sol1(36,:);sol1(37,:);sol1(38,:);sol1(39,:);sol1(40,:)];
Ymain_default = [sol.y(1,:); sol.y(2,:); sol.y(3,:); sol.y(4,:); sol.y(5,:); sol.y(6,:); sol.y(7,:); sol.y(8,:); sol.y(9,:); sol.y(10,:);sol.y(11,:); sol.y(12,:); sol.y(13,:); sol.y(14,:); sol.y(15,:); sol.y(16,:);sol.y(17,:);sol.y(18,:);sol.y(19,:);sol.y(20,:);sol.y(21,:);sol.y(22,:);sol.y(23,:);sol.y(24,:);sol.y(25,:);sol.y(26,:);sol.y(27,:);sol.y(28,:);sol.y(29,:);sol.y(30,:);sol.y(31,:);sol.y(32,:);sol.y(33,:);sol.y(34,:);sol.y(35,:);sol.y(36,:);sol.y(37,:);sol.y(38,:);sol.y(39,:);sol.y(40,:)];
[AT PH DT RP] = timescalculation(Ymain_default,sol.x);% default inflammation index values
clear sol

