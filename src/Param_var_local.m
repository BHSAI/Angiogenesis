% Function for calculating the logarithmic local sensitivity values for a given parameter set
function Gsen_local= Param_var_local(tmax,tspan,tspan_sens,Yinit,Param,delaytime,delaytime2)

% Function for calculation of local sensitivity indices by implementing finite difference formula
% Author: Sridevi Nagaraja
% July, 2017

%Inputs: tmax: maximum time of simulation
%        tspan: simulation time steps
%        Yinit: initial value vector of model variables for each iteration
%        Param: array of N default parameter sets created by Latin Hypercube Sampling
%        delaytime: lag in DDE for pro-inflammatory macrophages

%Outputs: Gsen_local: Matrix containing the logarthmic sensitivity values for all outputs for the given parameter set


delta_param = Param.*0.1; %Calculate 10% of every parameter value*** change here for different sensitivity range

% Calling DDE with default parameters
sol = dde23(@woundhealingequations2,[delaytime delaytime2],Yinit,tspan,[],Param);
sol1 = deval(sol,tspan_sens); % calculating solution at particular time points
Ymain_default = [sol1(1,:); sol1(2,:); sol1(3,:); sol1(4,:); sol1(5,:); sol1(6,:); sol1(7,:); sol1(8,:); sol1(9,:); sol1(10,:);sol1(11,:); sol1(12,:); sol1(13,:); sol1(14,:); sol1(15,:); sol1(16,:);sol1(17,:);sol1(18,:);sol1(19,:);sol1(20,:);sol1(21,:);sol1(22,:);sol1(23,:);sol1(24,:);sol1(25,:);sol1(26,:);sol1(27,:);sol1(28,:);sol1(29,:);sol1(30,:);sol1(31,:);sol1(32,:);sol1(33,:);sol1(34,:);sol1(35,:);sol1(36,:);sol1(37,:);sol1(38,:);sol1(39,:)];
clear sol sol1
clear woundhealingequations2

% Calculating outputs with kj+delta_kj for each parameter
for j = 1:size(Param,2)
   
Param_new_add = Param;
Param_new_add(j) = Param_new_add(j)+delta_param(j);
sol = dde23(@woundhealingequations2,[delaytime delaytime2],Yinit,tspan,[],Param_new_add);
sol1 = deval(sol,tspan_sens); % calculating solution at particular time points (per day)
Ymain_add{j} = [sol1(1,:); sol1(2,:); sol1(3,:); sol1(4,:); sol1(5,:); sol1(6,:); sol1(7,:); sol1(8,:); sol1(9,:); sol1(10,:);sol1(11,:); sol1(12,:); sol1(13,:); sol1(14,:); sol1(15,:); sol1(16,:);sol1(17,:);sol1(18,:);sol1(19,:);sol1(20,:);sol1(21,:);sol1(22,:);sol1(23,:);sol1(24,:);sol1(25,:);sol1(26,:);sol1(27,:);sol1(28,:);sol1(29,:);sol1(30,:);sol1(31,:);sol1(32,:);sol1(33,:);sol1(34,:);sol1(35,:);sol1(36,:);sol1(37,:);sol1(38,:);sol1(39,:);sol1(40,:)];
clear sol sol1 Param_new_add
clear woundhealingequations2

% Calculating outputs with kj-delta_kj for each parameter
Param_new_subs = Param;
Param_new_subs(j) = Param_new_subs(j)-delta_param(j);
sol = dde23(@woundhealingequations2,[delaytime delaytime2],Yinit,tspan,[],Param_new_subs);
sol1 = deval(sol,tspan_sens); % calculating solution at particular time points (per day)
Ymain_subs{j} = [sol1(1,:); sol1(2,:); sol1(3,:); sol1(4,:); sol1(5,:); sol1(6,:); sol1(7,:); sol1(8,:); sol1(9,:); sol1(10,:);sol1(11,:); sol1(12,:); sol1(13,:); sol1(14,:); sol1(15,:); sol1(16,:);sol1(17,:);sol1(18,:);sol1(19,:);sol1(20,:);sol1(21,:);sol1(22,:);sol1(23,:);sol1(24,:);sol1(25,:);sol1(26,:);sol1(27,:);sol1(28,:);sol1(29,:);sol1(30,:);sol1(31,:);sol1(32,:);sol1(33,:);sol1(34,:);sol1(35,:);sol1(36,:);sol1(37,:);sol1(38,:);sol1(39,:);sol1(40,:)];
clear sol sol1 Param_new_subs
clear woundhealingequations2
end

% Applying the central finite difference formula to calculate logarithmic sensitivites
for j = 1:size(Param,2)
    g1 = Ymain_add{j};
    g2 = Ymain_subs{j};
    G_num = (g1-g2);
    G_num(abs(G_num)<1e-10)=0;
    G_final{j}= (G_num./Ymain_default)./2e-1;%2e-2;% calculate the relative sensitivities by dividing with default values
    G_final{j}(isnan(G_final{j})) = 0; % remove any non-numbers
    clear g1 g2 G_num G_final1
end

Gsen_local= G_final;

