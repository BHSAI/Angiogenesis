% EXTENDED SENSITIVITY ANALYSIS

%Outputs to be analysed
 P = [26 29 30 35 36];
 Pname = {'Collagen fiber','VEGF','endothelial cells','blood vessel','oxygen'};
 
 % time points to be ana;ysed
 tps = [3 8 15 29 43];

% Ranking and sorting extended local sensitivity analysis outputs   
%removing erroneous simulations
 G1 = errorremovedmatrix(E1,G_global1);
 G2 = errorremovedmatrix(E2,G_global2);
 G3 = errorremovedmatrix(E3,G_global3);
 G4 = horzcat(G1,G2);
 G=horzcat(G4,G3); 

% sorting the top five ranked parameters for the five outputs
[P1 P2] = top5ranks(P,tps,G); %P1 stores the actual sensitivity value and P2 stores the parameter number that induces the sensitivity in a given output

rnk=5;% number of ranks to be anlysed
[H2 H3] = ranksorting(P2,rnk);
%% 
%PRCC ANALYSIS

%Load the Y_main_global and the Param1 variables saved at the end of the global sensitivity
%analysis

%Remove any erroneous simulations that did not converge
c = find(E==0);
for i=1:length(c)
Ymain_3{n}= Y_main_global{1,c(i)};
n=n+1;
end
%Remove any erroneous parameter sets that did not converge
Param = Param1(c,:);

% Extract the model output variable values for each of the 42 days post wounding
k=1;
for i=1:24:1009 %calculate the PRCCs for every day of the simulation
for j=1:size(Ymain_3,2)
S1 = Ymain_3{j}; %Extracting the output values at the specified time points for the 60,000 simulations
S2(:,j) = S1(:,i);
end
S3{k}=S2';
k=k+1;
end

% Calculating the PRCCs values between each parameter and each model
% output at 42 different days
for s=1:size(S3,2)
[prcc sign]=PRCC2(Param,S3{s});
M.prcc{s} = prcc;
M.sign{s} = sign;
end


