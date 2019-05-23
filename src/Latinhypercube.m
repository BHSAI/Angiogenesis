% Function that generates a user-defined number of random parameter sets
function P = Latinhypercube(Param,rangefactor,iter)

%Assigning parameter range
% Minimum values
k_min = Param./rangefactor;
% Maximum values
k_max = Param.*rangefactor;

%Creating Latin hypercube sampling matrix using Matlab function lhsdesign
X1 = lhsdesign(iter,size(Param,2));

%Creating randomized parameter sets using LHS
for i = 1:size(Param,2)
P(:,i) = (X1(:,i).*(k_max(i)-k_min(i)))+k_min(i);
end
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------























