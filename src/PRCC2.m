function [prcc sign]=PRCC2(Param1,Y);

%Defining Matrix size

[a k]=size(Param1); % Define the size of LHS matrix
[b out]=size(Y);

for i=1:k  % Loop for the whole submatrices
 for j =1:out
    F1 = Param1(:,i);
    F2 = Y(:,j);
    F3 = Param1;
    F3(:,i) = [];
    F4 = F3;
    [rho,p]=partialcorr(F1,F2,F4,'type','Spearman');
    prcc(i,j) = rho;
    sign(i,j) = p; 
 end
end