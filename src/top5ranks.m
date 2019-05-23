function [P1 P2] = top5ranks(P, tps ,G)

for i=1:size(P,2)
   for j=1:size(tps,2)
       S = extractsensitivity(P(i),tps(j),G);
       [d1 d2 D1 D2] = rankedsensitivity(S);    
       P1{i,j} = d1;
       P2{i,j} = d2;
   end
end