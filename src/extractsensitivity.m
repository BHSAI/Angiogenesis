function S = extractsensitivity(P,tps,G1)

for i=1:size(G1,2)
    F = G1{i};
    for j=1:size(F,2)
        F1 = F{j};
        S(j,i) = F1(P,tps);
    end      
end
