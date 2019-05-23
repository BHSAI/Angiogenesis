function [H2 H3] = ranksorting(P2,rnk)

for i=1:size(P2,1)
    for j=1:size(P2,2)   
H1 = P2{i,j};
for k=1:rnk
edges = unique(H1(k,:));
counts = histc(H1(k,:), edges);
[f1 f2] = sort(counts,'descend');
topcounts(k,:) = f1(1,1:3);
topparam(k,:) = edges(f2(1,1:3));
end

H2{i,j} = topcounts;
H3{i,j} = topparam;
    end
end