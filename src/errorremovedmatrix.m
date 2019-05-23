function G = errorremovedmatrix(E,Global)
c = find(E==0);
 for i=1:size(c,1)
    G{i} = Global{1,c(i)};
 end 
 
 