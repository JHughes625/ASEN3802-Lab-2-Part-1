function [g, Mexp] = M_exp(expData)
l = linspace(0,8*0.0127);
for i=1:5
    initialVec(i,:) = expData.values(1,2:9);
    T0(i) = initialVec(i,1);
    x = 1:length(initialVec);
    p = polyfit(x,initialVec(i,:),1);
    Mexp(i) = p(1);
    g(i) = Mexp(i)*l + T0(i);
end
end