function [g, Mexp, l] = M_exp(expData)
l = linspace(0,8*0.0127);
    initialVec = expData(1,2:9);
    T0 = initialVec(1);
    x = 1:length(initialVec);
    p = polyfit(x,initialVec,1);
    Mexp = p(1);
    g = Mexp*l + T0;
end