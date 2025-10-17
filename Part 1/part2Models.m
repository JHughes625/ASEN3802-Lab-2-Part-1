function [t, u] = part2Models(Han,k,rho,cp,T_0,L_Rod)
    x = .0127:0.0127:8*.0127;
    t = 1:.1:1000;
    u = 0;
    for n=1:50
        b = -8*Han*L_Rod*(((-1)^(n+1))/(pi*pi*(2*n-1)^2));
        lam = (2*n-1)*pi/(2*L_Rod);
        alph = k/rho/cp;
        for j=1:length(x)
            u1(j,:) = b*sin(lam*x(j))*exp(-1*lam^2 * alph * t);
        end
        u = u1 + u;
    end
end
