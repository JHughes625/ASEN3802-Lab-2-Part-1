function [t, u1] = part2Models(Han,alpha,T_0,L_Rod,tend)
    x = .0127:0.0127:8*.0127;
    t = 1:1:tend;
    u1 = 0;
    for n=1:length(t)
        for j=1:length(x)
            for i=1:10
                b = -8*Han*L_Rod*(((-1)^(i+1))/(pi*pi*(2*i-1)^2));
                lam = (2*i-1)*pi/(2*L_Rod);
                %alph = k/rho/cp;
                if i-1 == 0
                    u1(j,n,i) = b*sin(lam*x(j))*exp(-1*lam^2 * alpha * t(n));
                else
                    u1(j,n,i) = b*sin(lam*x(j))*exp(-1*lam^2 * alpha * t(n)) + u1(j,n,i-1);
                end
            end
            u1(j,n,end) = u1(j,n,end) + T_0 + Han*x(j);
        end
    end

    for i=1:length(x)
        u(i,:) = u1(i,:,end);
    end
end
