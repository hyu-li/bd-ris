function [phy_t_m,phy_r_m] = golden_section_search(B_t,B_r,c_t,c_r,phy_t,phy_r,m)
    err = 1e-4;
    x1 = 0;
    x2 = 1;
    a1 = x1 + 0.382*(x2 - x1);
    x_t = B_t(m,:)*phy_t - B_t(m,m)*phy_t(m) - c_t(m);
    x_r = B_r(m,:)*phy_r - B_r(m,m)*phy_r(m) - c_r(m);
    theta_t = angle(x_t) + pi;
    theta_r = angle(x_r) + pi;
    y1 = (B_t(m,m) - B_r(m,m))*a1^2 + B_r(m,m)...
        -2*abs(x_t)*a1 - 2*abs(x_r)*sqrt(1 - a1^2);
    a2 = x1 + 0.618*(x2 - x1);
    y2 = (B_t(m,m) - B_r(m,m))*a2^2 + B_r(m,m)...
        -2*abs(x_t)*a2 - 2*abs(x_r)*sqrt(1 - a2^2);
    while a2 - a1 > err
        if y1 <= y2
            x2 = a2;
            a2 = a1;
            y2 = y1;
            a1 = x1 + 0.382*(x2 - x1);
            y1 = (B_t(m,m) - B_r(m,m))*a1^2 + B_r(m,m)...
                    -2*abs(x_t)*a1 - 2*abs(x_r)*sqrt(1 - a1^2);
        else
            x1 = a1;
            a1 = a2;
            y1 = y2;
            a2 = x1 + 0.618*(x2 - x1);
            y2 = (B_t(m,m) - B_r(m,m))*a2^2 + B_r(m,m)...
                    -2*abs(x_t)*a2 - 2*abs(x_r)*sqrt(1 - a2^2);
        end
    end
    alpha_t = (a1 + a2)/2;
    alpha_r = sqrt(1 - alpha_t^2);
    phy_t_m = alpha_t*exp(1j*theta_t);
    phy_r_m = alpha_r*exp(1j*theta_r);

    %% test
%     x = [0 : 0.01 : 1];
%     y = zeros(length(x),1);
%     for i = 1 : length(x)
%         y(i) = (B_t(m,m) - B_r(m,m))*x(i)^2 + B_r(m,m)...
%         -2*abs(x_t)*x(i) - 2*abs(x_r)*sqrt(1 - x(i)^2);
%     end
%     y_alpha_t = (B_t(m,m) - B_r(m,m))*alpha_t^2 + B_r(m,m)...
%         -2*abs(x_t)*alpha_t - 2*abs(x_r)*sqrt(1 - alpha_t^2);
%     plot(x,y,'-k');
%     hold on
%     plot(alpha_t,y_alpha_t,'ro');
%     grid on
    
end

