function  [aux,sum_rate] = sum_rate_cal(H,obj)
%% this function caculates the sum-rate and auxillary variables r, t
    
    global sigma
    [~,K_t] = size(H.Ht);
    [~,K_r] = size(H.Hr);
    K = K_r + K_t;
    aux.r = zeros(K,1);
    aux.t = zeros(K,1);
    alpha_t = zeros(K_t,1);
    alpha_r = zeros(K_r,1);
    
    for i = 1 : K_t
        for j = 1 : K
            alpha_t(i) = alpha_t(i) + abs(H.Ht(:,i)'*obj.phy_t*H.G*obj.w(:,j))^2;
        end
        aux.r(i) = abs(H.Ht(:,i)'*obj.phy_t*H.G*obj.w(:,i))^2/(alpha_t(i)...
            - abs(H.Ht(:,i)'*obj.phy_t*H.G*obj.w(:,i))^2 + sigma^2);
        aux.t(i) = sqrt(1 + aux.r(i))*H.Ht(:,i)'*obj.phy_t*H.G*obj.w(:,i)...
            /(alpha_t(i) + sigma^2);
    end
    for i = 1 : K_r
        for j = 1 : K
            alpha_r(i) = alpha_r(i) + abs(H.Hr(:,i)'*obj.phy_r*H.G*obj.w(:,j))^2;
        end
        aux.r(K_t + i) = abs(H.Hr(:,i)'*obj.phy_r*H.G*obj.w(:,K_t + i))^2/...
            (alpha_r(i) - abs(H.Hr(:,i)'*obj.phy_r*H.G*obj.w(:,K_t + i))^2 + sigma^2);
        aux.t(K_t + i) = sqrt(1 + aux.r(K_t + i))*H.Hr(:,i)'*obj.phy_r*H.G*obj.w(:,K_t + i)...
            /(alpha_r(i) + sigma^2);
    end    
    sum_rate = real(sum(log2(1 + aux.r)));
end