function  [aux,sum_rate] = sum_rate_cal_tor(H,G,obj_tor)
%% this function caculates the sum-rate and auxillary variables r, t
    
    global sigma
    [~,K] = size(H);
    aux.r = zeros(K,1);
    aux.t = zeros(K,1);
    alpha = zeros(K,1);
    
    for i = 1 : K
        for j = 1 : K
            alpha(i) = alpha(i) + abs(H(:,i)'*obj_tor.phy*G*obj_tor.w(:,j))^2;
        end
        aux.r(i) = abs(H(:,i)'*obj_tor.phy*G*obj_tor.w(:,i))^2/(alpha(i)...
            - abs(H(:,i)'*obj_tor.phy*G*obj_tor.w(:,i))^2 + sigma^2);
        aux.t(i) = sqrt(1 + aux.r(i))*H(:,i)'*obj_tor.phy*G*obj_tor.w(:,i)...
            /(alpha(i) + sigma^2);
    end
    sum_rate = real(sum(log2(1 + aux.r)));
end