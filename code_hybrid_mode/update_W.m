function W = update_W(H,P,aux,phy_t,phy_r)
    [~,K_t] = size(H.Ht);
    [~,K_r] = size(H.Hr);
    [~,N_t] = size(H.G);
    K = K_t + K_r;
    global N_max
    global eta
    
    W = zeros(N_t,K);
    
    x1 = 0;
    x2 = 5000;
    mid = (x1 + x2)/2;  
    A = zeros(N_t,N_t);
    for i = 1 : K_t
        A = A + (aux.t(i)'*H.Ht(:,i)'*phy_t*H.G)'*...
            (aux.t(i)'*H.Ht(:,i)'*phy_t*H.G);
    end
    for i = 1 : K_r
        A = A + (aux.t(K_t + i)'*H.Hr(:,i)'*phy_r*H.G)'*...
            (aux.t(K_t + i)'*H.Hr(:,i)'*phy_r*H.G);        
    end
    for k = 1 : K_t
        W(:,k) = (A + mid*eye(N_t))\(sqrt(1 + aux.r(k))*(aux.t(k)'*...
            H.Ht(:,k)'*phy_t*H.G)');
    end
    for k = 1 : K_r
        W(:,K_t + k) = (A + mid*eye(N_t))\(sqrt(1 + aux.r(K_t + k))*(aux.t(K_t + k)'*...
            H.Hr(:,k)'*phy_r*H.G)');
    end       
    sum_p = norm(W,'fro')^2;
    
    nn = 1;           
    while abs(P - sum_p)/P >= eta && nn < N_max                            
        if P - sum_p > 0
            x2 = mid;
        else
            x1 = mid;
        end
        mid = (x1 + x2)/2;
        for k = 1 : K_t
            W(:,k) = (A + mid*eye(N_t))\(sqrt(1 + aux.r(k))*(aux.t(k)'*...
                H.Ht(:,k)'*phy_t*H.G)');
        end
        for k = 1 : K_r
            W(:,K_t + k) = (A + mid*eye(N_t))\(sqrt(1 + aux.r(K_t + k))*(aux.t(K_t + k)'*...
                H.Hr(:,k)'*phy_r*H.G)');
        end       
        sum_p = norm(W,'fro')^2;
        nn = nn + 1;
    end 
    if abs(P - sum_p)/P >= eta
        for k = 1 : K_t
            W(:,k) = A\(sqrt(1 + aux.r(k))*(aux.t(k)'*...
                H.Ht(:,k)'*phy_t*H.G)');
        end
        for k = 1 : K_r
            W(:,K_t + k) = A\(sqrt(1 + aux.r(K_t + k))*(aux.t(K_t + k)'*...
                H.Hr(:,k)'*phy_r*H.G)');
        end       
        W = sqrt(P)*W/norm(W,'fro');
    end
    

%% test
%     lamda_range = [0:0.1:100];
%     ff = zeros(1,length(lamda_range));
%     for nn = 1 : length(lamda_range)
%         lambda = lamda_range(nn);
%         for k = 1 : K_t
%             W(:,k) = pinv(A + lambda*eye(N_t))*(sqrt(1 + aux.r(k))*(aux.t(k)'*...
%                 H.Ht(:,k)'*phy_t*H.G)');
%         end
%         for k = 1 : K_r
%             W(:,K_t + k) = pinv(A + lambda*eye(N_t))*(sqrt(1 + aux.r(K_t + k))*(aux.t(K_t + k)'*...
%                 H.Hr(:,k)'*phy_r*H.G)');
%         end       
%         sum_p = norm(W,'fro')^2;
%         ff(nn) = P - sum_p;
%     end
%     plot(lamda_range,ff)
end

