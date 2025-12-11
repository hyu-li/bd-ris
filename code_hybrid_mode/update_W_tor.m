function W = update_W_tor(H,G,P,aux,phy)
    [~,K] = size(H);
    [~,N_t] = size(G);
    global N_max
    global eta
    
    W = zeros(N_t,K);
    
    x1 = 0;
    x2 = 5000;
    mid = (x1 + x2)/2;  
    A = zeros(N_t,N_t);
    for i = 1 : K
        A = A + (aux.t(i)'*H(:,i)'*phy*G)'*...
            (aux.t(i)'*H(:,i)'*phy*G);
    end
    for k = 1 : K
        W(:,k) = (A + mid*eye(N_t))\(sqrt(1 + aux.r(k))*(aux.t(k)'*...
            H(:,k)'*phy*G)');
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
        for k = 1 : K
            W(:,k) = (A + mid*eye(N_t))\(sqrt(1 + aux.r(k))*(aux.t(k)'*...
                H(:,k)'*phy*G)');
        end     
        sum_p = norm(W,'fro')^2;
        nn = nn + 1;
    end 
    if abs(P - sum_p)/P >= eta
        for k = 1 : K
            W(:,k) = A\(sqrt(1 + aux.r(k))*(aux.t(k)'*...
                H(:,k)'*phy*G)');
        end       
        W = sqrt(P)*W/norm(W,'fro');
    end
    

%% test
%     lamda_range = [0:0.1:100];
%     ff = zeros(1,length(lamda_range));
%     for nn = 1 : length(lamda_range)
%         lambda = lamda_range(nn);
%         for k = 1 : K
%             W(:,k) = pinv(A + lambda*eye(N_t))*(sqrt(1 + aux.r(k))*(aux.t(k)'*...
%                 H(:,k)'*phy*G)');
%         end      
%         sum_p = norm(W,'fro')^2;
%         ff(nn) = P - sum_p;
%     end
%     plot(lamda_range,ff)
end

