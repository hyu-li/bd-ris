function phy = update_phy_single_tor(H,G,W,aux,phy)
    [M,K] = size(H);
    
    global SIM
    global eta
    
    B = zeros(M,M);
    c = zeros(M,1);

    
    for i = 1 : K
        V = zeros(M,M);
        for j = 1 : K
            v = (H(:,i)'*diag(G*W(:,j)*(aux.t(i))'))';
            V = V + v*v';
        end
        B = B + V; 
        c = c + sqrt(1 + aux.r(i))*(H(:,i)'*diag(G*W(:,i)*(aux.t(i))'))';
    end
    
    %% update RIS    
    phy_old = zeros(M,1);
    phy_vec = zeros(M,1);
    for m = 1 : M
        phy_vec(m) = phy(m,m);
    end
    nn = 1;
    while norm(phy_vec - phy_old,2)^2 > eta && nn < SIM
        phy_old = phy_vec;
        for m = 1 : M
            x = c(m) - B(m,:)*phy_vec + B(m,m)*phy_vec(m);
            phy_vec(m) = x/abs(x);
        end
        nn = nn + 1;
    end 
    phy = diag(phy_vec);
end

