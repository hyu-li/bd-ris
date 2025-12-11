function [phy_t,phy_r] = update_phy_single(H,W,aux,phy_t,phy_r)
    [M,K_t] = size(H.Ht);
    [~,K_r] = size(H.Hr);
    
    K = K_t + K_r;
    global SIM
    global eta
    
    B_t = zeros(M,M);
    B_r = zeros(M,M);
    c_t = zeros(M,1);
    c_r = zeros(M,1);
    
    for i = 1 : K_t
        V = zeros(M,M);
        for j = 1 : K
            v = (H.Ht(:,i)'*diag(H.G*W(:,j)*(aux.t(i))'))';
            V = V + v*v';
        end
        B_t = B_t + V; 
        c_t = c_t + sqrt(1 + aux.r(i))*(H.Ht(:,i)'*diag(H.G*W(:,i)*(aux.t(i))'))';
    end
    for i = 1 : K_r
        V = zeros(M,M);
        for j = 1 : K
            v = (H.Hr(:,i)'*diag(H.G*W(:,j)*(aux.t(K_t + i))'))';
            V = V + v*v';
        end
        B_r = B_r + V; 
        c_r = c_r + sqrt(1 + aux.r(K_t + i))*(H.Hr(:,i)'*diag(H.G*W(:,K_t + i)*(aux.t(K_t + i))'))';
    end    

    
    %% update RIS    
    phy_t_old = zeros(M,1);
    phy_r_old = zeros(M,1);
    phy_t_vec = zeros(M,1);
    phy_r_vec = zeros(M,1);
    for m = 1 : M
        phy_t_vec(m) = phy_t(m,m);
        phy_r_vec(m) = phy_r(m,m);
    end
    nn = 1;
    while norm(phy_t_vec - phy_t_old,2)^2 + norm(phy_r_vec - phy_r_old,2)^2 > eta && nn < SIM
        phy_t_old = phy_t_vec;
        phy_r_old = phy_r_vec;
        for m = 1 : M
            [phy_t_vec(m),phy_r_vec(m)] = golden_section_search(B_t,B_r,c_t,c_r,phy_t_vec,phy_r_vec,m);
        end
        nn = nn + 1;
    end  
    phy_t = diag(phy_t_vec);
    phy_r = diag(phy_r_vec);
end

