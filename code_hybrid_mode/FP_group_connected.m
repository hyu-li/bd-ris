function [obj_group,R] = FP_group_connected(H,P,S)
 
    [M,N_t] = size(H.G);
    global SIM
    global sigma


    %% initialize RIS and digital beamforming
    obj_group.phy_r = 1/sqrt(2)*diag(exp(1j*2*pi*rand(M,1)));
    obj_group.phy_t = 1/sqrt(2)*diag(exp(1j*2*pi*rand(M,1)));
%     obj_group.phy_r = 1/sqrt(2)*diag(ones(M,1));
%     obj_group.phy_t = 1/sqrt(2)*diag(ones(M,1));
    Htr = blkdiag(H.Ht,H.Hr);
    Phy = [obj_group.phy_t;obj_group.phy_r];
    HH = Htr'*Phy*H.G;
    obj_group.w = inv(HH'*HH + sigma^2*eye(N_t))*HH';
    obj_group.w = sqrt(P)*obj_group.w/norm(obj_group.w,'fro');
       
    R_old = 0.1;
    [aux_group,R_new] = sum_rate_cal(H,obj_group);    
    sim = 0;
    
    %% BCD update
%     R = zeros(1,SIM);
    while sim < SIM && abs(R_new - R_old)/R_old > 1e-2
        R_old = R_new;
        sim = sim + 1;

        %% update precoder W
        obj_group.w = update_W(H,P,aux_group,obj_group.phy_t,obj_group.phy_r);

        %% update IRS matrix
        [obj_group.phy_t,obj_group.phy_r] = update_phy_group(H,obj_group.w,aux_group,S,obj_group.phy_t,obj_group.phy_r); 
        
        %% update auxiliary variables rho and beta
        [aux_group,R_new] = sum_rate_cal(H,obj_group); 
%         R_group = R_new;
%         R(sim) = R_new;
    end   
    R = R_new;
end


 
