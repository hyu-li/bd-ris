function [obj_group_tor,R] = FP_group_connected_tor(H,G,P,S)
  
    [M,N_t] = size(G);
    global SIM
    global eta
    global sigma


    %% initialize RIS and digital beamforming
    obj_group_tor.phy = diag(exp(1j*2*pi*rand(M,1)));
%     obj_group_tor.phy = diag(ones(M,1));
    HH = H'*obj_group_tor.phy*G;
    obj_group_tor.w = pinv(HH'*HH + sigma^2*eye(N_t))*HH';
    obj_group_tor.w = sqrt(P)*obj_group_tor.w/norm(obj_group_tor.w,'fro');
       
    R_old = 0.1;
    [aux,R_new] = sum_rate_cal_tor(H,G,obj_group_tor);    
    sim = 0;
    
    %% BCD update
    while sim < SIM && abs(R_new - R_old)/R_old > eta
        R_old = R_new;
        sim = sim + 1;

        %% update precoder W
        obj_group_tor.w = update_W_tor(H,G,P,aux,obj_group_tor.phy);

        %% update IRS matrix
        obj_group_tor.phy = update_phy_group_tor(H,G,obj_group_tor.w,aux,obj_group_tor.phy,S); 
        
        %% update auxiliary variables rho and beta
        [aux,R_new] = sum_rate_cal_tor(H,G,obj_group_tor);   
        R_full = R_new;
    end   
    R = R_new;
end


 
