function [obj_single_tor,R] = FP_single_connected_tor(H,G,P)
  
    [M,~] = size(H);
    global SIM
    global eta


    %% initialize RIS and digital beamforming
    obj_single_tor.phy = diag(exp(1j*2*pi*rand(M,1)));
    HH = H'*obj_single_tor.phy*G;
    obj_single_tor.w = HH'*pinv(HH*HH');
    obj_single_tor.w = sqrt(P)*obj_single_tor.w/norm(obj_single_tor.w,'fro');
       
    R_old = 0.1;
    [aux,R_new] = sum_rate_cal_tor(H,G,obj_single_tor);    
    sim = 0;
    
    %% BCD update
    while sim < SIM && abs(R_new - R_old)/R_old > eta
        R_old = R_new;
        sim = sim + 1;

        %% update precoder W
        obj_single_tor.w = update_W_tor(H,G,P,aux,obj_single_tor.phy);

        %% update IRS matrix
        obj_single_tor.phy = update_phy_single_tor(H,G,obj_single_tor.w,aux,obj_single_tor.phy); 
        
        %% update auxiliary variables rho and beta
        [aux,R_new] = sum_rate_cal_tor(H,G,obj_single_tor);        
    end   
    R = R_new;
end


 
