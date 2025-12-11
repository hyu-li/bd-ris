function [obj_full_tor,R] = FP_fully_connected_tor(H,G,P)
  
    [M,~] = size(H);
    global SIM
    global eta


    %% initialize RIS and digital beamforming
    obj_full_tor.phy = diag(exp(1j*2*pi*rand(M,1)));
    HH = H'*obj_full_tor.phy*G;
    obj_full_tor.w = HH'*pinv(HH*HH');
    obj_full_tor.w = sqrt(P)*obj_full_tor.w/norm(obj_full_tor.w,'fro');
       
    R_old = 0.1;
    [aux,R_new] = sum_rate_cal_tor(H,G,obj_full_tor);    
    sim = 0;
    
    %% BCD update
    while sim < SIM && abs(R_new - R_old)/R_old > eta
        R_old = R_new;
        sim = sim + 1;

        %% update precoder W
        obj_full_tor.w = update_W_tor(H,G,P,aux,obj_full_tor.phy);

        %% update IRS matrix
        obj_full_tor.phy = update_phy_full_tor(H,G,obj_full_tor.w,aux,obj_full_tor.phy); 
        
        %% update auxiliary variables rho and beta
        [aux,R_new] = sum_rate_cal_tor(H,G,obj_full_tor);    
%         R_new
    end 
    
    R = R_new;
end


 
