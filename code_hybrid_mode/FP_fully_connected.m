function [obj_full,R] = FP_fully_connected(H,P)
  
    [M,~] = size(H.Ht);
    global SIM
    global eta


    %% initialize RIS and digital beamforming
    obj_full.phy_r = 1/sqrt(2)*diag(exp(1j*2*pi*rand(M,1)));
    obj_full.phy_t = 1/sqrt(2)*diag(exp(1j*2*pi*rand(M,1)));
%     obj_full.phy_r = phy_r_ini;
%     obj_full.phy_t = phy_t_ini;
    Htr = blkdiag(H.Ht,H.Hr);
    Phy = [obj_full.phy_t;obj_full.phy_r];
    HH = Htr'*Phy*H.G;
    obj_full.w = HH'*inv(HH*HH');
    obj_full.w = sqrt(P)*obj_full.w/norm(obj_full.w,'fro');
       
    R_old = 0.1;
    [aux_full,R_new] = sum_rate_cal(H,obj_full);    
    sim = 0;
    
    %% BCD update
    R = zeros(1,SIM);
    while sim < SIM && abs(R_new - R_old)/R_old > eta
        R_old = R_new;
        sim = sim + 1;

        %% update precoder W        
        obj_full.w = update_W(H,P,aux_full,obj_full.phy_t,obj_full.phy_r);
        
        %% update IRS matrix
        [obj_full.phy_t,obj_full.phy_r] = update_phy_full(H,obj_full.w,aux_full,obj_full.phy_t,obj_full.phy_r); 
        
        %% update auxiliary variables rho and beta
        [aux_full,R_new] = sum_rate_cal(H,obj_full); 
%         R(sim) = R_new;
%         R_full = R_new
    end   
    R = R_new;
end


 
