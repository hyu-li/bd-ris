function [obj,R] = FP_single_connected(H,P)
  
    [M,~] = size(H.Ht);
    global SIM
    global eta


    %% initialize RIS and digital beamforming
    obj.phy_r = 1/sqrt(2)*diag(exp(1j*2*pi*rand(M,1)));
    obj.phy_t = 1/sqrt(2)*diag(exp(1j*2*pi*rand(M,1)));
%     obj.phy_r = 1/sqrt(2)*diag(ones(M,1));
%     obj.phy_t = 1/sqrt(2)*diag(ones(M,1));
    Htr = blkdiag(H.Ht,H.Hr);
    Phy = [obj.phy_t;obj.phy_r];
    HH = Htr'*Phy*H.G;
    obj.w = HH'*inv(HH*HH');
    obj.w = sqrt(P)*obj.w/norm(obj.w,'fro');
       
    R_old = 0.1;
    [aux,R_new] = sum_rate_cal(H,obj);    
    sim = 0;
    
    %% BCD update
    R = zeros(1,SIM);
    while sim < SIM && abs(R_new - R_old)/R_old > eta
        R_old = R_new;
        sim = sim + 1;

        %% update precoder W
        obj.w = update_W(H,P,aux,obj.phy_t,obj.phy_r);

        %% update IRS matrix
        [obj.phy_t,obj.phy_r] = update_phy_single(H,obj.w,aux,obj.phy_t,obj.phy_r); 
        
        %% update auxiliary variables rho and beta
        [aux,R_new] = sum_rate_cal(H,obj);  
%         R_single = R_new
%         R(sim) = R_new;
    end   
    R = R_new;
end


 
