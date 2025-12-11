clear all
close all

%%----------------------system parameters---------------------------


N_t = 4; %% number of antennas
K_t = 2; %% number of users
K_r = 2;
K = K_t + K_r;
M = 32; %% number of elements of RIS
S = 8;
global sigma  %% noise power
global SIM
global N_max
global eta
sigma = sqrt(1e-11); 
SIM = 20;
N_max = 50;
eta = 1e-3;

d_bi = 50;
d_iu = 2.5;





%%----------------------simulation settings-----------------------


N_sim = 10; %% number of channel realizations

SNR = 0 : 2 : 10; %% range of transmit power

C_hybrid_single = zeros(1,length(SNR));
C_hybrid_group = zeros(1,length(SNR));
C_hybrid_full = zeros(1,length(SNR));
C_t_single = zeros(1,length(SNR));
C_r_single = zeros(1,length(SNR));
C_t_group = zeros(1,length(SNR));
C_r_group = zeros(1,length(SNR));
C_t_full = zeros(1,length(SNR));
C_r_full = zeros(1,length(SNR));



for nn = 1 : N_sim
    
    nn
    
    %% channel   
    H.Ht = zeros(M,K_t);
    H.Hr = zeros(M,K_r);
    for k = 1 : K_t
        H.Ht(:,k) = channel_construction(1,M,d_iu);
    end
    for k = 1 : K_r
        H.Hr(:,k) = channel_construction(1,M,d_iu);
    end
%     H.Ht = channel_construction(K_t,M,d_iu);
%     H.Hr = channel_construction(K_r,M,d_iu);
    H.G = channel_construction(N_t,M,d_bi);

    for n = 1 : length(SNR)
        
        P = 10^((SNR(n) - 30)/10);

     
        %% fully-connected, hybrid 
        [obj_h_f,R_full] = FP_fully_connected(H,P);
        C_hybrid_full(n) = C_hybrid_full(n) + R_full; 
        
        %% group-connected, hybrid 
        [obj_h_g,R_group] = FP_group_connected(H,P,S);
        C_hybrid_group(n) = C_hybrid_group(n) + R_group; 
        
        %% single-connected, hybrid 
        [obj_h_s,R_single] = FP_single_connected(H,P);
        C_hybrid_single(n) = C_hybrid_single(n) + R_single;  

        %% fully-connected, transmissive 
        [obj_t_f,R_t_full] = FP_fully_connected_tor(H.Ht,H.G,P);
        C_t_full(n) = C_t_full(n) + R_t_full;
        
        %% group-connected, transmissive 
        [obj_t_g,R_t_group] = FP_group_connected_tor(H.Ht,H.G,P,S);
        C_t_group(n) = C_t_group(n) + R_t_group;   
        
        %% single-connected, transmissive 
        [obj_t_s,R_t_single] = FP_single_connected_tor(H.Ht,H.G,P);
        C_t_single(n) = C_t_single(n) + R_t_single;          

        %% fully-connected, reflective 
        [obj_r_f,R_r_full] = FP_fully_connected_tor(H.Hr,H.G,P);
        C_r_full(n) = C_r_full(n) + R_r_full; 
        
        %% group-connected, reflective 
        [obj_r_g,R_r_group] = FP_group_connected_tor(H.Hr,H.G,P,S);
        C_r_group(n) = C_r_group(n) + R_r_group;   
        
        %% single-connected, reflective 
        [obj_r_s,R_r_single] = FP_single_connected_tor(H.Hr,H.G,P);
        C_r_single(n) = C_r_single(n) + R_r_single;        
    end
    
end
 

C_hybrid_single = C_hybrid_single/N_sim;
C_hybrid_group = C_hybrid_group/N_sim;
C_hybrid_full = C_hybrid_full/N_sim;
C_t_single = C_t_single/N_sim;
C_r_single = C_r_single/N_sim;
C_t_group = C_t_group/N_sim;
C_r_group = C_r_group/N_sim;
C_t_full = C_t_full/N_sim;
C_r_full = C_r_full/N_sim;



figure;
plot(SNR,C_hybrid_full,'-rs','LineWidth',1.5);
hold on
plot(SNR,C_hybrid_group,'--rs','LineWidth',1.5);
hold on
plot(SNR,C_hybrid_single,'-.rs','LineWidth',1.5);
hold on
plot(SNR,C_t_full,'-b>','LineWidth',1.5);
hold on
plot(SNR,C_t_group,'--b>','LineWidth',1.5);
hold on
plot(SNR,C_t_single,'-.b>','LineWidth',1.5);
hold on
plot(SNR,C_r_full,'-ko','LineWidth',1.5);
hold on
plot(SNR,C_r_group,'--ko','LineWidth',1.5);
hold on
plot(SNR,C_r_single,'-.ko','LineWidth',1.5);
hold off


xlabel('{\it P}(dBm)');
ylabel('Sum-rate(b/s/Hz)');
legend('Hybrid, Full','Hybrid, Group','Hybrid, Single','Transmissive, Full',...
    'Transmissive, Group','Transmissive, Single','Reflective, Full',...
    'Reflective, Group','Reflective, Single');
grid on;
set(gca, 'GridLineStyle', ':','GridAlpha', 1);



