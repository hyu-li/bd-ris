function H = channel_construction(M,N,d)

    K_TR = 10^(5/10);
    PL0 = 1e-3;
    d0 = 1;
    beta = 2.2;
    PL = sqrt(PL0*(d/d0)^(-beta));

    % LoS channel
    ind_Tx = [0 : 1: M - 1]';
    ind_Rx = [0 : 1: N - 1]'; 
    AoD = pi*rand - pi/2;
    AoA = pi*rand - pi/2;
    alpha = 1;
%     alpha = sqrt(1/2)*(randn + 1j*randn);
    a_T = sqrt(1/M)*exp(1j*pi*ind_Tx*cos(AoD));
    a_R = sqrt(1/N)*exp(1j*pi*ind_Rx*cos(AoA));
    H_LoS = sqrt(M*N)*alpha*a_R*a_T';

    % NLos channel
    H_NLoS = sqrt(1/2)*(randn(N,M) + 1j*randn(N,M));

%     H = PL*(sqrt(K_TR/(K_TR + 1))*H_LoS + sqrt(1/(K_TR + 1))*H_NLoS);
    H = PL*H_NLoS;
%     H = PL*H_LoS;
end