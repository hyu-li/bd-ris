function [phy_t,phy_r] = update_phy_group(H,W,aux,S,phy_t,phy_r)
    [M,K_t] = size(H.Ht);
    [~,K_r] = size(H.Hr);
    global SIM
    global eta
    K = K_t + K_r;
    
    X_t = zeros(M,M);
    X_r = zeros(M,M);
    Y = zeros(M,M);
    Z_t = zeros(M,M);
    Z_r = zeros(M,M);
    for i = 1 : K_t
        X_t = X_t + sqrt(1 + aux.r(i))*H.G*W(:,i)*aux.t(i)'*H.Ht(:,i)';
        Z_t = Z_t + abs(aux.t(i))^2*H.Ht(:,i)*H.Ht(:,i)';
    end
    for i = 1 : K_r
        X_r = X_r + sqrt(1 + aux.r(K_t + i))*H.G*W(:,K_t + i)*...
            aux.t(K_t + i)'*H.Hr(:,i)';
        Z_r = Z_r + abs(aux.t(K_t + i))^2*H.Hr(:,i)*H.Hr(:,i)';
    end  
    for j = 1 : K
        Y = Y + H.G*W(:,j)*(H.G*W(:,j))';       
    end

    phy_t_old = sqrt(1/2)*eye(M);
    phy_r_old = sqrt(1/2)*eye(M);
    nn = 1;
    while norm(phy_t - phy_t_old,'fro')^2 + norm(phy_r - phy_r_old,'fro')^2 > eta && nn < SIM
        phy_t_old = phy_t;
        phy_r_old = phy_r;
       
        for s = 1 : S
            
            A_s_t = X_t((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S); 
%             - Y((s-1)*M/S+1:s*M/S,:)*...
%                   phy_t'*Z_t(:,(s-1)*M/S+1:s*M/S) + Y((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S)*...
%                   phy_t((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S)'*Z_t((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S)
            A_s_r = X_r((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S);
%             - Y((s-1)*M/S+1:s*M/S,:)*...
%                   phy_r'*Z_r(:,(s-1)*M/S+1:s*M/S) + Y((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S)*...
%                   phy_r((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S)'*Z_r((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S)
            A_s = [A_s_t,A_s_r];
            A_s = A_s./abs(X_t(1,1));
            Y_s = Y((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S);
            Y_s = Y_s./abs(X_t(1,1));
            Z_s = blkdiag(Z_t((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S),...
                Z_r((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S));
            
            
            % Create the problem structure
            manifold = stiefelcomplexfactory(2*M/S,M/S);
            problem.M = manifold;
            
            % Define the problem function and its Euclidean gradient
            problem.cost = @(x) trace(x*Y_s*x'*Z_s) - 2*real(trace(x*A_s));
            problem.egrad = @(x) 2*Z_s*x*Y_s - 2*A_s';
    
            % Settings
            options.tolgradnorm = 1e-6; %%% tolerance of gradnorm
            options.maxiter = 100; %%% maximum iterations
            options.minstepsize = 1e-6; %%% minimum stepsize
            
%             v_ini = 1/sqrt(2)*[diag(exp(1j*2*pi*rand(M/S,1)));diag(exp(1j*2*pi*rand(M/S,1)))];
            v_ini = [phy_t((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S);...
                phy_r((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S)];

            % Solve
            [v,~,~,~] = conjugategradient(problem,v_ini,options);
    
            phy_t((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S) = v(1:M/S,:);
            phy_r((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S) = v(M/S+1:2*M/S,:);
        end
        nn = nn + 1;
    end  
    

end

