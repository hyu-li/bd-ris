function [phy_t,phy_r] = update_phy_full(H,W,aux,phy_t,phy_r)
    [M,K_t] = size(H.Ht);
    [~,K_r] = size(H.Hr);
    
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
    X = [X_t,X_r];
    Z = blkdiag(Z_t,Z_r);
    for j = 1 : K
        Y = Y + H.G*W(:,j)*(H.G*W(:,j))';       
    end
    
    % Create the problem structure
    manifold = stiefelcomplexfactory(2*M,M);
    problem.M = manifold;
    
    % Define the problem function and its Euclidean gradient
    problem.cost = @(x) trace(x*Y*x'*Z) - 2*real(trace(x*X));
    problem.egrad = @(x) 2*Z*x*Y - 2*X';
    
    % Settings
    options.tolgradnorm = 1e-6; %%% tolerance of gradnorm
    options.maxiter = 100; %%% maximum iterations
    options.minstepsize = 1e-6; %%% minimum stepsize
    v_ini = [phy_t;phy_r];
%     v_ini = 1/sqrt(2)*[diag(exp(1j*2*pi*rand(M,1)));diag(exp(1j*2*pi*rand(M,1)))];
    
    % Solve
    [v,~,~,~] = conjugategradient(problem,v_ini,options);
%     [v,~,~,~] = rlbfgs(problem,v_ini,options);  
    phy_t = v(1:M,:);
    phy_r = v(M+1:2*M,:);
end

