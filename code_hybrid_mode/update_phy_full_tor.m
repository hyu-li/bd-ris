function phy = update_phy_full_tor(H,G,W,aux,phy)
    [M,K] = size(H);
   
    X = zeros(M,M);
    Y = zeros(M,M);
    Z = zeros(M,M);

    for i = 1 : K
        X = X + sqrt(1 + aux.r(i))*G*W(:,i)*aux.t(i)'*H(:,i)';
        Z = Z + abs(aux.t(i))^2*H(:,i)*H(:,i)';
        Y = Y + G*W(:,i)*(G*W(:,i))'; 
    end
    normalization_factor = abs(X(1,1));
    X = X./normalization_factor;
    Y = Y./normalization_factor;
    
    % Create the problem structure
    manifold = stiefelcomplexfactory(M,M);
    problem.M = manifold;
    
    % Define the problem function and its Euclidean gradient
    problem.cost = @(x) trace(x*Y*x'*Z) - 2*real(trace(x*X));
    problem.egrad = @(x) 2*Z*x*Y - 2*X';
    
    % Settings
    options.tolgradnorm = 1e-6; %%% tolerance of gradnorm
    options.maxiter = 200; %%% maximum iterations
    options.minstepsize = 1e-6; %%% minimum stepsize
    v_ini = phy;
%     v_ini = 1/sqrt(2)*[diag(exp(1j*2*pi*rand(M,1)));diag(exp(1j*2*pi*rand(M,1)))];
    
    % Solve
    [v,~,~,~] = conjugategradient(problem,v_ini,options);
      
    phy = v;

end

