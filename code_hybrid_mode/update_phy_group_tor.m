function phy = update_phy_group_tor(H,G,W,aux,phy,S)
    [M,K] = size(H);
    
    global SIM
    global eta
    
    X = zeros(M,M);
    Y = zeros(M,M);
    Z = zeros(M,M);
    for i = 1 : K
        X = X + sqrt(1 + aux.r(i))*G*W(:,i)*aux.t(i)'*H(:,i)';
        Z = Z + abs(aux.t(i))^2*H(:,i)*H(:,i)';
        Y = Y + G*W(:,i)*(G*W(:,i))';  
    end 
       
%     phy_old = eye(M);
%     nn = 1;
%     while norm(phy - phy_old,'fro')^2 > eta && nn < SIM
%         phy_old = phy;
        for s = 1 : S
            
            A_s = X((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S); %- Y((s-1)*M/S+1:s*M/S,:)*...
%                phy'*Z(:,(s-1)*M/S+1:s*M/S) + Y((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S)*...
%                phy((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S)'*Z((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S);
            Y_s = Y((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S);
            Z_s = Z((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S);
            A_s = A_s./abs(X(1,1));
            Y_s = Y_s./abs(X(1,1));
                        
            % Create the problem structure
            manifold = stiefelcomplexfactory(M/S,M/S);
            problem.M = manifold;
            
            % Define the problem function and its Euclidean gradient
            problem.cost = @(x) trace(x*Y_s*x'*Z_s) - 2*real(trace(x*A_s));
            problem.egrad = @(x) 2*Z_s*x*Y_s - 2*A_s';

            % Settings
            options.tolgradnorm = 1e-6; %%% tolerance of gradnorm
            options.maxiter = 200; %%% maximum iterations
            options.minstepsize = 1e-6; %%% minimum stepsize
            v_ini = phy((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S);
%             v_ini = diag(exp(1j*2*pi*rand(M/S,1)));
            % Solve
            [v,~,~,~] = conjugategradient(problem,v_ini,options);

            phy((s-1)*M/S+1:s*M/S,(s-1)*M/S+1:s*M/S) = v;
        end
%         nn = nn + 1;
%     end  
end

