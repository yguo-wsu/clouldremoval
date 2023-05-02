function [ L, C] = bcs_exact( D, lambda, tol, maxIter, method, verbose)
% Minimising the following:
%      |L|_* + \lambda|C|_1
% s.t. X = C + (1-C) .* L
%      C,L \in [0,1]
%
% Minimisation uses linearisation for nuclear norm with approximal and feasible set projection.
% J = |L|_* + \lambda|C|_1 +\mu/2 |X - C - (1-C) .* L + Y/\mu|_F^2
% dJ_Ck = \lambda sign(C) + \mu(L-1) .* (L-1) .* C + \mu (X - L + Y/\mu) .* (L-1)
% f(Lk) = \mu/2 |X - C - (1-C) .* Lk + Y/\mu|_F^2
% df_Lk = -\mu(X - C - (1-C) .* Lk + Y/\mu) .* (1-C)
MAXITERS=1000;
if ~exist('verbose', 'var')
    verbose = 0;
end

if ~exist('method', 'var')
    method = 1;
end

if ~exist('lambda', 'var')
    lambda = 1 / sqrt(m);
end

if ~exist('tol', 'var')
    tol = 1e-7;
end

if ~exist('maxIter', 'var')
    maxIter = MAXITERS;
end

switch method
        case 1
        %% Version 1
        [m n] = size(D);
      
        % initialize
        Y = D;
        norm_two = norm(Y, 2);
        norm_inf = norm( Y(:), inf) / lambda;
        dual_norm = max(norm_two, norm_inf);
        Y = Y / dual_norm;
        
        L = zeros( m, n);%D;
        C = zeros( m, n);
        mu = 1.25/norm_two; % this one can be tuned
        mu_bar = mu * 1e7;
%         mu = mu_bar/2; 
        rho = 1.5 ;         % this one can be tuned
        d_norm = norm(D, 'fro');
        rho_lineariser = nuclearnorm(D); %d_norm;
        
        iter = 0;
        total_svd = 0;
        converged = false;
        stopCriterion = 1;
        while ~converged
            iter = iter + 1;
            % Update C
            if 0
                df_Ck = mu*(D-C-(1-C).*L+Y/mu).*(L-1);
                tmp_H = C - df_Ck/rho_lineariser;
                tmp_A = ones(size(tmp_H))*lambda;
            else
                tmp_A = mu*(L-1).*(L-1);
                tmp_A(tmp_A==0) = eps;
                tmp_H = -mu*(D - L + Y/mu).*(L-1)./tmp_A;
                tmp_A = lambda./tmp_A;
            end
            C = max(tmp_H - tmp_A, 0) + min(tmp_H + tmp_A, 0);
            C(C<0) = 0;
            C(C>1) = 1;
            
            % Update L
            L = PAPGnuclearhadamard(D, C, mu, Y, L); 
%             L = PAPGavgnuclearhadamard(D, C, mu, Y, L); 
            L(L<0) = 0;
            L(L>1) = 1;
            total_svd = total_svd + 1; 
            Z = D - C - (1-C).*L ;
            
            Y = Y + mu*Z;
            mu = min(mu*rho, mu_bar);
            
            %% stop Criterion
            stopCriterion = norm(Z, 'fro') / d_norm;
            [U S V] = svd(L, 'econ');
            funval = sum(diag(S)) + lambda * sum(abs(C(:)));
            if stopCriterion < tol
                converged = true;
            end
            
            if verbose%mod(iter,5)==0
                disp(['#svd ' num2str(total_svd) ', r(L) ' num2str(rank(L))...
                    ', |C|_0 ' num2str(length(find(abs(C)>0))/length(C(:)))...
                    ', Obj ' num2str(funval), ', EQ: ',num2str(stopCriterion), ', mu: ',num2str(mu)]);
            end
            
            if ~converged && iter >= maxIter
                if verbose
                    disp('Maximum iterations reached');
                end
                converged = 1 ;
            end
        end

    case 2
        %% Version 2
        [m n] = size(D);
      
        % initialize
        Y = D;
        norm_two = norm(Y, 2);
        norm_inf = norm( Y(:), inf) / lambda;
        dual_norm = max(norm_two, norm_inf);
        Y = Y / dual_norm;
        
        L = D;
        C = zeros( m, n);
        mu = 1.25/norm_two; % this one can be tuned
        mu_bar = mu * 1e7;
%         mu = mu_bar/2; 
        rho = 1.5 ;         % this one can be tuned
        d_norm = norm(D, 'fro');
        rho_lineariser = nuclearnorm(D); %d_norm;
        
        iter = 0;
        total_svd = 0;
        converged = false;
        stopCriterion = 1;
        while ~converged
            iter = iter + 1;
            % Update C
            if 0
                df_Ck = mu*(D-C-(1-C).*L+Y/mu).*(L-1);
                tmp_H = C - df_Ck/rho_lineariser;
                tmp_A = ones(size(tmp_H))*lambda;
            else
                tmp_A = mu*(L-1).*(L-1);
                tmp_A(tmp_A==0) = eps;
                tmp_H = -mu*(D - L + Y/mu).*(L-1)./tmp_A;
                tmp_A = lambda./tmp_A;
            end
            C = max(tmp_H - tmp_A, 0) + min(tmp_H + tmp_A, 0);
            C(C<0) = 0;
            C(C>1) = 1;
            
            % Update L
            Df_Lk = - mu*(D-C-(1-C).*L+Y/mu).*(1-C);
            [U S V] = svd(L - Df_Lk/rho_lineariser, 'econ');
            diagS = diag(S);
            svp = length(find(diagS > 1/rho_lineariser));
            L = U(:, 1:svp) * diag(diagS(1:svp) - 1/rho_lineariser) * V(:, 1:svp)';
            L(L<0) = 0;
            L(L>1) = 1;
            total_svd = total_svd + 1;
            
            Z = D - C - (1-C).*L ;
            
            Y = Y + mu*Z;
            mu = min(mu*rho, mu_bar);
            
            %% stop Criterion
            stopCriterion = norm(Z, 'fro') / d_norm;
            [U S V] = svd(L, 'econ');
            funval = sum(diag(S)) + lambda * sum(abs(C(:)));
            if stopCriterion < tol
                converged = true;
            end
            
            if 1%mod(iter,5)==0
                disp(['#svd ' num2str(total_svd) ', r(L) ' num2str(rank(L))...
                    ', |C|_0 ' num2str(length(find(abs(C)>0))/length(C(:)))...
                    ', Obj ' num2str(funval), ', EQ: ',num2str(stopCriterion), ', mu: ',num2str(mu)]);
            end
            
            if ~converged && iter >= maxIter
                disp('Maximum iterations reached') ;
                converged = 1 ;
            end
        end
        
        
        %% Version 3
    case 3
        max_iterations = 1000;
        
        func_vals = zeros(max_iterations, 1);
        
        L = D;
        C = zeros(size(D));
        Y = D;
        norm_two = norm(Y, 2);
        norm_inf = norm( Y(:), inf) / lambda;
        dual_norm = max(norm_two, norm_inf);
        Y = Y / dual_norm;
        
        mu = 0.5;
        mu_max = 1e+6;
        rho_lineariser = log(norm(D,2)) * 1.2;
        gamma_0 = 1.1;
        
        tol_1 = 1*10^-6;
        tol_2 = 1*10^-6;
        
        normfX = norm(D,'fro');
        
        for k = 1 : max_iterations
            C_prev = C;
            L_prev = L;
            
            % Update C
            tmp_A = mu*(L-1).*(L-1);
            tmp_A(tmp_A==0) = eps;
            tmp_H = -mu*(D - L + Y/mu).*(L-1)./tmp_A;
            tmp_A = lambda./tmp_A;
            C = max(tmp_H - tmp_A, 0) + min(tmp_H + tmp_A, 0);
            C(C<0) = 0;
            C(C>1) = 1;
            
            % Update L
            Df_Lk = - mu*(D-C-(1-C).*L+Y/mu).*(1-C);
            [U S V] = svd(L - Df_Lk/rho_lineariser, 'econ');
            diagS = diag(S);
            svp = length(find(diagS > 1/rho_lineariser));
            L = U(:, 1:svp) * diag(diagS(1:svp) - 1/rho_lineariser) * V(:, 1:svp)';
            L(L<0) = 0;
            L(L>1) = 1;
            
            % Check convergence
            [U S V] = svd(L, 'econ');
            func_vals(k) = lambda*norm(C, 1) + sum(diag(S));
            
            stop_2 = mu * max([ sqrt(rho_lineariser)*norm(L - L_prev,'fro'), norm(C - C_prev, 'fro')]) / normfX;
            
            if ( (norm(D-C-(1-C).*L, 'fro')/normfX) < tol_1 ...
                    && stop_2 < tol_2)
%                 break;
            end
            
            % Update Y
            Y = Y + mu*(D-C-(1-C).*L);
            
            % Uppdate rate parameters
            if (stop_2 < tol_2)
                gamma = gamma_0;
            else
                gamma = 1;
            end
            
            mu = min(mu_max, gamma * mu);
             if mod(k,50)==0
            disp(['Iter: ', num2str(k), ', Obj: ',num2str(func_vals(k))])
             end
        end
        
end
