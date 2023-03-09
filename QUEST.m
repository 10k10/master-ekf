function [q,exe] = QUEST(a,b,n,method)
    n_vec = size(a,2);
    
    m = 0;
    B = zeros(3);
    
    for i=1:n_vec
        B = B + a(i)*b(:,i)*n(:,i)';
        m = m + a(i);
     
    end
    
    Z = [B(2,3)-B(3,2); B(3,1)-B(1,3); B(1,2)-B(2,1)];
    sigma = trace(B);
    
    % Normalize using weights
    sigma = (1/m).*sigma;
    B = (1/m).*(B);
    Z = (1/m).*(Z);    
    S = B + B';
        
    if strcmp(method,'QUEST-nr')
        tic
        
        s = sigma;      
        a = s^2 -trace(adjoint(S));
        b = s^2 + Z'*Z;
        c = det(S) + Z'*S*Z;
        d = Z'*S^2*Z;
    
        iter = 10;
        l = m; % Palpite inicial é o somatório dos pesos (se normalizado = 1)
        l_opt = 0;

        % Newton-Raphson method to find the max eigen value
        for k=1:iter
            l_opt = l - (l^4 -(a+b)*l^2 -c*l + (a*b + c*s -d))/(4*l^3 -2*(a+b)*l -c);
            l = l_opt;
        end
        
        lambda = l_opt;
       
        y = ((lambda + sigma).*eye(3)-S)^-1*Z;
        out = ((1/sqrt(1+norm(y)^2)).*[1;y])';
    end
    
    if strcmp(method,'Q-method')
        tic
        
        K = [sigma      Z';
                Z   S-sigma.*eye(3)];
      
        [V,D] = eig(K,'nobalance');   
        [~,id] = max(diag(D));
        
        out = V(:,id);
    end
        
    if strcmp(method,'QUEST-eig')
        tic
        
         K = [sigma      Z';
          Z   S-sigma.*eye(3) ];
      
        [~,D] = eig(K,'nobalance');   
        [lambda,~] = max(diag(D));
         
        y = ((lambda + sigma).*eye(3)-S)^-1*Z;
        out = ((1/sqrt(1+norm(y)^2)).*[1;y])';
    end   
     q = out;
     exe = toc;
end