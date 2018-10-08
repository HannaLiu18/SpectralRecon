% Transfer TV regularization to a convex problem and solve it by
% quasi-Newton
% norm(A*X_Vector - b) + tao * sqrt(Dx^2 + gamma) +  tao * sqrt(Dy^2 + gamma)

function [X,Obj,XTotal] = quasinewton_TV...
    (A,b,X0,X0_Square,NImage,gamma,tao,NIter,tol)
    
    k = 1;X = X0; X_Square = X0_Square;
    Diff = (A * X - b);
    G = 2 * A' * Diff;
    G_TV = computeGradient_TV(X_Square,NImage,gamma);
    G = tao * G + G_TV;
    
    N = size(A,2);
    H = opEye(N);
    XTotal = zeros(length(X0),NIter);
    
    while ((k <= NIter))% && (norm(G) > tol))
        % k
        % 1) Compute search direction
        Direction = -H * G;
        % 2) Compute step size  
        % better than fminunc, because steepest descent is using here.
%         temp = A * Direction;
%         alpha = (G' * H * G) / (2 * temp' * temp);
        alpha = 0.01;    % search step, back line search
        
        % Backtracking line search
        t = 1; alpha = 0.4; beta = 0.6;
        F_One = Func_Obj(X + t * Direction, A,b,tao,NImage,gamma);
        F_Two = Func_Obj(X, A,b,tao,NImage,gamma);
        F_Two_Inc = G' * Direction;
        while (F_One > F_Two + F_Two_Inc * alpha * t)
            t = beta * t;
            F_One = Func_Obj(X + t * Direction, A,b,tao,NImage,gamma);
        end
        alpha = t;
        
        % 3) Update X
        XOld = X;
        X = X + alpha * Direction;
        X_Square = reshape(X,NImage,NImage);
        % 4) Update Hessian-Inverse matrix
        GOld = G;
        Diff = (A * X - b);
        G = 2 * A' * Diff;
        G_TV = computeGradient_TV(X_Square,NImage,gamma);
        G = tao * G + G_TV;
        
        deltaGrad = opMatrix(G - GOld);
        % HdeltaGrad = (H*deltaGrad);
        deltaX = opMatrix(X - XOld);
        % deltaXDeltaGrad = double(deltaX'*deltaGrad);
        H = (opEye(N) - (deltaX * deltaGrad')/(deltaGrad' * deltaX)) ...
            * H * (opEye(N) - (deltaGrad * deltaX')/(deltaGrad' * deltaX))+...
            (deltaX * deltaX')/(deltaGrad' * deltaX);
        Obj(k) = Func_Obj(X, A,b,tao,NImage,gamma);
        % Func_Obj(X, A,b,tao,NImage,0)
        XTotal(:,k) = X;
        k = k + 1;    
    end
end

function F = Func_Obj(X, A,b,tao,NImage,gamma)
    Diff = A * X - b;
    F = Diff' * Diff * tao;
    
    X_Square = reshape(X,NImage,NImage);
    X_extend = zeros(NImage + 2, NImage + 2);
    X_extend(2:end-1, 2:end-1) = X_Square;
    X_extend(1,:) = X_extend(2,:);
    X_extend(end,:) = X_extend(end-1,:);
    X_extend(:,1) = X_extend(:,2);
    X_extend(:,end) = X_extend(:,end-1);  
    
    G = zeros(NImage, NImage);
    
    % ANISOTROPIC
    Temp = X_extend(3:end,2:end-1) - X_extend(2:end-1,2:end-1);
    G = G + sqrt(Temp.^2 + gamma);
    Temp = X_extend(2:end-1,3:end) - X_extend(2:end-1,2:end-1);
    G = G + sqrt(Temp.^2 + gamma);
    F = sum(G(:)) + F;

%     % ISOTROPIC
%     Temp = X_extend(3:end,2:end-1) - X_extend(2:end-1,2:end-1);
%     G = G + Temp.^2;
%     Temp = X_extend(2:end-1,3:end) - X_extend(2:end-1,2:end-1);
%     G = G + Temp.^2 + gamma;
%     F = sum(sqrt(G(:))) + F;
    
end

% Compute gradient for approximate TV
function G = computeGradient_TV(X,NImage,gamma)
    X_extend = zeros(NImage + 2, NImage + 2);
    X_extend(2:end-1, 2:end-1) = X;
    X_extend(1,:) = X_extend(2,:);
    X_extend(end,:) = X_extend(end-1,:);
    X_extend(:,1) = X_extend(:,2);
    X_extend(:,end) = X_extend(:,end-1);  
    
    G = zeros(NImage, NImage);
    
    % ISOTROPIC regu is smaller than ANISOTROPIC, therefore requires
    % smaller lamda for data fidity term.
    
%     % ISOTROPIC
%     Temp_One = X_extend(2:end-1,1:end-2) - X_extend(3:end,1:end-2);
%     Temp_Two = X_extend(2:end-1,2:end-1) - X_extend(2:end-1,1:end-2);    
% 
%     Temp_Three = X_extend(2:end-1,2:end-1) - X_extend(1:end-2,2:end-1);    
%     Temp_Four = X_extend(1:end-2,2:end-1) - X_extend(1:end-2,3:end);   
%     
%     Temp_Five = X_extend(2:end-1,2:end-1) - X_extend(3:end,2:end-1);  
%     Temp_Six = X_extend(2:end-1,2:end-1) - X_extend(2:end-1,3:end);  
%     
%     G = G + Temp_Two./sqrt(Temp_One.^2 + Temp_Two.^2 + gamma);
%     G = G + Temp_Three./sqrt(Temp_Three.^2 + Temp_Four.^2 + gamma);
%     G = G + (Temp_Five + Temp_Six)./sqrt(Temp_Five.^2 + Temp_Six.^2 + gamma);

    % ANISOTROPIC
    Temp = X_extend(2:end-1,2:end-1) - X_extend(1:end-2,2:end-1);
    G = G + Temp./sqrt(Temp.^2 + gamma);
    Temp = X_extend(2:end-1,2:end-1) - X_extend(3:end,2:end-1);
    G = G + Temp./sqrt(Temp.^2 + gamma);
    Temp = X_extend(2:end-1,2:end-1) - X_extend(2:end-1,1:end-2);
    G = G + Temp./sqrt(Temp.^2 + gamma);
    Temp = X_extend(2:end-1,2:end-1) - X_extend(2:end-1,3:end);
    G = G + Temp./sqrt(Temp.^2 + gamma);
    
    G = G(:);
end