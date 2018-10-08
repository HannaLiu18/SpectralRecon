function [X,Obj,XTotal] = quasinewton(A,b,X0,NIter,tol)
    % minimize norm(Ax-b)
    % quasi-newton
    
    k = 1;X = X0;
    Diff = (A * X - b);
    G = 2 * A' * Diff;
    M = size(A,1);
    N = size(A,2);
    H = opEye(N);
    XTotal = zeros(length(X0),NIter);
    
    while ((k <= NIter))% && (norm(G) > tol))
        % k
        % 1) Compute search direction
        Direction = -H * G;
        % 2) Compute step size  
        % better than fminunc, because steepest descent is using here.
        temp = A * Direction;
        alpha = (G' * H * G) / (2 * temp' * temp);
        % 3) Update X
        XOld = X;
        X = X + alpha * Direction;
        % 4) Update Hessian-Inverse matrix
        GOld = G;
        Diff = (A * X - b);
        G = 2 * A' * Diff;
        
        deltaGrad = opMatrix(G - GOld);
        % HdeltaGrad = (H*deltaGrad);
        deltaX = opMatrix(X - XOld);
        % deltaXDeltaGrad = double(deltaX'*deltaGrad);
        H = (opEye(N) - (deltaX * deltaGrad')/(deltaGrad' * deltaX)) ...
            * H * (opEye(N) - (deltaGrad * deltaX')/(deltaGrad' * deltaX))+...
            (deltaX * deltaX')/(deltaGrad' * deltaX);
        Obj(k) = Diff'*Diff;
        XTotal(:,k) = X;
        k = k + 1;    
    end
end