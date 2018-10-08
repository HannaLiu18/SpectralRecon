function [X,Obj] = gradientsteepest(A,b,X0,NIter,tol)
    % minimize norm(Ax-b)
    % gradient descent (Steepest descent)

    k = 1;X = X0;
    Diff = (A * X - b);
    G = 2 * A' * Diff;

    while ((k <= NIter)&& (norm(G) > tol))
        % k
        % Compute step size, steepest descent for quadratic problem
        temp = A * G;
        alpha = 0.5 * (G'*G)/(temp' * temp);
        X = X - alpha * G;
        Diff = (A * X - b);
        G = 2 * A' * Diff;
        Obj(k) = Diff'*Diff;
        k = k + 1;
    end
end