% A general limited-memory Quasi-Newton function...
% X0: Initial guess;
% NIter: Number of Iterations
% Func_Obj: objective function handle
% Func_Grad: gradient function handle
% NMemory: number of vectors to estimate hessian inverse (Limited memory)

% Reference: Representations of quasi-Newton matrices and their use in
% limited memory methods (Byrd, Nocedal, 1992)


function [X,Obj,XTotal] = Quasi_Newton_LBFGS...
    (X0, NIter, NMemory, Func_Obj, Func_Grad, Data,isQuadratic)
    k = 1; X = X0;
    N = length(X0);
    XTotal = zeros(N,NIter);
    G = Func_Grad(X,Data);
    Obj = zeros(1,NIter);

    % vectors for estimating hessian inverse
    Sk = [];
    Yk = [];
    SandY = [];
    Mk = [];
    Rk_inv = [];
    Dk = [];

    t_old = 1;
    while (k <= NIter)
        k
        % 1) Compute search direction
        % Direction = -H * G;
        if isempty(Sk)
            Direction = -G;
        else
            % Direction = - (G + SandY * (Mk * (SandY' * G)));
            
            Sg = Sk' * G;
            Yg = Yk' * G;
            Rk_inv_Sg = Rk_inv * Sg;
            Direction = - (G + Sk * (Rk_inv' * (Qk * Rk_inv_Sg)) ...
                - Sk * (Rk_inv' * Yg) - Yk * Rk_inv_Sg);
        end
        
        % 2) Compute step size 
        if isempty(isQuadratic)
            % Backtracking line search
            t = t_old * 2; alpha = 0.4; beta = 0.6;
            % t = 10^(-7);% 10^(-9);
            F_One = Func_Obj(X + t * Direction, Data);
            F_Two = Func_Obj(X, Data);
            F_Two_Inc = G' * Direction;
            while (F_One > F_Two + F_Two_Inc * alpha * t)
                t = beta * t;
                F_One = Func_Obj(X + t * Direction, Data);
            end
            alpha = t;
            t_old = t / beta;
            % t
        else
            temp = Data.A * Direction;
            alpha = (-G' * Direction) / (2 * temp' * temp);
            F_One = Func_Obj(X + alpha * Direction, Data);
        end

        
        % 3) Update X
        XOld = X;
        X = X + alpha * Direction;
        
%         if (G' * Direction > -10^(-7))
%             G' * Direction
%             pause; % need restart here~~~
%         end

        % 4) Update Hessian-Inverse matrix
        GOld = G;
        G = Func_Grad(X,Data);

        deltaX = X - XOld;  
        deltaGrad = G - GOld; 
        
        
%         Rouk = deltaX' * deltaGrad;
%         Rouk_inv = 1 / Rouk;
%         
%         if isempty(Rk_inv)
%             Rk_inv = Rouk_inv;
%         else
%             Rk_inv = [Rk_inv, -Rouk_inv * Rk_inv * (Sk' * deltaGrad)];
%             Rk_inv(k, k) = Rouk_inv;
%         end
        
        Sk = [Sk,deltaX];
        Yk = [Yk,deltaGrad];
        
        if (k > NMemory)
            Sk = Sk(:,2:end);
            Yk = Yk(:,2:end);
        end
        
        KSize = size(Sk,2);
        Temp = Sk' * Yk;
        Rk = triu(Temp);
        Rk_inv = inv(Rk);    
        Dk = eye(KSize).* Temp;
        
        
%         % SandY = [Sk, Yk];  % slow
%         SandY(:,1:KSize) = Sk;
%         SandY(:,KSize + 1:2 * KSize) = Yk;
%         Dk(KSize,KSize) = deltaX' * deltaGrad;
        
        Qk = Dk + Yk' * Yk;
        

%         Mk = [Rk_invT * (Dk + Yk' * Yk) * Rk_inv, -Rk_invT;  
%               -Rk_inv, zeros(KSize)]; % slow     
%         Mk = zeros(2 * KSize);  
%         Mk(1:KSize,1:KSize) = Rk_inv' * Qk * Rk_inv; % slow
%         Mk(1:KSize,KSize + 1:2 * KSize) = -Rk_inv';
%         Mk(KSize + 1:2 * KSize,1:KSize) = -Rk_inv;

  
%         H = (opEye(N) - (deltaX * deltaGrad')/(deltaGrad' * deltaX)) ...
%             * H * (opEye(N) - (deltaGrad * deltaX')/(deltaGrad' * deltaX))+...
%             (deltaX * deltaX')/(deltaGrad' * deltaX);
        
        Obj(k) = Func_Obj(X, Data);% F_One; % Func_Obj(X, Data);
        XTotal(:,k) = X;
        k = k + 1;    
    end
end

