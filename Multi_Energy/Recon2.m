% Recon Experiment 2
% Date: 02/13/2018
clc;
clear;load('ATotal.mat');
% [structGeo,structProj] = setGeometryFwd(structImg);

% Compute the covariance matrix
Nb = length(structEBas);
Ne = structDet.nEnergyBin;
Np = structProj.nAngle * structProj.nChannel;
Nim = structImg.nPixel;
MatrixF = zeros(Ne,Nb);
for m = 1:Ne
    for n = 1:Nb
        nSrcLoc = find(structEBas(n).pfEnergy >= structDet.pfEnergyLow(m)&...
                       structEBas(n).pfEnergy < structDet.pfEnergyHigh(m));
        MatrixF(m,n) = mean(structEBas(n).pfBase(nSrcLoc));
    end
end

temp = cell(structProj.nChannel,structProj.nAngle);
for m = 1:structProj.nAngle 
    for n = 1:structProj.nChannel
        cc = structProj.pfPhotonrnd(n,m,:);cc = cc(:);
        cc = (MatrixF' * diag(cc) * MatrixF) * 10^(-6);
%         cc = chol(MatrixF' * diag(cc) * MatrixF)'*10^(-4);
%         cc = opMatrix(cc^(-1));
%         temp(n,m) = {eye(Nb).* cc};
        temp(n,m) = {(cc)};
    end
end
temp = temp';temp = temp(:);
 Cov = blkdiag(temp{:});


% Reconstruction parameters setup
[structGeo,structProj] = setGeometryFwd(structImg);
tao = 3600;% 60^2; % 60^2;% [10:10:100];
y = permute(ATotal,[3,2,1]); y = y(:);
WNew = opKron(structGeo.W,opEye(Nb));

% % ADMM -- should have faster convergence
% X = zeros(Nim * Nim * Nb,1);
% Z = zeros(Np * Nb,1);% y;% zeros(Np * Nb,1);
% U = zeros(Np * Nb,1);% -Z;% zeros(Np * Nb,1);
% m = 1;rou = 4;
% 
% 
% % Compute (2 * Cov + rou * I)^(-1), for Z update,closed form
% for nn = 1: length(temp)
%     tempinv{nn} = (temp{nn} * 2 + rou * eye(Nb))^(-1);
% end

% % Main iteration
% while (m < 100)
%     m
%     % X update
%     k = 1;
%     Proj = WNew * X;
%     Diff = Proj - Z + U;
%     G = rou * WNew' * Diff + 2 * tao * X;
%     H = opEye(Nim * Nim * Nb);
%     while ((k <= 5) && (norm(G) > 0.0001))
%         % k
%         % 1) Compute search direction
%         Direction = H * G;
%         % 2) Compute step size    
%         alpha = (G' * Direction) / ...
%             (Direction' * (2 * tao * Direction + ...
%             rou * WNew' * WNew * Direction));
%         % 3) Update X
%         XOld = X;
%         X = XOld - alpha * Direction;
%         % 4) Update Hessian-Inverse matrix    
%         GOld = G;
%         Proj = WNew * X;
%         Diff = Proj - Z + U;
%         G = rou * WNew' * Diff + 2 * tao * X;
%     
%         deltaGrad = opMatrix(G - GOld);
%         HdeltaGrad = (H*deltaGrad);
%         deltaX = opMatrix(X - XOld);
%         deltaXDeltaGrad = double(deltaX'*deltaGrad);
%             
%         % BFGS update
%         H = (opEye(Nim * Nim * Nb) - (deltaX * deltaGrad')/(deltaGrad' * deltaX)) ...
%             * H * (opEye(Nim * Nim * Nb) - (deltaGrad * deltaX')/(deltaGrad' * deltaX))+...
%             (deltaX * deltaX')/(deltaGrad' * deltaX);
%     
%         % 5) Display Result
%         % 0.5 * rou * Diff' * Diff + tao * norm(X)^2
%         k = k + 1;
%     end
%     % Z update
%     Proj = WNew * X;
%     Z = CovxData(tempinv,rou * (Proj + U) + 2 * CovxData(temp,y,Nb,Np),Nb,Np);
%     % Z = (rou * (Proj + U)+ 2 * y)/(2+rou); %Not use this
%     % U update
%     U = U + Proj - Z;
%     % Record objective function
%     Diff = Proj - y;    
%     DiffCov = CovxData(temp,Diff,Nb,Np);
%     Obj(m) = Diff' * DiffCov + tao * norm(X)^2;
%     [Diff' * Diff,Obj(m)],
%     norm(U),
%     m = m + 1;
% end
% X1 = reshape(X,[Nb,Nim,Nim]);
% X1 = permute(X1,[2,3,1]);ArtRecon2 = X1(:,:,3);
% figure;imagesc(X1(:,:,3));colorbar; 

% % fminunc
% X0 = zeros(Nim * Nim * Nb,1);
% options = optimoptions(@fminunc,'Algorithm','quasi-newton','MaxIterations',2,...
%     'SpecifyObjectiveGradient',true);
% [XOpt,FVAL,EXITFLAG,OUTPUT] = fminunc(@(X) myfun(X,tao,WNew,Nb,Np,y,temp),X0,options);
% pause;

% BFGS(quasi-newton)
X = zeros(Nim * Nim * Nb,1);
k = 1;
Diff = WNew * X - y;
DiffCov = CovxData(temp,Diff,Nb,Np);
G = 2 * WNew' * DiffCov + 2 * tao * X;
H = opEye(Nim * Nim * Nb);

while ((k < 10000) &&(norm(G) > 0.001))
    k
    % 1) Compute search direction
    Direction = H * G;
    % 2) Compute step size    
    alpha = (G' * Direction) / ...
        (Direction' * (2 * tao * Direction + ...
        2 * WNew' * CovxData(temp,WNew * Direction,Nb,Np)))
    % 3) Update X
    XOld = X;
    X = XOld - alpha * Direction;
    % 4) Update Hessian-Inverse matrix    
    GOld = G;
    Diff = WNew * X - y;
    DiffCov = CovxData(temp,Diff,Nb,Np);
    G = 2 * WNew' * DiffCov + 2 * tao * X;
    
    deltaGrad = opMatrix(G - GOld);
    HdeltaGrad = (H*deltaGrad);
    deltaX = opMatrix(X - XOld);
    deltaXDeltaGrad = double(deltaX'*deltaGrad);
    
    % Should consider a limited-memory version!!
        
    % BFGS update
    % NOTE: the outer product (deltaX*deltaX') is not parenthesized on
    % purpose. This is how the computation was originally done. Adding the
    % parentheses changes the results subtly, creating an incompatibility.
    % Since it doesn't improve the solver, this is undesirable.
%     H = H + (1 + deltaGrad'*HdeltaGrad/deltaXDeltaGrad) * ...
%         deltaX*deltaX'/deltaXDeltaGrad - (deltaX*HdeltaGrad' + ... 
%         HdeltaGrad*deltaX')/deltaXDeltaGrad; %#ok
    H = (opEye(Nim * Nim * Nb) - (deltaX * deltaGrad')/(deltaGrad' * deltaX)) ...
        * H * (opEye(Nim * Nim * Nb) - (deltaGrad * deltaX')/(deltaGrad' * deltaX))+...
        (deltaX * deltaX')/(deltaGrad' * deltaX);
    

    % 5) Record Result
    Obj(k) = Diff' * DiffCov + tao * norm(X)^2;
    Diff' * Diff,% + tao * norm(X)^2
    Obj(k)
    k = k + 1;
end
X1 = reshape(X,[Nb,Nim,Nim]);
X1 = permute(X1,[2,3,1]);ArtRecon2 = X1(:,:,3);
figure;imagesc(X1(:,:,3));colorbar;  
% 
% 
% % Gradient descent method
% ArtRecon = zeros(Nim * Nim * Nb,1);
% k = 1; Diff = WNew * ArtRecon - y;
% DiffCov = CovxData(temp,Diff,Nb,Np);
% G = 2 * WNew' * DiffCov + 2 * tao * ArtRecon;
% Diff'*Diff
% Diff'*DiffCov
% while ((k < 10000) &&(norm(G) > 0.01))
%     k
%     Stepsize =(G' * G) /(G' * (2 * tao * G + 2 * WNew' * CovxData(temp,WNew * G,Nb,Np)))
%     ArtRecon = ArtRecon - Stepsize * G;    
%     Diff = WNew * ArtRecon - y;
%     DiffCov = CovxData(temp,Diff,Nb,Np);
%     G = 2 * WNew' * DiffCov + 2 * tao * ArtRecon;
%     Obj(k) = Diff' * DiffCov + tao * norm(ArtRecon)^2;
%     Diff' * Diff + tao * norm(ArtRecon)^2
%     Obj(k)
%     k = k + 1;
% end
% ArtRecon1 = reshape(ArtRecon,[Nb,Nim,Nim]);
% ArtRecon1 = permute(ArtRecon1,[2,3,1]);ArtRecon2 = ArtRecon1(:,:,3);
% figure;imagesc(ArtRecon1(:,:,3));colorbar;  

function z = CovxData(Cov,y,Nb,Np)
    z = zeros(Nb * Np,1);
    for m = 1:Np
        z((m-1)*Nb + 1:m * Nb) = Cov{m} * y((m-1)*Nb + 1:m * Nb);
        % z((m-1)*Nb + 1:m * Nb) = y((m-1)*Nb + 1:m * Nb);
    end
%    z = y;
end

function [f,g] = myfun(X,tao,WNew,Nb,Np,y,temp)
    Diff = WNew * X - y;
    DiffCov = CovxData(temp,Diff,Nb,Np);
    f = Diff' * DiffCov + tao * X'*X;
    g = 2 * (WNew' * DiffCov + tao * X);
end



