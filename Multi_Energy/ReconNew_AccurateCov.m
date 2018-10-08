% ReconNew
% Date: 2018/02/20
% Compare Gradient, Conjugate gradient, and Quasi-Newton
% Weighted or unweighted objective function

close all;
clear;load('ATotal.mat');
set(0, 'DefaultLineLineWidth', 2);

% Define system size
Nb = length(structEBas);
Ne = structDet.nEnergyBin;
Np = structProj.nAngle * structProj.nChannel;
Nim = structImg.nPixel;

%------------------ Compute the covariance matrix ------------------------%
MatrixF = zeros(Ne,Nb);
for m = 1:Ne
    for n = 1:Nb
        nSrcLoc = find(structEBas(n).pfEnergy >= structDet.pfEnergyLow(m)&...
                       structEBas(n).pfEnergy < structDet.pfEnergyHigh(m));
        MatrixF(m,n) = mean(structEBas(n).pfBase(nSrcLoc));
    end
end

% covariance matrix is sparse and block diagonal
% be careful about sparse matrix generation
tic
X = zeros(Np * Nb * Nb,1);
Y = zeros(Np * Nb * Nb,1);
Data = zeros(Np * Nb * Nb,1);

% %%%%%%%%%% An accurate method for computing covariance matrix %%%%%
% for n = 1:Nb
%     nBaseLoc = (structEBas(n).pfEnergy >=structDet.pfEnergyLow(1)&...
%         structEBas(n).pfEnergy < structDet.pfEnergyHigh(end));
%     F3dimention(:,:,n) = reshape(structEBas(n).pfBase(nBaseLoc),...
%         structDet.fEnerghWidth,structDet.nEnergyBin);
% end
% F3dimention = permute(F3dimention,[2,1,3]);

for n = 1:structProj.nChannel 
    for m = 1:structProj.nAngle
        nPos = (n-1) * structProj.nAngle + m;
        
%         %%%%%%%%%% An accurate method for computing covariance matrix %%%%%
%         % Doesn't seem to improve...
%         pfPhoton = PhotonFromBaseContinuous(ATotal(n,m,:),structEBas,structDet,...
%                                             structSrc,structProj);
%         pfPhotonDev0 = sum(pfPhoton,2);
%                                         
%         for mm = 1:Nb                                
%             pfPhotonDev1(mm,:) = -sum(F3dimention(:,:,mm).* pfPhoton,2);
%         end
%         for mm = 1:Nb  
%             for nn = 1:Nb
%                 pfPhotonDev2(mm,nn,:) = sum(F3dimention(:,:,nn).* F3dimention(:,:,mm).* pfPhoton,2);
%             end
%         end
% 
%         % Compute Nb * Nb matrix (mm,nn)
%         CSmall = zeros(Nb,Nb);
%         for mm = 1:Nb
%             for nn = mm:Nb
%                 temp = 0;
%                 for ii = 1:structDet.nEnergyBin
%                     temp = temp + structProj.pfPhotonrnd(n,m,ii) * ...
%                         (pfPhotonDev1(mm,ii) * pfPhotonDev1(nn,ii) - pfPhotonDev0(ii) * pfPhotonDev2(mm,nn,ii))/...
%                         pfPhotonDev0(ii)^2 - pfPhotonDev2(mm,nn,ii);
%                 end
%                 CSmall(mm,nn) = temp;
%             end
%         end
%         CSmall = CSmall + CSmall' - diag(diag(CSmall));
%         CSmall = -CSmall * 10^(-6);  
%         
%         [xx,yy] = meshgrid((nPos-1) * Nb + 1: nPos * Nb, (nPos-1) * Nb + 1: nPos * Nb);
%         CSmall = chol(CSmall)';
%         X((nPos - 1) * Nb * Nb + 1: nPos * Nb * Nb) = xx(:);
%         Y((nPos - 1) * Nb * Nb + 1: nPos * Nb * Nb) = yy(:);
%         Data((nPos - 1) * Nb * Nb + 1: nPos * Nb * Nb) = CSmall(:);
%         
%         %%%%%%%%%% An accurate method for computing covariance matrix %%%%%
        
        %%%%%%%%%% An original method for covariance matrix %%%%%%%%%%%%%%%
        [xx,yy] = meshgrid((nPos-1) * Nb + 1: nPos * Nb, (nPos-1) * Nb + 1: nPos * Nb);
        cc = structProj.pfPhotonrnd(n,m,:);cc = cc(:);
        cc = (MatrixF' * diag(cc) * MatrixF) * 10^(-6);
        % cc = cc .* eye(Nb); % remove off-diagonal covariance?
        cc = chol(cc)';
        X((nPos - 1) * Nb * Nb + 1: nPos * Nb * Nb) = xx(:);
        Y((nPos - 1) * Nb * Nb + 1: nPos * Nb * Nb) = yy(:);
        Data((nPos - 1) * Nb * Nb + 1: nPos * Nb * Nb) = cc(:);
        
        % CSmall,cc,pause;
        %%%%%%%%%% An original method for covariance matrix %%%%%%%%%%%%%%%

    end
end
toc

Cov = opMatrix(sparse(X,Y,Data));
% full(Cov(1:6,1:6)'*Cov(1:6,1:6)) % validate the covariance matrix

%------------------ Compute the covariance matrix ------------------------%

% Define the problem
[structGeo,structProj] = setGeometryFwd(structImg);
tao = 3600;% 60^2; % 60^2;% [10:10:100];
y = permute(ATotal,[3,2,1]); y = y(:);
WNew = opKron(structGeo.W,opEye(Nb));
NIter = 20;
tol = 10^(-10);
% A good initial guess is essential
X0 = zeros(Nim * Nim * Nb,1);

Warmup = 20;
[X0,Obj1] = gradientsteepest(WNew,y,X0,Warmup,tol);
WNew = Cov * WNew;
y = Cov * y;



% Method 1--Steepest descent
tic
[X1,Obj1] = gradientsteepest(WNew,y,X0,NIter,tol);
toc
X1 = reshape(X1,[Nb,Nim,Nim]);
X1 = permute(X1,[2,3,1]);
figure;imagesc(X1(:,:,3));colorbar; 

% Method 2-- Quasi-Newton
tic 
[X2,Obj2] = quasinewton(WNew,y,X0,NIter,tol);
toc
X2 = reshape(X2,[Nb,Nim,Nim]);
X2 = permute(X2,[2,3,1]);
figure;imagesc(X2(:,:,3));colorbar; 

% Method 3-- pcg(Matlab build-in)
tic
[X3,xtotal] = pcg(WNew' * WNew, WNew' * y,tol,NIter,[],[],X0); 
toc
for k = 1:size(xtotal,2)
    Diff = WNew * xtotal(:,k) - y;
    Obj3(k) = Diff' * Diff;
end
X3 = reshape(X3,[Nb,Nim,Nim]);
X3 = permute(X3,[2,3,1]);
figure;imagesc(X3(:,:,3));colorbar; 

figure;plot(Obj1);hold on;plot(Obj2);plot(Obj3);
legend('GradientDescent','Quasi-Newton','CG');

