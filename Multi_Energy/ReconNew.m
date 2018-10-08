% ReconNew
% Date: 2018/02/20
% Compare Gradient, Conjugate gradient, and Quasi-Newton
% Weighted or unweighted objective function

close all;
clear;load('ATotal.mat');
load('ATotal_Update.mat');
set(0, 'DefaultLineLineWidth', 2);

% Define system size
Nb = length(structEBas);
Ne = structDet.nEnergyBin;
Np = structProj.nAngle * structProj.nChannel;
Nim = structImg.nPixel;

XLoc = [];
YLoc = [];
Data = [];
nn = 0;Nx = Nim; Ny = Nim;
for yy = 1:(Ny-1)
    for xx = 1:(Nx-1)
        nLoc = (yy-1) * Nx + xx;
        XLoc(nn * 4 + 1:nn * 4 + 4) = [nn *2 + 1,nn*2 + 1,nn*2+2,nn*2+2];
        YLoc(nn * 4 + 1:nn * 4 + 4) = [nLoc, nLoc + 1, nLoc, nLoc + Nx];
        Data(nn * 4 + 1:nn * 4 + 4) = [1,-1,1,-1];
        nn = nn + 1;
    end
end
XLine = nn * 2;
nn = nn * 4;

yy = Ny;
for xx = 1:(Nx-1)
    nLoc = (yy-1) * Nx + xx;
    XLoc(nn + 1:nn + 2) = [XLine + 1, XLine + 1];
    YLoc(nn + 1:nn + 2) = [nLoc, nLoc + 1];
    Data(nn + 1:nn + 2) = [1,-1];
    nn = nn + 2;
    XLine = XLine + 1;
end

xx = Nx;
for yy = 1:(Ny-1)
    nLoc = (yy-1) * Nx + xx;
    XLoc(nn + 1:nn + 2) = [XLine + 1, XLine + 1];
    YLoc(nn + 1:nn + 2) = [nLoc, nLoc + Nx];
    Data(nn + 1:nn + 2) = [1,-1];
    nn = nn + 2;
    XLine = XLine + 1;
end

MaskM = (sparse(XLoc,YLoc,Data,XLine, Nx * Ny, nn));
MaskM = kron(MaskM,speye(Nb,Nb));


%------------------ Compute the covariance matrix ------------------------%
PhotonThres = 1000;
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

% CovCond = zeros(structProj.nChannel,structProj.nAngle);
% CovEigMax = zeros(structProj.nChannel,structProj.nAngle);
% CovEigMin = zeros(structProj.nChannel,structProj.nAngle);

CovComp = zeros(structProj.nChannel,structProj.nAngle,6);

for n = 1:structProj.nChannel 
    for m = 1:structProj.nAngle
        nPos = (n-1) * structProj.nAngle + m;        
        %%%%%%%%%% An original method for covariance matrix %%%%%%%%%%%%%%%
        [xx,yy] = meshgrid((nPos-1) * Nb + 1: nPos * Nb, (nPos-1) * Nb + 1: nPos * Nb);
        cc = structProj.pfPhotonrnd(n,m,:);cc = cc(:);
        cc = min(cc,PhotonThres);
        cc = (MatrixF' * diag(cc) * MatrixF) * 10^(-6);
        
%         CovCond(n,m) = cond(cc);
%         CovEigMax(n,m) = max(eig(cc));
%         CovEigMin(n,m) = min(eig(cc));

        CovComp(n,m,1) = cc(1,1);
        CovComp(n,m,2) = cc(2,2);
        CovComp(n,m,3) = cc(3,3);
        CovComp(n,m,4) = cc(1,2);
        CovComp(n,m,5) = cc(1,3);
        CovComp(n,m,6) = cc(2,3);
        
        % cc = cc .* eye(Nb); % remove off-diagonal covariance?
        DataFull((nPos - 1) * Nb * Nb + 1: nPos * Nb * Nb) = cc(:);
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
CovFull = opMatrix(sparse(X,Y,DataFull));
% full(Cov(1:6,1:6)'*Cov(1:6,1:6)) % validate the covariance matrix

%------------------ Compute the covariance matrix ------------------------%

% Define the problem
[structGeo,structProj] = setGeometryFwd(structImg);
tao = 60;% 60^2; % 60^2;% [10:10:100];
y = permute(ATotal,[3,2,1]); y = y(:);
WNew = opKron(structGeo.W,opEye(Nb));
NIter = 20;
tol = 10^(-10);
% A good initial guess is essential
X0 = zeros(Nim * Nim * Nb,1);

% %%
% clear;
% NX = 150;
% WNew = rand(NX,100);  % Only when WNew is an over-determined system,but why.....
% y = WNew * ones(100,1);
% CovFull =  diag([1:NX]) + 0 * eye(NX); % eye(100);%
% % CovFull(1) = 100;
% Cov = sqrt(CovFull);
% X0 = zeros(100,1);

% % % Method 4--ADMM(Is it really faster than quasi-newton?)
% X4 = zeros(Nim * Nim * Nb,1);
% Z =  zeros(Np * Nb,1);% y;% zeros(Np * Nb,1);
% U =  zeros(Np * Nb,1);% -Z;% zeros(Np * Nb,1);
% rou = 20; NIterx = 5;k = 1;
% alpha = 1.0;
% mu = 10;
% tic
% while (k <= 300/NIterx)
%     % k
% 
%     % X update
%     X4 = quasinewton(WNew,Z-U,X4,NIterx,tol);
% %   X4 = alpha* X4 + (1 - alpha)*Z;
%     % Z update
%     Proj = WNew * X4;
%     ZOld = Z;
%     Z = (rou * opEye(Np * Nb) + 2 * CovFull)\(2 * CovFull * y + rou * (Proj + U));
%     % U update
%     U = U + Proj - Z;
%     % S update
%     S = rou * WNew' * (Z - ZOld);
%     PrimalRes(k) = U' * U;
%     DualRes(k) = S' * S;
%     
%     Diff = Cov * (Proj - y);
%     Obj4(k) = Diff' * Diff;
%     PrimalRes,DualRes,Obj4
%     
% %     % Adjusting rou
% %     if (PrimalRes(k) > mu * DualRes(k))
% %         rou = rou * 2
% %     end
% %     if (DualRes(k) > mu * PrimalRes(k))
% %         rou = rou / 2
% %     end
%     
%     k = k + 1;
% end
% %%
% toc
% X4 = reshape(X4,[Nb,Nim,Nim]);
% X4 = permute(X4,[2,3,1]);
% figure;imagesc(X4(:,:,3));colorbar; 

Warmup = 20;
[Xwarm,Obj1] = gradientsteepest(WNew,y,X0,Warmup,tol);
XwarmImg = reshape(Xwarm,[Nb,Nim,Nim]);
XwarmImg = permute(XwarmImg,[2,3,1]);
figure;imagesc(XwarmImg(:,:,3));colorbar; 

WNew = Cov * WNew;
y = Cov * y;


% Quasi-Newton -- Preconditioning & L2-Norm
size([WNew;tao * opEye(Nim * Nim * Nb)])
size([y;zeros(Nim * Nim * Nb,1)])

tao = 60;
NIter = 20;
[X2,Obj,XTotal] = quasinewton([WNew;tao * opEye(Nim * Nim * Nb)],...
    [y;zeros(Nim * Nim * Nb,1)],Xwarm,NIter,tol);


xx = sqrt(max(Xwarm,10^(-8)));
PreCon = opMatrix(spdiags(xx, 0, Nim * Nim * Nb, Nim * Nim * Nb));
XwarmP = Xwarm ./xx;
% [X2,Obj,XTotal] = quasinewton([WNew;tao * opEye(Nim * Nim * Nb)] *PreCon,...
%     [y;zeros(Nim * Nim * Nb,1)],XwarmP,NIter,tol);

[X2,Obj,XTotal] = quasinewton([WNew;tao * MaskM] *PreCon,...
    [y;zeros(size(MaskM,1),1)],XwarmP,NIter,tol);
X2 = PreCon * X2;
XTotal = PreCon * XTotal;

XX = reshape(X2,[Nb,Nim,Nim]);
XX = permute(XX,[2,3,1]);
figure;imagesc(XX(:,:,3));colorbar; 

pause;




% % Method 1--Steepest descent
% tic
% [X1,Obj1] = gradientsteepest(WNew,y,Xwarm,NIter,tol);
% toc
% X1 = reshape(X1,[Nb,Nim,Nim]);
% X1 = permute(X1,[2,3,1]);
% figure;imagesc(X1(:,:,3));colorbar; 




% % Method 2-- Quasi-Newton
% [X2,Obj,XTotal] = quasinewton(WNew,y,zeros(Nim * Nim * Nb,1),NIter,tol);
[X2,Obj,XTotal] = quasinewton(WNew,y,Xwarm,NIter,tol);
% [X2,Obj2] = quasinewton([WNew;tao * MaskM],...
%     [y;zeros(size(MaskM,1)* Nb,1)],zeros(Nim * Nim * Nb,1),NIter,tol);
X2 = reshape(X2,[Nb,Nim,Nim]);
X2 = permute(X2,[2,3,1]);
figure;imagesc(X2(:,:,3));colorbar; 

XTotal = reshape(XTotal,[Nb,Nim,Nim,NIter]);
XTotal = permute(XTotal,[2,3,4,1]);
for n = 1:NIter  
    gca = imagesc(XTotal(:,:,n,1));colorbar;
    saveas(gca,['X1_',num2str(n),'.png']);
    close all;
    
    gca = imagesc(XTotal(:,:,n,2));colorbar;
    saveas(gca,['X2_',num2str(n),'.png']);
    close all;
    
    gca = imagesc(XTotal(:,:,n,3));colorbar;
    saveas(gca,['X3_',num2str(n),'.png']);
    close all;
end
pause;



% % Method 5-- Preconditioned Quasi-Newton???
% tic 
% xx = sqrt(max(Xwarm,10^(-8)));
% PreCon = opMatrix(spdiags(xx, 0, Nim * Nim * Nb, Nim * Nim * Nb));
% % xx = 1./(WNew' * WNew * ones(Nim * Nim * Nb,1)) * 10^(5);
% % PreCon = opMatrix(spdiags(xx, 0, Nim * Nim * Nb, Nim * Nim * Nb));
% WNewP = WNew * PreCon;
% [X5,Obj5] = quasinewton(WNewP,y,zeros(Nim * Nim * Nb,1),NIter,tol);
% X5 = PreCon * X5;
% Diff = (WNew * X5 - y);
% Diff' * Diff
% toc
% X5 = reshape(X5,[Nb,Nim,Nim]);
% X5 = permute(X5,[2,3,1]);
% figure;imagesc(X5(:,:,3));colorbar; 


%% % Algorithm One(Page 4)--- Regu term |x|_1---Primal Dual
lamda = 1/36;
gamma = 0.01;
tao = 9 * 10 ^(-5) ;
sigma = 9 * 10 ^(-5);
theta = 1;
X = Xwarm;% zeros(Nb * Nim * Nim,1);
XBar = X;
Y = zeros(Nb * Np,1);
k = 1;
while (k < NIter)
    k
    Temp = (Y + sigma * WNew * XBar);
    Y = (Temp - sigma * y)/(1 + sigma / (2 * lamda));

    Temp = X - tao * WNew' * Y;
    % XNew = max(1 - tao ./ abs(Temp),0).* Temp;  % Regu term |x|_1
    XNew = (speye(Nb * Nim * Nim) + 2 * tao * MaskM' * MaskM)\Temp; % Regu term |Dx|_2
    
%     theta = 1/sqrt(1 + 2 * gamma * tao);
%     tao = tao / theta; sigma = sigma * theta;
    
    XBar = XNew + theta * (XNew - X);
    X = XNew;
    Obj(k) = lamda * norm(WNew * X - y)^2 + norm(MaskM * X)^2
    
%     if (k> 2)&&(Obj(k) > Obj(k-1))
%         tao = tao / 1.1;
%         sigma = sigma / 1.1;
%     end
    k = k + 1;
end
lamda * norm(WNew * X - y)^2 + norm(MaskM * X)^2
lamda * norm(WNew * Xwarm - y)^2 + norm(MaskM * Xwarm)^2

figure;semilogy(Obj);
result = reshape(X,[Nb,Nim,Nim]);
result = permute(result,[2,3,1]);
figure;imagesc(result(:,:,3));colorbar; 

% 
% % Method 3-- pcg(Matlab build-in)
% tic
% [X3,xtotal] = pcg(WNew' * WNew, WNew' * y,tol,NIter,[],[],Xwarm); 
% toc
% for k = 1:size(xtotal,2)
%     Diff = WNew * xtotal(:,k) - y;
%     Obj3(k) = Diff' * Diff;
% end
% X3 = reshape(X3,[Nb,Nim,Nim]);
% X3 = permute(X3,[2,3,1]);
% figure;imagesc(X3(:,:,3));colorbar; 
% 
% figure;plot(Obj1);hold on;plot(Obj2);plot(Obj3);
% plot([1:length(Obj4)] * NIterx, Obj4);
% legend('GradientDescent','Quasi-Newton','CG','ADMM');
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Generate a true composition signal for comparison
Contrast = zeros(structImg.nPixel,structImg.nPixel);
% Region 5
fTemp = (structImg.pfImgLabel(:,:) == 5);
fMAC = 0.03;
fTemp = fTemp * fMAC * structImg.pfDensity(5);
Contrast = Contrast + fTemp;
% Region 6
fTemp = (structImg.pfImgLabel(:,:) == 6);
fMAC = 0.01;
fTemp = fTemp * fMAC * structImg.pfDensity(6);
Contrast = Contrast + fTemp;
% Region 7
fTemp = (structImg.pfImgLabel(:,:) == 7);
fMAC = 0.01;
fTemp = fTemp * fMAC * structImg.pfDensity(7);
Contrast = Contrast + fTemp;
Contrast = Contrast * 0.1;
[true_id, true] = astra_create_sino(...
            Contrast, structGeo.proj_id);
plot(ATotal(:,1,3));hold on;plot(true(1,:)*0.1);

%%%% ASTRA-FBP Comparison
FBP1 = ReconDecomp(structGeo,ATotal(:,:,1)');
FBP2 = ReconDecomp(structGeo,ATotal(:,:,2)');
FBP3 = ReconDecomp(structGeo,ATotal(:,:,3)');

figure;subplot(2,3,1);imagesc(FBP1);colorbar;
subplot(2,3,2);imagesc(FBP2);colorbar;
subplot(2,3,3);imagesc(FBP3);colorbar;

subplot(2,3,4);imagesc(X2(:,:,1));colorbar;
subplot(2,3,5);imagesc(X2(:,:,2));colorbar;
subplot(2,3,6);imagesc(X2(:,:,3));colorbar;


%%%%%% Error Component3
for n = 1:NIter
    Error3(n) = norm(XTotal(:,:,n,3)-Contrast);
end
figure;plot(Error3,'*-');
title('Error Norm for recon component3--No \Sigma matrix');grid on;

%%%%%% Error energy images
% generate energy image at 60keV
Im_Energy = structEBas(1).pfBase(51) * X2(:,:,1) + ...
            structEBas(2).pfBase(51) * X2(:,:,2) + ...
            structEBas(3).pfBase(51) * X2(:,:,3);
        
figure;subplot(2,2,1);imagesc(Im_Energy);colorbar;title('Recon-60Kev');
subplot(2,2,2);imagesc(structImg.pfImgLac(:,:,46));colorbar;title('True-60Kev');


% generate energy image at 40keV
Im_Energy = structEBas(1).pfBase(31) * X2(:,:,1) + ...
            structEBas(2).pfBase(31) * X2(:,:,2) + ...
            structEBas(3).pfBase(31) * X2(:,:,3);
        
subplot(2,2,3);imagesc(Im_Energy);colorbar;title('Recon-40Kev');
subplot(2,2,4);imagesc(structImg.pfImgLac(:,:,26));colorbar;title('True-40Kev');


% generate energy image at 60keV
Im_Energy = structEBas(1).pfBase(51) * XTotal(:,:,:,1) + ...
            structEBas(2).pfBase(51) * XTotal(:,:,:,2) + ...
            structEBas(3).pfBase(51) * XTotal(:,:,:,3);
        
% Error for energy image
for n = 1:NIter
    ErrorEnergy(n) = norm(Im_Energy(:,:,n)-structImg.pfImgLac(:,:,46));
end
figure;plot(ErrorEnergy,'*-');
title('Error Norm for energy image at 60kev--No \Sigma matrix');grid on;

