
% close all;
clear;load('ATotal_NewG.mat');
load('ATotal_NewG_Update.mat');
set(0, 'DefaultLineLineWidth', 2);
load Proj.mat;

% Define system size
Nb = length(structEBas);
Ne = structDet.nEnergyBin;
Np = structProj.nAngle * structProj.nChannel;
Nim = structImg.nPixel;

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
DataDiag = zeros(Np * Nb * Nb,1);

CovCond = zeros(structProj.nChannel,structProj.nAngle);
CovEigMax = zeros(structProj.nChannel,structProj.nAngle);
CovEigMin = zeros(structProj.nChannel,structProj.nAngle);

CovComp = zeros(structProj.nChannel,structProj.nAngle,6);

for n = 1:structProj.nChannel 
    for m = 1:structProj.nAngle
        nPos = (n-1) * structProj.nAngle + m;        
        %%%%%%%%%% An original method for covariance matrix %%%%%%%%%%%%%%%
        [xx,yy] = meshgrid((nPos-1) * Nb + 1: nPos * Nb, (nPos-1) * Nb + 1: nPos * Nb);
        cc = structProj.pfPhotonrnd(n,m,:);cc = cc(:);
        cc = min(cc,PhotonThres);
        cc = (MatrixF' * diag(cc) * MatrixF) * 10^(-6);
        
        % cc = cc .* eye(Nb); % remove off-diagonal covariance?
        CovCond(n,m) = cond(cc);
        CovEigMax(n,m) = max(eig(cc));
        CovEigMin(n,m) = min(eig(cc));

        CovComp(n,m,1) = cc(1,1);
        CovComp(n,m,2) = cc(2,2);
        CovComp(n,m,3) = cc(3,3);
        CovComp(n,m,4) = cc(1,2);
        CovComp(n,m,5) = cc(1,3);
        CovComp(n,m,6) = cc(2,3);
        
       
        DataFull((nPos - 1) * Nb * Nb + 1: nPos * Nb * Nb) = cc(:);
        cc = chol(cc)';
        X((nPos - 1) * Nb * Nb + 1: nPos * Nb * Nb) = xx(:);
        Y((nPos - 1) * Nb * Nb + 1: nPos * Nb * Nb) = yy(:);
        Data((nPos - 1) * Nb * Nb + 1: nPos * Nb * Nb) = cc(:);
        ccdiag = diag(diag(cc));
        DataDiag((nPos - 1) * Nb * Nb + 1: nPos * Nb * Nb) = ccdiag(:);
        
        % CSmall,cc,pause;
        %%%%%%%%%% An original method for covariance matrix %%%%%%%%%%%%%%%
        
    end
end
toc

Cov = opMatrix(sparse(X,Y,Data));
CovDiag = opMatrix(sparse(X,Y,DataDiag));



% full(Cov(1:6,1:6)'*Cov(1:6,1:6)) % validate the covariance matrix

%------------------ Compute the covariance matrix ------------------------%

% % %% Generate a true composition signal for comparison
% % Contrast = zeros(structImg.nPixel,structImg.nPixel);
% % % Region 5
% % fTemp = (structImg.pfImgLabel(:,:) == 5);
% % fMAC = 0.03;
% % fTemp = fTemp * fMAC * structImg.pfDensity(5);
% % Contrast = Contrast + fTemp;
% % % Region 6
% % fTemp = (structImg.pfImgLabel(:,:) == 6);
% % fMAC = 0.01;
% % fTemp = fTemp * fMAC * structImg.pfDensity(6);
% % Contrast = Contrast + fTemp;
% % % Region 7
% % fTemp = (structImg.pfImgLabel(:,:) == 7);
% % fMAC = 0.01;
% % fTemp = fTemp * fMAC * structImg.pfDensity(7);
% % Contrast = Contrast + fTemp;
% % Contrast = Contrast * 0.1;

%% ASTRA-FBP Comparison
Image_FBP(:,:,1) = ReconDecomp(structGeo,ATotal(:,:,1)');
Image_FBP(:,:,2) = ReconDecomp(structGeo,ATotal(:,:,2)');
Image_FBP(:,:,3) = ReconDecomp(structGeo,ATotal(:,:,3)');

Image_FBP = Image_FBP * 0.1860; % Experience....

figure;imagesc(Image_FBP(:,:,1));colorbar;
figure;imagesc(Image_FBP(:,:,2));colorbar;
figure;imagesc(Image_FBP(:,:,3));colorbar;


% Define the problem
[structGeo,structProj] = setGeometryFwd(structImg);
clear Data;
pause;


%% 0) One-step algorithm?!
NeStart = 4;  %% ultra-low energy data is thrown away...
NeP = Ne - NeStart + 1;
load Proj.mat
y = structProj.pfLACrnd(:,:,NeStart:end);
y = permute(y,[2,1,3]);y = y(:);

WNew = opKron(MatrixF(NeStart:end,:),structGeo.W);
Xwarm = Image_FBP(:);

% Xwarm = X(:); % Continue...

Cov = opDiag(ones(Np * NeP,1)); 
Cov = structProj.pfPhotonrnd(:,:,NeStart:end);
Cov = permute(Cov,[2,1,3]);
Cov = opDiag(sqrt(Cov(:)))/1000;

% Cov = opDiag(ones(Np * (Ne-NeStart+1),1));

% tao = [0.08,0.6,0.2] * 0.5; % AlphaPara1.*AlphaPara2 * 1;
tao = [0.08 0.02 0.1];


% tao = [0.00 0.005 0.0] * 0.5;
% tao = [0.00 0.000 0.05] * 0.5;
% tao = [0.005 0.000 0.0] * 0.5;
tao = [0.006 0.005 0.05] * 0.5;

% 2.418616642147577e+04? 4.319350261425175e+03? 3.108056750899009e+04



% Preconditioning with a unweighted guess improves the convergence rate!
Data.A = Cov * WNew; 
Data.b = Cov * y; 
Data.NImage = structImg.nPixel;
Data.Nb = Nb;
Data.gamma = (1*10^(-5))^2;
Data.tao = tao;

%  Preconditioning
PreCon = ones(Nim * Nim * Nb,1); % 
PreConM = opDiag(PreCon); 
Data.PreCon = PreCon;

NIter = ones(1,50)*20;
NMemory = 20;
Func_Obj = @TV_denoising_Obj_OS;
Func_Grad = @TV_denoising_Grad_OS;

tic
X_LBFGS = Xwarm;
Obj_BFGS = [];
XTotal = [];

for mm = 1:50
    
    Data.PreCon = ones(Nb * Nim * Nim,1); 
    PreCon = sqrt(10^6./TV_denoising_Grad2_E(X_LBFGS(:),Data));

    PreConM = opDiag(PreCon); 
    Data.PreCon = PreCon;
% % 
    [X_LBFGS,Obj,XTotal_LBFGS] = Quasi_Newton_LBFGS...
        (X_LBFGS(:)./PreCon, NIter(mm),NMemory, Func_Obj, Func_Grad, Data,[]); % Image_FBP,X_LBFGS
    
    
    for m = 1:NIter(mm)
        XTotal(:,m + (mm - 1) * NMemory) = PreConM * XTotal_LBFGS(:,m);
    end
    
    X_LBFGS = PreConM * X_LBFGS;

    Obj_BFGS = [Obj_BFGS,Obj];
    Obj_BFGS(end)

    X = reshape(X_LBFGS,Nim,Nim,Nb);

    figure(8);imagesc(X(:,:,1),[0 0.002]);colorbar;
    figure(9);imagesc(X(:,:,2),[0.015 0.020]);colorbar;
    figure(10);imagesc(X(:,:,3),[0 0.002]);colorbar;
        
end

figure;loglog(Obj_BFGS);hold on;
pause;

%% 1) Tikhonov smoothing regularization (No Weight)
for k = 1:Nb
    Xwarm = Image_FBP(:,:,k); Xwarm = Xwarm(:);
    y = ATotal(:,:,k); y = y';y = y(:);
    NIter = 50;
    NMemory = 50;
    tao = [0,0,0];
    
    % Statistical weight matrix
    ErrorW = opDiag(ones(structProj.nChannel * structProj.nAngle,1)); % Cov;
    % Renewed Preconditioning matrix (After 20180903)
    CovIm_Row_Cov = structGeo.W' * ErrorW' * ErrorW * structGeo.W * ones(Nim * Nim,1);
    RTemp = FirstDerivativeM(Nim);
    CovIm_Row_Cov = CovIm_Row_Cov + tao(k)^2 * (RTemp' * RTemp * ones(Nim * Nim,1));
    CovIm_Row_Cov = 1./CovIm_Row_Cov;
    CovIm_Row_Cov = CovIm_Row_Cov * 10^10;
    
    % Preconditioning matrix
    PreCon =  CovIm_Row_Cov; % ones(Nim * Nim,1);% CovIm_Row_Cov;
    PreConM = opDiag(PreCon); 
    
    A  = ErrorW * structGeo.W * PreConM;
    R = FirstDerivativeM(Nim) * PreConM;
    y0 =  ErrorW * y;
    yR = zeros(size(R,1),1);
    
    
    XwarmP = Xwarm ./ PreCon;

    Func_Diff = @(X,Data)(Data.A * X - Data.y);
    Func_Obj = @(X,Data) norm(Func_Diff(X,Data))^2;
    Func_Grad =  @(X,Data) (2 * Data.A' * Func_Diff(X,Data));
    isQuadratic = 1;

    Data.A = [A;tao(k)* R];
    Data.y = [y0;yR];
    
    
    [X0,Obj_BFGS,XTotal] = Quasi_Newton_LBFGS...
        (XwarmP, NIter, NMemory,Func_Obj, Func_Grad, Data,isQuadratic);
    toc

    X = PreConM * X0;
    X_Image(:,:,k) = reshape(X,structImg.nPixel,structImg.nPixel);
    for m = 1:NIter
        XTotal(:,m) = PreConM * XTotal(:,m);
    end

    figure;loglog(Obj_BFGS);hold on;
    figure;imagesc(X_Image(:,:,k));
end

%% 2) Tikhonov smoothing regularization (Diagonal weight)
for k = 3:3 %  1:Nb
    Xwarm = Image_FBP(:,:,k); Xwarm = Xwarm(:);
    y = ATotal(:,:,k); y = y';y = y(:);
    NIter = 100;
    NMemory = 100;
    tao = [0,0,3];
    % Statistical weight matrix
    ErrorW = CovDiag * ones(Np * Nb,1); 
    ErrorW = reshape(ErrorW,Nb,structProj.nAngle,structProj.nChannel);
    ErrorW = squeeze(ErrorW(k,:,:));
    ErrorW = opDiag(ErrorW(:));
    
    % ErrorW = opDiag(ones(Np,1)); %%% Identity covariance
    
    % Renewed Preconditioning matrix (After 20180903)
    CovIm_Row_Cov = structGeo.W' * ErrorW' * ErrorW * structGeo.W * ones(Nim * Nim,1);
    RTemp = FirstDerivativeM(Nim);
    CovIm_Row_Cov = CovIm_Row_Cov + tao(k)^2 * (RTemp' * RTemp * ones(Nim * Nim,1));
    CovIm_Row_Cov = 1./CovIm_Row_Cov;
    CovIm_Row_Cov = CovIm_Row_Cov * 10^10;
    CovIm_Row_Cov = sqrt(CovIm_Row_Cov);
    
    % Preconditioning matrix
    PreCon = CovIm_Row_Cov; % ones(Nim * Nim,1);% CovIm_Row_Cov;
    PreConM = opDiag(PreCon); 
    
    A  = ErrorW * structGeo.W * PreConM;
    R = FirstDerivativeM(Nim) * PreConM;
    y0 =  ErrorW * y;
    yR = zeros(size(R,1),1);
    
    
    XwarmP = Xwarm ./ PreCon;

    Func_Diff = @(X,Data)(Data.A * X - Data.y);
    Func_Obj = @(X,Data) norm(Func_Diff(X,Data))^2;
    Func_Grad =  @(X,Data) (2 * Data.A' * Func_Diff(X,Data));
    isQuadratic = 1;

    Data.A = [A;tao(k)* R];
    Data.y = [y0;yR];
    
    
    [X0,Obj_BFGS,XTotal] = Quasi_Newton_LBFGS...
        (XwarmP, NIter, NMemory,Func_Obj, Func_Grad, Data,isQuadratic);
    toc

    X = PreConM * X0;
    X_Image(:,:,k) = reshape(X,structImg.nPixel,structImg.nPixel);
    for m = 1:NIter
        XTotal(:,m) = PreConM * XTotal(:,m);
    end

    figure;loglog(Obj_BFGS);hold on;
    figure;imagesc(X_Image(:,:,k));
    
    save(['TKCovDiag_Base',num2str(k),'_tao_',num2str(tao(k)),'.mat'],'X','Obj_BFGS','XTotal');
end

%% Noise-Resolution trade-off
X = squeeze(X_Image(:,:,3));
Area = X(60:70,60:70);
std(Area(:))
% Resolution
Loc = [97,97;35,45]; % Strong
% Loc = [45,45;27,36]; % Week
Line = X(Loc(1,1):Loc(1,2),Loc(2,1):Loc(2,2));
Line_True = Contrast(Loc(1,1):Loc(1,2),Loc(2,1):Loc(2,2));
Padding = 30;
Line_pad = padarray(Line_True,[0,Padding],'replicate','both');

gauss_sigma = [0.1:0.01:5.0];% 5.0];
for kk = 1:length(gauss_sigma)
    gauss_func  = normpdf([-15:1:15],0,gauss_sigma(kk));%  * sqrt(2*pi) * gauss_sigma;
    gauss_func = gauss_func/sum(gauss_func);
    % figure;plot(gauss_func);

    Line_conv = conv(Line_pad,gauss_func,'same');
    Line_conv = Line_conv(1+Padding:end-Padding);
    % figure;plot(Line);hold on; plot(Line_True);plot(Line_conv);
    Res(kk) = norm(Line_conv - Line);
end
figure;plot(gauss_sigma,Res,'-');

[Min_res,pos] = min(Res);
gauss_sigma(pos)
gauss_func  = normpdf([-15:1:15],0,gauss_sigma(pos));%  * sqrt(2*pi) * gauss_sigma;
gauss_func = gauss_func/sum(gauss_func);
Line_conv = conv(Line_pad,gauss_func,'same');
Line_conv = Line_conv(1+Padding:end-Padding);
figure;plot(Line);hold on; plot(Line_True);plot(Line_conv);

std(X(:) - Contrast(:))


%% 3) Tikhonov smoothing regularization (Full weight)

y = permute(ATotal,[3,2,1]); y = y(:);
WNew = opKron(structGeo.W,opEye(Nb));
Xwarm = permute(Image_FBP,[3,1,2]); Xwarm = Xwarm(:);
% Xwarm = XTotal(:,end); % Continue iteration...


ErrorW = Cov * ones(Np * Nb,1);
ErrorW = reshape(ErrorW,Nb,structProj.nAngle,structProj.nChannel);
AlphaPara1 = mean(mean(ErrorW.^2,2),3);
AlphaPara2 = squeeze(mean(mean(ATotal,2),1));
    
NIter = 50;
NMemory = 50;
tao = [1,1,2] * 2; % AlphaPara1.*AlphaPara2 * 1;


% Renewed Preconditioning matrix (After 20180903)
CovIm_Row_Cov = WNew' * Cov' * Cov * WNew * ones(Nim * Nim * Nb,1);
RTemp = FirstDerivativeM(Nim);
RTempNew = opKron (RTemp, opEye(Nb) * opDiag(tao));
    
CovIm_Row_Cov = CovIm_Row_Cov + (RTempNew' * RTempNew * ones(Nim * Nim * Nb,1));
CovIm_Row_Cov = 1./CovIm_Row_Cov;
CovIm_Row_Cov = CovIm_Row_Cov * 10^10;
CovIm_Row_Cov = sqrt(CovIm_Row_Cov);

% Statistical weight matrix
ErrorW = Cov;
%Cov; opDiag(ones(Np * Nb,1));
% Preconditioning matrix
PreCon = ones(Nim * Nim * Nb,1); % ones(Nim * Nim * Nb,1);% CovIm_Row_Cov;
PreConM = opDiag(PreCon); 
A  = ErrorW * WNew * PreConM;
R = RTempNew * PreConM;

y0 =  ErrorW * y;
yR = zeros(size(R,1),1);
    
    
XwarmP = Xwarm ./ PreCon;

Func_Diff = @(X,Data)(Data.A * X - Data.y);
Func_Obj = @(X,Data) norm(Func_Diff(X,Data))^2;
Func_Grad =  @(X,Data) (2 * Data.A' * Func_Diff(X,Data));
isQuadratic = 1;

Data.A = [A;R];
Data.y = [y0;yR];
    
    
[X0,Obj_BFGS,XTotal] = Quasi_Newton_LBFGS...
    (XwarmP, NIter, NMemory,Func_Obj, Func_Grad, Data,isQuadratic);

X = PreConM * X0; 
X_Image = reshape(X,Nb,Nim,Nim);
for m = 1:NIter
    XTotal(:,m) = PreConM * XTotal(:,m);
end

figure;loglog(Obj_BFGS);hold on;
X_Image = permute (X_Image,[2,3,1]);
figure;imagesc(X_Image(:,:,1));
figure;imagesc(X_Image(:,:,2));
figure;imagesc(X_Image(:,:,3),[0 0.004]);
figure;imagesc(Image_FBP(:,:,3),[0 0.004]);

%%
X =squeeze(X_Image(:,:,3));
% Noise-Resolution trade-off
Area = X(60:70,60:70);
std(Area(:))
% Resolution
 Loc = [97,97;35,45];
% Loc = [45,45;27,36];
Line = X(Loc(1,1):Loc(1,2),Loc(2,1):Loc(2,2));
Line_True = Contrast(Loc(1,1):Loc(1,2),Loc(2,1):Loc(2,2));
Padding = 30;
Line_pad = padarray(Line_True,[0,Padding],'replicate','both');

gauss_sigma = [0.1:0.01:5.0];% 5.0];
for kk = 1:length(gauss_sigma)
    gauss_func  = normpdf([-15:1:15],0,gauss_sigma(kk));%  * sqrt(2*pi) * gauss_sigma;
    gauss_func = gauss_func/sum(gauss_func);
    % figure;plot(gauss_func);

    Line_conv = conv(Line_pad,gauss_func,'same');
    Line_conv = Line_conv(1+Padding:end-Padding);
    % figure;plot(Line);hold on; plot(Line_True);plot(Line_conv);
    Res(kk) = norm(Line_conv - Line);
end
figure;plot(gauss_sigma,Res,'-');

[Min_res,pos] = min(Res);
gauss_sigma(pos)
gauss_func  = normpdf([-15:1:15],0,gauss_sigma(pos));%  * sqrt(2*pi) * gauss_sigma;
gauss_func = gauss_func/sum(gauss_func);
Line_conv = conv(Line_pad,gauss_func,'same');
Line_conv = Line_conv(1+Padding:end-Padding);
figure;plot(Line);hold on; plot(Line_True);plot(Line_conv);

std(X(:) - Contrast(:))
k = 3;
save(['TKCovFull_Base',num2str(k),'_tao_',num2str(tao(k)),'.mat'],'X','Obj_BFGS','XTotal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4) Smoothed TV (Diagonal weight)

NIter = 200;
NMemory = 200;
tao = [0,0.005,0.0];
% 0,0,0.02
% 0,0.02,0--> 0.01~0.02

for k = 2:2 % 1:Nb
    Xwarm = Image_FBP(:,:,k); Xwarm = Xwarm(:);
    Xwarm = X; Xwarm = Xwarm(:);
    y = ATotal(:,:,k); y = y';y = y(:);
    
    % Statistical weight matrix
    ErrorW = CovDiag * ones(Np * Nb,1); %%%% Debug: Oh this is wrong my god!
    ErrorW = reshape(ErrorW,Nb,structProj.nAngle,structProj.nChannel);
    ErrorW = squeeze(ErrorW(k,:,:));
    ErrorW = opDiag(ErrorW(:));
    
    % ErrorW = opDiag(ones(Np,1)); % Identity covariance
    
    % Preconditioning with a unweighted guess improves the convergence rate!
    Data.A = ErrorW * structGeo.W; 
    Data.b = ErrorW * y; 
    Data.NImage = structImg.nPixel;
    Data.Nb = 1;
    Data.gamma = (1*10^(-5))^2;
    Data.tao = tao(k);

    %  Preconditioning
    PreCon = ones(Nim * Nim,1); % 
    PreConM = opDiag(PreCon); 
    Data.PreCon = PreCon;

    NIter = ones(1,20)*20;
    NMemory = 20;
    Func_Obj = @TV_denoising_Obj_E;
    Func_Grad = @TV_denoising_Grad_E;

    tic
    X_LBFGS = Xwarm;
    Obj_BFGS = [];
    XTotal = [];figure(10);
    for mm = 1:10
        Data.PreCon = ones(Nim * Nim,1); 
        PreCon = sqrt(10^6./TV_denoising_Grad2_E(X_LBFGS(:),Data));

        PreConM = opDiag(PreCon); 
        Data.PreCon = PreCon;

        [X_LBFGS,Obj,XTotal_LBFGS] = Quasi_Newton_LBFGS...
            (X_LBFGS(:)./PreCon, NIter(mm),NMemory, Func_Obj, Func_Grad, Data,[]); % Image_FBP,X_LBFGS
    
        for m = 1:NIter(mm)
            XTotal(:,m + (mm - 1) * NMemory) = PreConM * XTotal_LBFGS(:,m);
        end
    
        X_LBFGS = PreConM * X_LBFGS;

        Obj_BFGS = [Obj_BFGS,Obj];
        Obj_BFGS(end)

        X = reshape(X_LBFGS,Nim,Nim);
        figure(10);imagesc(X);colorbar;
    end
    XD(:,:,k) = X;
    figure;loglog(Obj_BFGS);hold on;
    % save(['TVCovDiag_Base',num2str(k),'_tao_',num2str(tao(k)),'.mat'],'X','Obj_BFGS','XTotal');
end

%% Noise-Resolution trade-off
Area = X(60:70,60:70);
std(Area(:))
% Resolution
 Loc = [97,97;35,45];
% Loc = [45,45;27,36];
Line = X(Loc(1,1):Loc(1,2),Loc(2,1):Loc(2,2));
Line_True = Contrast(Loc(1,1):Loc(1,2),Loc(2,1):Loc(2,2));
Padding = 30;
Line_pad = padarray(Line_True,[0,Padding],'replicate','both');

gauss_sigma = [0.1:0.01:5.0];% 5.0];
for kk = 1:length(gauss_sigma)
    gauss_func  = normpdf([-15:1:15],0,gauss_sigma(kk));%  * sqrt(2*pi) * gauss_sigma;
    gauss_func = gauss_func/sum(gauss_func);
    % figure;plot(gauss_func);

    Line_conv = conv(Line_pad,gauss_func,'same');
    Line_conv = Line_conv(1+Padding:end-Padding);
    % figure;plot(Line);hold on; plot(Line_True);plot(Line_conv);
    Res(kk) = norm(Line_conv - Line);
end
figure;plot(gauss_sigma,Res,'-');

[Min_res,pos] = min(Res);
gauss_sigma(pos)
gauss_func  = normpdf([-15:1:15],0,gauss_sigma(pos));%  * sqrt(2*pi) * gauss_sigma;
gauss_func = gauss_func/sum(gauss_func);
Line_conv = conv(Line_pad,gauss_func,'same');
Line_conv = Line_conv(1+Padding:end-Padding);
figure;plot(Line);hold on; plot(Line_True);plot(Line_conv);

std(X(:) - Contrast(:))
%% 5) smoothed TV (full weight)
y = permute(ATotal,[3,2,1]); y = y(:);
WNew = opKron(structGeo.W,opEye(Nb));
Xwarm = permute(Image_FBP,[3,1,2]); Xwarm = Xwarm(:);

Xwarm = X_LBFGS;

ErrorW = Cov * ones(Np * Nb,1);
ErrorW = reshape(ErrorW,Nb,structProj.nAngle,structProj.nChannel);
AlphaPara1 = mean(mean(ErrorW.^2,2),3);
AlphaPara2 = squeeze(mean(mean(ATotal,2),1));
tao = [0,0,0.01] * 0; % AlphaPara1.*AlphaPara2 * 1;

% Preconditioning with a unweighted guess improves the convergence rate!
Data.A = Cov * WNew; 
Data.b = Cov * y; 
Data.NImage = structImg.nPixel;
Data.Nb = Nb;
Data.gamma = (1*10^(-5))^2;
Data.tao = tao;

%  Preconditioning
PreCon = ones(Nim * Nim * Nb,1); % 
PreConM = opDiag(PreCon); 
Data.PreCon = PreCon;

NIter = ones(1,20)*20;
NMemory = 20;
Func_Obj = @TV_denoising_Obj_E;
Func_Grad = @TV_denoising_Grad_E;

tic
X_LBFGS = Xwarm;
Obj_BFGS = [];
XTotal = [];figure(10);

for mm = 1:30
    
    Data.PreCon = ones(Nb * Nim * Nim,1); 
    PreCon = sqrt(10^6./TV_denoising_Grad2_E(X_LBFGS(:),Data));

    PreConM = opDiag(PreCon); 
    Data.PreCon = PreCon;
% % 
    [X_LBFGS,Obj,XTotal_LBFGS] = Quasi_Newton_LBFGS...
        (X_LBFGS(:)./PreCon, NIter(mm),NMemory, Func_Obj, Func_Grad, Data,[]); % Image_FBP,X_LBFGS
    
    
    for m = 1:NIter(mm)
        XTotal(:,m + (mm - 1) * NMemory) = PreConM * XTotal_LBFGS(:,m);
    end
    
    X_LBFGS = PreConM * X_LBFGS;


    Obj_BFGS = [Obj_BFGS,Obj];
    Obj_BFGS(end)

    X = reshape(X_LBFGS,Nb,Nim,Nim);
    X = permute(X,[2,3,1]);

    figure(8);imagesc(X(:,:,1),[0 0.002]);colorbar;
    figure(9);imagesc(X(:,:,2),[0.015 0.020]);colorbar;
    figure(10);imagesc(X(:,:,3),[0 0.002]);colorbar;
        
        
end
save(['TVCovFull_tao_',num2str(tao(3)),'.mat'],'X','Obj_BFGS','XTotal');

figure;loglog(Obj_BFGS);hold on;

%%
X = squeeze(X(:,:,3));
% Noise-Resolution trade-off
Area = X(60:70,60:70);
std(Area(:))
% Resolution
 Loc = [97,97;35,45];
% Loc = [45,45;27,36];
Line = X(Loc(1,1):Loc(1,2),Loc(2,1):Loc(2,2));
Line_True = Contrast(Loc(1,1):Loc(1,2),Loc(2,1):Loc(2,2));
Padding = 30;
Line_pad = padarray(Line_True,[0,Padding],'replicate','both');

gauss_sigma = [0.1:0.01:5.0];% 5.0];
for kk = 1:length(gauss_sigma)
    gauss_func  = normpdf([-15:1:15],0,gauss_sigma(kk));%  * sqrt(2*pi) * gauss_sigma;
    gauss_func = gauss_func/sum(gauss_func);
    % figure;plot(gauss_func);

    Line_conv = conv(Line_pad,gauss_func,'same');
    Line_conv = Line_conv(1+Padding:end-Padding);
    % figure;plot(Line);hold on; plot(Line_True);plot(Line_conv);
    Res(kk) = norm(Line_conv - Line);
end
figure;plot(gauss_sigma,Res,'-');

[Min_res,pos] = min(Res);
gauss_sigma(pos)
gauss_func  = normpdf([-15:1:15],0,gauss_sigma(pos));%  * sqrt(2*pi) * gauss_sigma;
gauss_func = gauss_func/sum(gauss_func);
Line_conv = conv(Line_pad,gauss_func,'same');
Line_conv = Line_conv(1+Padding:end-Padding);
figure;plot(Line);hold on; plot(Line_True);plot(Line_conv);

std(X(:) - Contrast(:))

%% Generate a true composition signal for comparison
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
plot(ATotal(:,1,3));hold on;plot(true(1,:));  

% true, Local concentration  (g/cm^3)

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
