% Thesis Chapter Two
% Single-Energy Experiment
clc; 
clear; close all;
set(0, 'DefaultLineLineWidth', 2);
% set(0, 'defaultFigureUnits', 'normalized')
% set(0, 'defaultFigurePosition', [0 0 1 1])
rng('default');rng(1); % Make sure the experiments are repeatable!

%% Set Image
structImg.nPixel = 256;% 256;% 512;
structImg.nPixel_True = structImg.nPixel * 1;  % Analytical simulation!?
structImg.fFOV = 300;
% Set attenuation coefficients for different tissues at 80keV;
% See Table 1.1 in the Chapter2;
MuBone = 0.0428;        % Cortical Bone
MuWater = 0.0184;       % Water
MuBrain = 0.0190;       % Brain
MuBlood = 0.0194;       % Blood
toft = [ MuBone              .69   .92    0     0     0   
        -MuBone + MuBrain   .6624 .8740   0  -.0184   0        
        -MuBrain + MuWater  .1100 .3100  .22    0    -18
        -MuBrain + MuWater  .1600 .4100 -.22    0     18
         MuBlood - MuBrain  .2100 .2500   0    .35    0
         MuBlood - MuBrain  .0460 .0460   0    .1     0
         MuBlood - MuBrain  .0460 .0460   0   -.1     0
         MuBlood - MuBrain  .0460 .0230 -.08  -.605   0 
         MuBlood - MuBrain  .0230 .0230   0   -.606   0
         MuBlood - MuBrain  .0230 .0460  .06  -.605   0];

% toft = [ MuWater 1 1 0 0 0];
% Generate a modified Shepp-Logan phantom;

structImg.ITrue = phantom(toft,structImg.nPixel_True);
structImg.I     = phantom(toft,structImg.nPixel);
Nim = structImg.nPixel;
figure;imshow(structImg.I,[MuBrain - 0.001, MuBrain + 0.001]);% colorbar;
% export_fig Figure1a_TruePhantom

%% Set geometry
structProj.nChannel = 672;
structProj.nChnCenter = structProj.nChannel / 2 + 0.5;
structProj.nAngle = 290/5; % 290/5; % 50;% 290; 
structProj.fFOV = 300;
structProj.fDetWidth = structProj.fFOV / structProj.nChannel;
% Set geometries for astra
% projection geometry
structGeo.proj_geom = astra_create_proj_geom('parallel', ...
        structProj.fDetWidth, structProj.nChannel, ...
        linspace(0,2*pi,structProj.nAngle));
% reconstruction image geometry    
structGeo.vol_geom = astra_create_vol_geom(...
         structImg.nPixel, structImg.nPixel,...
        -structImg.fFOV/2, structImg.fFOV/2,...
        -structImg.fFOV/2, structImg.fFOV/2);
structGeo.W = opTomo('linear', structGeo.proj_geom, structGeo.vol_geom); 

% structGeo.proj_id = astra_create_projector('linear', ...
%     structGeo.proj_geom, structGeo.vol_geom);
% matrix_id = astra_mex_projector('matrix', structGeo.proj_id);
% structGeo.Wf = astra_mex_matrix('get', matrix_id);

% true image geometry    
structGeo.vol_geom_true = astra_create_vol_geom(...
         structImg.nPixel_True, structImg.nPixel_True,...
        -structImg.fFOV/2, structImg.fFOV/2,...
        -structImg.fFOV/2, structImg.fFOV/2);
structGeo.proj_id_true = astra_create_projector('linear', ...
    structGeo.proj_geom, structGeo.vol_geom_true);
%% Compute Sinogram
[sinogram_id, sinogram] = astra_create_sino(...
    structImg.ITrue, structGeo.proj_id_true); % Discrete

% Analytic (They seem to be very different....)
% ProjData = EllipseProj(toft, structProj);     
% sinogram = ProjData; 
% figure;plot(sinogram(1,:));hold on;plot(ProjData(1,:));
% pause;

%% Compute noisy sinogram
% Photon statistics, Constant Photon_In
Photon_In = 5 * 2.5 * 10^5; % 1.0 * 10^6;
% 2 * 5 * 2 * 10^5; % 5 * 2 * 10^5; % 8 * 10^6; % 2 * 10^5;
Photon = Photon_In .* exp(- sinogram);
Photon_rnd = poissrnd(Photon);
disp(['Minimum detected photon number is ',num2str(min(Photon_rnd(:))),'.']);
sino_noisy = log(Photon_In./Photon_rnd);
figure;plot(sino_noisy(1,:));hold on;plot(sinogram(1,:));

%% Compute Weight matrix for error term
Cov = opDiag(sqrt(Photon_rnd(:)));

%% FBP reconstruction
theta = [0:360/structProj.nAngle: 360-360/structProj.nAngle];
Image_FBP = iradon(sino_noisy',theta,'spline','Hann',0.20,structProj.nChannel);
% 290 projections: Hanning 0.2; % 58 projections: Hanning 0.15;

Image_FBP = Image_FBP * structProj.nChannel / structImg.fFOV;
% figure;imshow(Image_FBP,[MuBrain - 0.001, MuBrain + 0.001]);
Image_FBP = imresize(Image_FBP,[structImg.nPixel,structImg.nPixel]);
% Circle mask for image
N = Nim;
I_Mask = zeros(N, N);
for m = 1:N
    for n = 1:N
        if ((m-(N/2+0.5))^2 + (n-(N/2+0.5))^2<(N/2)^2)
            I_Mask(m,n) = 1;
        end
    end
end
Circle_Label = I_Mask>0;
Circle_Label = ones(Nim,Nim);
Circle_Label  = Circle_Label>0;%%%
structGeo.Wf = structGeo.W(:,Circle_Label(:));
% structGeo.Wf = opMatrix(structGeo.Wf);

Image_FBP = Image_FBP .* I_Mask;
figure;imshow(Image_FBP,[MuBrain - 0.001, MuBrain + 0.001]);
% colorbar;
% Colorbar [0.018 0.02] mm^(-1).
% export_fig Figure1b_FBP_290Views
% export_fig Figure1c_FBP_58Views
% pause;

%% preconditioning (Method One) for Tikhonov-regularized problem
%----Result: Convergence slows down after 20 iterations (Obselete!)
% N_Coarse_grid = 64;
% structGeoC.vol_geom = astra_create_vol_geom(...
%          N_Coarse_grid, N_Coarse_grid,...
%         -structImg.fFOV/2, structImg.fFOV/2,...
%         -structImg.fFOV/2, structImg.fFOV/2);
% structGeoC.proj_geom = astra_create_proj_geom('parallel', ...
%         structProj.fDetWidth, structProj.nChannel, ...
%         linspace(0,2*pi,structProj.nAngle));
% structGeoC.proj_id = astra_create_projector('linear', ...
%     structGeoC.proj_geom, structGeoC.vol_geom);
% structGeoC.W = opTomo('linear', structGeoC.proj_geom, structGeoC.vol_geom); 
% matrix_id = astra_mex_projector('matrix', structGeoC.proj_id);
% W = astra_mex_matrix('get', matrix_id);
% 
% CovIm_Jac = full( W' * Cov' * Cov * W); %%
% CovIm_Jac = diag(CovIm_Jac);
% CovIm_Jac = sqrt(1./ CovIm_Jac);
% CovIm_Jac = reshape(CovIm_Jac,[N_Coarse_grid,N_Coarse_grid]);
% CovIm_Jac = imresize(CovIm_Jac,[structImg.nPixel,structImg.nPixel]);
% CovIm_Jac = CovIm_Jac';
% clear W

%% preconditioning (Method Two) for Tikhonov-regularized problem
%----Result:Convergence is much faster! (Choose this one currently!)
% Convergence is only faster with an over-smoothed image;

CovIm_Row = structGeo.W' * Cov * Cov * structGeo.W * ones(Nim * Nim,1);
CovIm_Row = 1./CovIm_Row;
CovIm_Row = CovIm_Row * 10^10;
CovIm_Row = sqrt(CovIm_Row);

CovIm_Row_NoCov = structGeo.W' * structGeo.W * ones(Nim * Nim,1);
CovIm_Row_NoCov = 1./CovIm_Row_NoCov;
CovIm_Row_NoCov = CovIm_Row_NoCov * 10^10;
CovIm_Row_NoCov = sqrt(CovIm_Row_NoCov);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% Tikhonov smoothing regularization 
% tol = 1;
% Xwarm = Image_FBP(:);
% y = sino_noisy(:);
% NIter = 50;
% tao = 5 * 2500 * 2/sqrt(5);% 2500 * 2;% [500:500:3000];
% 
% % Statistical weight matrix
% ErrorW = Cov; % opDiag(ones(structProj.nChannel * structProj.nAngle,1)); % Cov;
% % Preconditioning matrix
% PreCon = ones(Nim * Nim,1); % CovIm_Row(:); % CovIm_Jac(:); % ones(Nim * Nim,1);
% PreConM = opDiag(PreCon); 
% 
% A  = ErrorW * structGeo.W * PreConM;
% R = FirstDerivativeM(Nim) * PreConM;
% y0 =  ErrorW * y;
% yR = zeros(size(R,1),1);
% 
% XwarmP = Xwarm ./ PreCon;
% 
% for k = 1:length(tao)
%     tic
%     [X0,Obj,XTotal] = quasinewton([A;tao(k) * R],...
%         [y0;yR],XwarmP,NIter,tol);
%     toc
%     
%     X = PreConM * X0;
%     X_Image = reshape(X,structImg.nPixel,structImg.nPixel);
%     for m = 1:NIter
%         XTotal(:,m) = PreConM * XTotal(:,m);
%     end
%     save(['Tikhonov_',num2str(tao(k)),'.mat'],'X_Image','Obj');
%     Error(k) = norm(A * X0 - y0);
%     Regu(k) = norm(R * X);
%     
%     figure;imshow(X_Image,[MuBrain - 0.001, MuBrain + 0.001]);colorbar;
%     title(['\alpha=',num2str(tao(k))]);
% end
% 
% Error,Regu
% % figure;plot(Error,Regu,'*');

% %% Tikhonov smoothing regularization (New--Acceleration)
% Xwarm = Image_FBP(:); % X_Image(:); % Image_FBP(:);
% y = sino_noisy(:);
% NIter = 50;
% NMemory = 50;
% tao = 2000; % 500; % 50; % 5 * 2500 * 2/sqrt(5); % 2500 * 2;% [500:500:3000];
% % TK290:[1000?5000];
% % TK58: [];
% 
% % Statistical weight matrix
% ErrorW = Cov; % Cov;
% % Cov; % opDiag(ones(structProj.nChannel * structProj.nAngle,1)); % Cov;
% 
% % % Renewed Preconditioning matrix (Before 20180903)
% % CovIm_Row_Cov = structGeo.W' * ErrorW' * ErrorW * structGeo.W * ones(Nim * Nim,1);
% % RTemp = FirstDerivativeM(Nim);
% % CovIm_Row_Cov = CovIm_Row_Cov + tao^2 * diag(RTemp' * RTemp);
% % CovIm_Row_Cov = 1./CovIm_Row_Cov;
% % CovIm_Row_Cov = CovIm_Row_Cov * 10^10;
% % CovIm_Row_Cov = sqrt(CovIm_Row_Cov);
% 
% % Renewed Preconditioning matrix (After 20180903)
% CovIm_Row_Cov = structGeo.W' * ErrorW' * ErrorW * structGeo.W * ones(Nim * Nim,1);
% RTemp = FirstDerivativeM(Nim);
% CovIm_Row_Cov = CovIm_Row_Cov + tao^2 * (RTemp' * RTemp * ones(Nim * Nim,1));
% CovIm_Row_Cov = 1./CovIm_Row_Cov;
% CovIm_Row_Cov = CovIm_Row_Cov * 10^10;
% CovIm_Row_Cov = sqrt(CovIm_Row_Cov);
% 
% % Preconditioning matrix
% PreCon = CovIm_Row_Cov; % ones(Nim * Nim,1);% CovIm_Row_Cov; % ones(Nim * Nim,1);% CovIm_Row_Cov;
% % CovIm_Row_Cov(:); CovIm_Row(:); % CovIm_Jac(:); % ones(Nim * Nim,1);
% PreConM = opDiag(PreCon); 
% 
% A  = ErrorW * structGeo.W * PreConM;
% R = FirstDerivativeM(Nim) * PreConM;
% y0 =  ErrorW * y;
% yR = zeros(size(R,1),1);
% 
% XwarmP = Xwarm ./ PreCon;
% 
% Func_Diff = @(X,Data)(Data.A * X - Data.y);
% Func_Obj = @(X,Data) norm(Func_Diff(X,Data))^2;
% Func_Grad =  @(X,Data) (2 * Data.A' * Func_Diff(X,Data));
% isQuadratic = 1;
% 
% for k = 1:length(tao) 
%     Data.A = [A;tao(k) * R];
%     Data.y = [y0;yR];
%     
%     tic
%     [X0,Obj_BFGS,XTotal] = Quasi_Newton_LBFGS...
%         (XwarmP, NIter, NMemory,Func_Obj, Func_Grad, Data,isQuadratic);
%     toc
% 
%     X = PreConM * X0;
%     X_Image = reshape(X,structImg.nPixel,structImg.nPixel);
%     for m = 1:NIter
%         XTotal(:,m) = PreConM * XTotal(:,m);
%     end
%     save(['Tikhonov290_',num2str(tao(k)),'.mat'],'X_Image','Obj_BFGS');
%     Error(k) = norm(A * X0 - y0)^2;
%     Regu(k) = (norm(R * X))^2;
%     
%     figure;imshow(X_Image,[MuBrain - 0.001, MuBrain + 0.001]);colorbar;
%     title(['\alpha=',num2str(tao(k))]);
% end
% 
% Error, Regu, Obj_BFGS(end)
% % figure;plot(Error,Regu,'*');
% figure;loglog(Obj_BFGS);hold on;
% % plot(log(Obj));
% 
% pause;
% 
% %%%%% SNR in ROI %%%%%%%%%%%%%%%%%%%%%
% 
% I_Mask_SNR = (structImg.I < 0.04).*(structImg.I > 0);
% se = strel('square',11);
% % Erode the image with the structuring element.
% I_Mask_SNR_New = imerode(I_Mask_SNR,se);
% I_Mask_SNR_New = I_Mask_SNR;
% figure;imshow(I_Mask_SNR - I_Mask_SNR_New,[]);
% figure;imshow(X_Image.* I_Mask_SNR_New,[0.018 0.020]);
% 
% sum(norm(X_Image(:) - structImg.I(:)))
% sum(norm((X_Image(:) - structImg.I(:)).*I_Mask_SNR_New(:)))
% 
% %%%%% Enlarged regions %%%%%
% X = X_Image;
% Region = [191,220;101,150];
% Grid_X = ([Nim-Region(1,1):-1:Nim-Region(1,2)] - Nim/2 - 0.5) * (structImg.fFOV / Nim);
% Grid_Y = ([Region(2,1):Region(2,2)] - Nim/2 - 0.5) * (structImg.fFOV / Nim);
% 
% figure;imagesc(Grid_Y,Grid_X(end:-1:1),X(Region(1,1):Region(1,2),Region(2,1):Region(2,2)),...
%     [0.0190 0.0195]);colormap(gray);hold on;
% 
% % [MuBrain - 0.001*0.5, MuBrain + 0.001*0.5]
% 
% ThetaPlot = [0:pi * 2 / 300:pi*2];
% for ii = 0:2
%     aa = toft(end-ii,2) * structImg.fFOV/2;
%     bb = toft(end-ii,3) * structImg.fFOV/2;
%     x0 = toft(end-ii,4) * structImg.fFOV/2;
%     y0 = toft(end-ii,5) * structImg.fFOV/2;
%     plot(x0 + aa * cos(ThetaPlot),y0 + bb * sin(ThetaPlot),'r--');hold on;
% end
% 
% %%%% Resolution %%%%%
% Loc = [85,85;145,165];
% %%%% Loc = [67,67;96,115];
% %%%% Loc = [132,132;160,179];
% Loc = [133,133;71,90];
% Line = X(Loc(1,1):Loc(1,2),Loc(2,1):Loc(2,2));
% Line_True = structImg.I(Loc(1,1):Loc(1,2),Loc(2,1):Loc(2,2));
% Padding = 30;
% Line_pad = padarray(Line_True,[0,Padding],'replicate','both');
% 
% gauss_sigma = [0.1:0.01:5.0];% 5.0];
% for kk = 1:length(gauss_sigma)
%     gauss_func  = normpdf([-15:1:15],0,gauss_sigma(kk));%  * sqrt(2*pi) * gauss_sigma;
%     gauss_func = gauss_func/sum(gauss_func);
%     % figure;plot(gauss_func);
% 
%     Line_conv = conv(Line_pad,gauss_func,'same');
%     Line_conv = Line_conv(1+Padding:end-Padding);
%     % figure;plot(Line);hold on; plot(Line_True);plot(Line_conv);
%     Res(kk) = norm(Line_conv - Line);
% end
% figure;plot(gauss_sigma,Res,'-');
% figure;plot(Line);hold on; plot(Line_True);plot(Line_conv);
% [Min_res,pos] = min(Res);
% gauss_sigma(pos)
% 
% %%%% Noise %%%%%
% Flat_Region = [76,85;171,180];
% X_Flat = X(Flat_Region(1,1):Flat_Region(1,2),Flat_Region(2,1):Flat_Region(2,2));
% figure;imagesc(X_Flat);colormap(gray);colorbar;
% std(X_Flat(:))
% 
% pause;


%% TV Regularization (Corner Smoothing) --> transfer into a convex problem
% % 1) Periodic restart is required for fast convergence when solving
% % non-quadratic problems!!!!
% % 2) Should consider unweighted results as initial guess for weighted
% % problem, and weighted problem should give better results, but converge
% % very slowly. With FBP initial, weighted problem requires more than 500
% % iters to converge!
% % -------------------------------------------------------------------------
% 
% 1) Compute an initial guess (unweighted guess)
Data.A = structGeo.W; % Cov * structGeo.W;
Data.b = sino_noisy(:); % Cov * sino_noisy(:);
Data.NImage = structImg.nPixel;
Data.gamma = (1*10^(-5))^2;
Data.tao = 2; % 1; % 10^0 * 8;% 10^0 * 8;
% tao=1 for TK290; tao=? for TK58;

%  Preconditioning
PreCon = ones(Nim * Nim,1);% sqrt(10^10./TV_denoising_Grad2(Image_FBP(:),Data));
% % CovIm_Row_NoCov(:); % CovIm_Row(:); % CovIm_Jac(:); % ones(Nim * Nim,1);
PreConM = opDiag(PreCon); 
Data.PreCon = PreCon;

NIter = 1;% ones(1,10)*20;
NMemory = 1;% 20;
Func_Obj = @TV_denoising_Obj;
Func_Grad = @TV_denoising_Grad;

tic
X_LBFGS = Image_FBP;
Obj_BFGS = [];figure(10);
for mm = 1:1
    
%     PreCon =  sqrt(10^10./TV_denoising_Grad2(X_LBFGS(:),Data));
%     % % CovIm_Row_NoCov(:); % CovIm_Row(:); % CovIm_Jac(:); % ones(Nim * Nim,1);
%     PreConM = opDiag(PreCon); 
%     Data.PreCon = PreCon;

    [X_LBFGS,Obj,XTotal_LBFGS] = Quasi_Newton_LBFGS...
        (X_LBFGS(:)./PreCon, NIter(mm),NMemory, Func_Obj, Func_Grad, Data,[]); % Image_FBP,X_LBFGS
    
    
    
%     [X_LBFGS,Obj,XTotal_LBFGS] = Quasi_Newton_BFGS...
%         (X_LBFGS(:)./PreCon, NIter(mm), Func_Obj, Func_Grad, Data,[]); % Image_FBP,X_LBFGS   
    
    X_LBFGS = PreConM * X_LBFGS;
    X = X_LBFGS;
    Obj_BFGS = [Obj_BFGS,Obj];
    Obj_BFGS(end)
%     figure;plot(log(Obj_BFGS));
    X = reshape(X,structImg.nPixel,structImg.nPixel);
    X = X .* I_Mask;
    imshow(X,[MuBrain - 0.001, MuBrain + 0.001]);colorbar;
end
toc

Data.gamma = 0;
TV_denoising_Obj(X(:),Data)

figure;plot(log(Obj_BFGS));
X = reshape(X,structImg.nPixel,structImg.nPixel);
figure;imshow(X,[MuBrain - 0.001, MuBrain + 0.001]);colorbar;
X_Warm_TV = Image_FBP; % X;
 
%% 1) Initial guess (unweighted)(ROI region!!!)~~~
% Preconditioning with a unweighted guess improves the convergence rate!
Data.A = structGeo.Wf; % Cov * structGeo.W;
Data.b = sino_noisy(:); % Cov * sino_noisy(:);
Data.NImage = structImg.nPixel;
Data.Circle_Label = Circle_Label;
Data.gamma = (1*10^(-5))^2;
Data.tao = 2.5;% 3.125000000000000e-04;% 2 * 2 * 10^(-4);
% 4 * 10^(-4);% 10^0 * 8/20000;% 10^0 * 8;

%  Preconditioning
PreCon = ones(sum(Data.Circle_Label(:)),1); 
% sqrt(10^6./TV_denoising_Grad2_Circle(X_LBFGS(:),Data));
PreConM = opDiag(PreCon); 
Data.PreCon = PreCon;

NIter = ones(1,20)*20;
NMemory = 20;
Func_Obj = @TV_denoising_Obj_Circle;
Func_Grad = @TV_denoising_Grad_Circle;

tic
X_LBFGS = Image_FBP(Circle_Label);
Obj_BFGS = [];figure(10);
XTotal = [];

for mm = 1:2
    
%     Data.PreCon = ones(sum(Data.Circle_Label(:)),1); 
%     PreCon = sqrt(10^6./TV_denoising_Grad2_Circle(X_LBFGS(:),Data));
%     % sqrt(10^6./TV_denoising_Grad2_Circle(X_LBFGS(:),Data));
%     % ones(length(Data.Circle_Label),1); 
%     PreConM = opDiag(PreCon); 
%     Data.PreCon = PreCon;

    [X_LBFGS,Obj,XTotal_LBFGS] = Quasi_Newton_LBFGS...
        (X_LBFGS(:)./PreCon, NIter(mm),NMemory, Func_Obj, Func_Grad, Data,[]); % Image_FBP,X_LBFGS
    
    for m = 1:NIter(mm)
        XTotal(:,m) = PreConM * XTotal_LBFGS(:,m);
        MSE(m + (mm - 1) * NMemory) = sum(norm(XTotal(:,m) - structImg.I(Circle_Label)));
    end
    
    Obj_BFGS = [Obj_BFGS,Obj];
    Obj_BFGS(end)
    
    X_LBFGS = PreConM * X_LBFGS;
    
    X = zeros(structImg.nPixel,structImg.nPixel);
    X(Circle_Label) = X_LBFGS;
    figure(10);imshow(X,[MuBrain - 0.001, MuBrain + 0.001]);colorbar;
    
end
toc

% save(['TV58_',num2str(1/Data.tao),'.mat'],'X','Obj_BFGS');

% Data.gamma = 0;

Data1 = Data;
Data1.PreCon = ones(sum(Data.Circle_Label(:)),1);
Error = norm(Data1.A * X(Circle_Label) - Data1.b)
Obj = Error * Error * Data1.tao;
Data1.tao = 0;
Regu = TV_denoising_Obj_Circle(X(Circle_Label(:)),Data1)
Obj = Obj+ Regu
Obj_BFGS(end)
1/Data.tao

figure;plot(log(Obj_BFGS));
X = reshape(X,structImg.nPixel,structImg.nPixel);
figure;imshow(X,[MuBrain - 0.001, MuBrain + 0.001]);colorbar;
sum(norm(X(:) - structImg.I(:)))
X_Warm_TV = X;

% 2) Refine the result now (ROI region!!!)~~
% Preconditioning with a unweighted guess improves the convergence rate!
Data.A = Cov * structGeo.Wf; % Cov * structGeo.W;
Data.b = Cov * sino_noisy(:); % Cov * sino_noisy(:);
Data.NImage = structImg.nPixel;
Data.Circle_Label = Circle_Label;
Data.gamma = (1*10^(-5))^2;
Data.tao = 1 * 1/4000;% 3.125000000000000e-04;% 2 * 2 * 10^(-4);
% 4 * 10^(-4);% 10^0 * 8/20000;% 10^0 * 8;

%  Preconditioning
PreCon = ones(sum(Data.Circle_Label(:)),1); 
% sqrt(10^6./TV_denoising_Grad2_Circle(X_LBFGS(:),Data));
PreConM = opDiag(PreCon); 
Data.PreCon = PreCon;

NIter = ones(1,20)*20;
NMemory = 20;
Func_Obj = @TV_denoising_Obj_Circle;
Func_Grad = @TV_denoising_Grad_Circle;

tic
X_LBFGS = X_Warm_TV(Data.Circle_Label);
% X_Warm_TV(Data.Circle_Label);% Image_FBP(Data.Circle_Label);
Obj_BFGS = [];figure(10);
XTotal = [];

for mm = 1:5
    
    
%     Data.PreCon = ones(sum(Data.Circle_Label(:)),1); 
%     PreCon = sqrt(10^6./TV_denoising_Grad2_Circle(X_LBFGS(:),Data));
% %     if (mm < 3)
% %       PreCon = sqrt(10^6./TV_denoising_Grad2_Circle(X_LBFGS(:),Data));
% %     else
% %       PreCon = ones(sum(Data.Circle_Label(:)),1);
% %     end
%     % sqrt(10^6./TV_denoising_Grad2_Circle(X_LBFGS(:),Data));
%     % ones(length(Data.Circle_Label),1); 
%     PreConM = opDiag(PreCon); 
%     Data.PreCon = PreCon;

    [X_LBFGS,Obj,XTotal_LBFGS] = Quasi_Newton_LBFGS...
        (X_LBFGS(:)./PreCon, NIter(mm),NMemory, Func_Obj, Func_Grad, Data,[]); % Image_FBP,X_LBFGS
    
    for m = 1:NIter(mm)
        XTotal(:,m + (mm - 1) * NMemory) = PreConM * XTotal_LBFGS(:,m);
        MSE(m + (mm - 1) * NMemory) = sum(norm(XTotal(:,m + (mm - 1) * NMemory)...
            - structImg.I(Circle_Label)));
    end
    
    Obj_BFGS = [Obj_BFGS,Obj];
    Obj_BFGS(end)
    
    X_LBFGS = PreConM * X_LBFGS;
    
    X = zeros(structImg.nPixel,structImg.nPixel);
    X(Circle_Label) = X_LBFGS;
    figure(10);imshow(X,[MuBrain - 0.001, MuBrain + 0.001]);colorbar;
    
end
toc

% save(['TV58_',num2str(1/Data.tao),'.mat'],'X','Obj_BFGS');

% Data.gamma = 0;

Data1 = Data;
Data1.PreCon = ones(sum(Data.Circle_Label(:)),1);
Error = norm(Data1.A * X(Circle_Label) - Data1.b)
Obj = Error * Error * Data1.tao;
Data1.tao = 0;
Regu = TV_denoising_Obj_Circle(X(Circle_Label(:)),Data1)
Obj = Obj+ Regu
Obj_BFGS(end)
1/Data.tao

figure;loglog(Obj_BFGS);
X = reshape(X,structImg.nPixel,structImg.nPixel);
figure;imshow(X,[MuBrain - 0.001, MuBrain + 0.001]);colorbar;
sum(norm(X(:) - structImg.I(:)))
pause;
%% 2) Refine the result now~~
% Preconditioning with a unweighted guess improves the convergence rate!
Data.A = Cov * structGeo.W; % Cov * structGeo.W;
Data.b = Cov * sino_noisy(:); % Cov * sino_noisy(:);
Data.NImage = structImg.nPixel;
Data.gamma = (1*10^(-5))^2;
Data.tao = 1 * 1/4000;% 3.125000000000000e-04;% 2 * 2 * 10^(-4);
% 4 * 10^(-4);% 10^0 * 8/20000;% 10^0 * 8;

%  Preconditioning
PreCon = ones(Nim * Nim,1); % sqrt(10^6./TV_denoising_Grad2(X_Warm_TV(:),Data));
% % CovIm_Row_NoCov(:); % CovIm_Row(:); % CovIm_Jac(:); % ones(Nim * Nim,1);
PreConM = opDiag(PreCon); 
Data.PreCon = PreCon;

NIter = ones(1,20)*20;
NMemory = 20;
Func_Obj = @TV_denoising_Obj;
Func_Grad = @TV_denoising_Grad;

tic
X_LBFGS = Image_FBP;% X_Warm_TV;
X = X_LBFGS;
Obj_BFGS = [];
XTotal = [];figure(10);

for mm = 1:5
    
    PreCon = ones(Nim * Nim,1); % sqrt(10^6./TV_denoising_Grad2(X_LBFGS(:),Data));
    % % CovIm_Row_NoCov(:); % CovIm_Row(:); % CovIm_Jac(:); % ones(Nim * Nim,1);
    PreConM = opDiag(PreCon); 
    Data.PreCon = PreCon;

    [X_LBFGS,Obj,XTotal_LBFGS] = Quasi_Newton_LBFGS...
        (X_LBFGS(:)./PreCon, NIter(mm),NMemory, Func_Obj, Func_Grad, Data,[]); % Image_FBP,X_LBFGS
    
%     [X_LBFGS,Obj,XTotal_LBFGS] = Quasi_Newton_BFGS...
%         (X_LBFGS(:)./PreCon, NIter(mm), Func_Obj, Func_Grad, Data,[]); % Image_FBP,X_LBFGS   
    
    for m = 1:NIter(mm)
        XTotal(:,m) = PreConM * XTotal_LBFGS(:,m);
        MSE(m + (mm - 1) * NMemory) = sum(norm(XTotal(:,m) - structImg.I(:)));
    end
    
    X_LBFGS = PreConM * X_LBFGS;
    figure(7);imshow(X-reshape(X_LBFGS,Nim,Nim),[]);colorbar;
    X = X_LBFGS;
    Obj_BFGS = [Obj_BFGS,Obj];
    Obj_BFGS(end)
%     figure;plot(log(Obj_BFGS));
    X = reshape(X,structImg.nPixel,structImg.nPixel);
    figure(10);imshow(X,[MuBrain - 0.001, MuBrain + 0.001]);colorbar;
    
%     X_LBFGS(X_LBFGS<0.005) = 0; sum(norm(X_LBFGS(:) - structImg.I(:)))
%     X_LBFGS = X_LBFGS .* I_Mask(:); 
end
toc

% save(['TV58_',num2str(1/Data.tao),'.mat'],'X','Obj_BFGS');

% Data.gamma = 0;

Data1 = Data;
Data1.PreCon = ones(Nim * Nim,1);
Error = norm(Data1.A * X(:) - Data1.b)
Obj = Error * Error * Data1.tao;
Data1.tao = 0;
Regu = TV_denoising_Obj(X(:),Data1)
Obj = Obj+ Regu
Obj_BFGS(end)
1/Data.tao


    
figure;plot(log(Obj_BFGS));
X = reshape(X,structImg.nPixel,structImg.nPixel);
figure;imshow(X,[MuBrain - 0.001, MuBrain + 0.001]);colorbar;

X(X<0.005) = 0; sum(norm(X(:) - structImg.I(:)))

pause;

%% TV Regularization (Primal-dual algorithm)
% Almost there after about 500 iterations...

clear Regu,clear Obj,clear Error
I_noisy = Image_FBP;
% Function handles
% Amplitude = @(u)sqrt(sum(u.^2,3)); % ISOTROPIC
Amplitude = @(u)(sum(abs(u),3)); % ANISOTROPIC

F = @(u)sum(sum(Amplitude(u)));
ProxF = @(u,lambda)max(0,1-lambda./repmat(Amplitude(u), [1 1 2])).*u;
ProxFS = @(y,sigma)y-sigma*ProxF(y/sigma,1/sigma);
ProxG = @(f,tau,Lamda)  (f + tau * Lamda * I_noisy)/(1 + tau * Lamda);
                        %f + Phi(y - Phi(f));
% parameters
L = 8;
% % Algorithm 1~2
% sigma = 10;% 160;   % This significantly influence the convergence rate!
% tau = .9/(L*sigma);

% AHMOD
tau = 4.824789442826363e-04 * 8;% 4.824789442826363e-04 * 8;% 80;
sigma = 4 / (L * tau);

theta = 1;
Lamda = 10^0 * 8 * 2; % 10^0 * 8;
% initial values
f = I_noisy; % zeros(Nim,Nim);% I_noisy;
g = grad(I_noisy)*0;
f1 = f;
NIter = 2000;
gamma = 0.7 * Lamda * 36;% 1.2;

% Tikhonov sub-step parameters
tol = 1;
y = sino_noisy(:);
NIter_inner = 2;
% Statistical weight matrix
ErrorW = opDiag(ones(structProj.nChannel * structProj.nAngle,1));
       % opDiag(ones(structProj.nChannel * structProj.nAngle,1)); % Cov;
% Preconditioning matrix
PreCon = CovIm_Row_NoCov(:); % CovIm_Row(:); % CovIm_Jac(:); % ones(Nim * Nim,1);
PreConM = opDiag(PreCon); 
A  = ErrorW * structGeo.W * PreConM;
R = opEye(Nim * Nim) * PreConM;
y0 =  ErrorW * y;

figure(1);
p = 1;
for k = 1:NIter
    fold = f;
    
%     g = ProxFS(g + sigma * grad(f1), sigma);  % An alternative way
    Temp = g + sigma * grad(f1);  % This seems to converge faster
    g = Temp./max(ones(size(Temp)),abs(Temp));
    
    % f = ProxG(f + tau * div(g(:,:,1),g(:,:,2)), tau, Lamda); % denoising
    % reconstruction
    
    % Tikhonov smoothing regularization 
    Xwarm = f(:);
    yR = f + tau * div(g(:,:,1),g(:,:,2));
    yR = yR(:);

    XwarmP = Xwarm ./ PreCon;
    para = 1/sqrt(tau * Lamda);
    [X0,Obj_inner,XTotal] = quasinewton([A;para * R],...
        [y0;para * yR],XwarmP,NIter_inner,tol);
    X = PreConM * X0;
    X_Image = reshape(X,structImg.nPixel,structImg.nPixel);
%     for m = 1:NIter_inner
%         XTotal(:,m) = PreConM * XTotal(:,m);
%     end
    f = X_Image;
    
%     % Algorithm 1~2
%     theta = 0.5;% 1/sqrt(1+2*gamma*tau);
%     if (mod(k,2^p) == 0)
%         tau = tau * theta;
%         sigma = sigma / theta;
%         p = p + 1;
%     end
%     
%     theta = 1;
%     f1 = f + theta * (f-fold);
    
    % AHMOD
    theta = 1;%1/sqrt(1+2*gamma*tau);
    tau = tau * theta;
    sigma = sigma / theta;
    f1 = f;  
    
    figure(1);
    imshow(f1,[MuBrain - 0.001, MuBrain + 0.001]);pause(0.05);
    % Obj(k) = Lamda * norm(f1 - I_noisy)^2 + F(grad(f1)); % denoising
    % reconstruction
    Error(k) = norm(ErrorW * structGeo.W * f1(:) - ErrorW * y)^2;
    Regu(k) = F(grad(f1));
    Obj(k) =  0.5 * Lamda * Error(k) + Regu(k);
    Error,Regu,Obj
end
figure(2);imshow(f1,[MuBrain - 0.001, MuBrain + 0.001]);colorbar;
figure(3);
loglog([1:k],Obj(1:k));hold on;
loglog([1:k],20000./[1:k].^2);

%% Functions used in different algorithms
% Be really careful about the discrete gradient and divergence definition!

% Gradient
function fM = grad(M)
    fx = M([2:end end],:)-M;
    fy = M(:,[2:end end])-M;
    fM(:,:,1) = fx;
    fM(:,:,2) = fy;
end

% Divergence
function fd = div(Px,Py)
    fx = Px-Px([1 1:end-1],:,:);         
    fx(1,:,:)   = Px(1,:,:);        % boundary
    fx(end,:,:) = -Px(end-1,:,:);  
    
    fy = Py-Py(:,[1 1:end-1],:);
    fy(:,1,:)   = Py(:,1,:);        % boundary
    fy(:,end,:) = -Py(:,end-1,:);
    fd = fx+fy;
end




