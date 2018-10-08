% Thesis Chapter Two
% Single-Energy Experiment
clc; clear; close all;
set(0, 'DefaultLineLineWidth', 2);
rand('seed',1);  % Make sure the experiments are repeatable!

%% Set Image
structImg.nPixel = 256;% 512;
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
structImg.ITrue = phantom(toft,structImg.nPixel);
Nim = structImg.nPixel;
figure;imshow(structImg.ITrue,[MuBrain - 0.001, MuBrain + 0.001]);colorbar;

%% Set geometry
structProj.nChannel = 672;
structProj.nChnCenter = structProj.nChannel / 2 + 0.5;
structProj.nAngle = 50;% 290; 
structProj.fFOV = 300;
structProj.fDetWidth = structProj.fFOV / structProj.nChannel;
% Set geometries for astra
structGeo.vol_geom = astra_create_vol_geom(...
         structImg.nPixel, structImg.nPixel,...
        -structImg.fFOV/2, structImg.fFOV/2,...
        -structImg.fFOV/2, structImg.fFOV/2);
structGeo.proj_geom = astra_create_proj_geom('parallel', ...
        structProj.fDetWidth, structProj.nChannel, ...
        linspace(0,2*pi,structProj.nAngle));
structGeo.proj_id = astra_create_projector('linear', ...
    structGeo.proj_geom, structGeo.vol_geom);
structGeo.W = opTomo('linear', structGeo.proj_geom, structGeo.vol_geom); 

%% Compute Sinogram
[sinogram_id, sinogram] = astra_create_sino(...
    structImg.ITrue, structGeo.proj_id);

%% Compute noisy sinogram
% Photon statistics, Constant Photon_In
Photon_In = 8 * 10^6; % 2 * 10^5;
Photon = Photon_In .* exp(- sinogram);
Photon_rnd = poissrnd(Photon);
disp(['Minimum detected photon number is ',num2str(min(Photon_rnd(:))),'.']);
sino_noisy = log(Photon_In./Photon_rnd);
figure;plot(sino_noisy(1,:));hold on;plot(sinogram(1,:));

%% Compute Weight matrix for error term
Cov = opDiag(sqrt(Photon_rnd(:)));

%% FBP reconstruction
theta = [0:360/structProj.nAngle: 360-360/structProj.nAngle];
Image_FBP = iradon(sino_noisy',theta,'spline','Hann',0.5,structProj.nChannel);
Image_FBP = Image_FBP * structProj.nChannel / structImg.fFOV;
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
Image_FBP = Image_FBP .* I_Mask;
figure;imshow(Image_FBP,[MuBrain - 0.001, MuBrain + 0.001]);colorbar;

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
CovIm_Row = structGeo.W' * Cov * Cov * structGeo.W * ones(Nim * Nim,1);
CovIm_Row = 1./CovIm_Row;
CovIm_Row = CovIm_Row * 10^10;
CovIm_Row = sqrt(CovIm_Row);

% %% Tikhonov smoothing regularization 
% tol = 1;
% Xwarm = Image_FBP(:);
% y = sino_noisy(:);
% NIter = 50;
% % Statistical weight matrix
% ErrorW = Cov; % opDiag(ones(structProj.nChannel * structProj.nAngle,1)); % Cov;
% % Preconditioning matrix
% PreCon = CovIm_Row(:); % CovIm_Jac(:); % ones(Nim * Nim,1);
% PreConM = opDiag(PreCon); 
% 
% A  = ErrorW * structGeo.W * PreConM;
% R = FirstDerivativeM(Nim) * PreConM;
% y0 =  ErrorW * y;
% yR = zeros(size(R,1),1);
% 
% XwarmP = Xwarm ./ PreCon;
% tao = [500:500:3000];
% for k = 1:length(tao)
%     [X0,Obj,XTotal] = quasinewton([A;tao(k) * R],...
%         [y0;yR],XwarmP,NIter,tol);
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
% figure;plot(Error,Regu,'*');


%% TV Regularization (Primal-dual algorithm)
clear Regu,clear Obj,clear Error
I_noisy = Image_FBP;
% Function handles
Amplitude = @(u)sqrt(sum(u.^2,3));
F = @(u)sum(sum(Amplitude(u)));
ProxF = @(u,lambda)max(0,1-lambda./repmat(Amplitude(u), [1 1 2])).*u;
ProxFS = @(y,sigma)y-sigma*ProxF(y/sigma,1/sigma);
ProxG = @(f,tau,Lamda)  (f + tau * Lamda * I_noisy)/(1 + tau * Lamda);
                        %f + Phi(y - Phi(f));
% parameters
L = 8;
% % Algorithm 1~2
% sigma = 160;% 40;   % This significantly influence the convergence rate!
% tau = .9/(L*sigma);
% AHMOD
tau = 4.824789442826363e-04 * 8;% 80;
sigma = 4 / (L * tau);

theta = 1;
Lamda = 10^0 * 8;
% initial valius
f = f1;% zeros(Nim,Nim);% I_noisy;
g = g;% grad(I_noisy)*0;
f1 = f;
NIter = 120;
gamma = 0.7 * Lamda * 36;% 1.2;

% Tikhonov sub-step parameters
tol = 1;
y = sino_noisy(:);
NIter_inner = 2;
% Statistical weight matrix
ErrorW = opDiag(ones(structProj.nChannel * structProj.nAngle,1));
       % opDiag(ones(structProj.nChannel * structProj.nAngle,1)); % Cov;
% Preconditioning matrix
PreCon = CovIm_Row(:); % CovIm_Row(:); % CovIm_Jac(:); % ones(Nim * Nim,1);
PreConM = opDiag(PreCon); 
A  = ErrorW * structGeo.W * PreConM;
R = opEye(Nim * Nim) * PreConM;
y0 =  ErrorW * y;

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
    para = 1/sqrt(2 * tau * Lamda);
    [X0,Obj_inner,XTotal] = quasinewton([A;para * R],...
        [y0;para * yR],XwarmP,NIter_inner,tol);
    X = PreConM * X0;
    X_Image = reshape(X,structImg.nPixel,structImg.nPixel);
%     for m = 1:NIter_inner
%         XTotal(:,m) = PreConM * XTotal(:,m);
%     end
    f = X_Image;
    
%     % Algorithm 1~2
%     theta = 1/sqrt(1+2*gamma*tau);
%     tau = tau * theta;
%     sigma = sigma / theta;
%     f1 = f + theta * (f-fold);
    
    % AHMOD
    theta = 1;%1/sqrt(1+2*gamma*tau);
    tau = tau * theta;
    sigma = sigma / theta;
    f1 = f;  
    
    
    imshow(f1,[MuBrain - 0.001, MuBrain + 0.001]);pause(0.05);
    % Obj(k) = Lamda * norm(f1 - I_noisy)^2 + F(grad(f1)); % denoising
    % reconstruction
    Error(k) = norm(ErrorW * structGeo.W * f1(:) - ErrorW * y)^2;
    Regu(k) = F(grad(f1));
    Obj(k) =  Lamda * Error(k) + Regu(k);
    Error,Regu,Obj
end
figure;imshow(f1,[MuBrain - 0.001, MuBrain + 0.001]);colorbar;
figure;
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




