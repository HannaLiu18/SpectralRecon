%% Single-Energy Experiment
% clc;
clear;
rand('seed',1);  % Make sure the experiments are repeatable!

%% Set Image
structImg.nPixel = 256;% 256;%512
structImg.fFOV = 500;
%  At 100keV, Need this because we want the result to have physical
%  meanings;
MuB = 0.01855 *  1.920; % Cortical Bone
MuW = 0.01707;
toft = [  MuB   .69   .92    0     0     0   
        -MuB + MuW  .6624 .8740   0  -.0184   0
        
%         -MuW  .1100 .3100  .22    0    -18
%         -MuW  .1600 .4100 -.22    0     18
        
        -MuW * 0.04  .1100 .3100  .22    0    -18
        -MuW * 0.04  .1600 .4100 -.22    0     18
        
         MuW * 0.04  .2100 .2500   0    .35    0
         MuW * 0.04  .0460 .0460   0    .1     0
         MuW * 0.04  .0460 .0460   0   -.1     0
         MuW * 0.04 .0460 .0230 -.08  -.605   0 
         MuW * 0.04 .0230 .0230   0   -.606   0
         MuW * 0.04  .0230 .0460  .06  -.605   0   ];

structImg.ITrue = phantom(toft,structImg.nPixel);
% structImg.ITrue = phantom(structImg.nPixel);
Nim = structImg.nPixel;


%% Set geometry
structProj.nChannel = 672;
structProj.nChnCenter = structProj.nChannel / 2 + 0.5;
structProj.nAngle = 1160/4; % 1160;
structProj.fFOV = 500;
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

% matrix_id = astra_mex_projector('matrix', structGeo.proj_id);
% W = astra_mex_matrix('get', matrix_id);
%% Compute Sinogram
[sinogram_id, sinogram] = astra_create_sino(...
    structImg.ITrue, structGeo.proj_id);

%% Compute noisy sinogram and weight matrix
% Photon statistics
%  Should consider photon filtering for different channels;
%  Otherwise difficult to simulate low-dose condition with few starvations;
Photon0 =  5 * 10^3; % 1 * 10^6;
PhotonBase = 0;
R_Filter = 250;
Channel_Loc = (R_Filter)^2 - ....
(([1:structProj.nChannel] - structProj.nChnCenter) * structImg.fFOV / structProj.nChannel).^2;
Channel_Loc(Channel_Loc < 0) = 0;
Channel_Loc = sqrt(Channel_Loc);
Channel_Loc(Channel_Loc < 0) = 0;
% 0.01928 water attentuation at 100keV
Photon_In = Photon0 * exp(0.01707 * Channel_Loc) + PhotonBase;  % Photon;
figure;plot(Photon_In);
for k = 1:structProj.nAngle
    Photon(k,:) = Photon_In .* exp(-sinogram(k,:));
end
figure;plot(log(Photon(1160/4,:)));hold on;
sum(Photon_In),min(Photon(:))
% %Constant Photon In
Photon_In = 2 * 10^5 * 2; % 2 * 10^5;
Photon = Photon_In .* exp(- sinogram);
plot(log(Photon(1160/4,:)));sum(Photon_In) * structProj.nChannel,min(Photon(:))


Photon_rnd = poissrnd(Photon);
Photon_rnd(Photon_rnd < 5) = 5;
sino_noisy = log(Photon_In./Photon_rnd);



%% Compute Weight matrix for error and regularization (Original one and actually wrong!)
Cov = opDiag(sqrt(Photon_rnd(:)));
CovIm = structGeo.W' * Cov * Cov * structGeo.W * ...
    opEye(Nim * Nim) *  ones(Nim * Nim,1);

CovIm = 1./CovIm;
CovIm = CovIm * 10^10;
CovIm = sqrt(CovIm);
figure;imshow(reshape(CovIm,Nim,Nim),[]);

CovIm1 = structGeo.W' * Cov * Cov * structGeo.W * ...
    opEye(Nim * Nim) *  ones(Nim * Nim,1);
CovIm2 = structGeo.W' * structGeo.W * ...
    opEye(Nim * Nim) *  ones(Nim * Nim,1);

CovIm1 = CovIm1 ./CovIm2;
CovIm1 = sqrt(CovIm1);
figure;imshow(reshape(CovIm1,Nim,Nim),[]);

%% Estimate inv(W' * Cov * W) ===> Non-uniform spatial parameter
tic;
N_Coarse_grid = 64;
structGeoC.vol_geom = astra_create_vol_geom(...
         N_Coarse_grid, N_Coarse_grid,...
        -structImg.fFOV/2, structImg.fFOV/2,...
        -structImg.fFOV/2, structImg.fFOV/2);
structGeoC.proj_geom = astra_create_proj_geom('parallel', ...
        structProj.fDetWidth, structProj.nChannel, ...
        linspace(0,2*pi,structProj.nAngle));
structGeoC.proj_id = astra_create_projector('linear', ...
    structGeoC.proj_geom, structGeoC.vol_geom);
structGeoC.W = opTomo('linear', structGeoC.proj_geom, structGeoC.vol_geom); 
matrix_id = astra_mex_projector('matrix', structGeoC.proj_id);
W = astra_mex_matrix('get', matrix_id);

% Circle mask for image
N = N_Coarse_grid;
NTotal = 0;
I_MaskMatrix = zeros(N * N,1);
for m = 1:N
    for n = 1:N
        if ((m-(N/2+0.5))^2 + (n-(N/2+0.5))^2<(N/2)^2)
            I_Mask(m,n) = 1;
            NTotal = NTotal + 1;
            I_MaskMatrix((m-1)*N+n,NTotal) = 1;
        end
    end
end

I_MaskMatrix = sparse(I_MaskMatrix);

Cov1 = Photon_rnd';
Cov1 = opDiag(sqrt(Cov1(:)));

% CovIm = full(I_MaskMatrix' * W' * Cov1' * Cov1 * W * I_MaskMatrix);
% DD = full(diag(inv(CovIm))) * 10^10;
% 
% NTotal = 0;
% for m = 1:N
%     for n = 1:N
%         if ((m-(N/2+0.5))^2 + (n-(N/2+0.5))^2<(N/2)^2)
%             NTotal = NTotal + 1;
%             DDIm(m,n) = DD(NTotal);
%         end
%     end
% end

CovIm_New = full( W' * Cov1' * Cov1 * W );
DD = full(diag(inv(CovIm_New))) * 10^10;
DDIm = reshape(DD,[N_Coarse_grid,N_Coarse_grid]);
DDIm = DDIm';
CovIm_New = reshape(diag(CovIm_New),[N_Coarse_grid,N_Coarse_grid]);
CovIm_New = CovIm_New';
figure;imshow(DDIm,[]);
figure;imshow(CovIm_New,[]);
DDIm = imresize(DDIm,[structImg.nPixel,structImg.nPixel]);
toc;


%% Sparse inverse approximate ---> Precondition 1 (Seems not good here)
CovIm_New = full( W' * Cov1' * Cov1 * W );
for m = 1:size(CovIm_New,1)
    PreC_SIA(m) = CovIm_New(m,m)/norm(CovIm_New(:,m))^2;
end
PreC_SIA = reshape(PreC_SIA,[N_Coarse_grid,N_Coarse_grid]);
PreC_SIA = imresize(PreC_SIA,[structImg.nPixel,structImg.nPixel]);
PreC_SIA = PreC_SIA';
figure;imshow(PreC_SIA,[]);

%% Jacobian precondition  ----> Precondition 2
CovIm_Jac = full( W' * Cov1' * Cov1 * W );
CovIm_Jac = diag(CovIm_Jac);
CovIm_Jac = sqrt(1./ CovIm_Jac);
CovIm_Jac = reshape(CovIm_Jac,[N_Coarse_grid,N_Coarse_grid]);
CovIm_Jac = imresize(CovIm_Jac,[structImg.nPixel,structImg.nPixel]);
CovIm_Jac = CovIm_Jac';

%% Estimate Image variance by repeated simulation
for k = 1:10
    k
    Photon_rnd = poissrnd(Photon);
    Photon_rnd(Photon_rnd < 5) = 5;
    sino_noisy = log(Photon_In./Photon_rnd);
    theta = [0:360/structProj.nAngle: 360-360/structProj.nAngle];
    Image_FBP = iradon(sino_noisy',theta,'spline','Ram-Lak',structProj.nChannel);
    % Image_FBP = imresize(Image_FBP * structProj.nChannel / structImg.fFOV,[structImg.nPixel,structImg.nPixel]);
    I(k,:) = Image_FBP(:);
end

I = reshape(var(I),[structProj.nChannel,structProj.nChannel]);
I = imresize(I,[structImg.nPixel,structImg.nPixel]);

I1 = imresize(I,N_Coarse_grid/structImg.nPixel);
figure;imshow(I1,[]);

%% FBP reconstruction
theta = [0:360/structProj.nAngle: 360-360/structProj.nAngle];
Image_FBP = iradon(sino_noisy',theta,'spline','Ram-Lak',structProj.nChannel);
Image_FBP = imresize(Image_FBP * structProj.nChannel / structImg.fFOV,[structImg.nPixel,structImg.nPixel]);

% Image_FBP(241:270,371:400) = 1000;
% Image_FBP(441:470,241:270) = 1000;
% Image_FBP(291:320,111:140) = 1000;

figure;imshow(Image_FBP,[MuW - MuW * 0.1 MuW + MuW * 0.1]);colorbar;

% For Non-local filter
theta = [0:360/structProj.nAngle: 360-360/structProj.nAngle];
Image_FBP_Filter = iradon(sino_noisy',theta,'spline','Hann',0.2,structProj.nChannel);
Image_FBP_Filter = imresize(Image_FBP_Filter * structProj.nChannel / structImg.fFOV,[structImg.nPixel,structImg.nPixel]);


% Region1 = Image_FBP(241:270,371:400);
% mean(Region1(:)),var(Region1(:))
% Region2 = Image_FBP(441:470,241:270);
% mean(Region2(:)),var(Region2(:))
% Region3 = Image_FBP(291:320,111:140);
% mean(Region3(:)),var(Region3(:))

Region1 = Image_FBP(241/2:270/2,371/2:400/2);
mean(Region1(:)),var(Region1(:))
Region2 = Image_FBP(441/2:470/2,241/2:270/2);
mean(Region2(:)),var(Region2(:))
Region3 = Image_FBP(291/2:320/2,111/2:140/2);
mean(Region3(:)),var(Region3(:))



tao = 60;
NIter = 100;
tol = 1;
Xwarm = Image_FBP(:);
y = sino_noisy(:);

% %L2_NoWeight
% A  = structGeo.W; % Cov * structGeo.W
% R = opEye(Nim * Nim); 
% y0 = y;% Cov * y
% yR = zeros(Nim * Nim,1);
% tao = [5:5:50];
% NIter = 50;
% for k = 1:length(tao)
%     [X,Obj,XTotal] = quasinewton([A;tao(k) * R],...
%         [y0;yR],Xwarm,NIter,tol);
%     X_Image = reshape(X,structImg.nPixel,structImg.nPixel);
%     save(['Image_L2_NoWeight_',num2str(tao(k)),'.mat'],'X_Image','Obj');
%     Error(k) = norm(A * X - y0);
%     Regu(k) = norm(R * X);
% end
% figure;plot(Error,Regu);

% %L2_Weight
% A  =  Cov * structGeo.W;
% R =  opEye(Nim * Nim); 
% y0 = Cov * y;
% yR = zeros(Nim * Nim,1);
% tao = [100:100:500];% 10.^(1:5);
% NIter = 50;
% for k = 1:length(tao)
%     [X,Obj,XTotal] = quasinewton([A;tao(k) * R],...
%         [y0;yR],Xwarm,NIter,tol);
%     X_Image = reshape(X,structImg.nPixel,structImg.nPixel);
%     save(['Image_L2_Weight_',num2str(tao(k)),'.mat'],'X_Image','Obj');
%     Error(k) = norm(A * X - y0);
%     Regu(k) = norm(R * X);
% end
% figure;plot(Error,Regu);

% %L2_Smooth_NoWeight
% A  = structGeo.W; % Cov * structGeo.W
% R = FirstDerivativeM(Nim);
% y0 = y;% Cov * y
% yR = zeros(size(R,1),1);
% tao = [25:5:50];
% NIter = 50;
% for k = 1:length(tao)
%     [X,Obj,XTotal] = quasinewton([A;tao(k) * R],...
%         [y0;yR],Xwarm,NIter,tol);
%     X_Image = reshape(X,structImg.nPixel,structImg.nPixel);
%     save(['Image_L2D_NoWeight_',num2str(tao(k)),'.mat'],'X_Image','Obj');
%     Error(k) = norm(A * X - y0);
%     Regu(k) = norm(R * X);
% end
% figure;plot(Error,Regu);

% %L2_Smooth_Weight
% A  = Cov * structGeo.W;
% R = FirstDerivativeM(Nim);
% y0 =  Cov * y;
% yR = zeros(size(R,1),1);
% tao = 2000; % [600:200:1600]; % 10.^(0:5);
% NIter = 200;
% for k = 1:length(tao)
%     [X,Obj,XTotal] = quasinewton([A;tao(k) * R],...
%         [y0;yR],Xwarm,NIter,tol);
%     X_Image = reshape(X,structImg.nPixel,structImg.nPixel);
%     save(['Image_L2D_Weight_',num2str(tao(k)),'.mat'],'X_Image','Obj');
%     Error(k) = norm(A * X - y0);
%     Regu(k) = norm(R * X);
% end
% figure;plot(Error,Regu);

% L2_Smooth_Weight_ReguW
% A  = Cov * structGeo.W;
% R = FirstDerivativeMW(Nim,CovIm);
% y0 =  Cov * y;
% yR = zeros(size(R,1),1);
% tao = 300;% [100:100:300]; % 10.^(0:5);
% NIter = 200; % 50;
% for k = 1:length(tao)
%     [X,Obj,XTotal] = quasinewton([A;tao(k) * R],...
%         [y0;yR],Xwarm,NIter,tol);
%     X_Image = reshape(X,structImg.nPixel,structImg.nPixel);
%     save(['Image_L2D_Weight_ReguW_',num2str(tao(k)),'.mat'],'X_Image','Obj');
%     Error(k) = norm(A * X - y0);
%     Regu(k) = norm(R * X);
% end
% figure;plot(Error,Regu);


% %L2_Smooth_Weight_PreC1
% PreW = opDiag(Xwarm);
% A  = Cov * structGeo.W * PreW;
% R = FirstDerivativeM(Nim) * PreW;
% y0 =  Cov * y;
% yR = zeros(size(R,1),1);
% tao = 2000; % [600:200:1600]; % 10.^(0:5);
% XwarmP = Xwarm ./ Xwarm ;
% NIter = 50;
% for k = 1:length(tao)
%     [X,Obj,XTotal] = quasinewton([A;tao(k) * R],...
%         [y0;yR],XwarmP,NIter,tol);
%     X = PreW * X;
%     X_Image = reshape(X,structImg.nPixel,structImg.nPixel);
%     for m = 1:NIter
%         XTotal(:,m) = PreW * XTotal(:,m);
%     end
%     save(['Image_L2D_Weight_PreC1_',num2str(tao(k)),'.mat'],'X_Image','Obj');
%     Error(k) = norm(A * X - y0);
%     Regu(k) = norm(R * X);
% end
% figure;plot(Error,Regu);


% %L2_Smooth_Weight_PreC2
% XwarmP = Xwarm;
% XwarmP(Xwarm < 0) = 0;
% XwarmP = sqrt(XwarmP);
% PreW = opDiag(XwarmP);
% A  = Cov * structGeo.W * PreW;
% R = FirstDerivativeM(Nim) * PreW;
% y0 =  Cov * y;
% yR = zeros(size(R,1),1);
% tao = 2000; % [600:200:1600]; % 10.^(0:5);
% % XwarmP = Xwarm ./ Xwarm ;
% NIter = 50;
% for k = 1:length(tao)
%     [X,Obj,XTotal] = quasinewton([A;tao(k) * R],...
%         [y0;yR],XwarmP,NIter,tol);
%     X = PreW * X;
%     X_Image = reshape(X,structImg.nPixel,structImg.nPixel);
%     for m = 1:NIter
%         XTotal(:,m) = PreW * XTotal(:,m);
%     end
%     save(['Image_L2D_Weight_PreC2_',num2str(tao(k)),'.mat'],'X_Image','Obj');
%     Error(k) = norm(A * X - y0);
%     Regu(k) = norm(R * X);
% end
% figure;plot(Error,Regu);

% %L2_Smooth_Weight_PreC3
% PreW = opDiag(CovIm);
% A  = Cov * structGeo.W * PreW;
% R = FirstDerivativeM(Nim) * PreW;
% y0 =  Cov * y;
% yR = zeros(size(R,1),1);
% tao = 2000; % [600:200:1600]; % 10.^(0:5);
% XwarmP = Xwarm ./ CovIm ;
% NIter = 50;
% for k = 1:length(tao)
%     [X,Obj,XTotal] = quasinewton([A;tao(k) * R],...
%         [y0;yR],XwarmP,NIter,tol);
%     X = PreW * X;
%     X_Image = reshape(X,structImg.nPixel,structImg.nPixel);
%     for m = 1:NIter
%         XTotal(:,m) = PreW * XTotal(:,m);
%     end
%     save(['Image_L2D_Weight_PreC3_',num2str(tao(k)),'.mat'],'X_Image','Obj');
%     Error(k) = norm(A * X - y0);
%     Regu(k) = norm(R * X);
% end
% figure;plot(Error,Regu);
% 
% %L2_Smooth_Weight_PreC5
PreW = opDiag(sqrt(PreC_SIA(:)));
A  = Cov * structGeo.W * PreW;
R = FirstDerivativeM(Nim) * PreW;
y0 =  Cov * y;
yR = zeros(size(R,1),1);
tao = 2000; % [600:200:1600]; % 10.^(0:5);
XwarmP = Xwarm ./ sqrt(PreC_SIA(:)) ;
NIter = 50;
for k = 1:length(tao)
    [X,Obj,XTotal] = quasinewton([A;tao(k) * R],...
        [y0;yR],XwarmP,NIter,tol);
    X = PreW * X;
    X_Image = reshape(X,structImg.nPixel,structImg.nPixel);
    for m = 1:NIter
        XTotal(:,m) = PreW * XTotal(:,m);
    end
    save(['Image_L2D_Weight_PreC5_',num2str(tao(k)),'.mat'],'X_Image','Obj');
    Error(k) = norm(A * X - y0);
    Regu(k) = norm(R * X);
end
% figure;plot(Error,Regu);

% %L2_Smooth_Weight_PreC4
% PreW = opDiag((CovIm).^2);
% A  = Cov * structGeo.W * PreW;
% R = FirstDerivativeM(Nim) * PreW;
% y0 =  Cov * y;
% yR = zeros(size(R,1),1);
% tao = 2000; % [600:200:1600]; % 10.^(0:5);
% XwarmP = Xwarm ./ (CovIm).^2;
% NIter = 50;
% for k = 1:length(tao)
%     [X,Obj,XTotal] = quasinewton([A;tao(k) * R],...
%         [y0;yR],XwarmP,NIter,tol);
%     X = PreW * X;
%     X_Image = reshape(X,structImg.nPixel,structImg.nPixel);
%     for m = 1:NIter
%         XTotal(:,m) = PreW * XTotal(:,m);
%     end
%     save(['Image_L2D_Weight_PreC4_',num2str(tao(k)),'.mat'],'X_Image','Obj');
%     Error(k) = norm(A * X - y0);
%     Regu(k) = norm(R * X);
% end
% figure;plot(Error,Regu);

% %L2_NLM_NoWeight
% A  = structGeo.W; % Cov * structGeo.W;
% sigma = 10^(-6); % 0.01; %20;
% NPatch = 2;
% NNeighbour = 2;
% R = NLM_Matrix(Nim,Image_FBP_Filter,sigma,NPatch,NNeighbour);
% y0 =  y;% Cov * y;
% yR = zeros(size(R,1),1);
% tao = 300; % [600:200:1600]; % 10.^(0:5);
% NIter = 20;
% 
% for k = 1:length(tao)
%     [X,Obj,XTotal] = quasinewton([A;tao(k) * R],...
%         [y0;yR],Xwarm,NIter,tol);
%     X_Image = reshape(X,structImg.nPixel,structImg.nPixel);
%     save(['Image_L2NLM_NoWeight_',num2str(tao(k)),'.mat'],'X_Image','Obj');
%     Error(k) = norm(A * X - y0);
%     Regu(k) = norm(R * X);
% end
% figure;plot(Error,Regu);


% %L2_NLM_Weight
% PreW = opDiag(CovIm);
% A  = Cov * structGeo.W * PreW;
% sigma = 10^(-6); % very sensitive!
% NPatch = 2;
% NNeighbour = 2;  % Larger search neighbourhood means larger matrix!
% R = NLM_Matrix(Nim,Image_FBP_Filter,sigma,NPatch,NNeighbour) * PreW;
% % Choose the proper image for pre-computed weights!
% sum(isnan(R))
% y0 =  Cov * y;
% yR = zeros(size(R,1),1);
% tao = 6000; % [600:200:1600]; % 10.^(0:5);
% XwarmP = Xwarm ./ CovIm ;
% NIter = 50;
% 
% for k = 1:length(tao)
%     [X,Obj,XTotal] = quasinewton([A;tao(k) * R],...
%         [y0;yR],XwarmP,NIter,tol);
%     X = PreW * X;
%     X_Image = reshape(X,structImg.nPixel,structImg.nPixel);
%     for m = 1:NIter
%         XTotal(:,m) = PreW * XTotal(:,m);
%     end
%     save(['Image_L2NLM_Weight_FBPFilter_tao',num2str(tao(k)),'sigma_',num2str(sigma),'.mat'],'X_Image','Obj');
%     Error(k) = norm(A * X - y0);
%     Regu(k) = norm(R * X);
% end



% 
A = structGeo.W' * Cov * Cov * structGeo.W + tao * tao * FirstDerivativeM(Nim)' * FirstDerivativeM(Nim);
b = structGeo.W' * Cov * Cov * y;
PreC = CovIm_Jac(:);% ones(65536,1);% PreC_SIA(:);
[x,xtotal] = pcg(A * opDiag(PreC) ,b,10^(-10),NIter,[],[],Xwarm./PreC);
for k = 1:NIter
    Obj_pcg(k) = norm(Cov * (structGeo.W * (xtotal(:,k).* PreC) - y))^2 + ...
        tao * tao * norm(FirstDerivativeM(Nim) * (xtotal(:,k).* PreC))^2;
end

