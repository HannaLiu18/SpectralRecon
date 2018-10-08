% Main Run for material decomposition
% Author:   Hanna Liu
% Date:     01/12/2017
clc;clear;close all;
rng(100); % generate repeatable simulation data

%% Data structure
% structDet: energy discreminating detector properties
% structMas: mass attenuation coefficients
% structSrc: X-ray tube spectral
% structImg: true spectral image 
% structContrast: basic properties of contrast agent
% structGeo: geometry for forward and backward projections using astra

%% Data generation and model setup
structMas = setMasCurve();  % MAC curves for generating simulated data
structSrc = setSource();    % X-ray source setup
structDet = setDetector();  % Photon discreminating detector setup
[structImg,structContrast] = setTrueImg(structMas);
                            % Generate true image
[structGeo,structProj] = setGeometryFwd(structImg);
                            % Geometry setup for foward and backward proj
structProj = computeProjME(structImg,structProj,structSrc,...
                           structGeo,structMas,structDet);
                            % Generate simulated data
structProj.pfPhotonrnd(336/2,1,9),pause;
%% Material decomposition 
structEBas = setEnergyBasis();
                            % Generate MAC curves for decomposition                   
tic
ATotal = zeros(structProj.nChannel,structProj.nAngle,size(structEBas,2));

A = zeros(1,size(structEBas,2));% * length(nChannel) * length(nAngle));
options = optimset('TolFun',1e-3);
for m = 1: structProj.nAngle % for every projection
    m,tic
    for n =  1: structProj.nChannel % for every channel
        Astart = A;
        A = fminsearch(@(A) decomplikelihood(A,...
            structEBas,structDet,structSrc,structProj,n,m),Astart,options);
        ATotal(n,m,:) = A;
    end % end channel
    toc
end % end projection

figure; imagesc(abs(ATotal(:,:,3)-true'*0.1),[0 0.05]);colorbar;
xlabel('Angle');ylabel('Channel');axis square;
toc

%% Material decomposition (Weighted least-square)
structEBas = setEnergyBasis();
Nb = length(structEBas);
Ne = structDet.nEnergyBin;
MatrixF = zeros(Ne,Nb);
for m = 1:Ne
    for n = 1:Nb
        nSrcLoc = find(structEBas(n).pfEnergy >= structDet.pfEnergyLow(m)&...
                       structEBas(n).pfEnergy < structDet.pfEnergyHigh(m));
        MatrixF(m,n) = mean(structEBas(n).pfBase(nSrcLoc));
    end
end

tic
Ne_Start = 4;
A = zeros(1,size(structEBas,2));% * length(nChannel) * length(nAngle));
ATotal_LS = zeros(structProj.nChannel,structProj.nAngle,size(structEBas,2));
for m = 1: structProj.nAngle % for every projection
   %  m,tic
    for n =  1: structProj.nChannel % for every channel
        temp = max(squeeze(structProj.pfPhotonrnd(n,m,:)),ones(Ne,1));
        Cov = diag(temp);
        A = (MatrixF(Ne_Start:end,:)' * Cov(Ne_Start:end,Ne_Start:end) * MatrixF(Ne_Start:end,:))^(-1) * ...
            MatrixF(Ne_Start:end,:)' * Cov(Ne_Start:end,Ne_Start:end) * ...
            squeeze(structProj.pfLACrnd(n,m,Ne_Start:end));
        ATotal_LS(n,m,:) = A;
    end % end channel
   %  toc
end % end projection
toc
figure; imagesc(ATotal_LS(:,:,3),[0 0.10]);colorbar;
xlabel('Angle');ylabel('Channel');axis square;

pause;

%%%%%%%%%%%%%%%%%% Result validation (Temporary)
R = zeros(structProj.nChannel,structProj.nAngle,structDet.nEnergyBin);
for n = 1: structProj.nChannel
    for m = 1:structProj.nAngle
    A = squeeze(ATotal(n,m,1:3));
    R(n,m,:) = computePhotonFromBase(A,structEBas,structDet,structSrc,structProj);
    end
end

figure;
subplot(3,1,1);plot(ATotal(:,1,1));
subplot(3,1,2);plot(ATotal(:,1,2));
subplot(3,1,3);plot(ATotal(:,1,3));

figure;
plot(ATotal(:,1,2));
hold on;plot(ATotal(:,1,3)*100);

% reconstruct
Image = ReconDecomp(structGeo,ATotal(:,:,3)');
imshow(Image, []);colorbar;

mean(Image(structImg.pfImgLabel == 5))
mean(Image(structImg.pfImgLabel == 6))
mean(Image(structImg.pfImgLabel == 7))

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
[true_id, true] = astra_create_sino(...
            Contrast, structGeo.proj_id);
plot(ATotal(:,1,3));hold on;plot(true(1,:)*0.1);
        


%%% Iterative reconstruction (Comparison)
y = ATotal(:,:,3)';
% y = structProj.pfLac(:,:,91) + 0.2 * rand(672,90);y = y';
% 1-- Min-square 
ArtRecon1 = (structGeo.W' * structGeo.W)^(-1)*structGeo.W' * y(:);
Obj1 = structGeo.W * ArtRecon1 - y(:);
Obj1 = Obj1' * Obj1
ArtRecon1 = reshape(ArtRecon1,[128,128]);
figure;imagesc(ArtRecon1);colorbar;
% 2-- pseudo-inverse
ArtRecon2 = structGeo.W \ y(:);
Obj2 = structGeo.W * ArtRecon2 - y(:);
Obj2 = Obj2' * Obj2
ArtRecon2 = reshape(ArtRecon2,[128,128]);
figure;imagesc(ArtRecon2);colorbar;
% 3-- astra-fbp
Image = ReconDecomp(structGeo,y);
Obj3 = structGeo.W * Image(:) - y(:);
Obj3 = Obj3' * Obj3
[Obj3_id, Obj3new] = astra_create_sino(Image, structGeo.proj_id);
% Obj3new = Obj3new(:)' * Obj3new(:)
figure;imagesc(Image);colorbar;


% 4-- build-in iradon
y1 = imresize(y,[90,128]);
sino = radon(structImg.pfImgLabel,linspace2(0,360-4,90));
sino(29:end-29,:) = y1';
Vbuild = iradon(sino,linspace2(0,360-4,90));
figure;imagesc(Vbuild);colorbar;

bound = floor((size(y,2)*sqrt(2) - size(y,2)-2)/2);
y1(:,bound+1:bound+672) = y;
y1(:,bound+673:bound + 673+bound) = 0;
Image2 = iradon(y1',[0:4:360-4]);
figure;imagesc(Image2);colorbar;


% 5-- lsqr (Comparing to 1 and 2)--> equavalent to 2 arithmetic
% RELRES: relative residual NORM(B-A*X)/NORM(B)
% RESVEC: NORM(B-A*X0) at every iteration
[ArtRecon5,FLAG,RELRES,ITER,RESVEC] = lsqr(structGeo.W,y(:));
ArtRecon5 = reshape(ArtRecon5,[128,128]);
figure;imagesc(ArtRecon5);colorbar;

% 6-- lsqr with regularization
tao = 60;% [10:10:100];
yTilt = [y(:);zeros(128*128,1)];
for k = 1:length(tao)
    k
    WTilt = [structGeo.W;tao(k) * opEye(128*128)];
    [ArtRecon6,FLAG,RELRES,ITER,RESVEC] = lsqr(WTilt,yTilt);
    xnorm(k) = norm(ArtRecon6);
    resnorm(k) = norm(structGeo.W * ArtRecon6(:) - y(:));
    ArtRecon6 = reshape(ArtRecon6,[128,128]);
    figure;imagesc(ArtRecon6);colorbar;  
end
% figure;plot(xnorm,resnorm);
mean(ArtRecon6(structImg.pfImgLabel == 5))
mean(ArtRecon2(structImg.pfImgLabel == 5))
mean(Contrast(structImg.pfImgLabel == 5))

