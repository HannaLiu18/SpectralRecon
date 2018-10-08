% Recon Experiment 1
% Date: 02/12/2018
clc;clear;load('ATotal.mat');
[structGeo,structProj] = setGeometryFwd(structImg);

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


% % 4-- build-in iradon
% y1 = imresize(y,[90,128]);
% sino = radon(structImg.pfImgLabel,linspace2(0,360-4,90));
% sino(29:end-29,:) = y1';
% Vbuild = iradon(sino,linspace2(0,360-4,90));
% figure;imagesc(Vbuild);colorbar;
% 
% bound = floor((size(y,2)*sqrt(2) - size(y,2)-2)/2);
% y1(:,bound+1:bound+672) = y;
% y1(:,bound+673:bound + 673+bound) = 0;
% Image2 = iradon(y1',[0:4:360-4]);
% figure;imagesc(Image2);colorbar;


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

