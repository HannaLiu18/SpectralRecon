clc;clear;
MuB = 0.01855 *  1.920; % Cortical Bone
MuW = 0.01707;

% Define regions and lines for evaluating image qualities
Line1 = [412,210:293];
Line2 = [145,280:335];
Line3 = [263,321:370];

% L2_NoWeight
tao = [5:5:50];
for k = 1:length(tao)
    load(['Image_L2_NoWeight_',num2str(tao(k)),'.mat']);
    figure;imshow(X_Image,[MuW - MuW * 0.1 MuW + MuW * 0.1]);
    Norm(k) = norm(X_Image(:));
end

% L2_Weight  --- Regularization is too small for weighted condition...
tao = [100:100:500]; % 10.^(1:5); % Should be some value between 2~3
for k = 1:length(tao)
    load(['Image_L2_Weight_',num2str(tao(k)),'.mat']);
    figure;imshow(X_Image,[MuW - MuW * 0.1 MuW + MuW * 0.1]);
    Norm(k) = norm(X_Image(:));
end

% L2D_NoWeight
tao = [5:5:50];
for k = 1:length(tao)
    load(['Image_L2D_NoWeight_',num2str(tao(k)),'.mat']);
    figure;imshow(X_Image,[MuW - MuW * 0.1 MuW + MuW * 0.1]);
    Norm(k) = norm(X_Image(:));
end

% L2D_Weight
tao = 2000; % [600:200:1600];
for k = 1:length(tao)
    load(['Image_L2D_Weight_',num2str(tao(k)),'.mat']);
    figure;imshow(X_Image,[MuW - MuW * 0.1 MuW + MuW * 0.1]);
    Norm(k) = norm(X_Image(:));
end
X_Image1 = X_Image;

% L2D_Weight_ImgRegu
tao = 600;% [100:100:300];% 600;
for k = 1:length(tao)
    load(['Image_L2D_Weight_ReguW_',num2str(tao(k)),'.mat']);
    figure;imshow(X_Image,[MuW - MuW * 0.1 MuW + MuW * 0.1]);
    Norm(k) = norm(X_Image(:));
end
X_Image2 = X_Image;
