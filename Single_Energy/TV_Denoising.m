% TV denoising tests
clc;clear;close all;
I = double(imread('lena512.bmp'));
% I = zeros(512,512);
% I(200:300,300:380) = 50;
[M,N] = size(I);
I_noisy = I + 50 * (rand(M,N) - 0.5);
figure;imshow(I,[]);
figure;imshow(I_noisy,[]);

%% Approximate TV--LBFGS
Data.A = opEye(M * N);
Data.b = I_noisy(:);
Data.PreConM = opEye(M * N);
Data.NImage = M;
Data.gamma = 0.01;
Data.tao = 0.16/2;

NIter = 100;
NMemory = 100;
Func_Obj = @TV_denoising_Obj;
Func_Grad = @TV_denoising_Grad;

tic
[X_LBFGS,Obj_BFGS,XTotal_LBFGS] = Quasi_Newton_LBFGS...
    (I_noisy(:), NIter, Func_Obj, Func_Grad, Data);
toc
X_LBFGS = reshape(X_LBFGS,M,N);
figure;imshow(X_LBFGS,[]);

pause;

%% Approximate TV
Ident = opEye(M * N);
X0_Square = I_noisy;
X0 = X0_Square(:);
gamma = 0.01;
tao = 0.16/2;
NIter = 100;
tol = 1;
tic
[X,Obj_appro,XTotal] = quasinewton_TV(Ident,I_noisy(:),X0,X0_Square,M,gamma,tao,NIter,tol);
toc
X = reshape(X,M,N);
figure;imshow(X,[]);
figure;plot(log(Obj_appro));hold on;
plot(log(Obj_BFGS),'*--');



% %% ADMM algorithm
% rou = 1;
% X = zeros(M,N);
% Z = zeros(M,N,2);
% U = zeros(M,N);
% Lamda = 0.16;
% f = I_noisy;
% NIter = 100;
% for k = 1:NIter
%     
% end
% 
% %% Nesterov algorithm (FISTA) --> Solving TV is much harder than simple L1!

%% Primal-Dual Algorithm (ANISOTROPIC regu using now...)
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
% Result: It seems that the regularization parameter doesn't influence the
% selection of sigma here?!
% % Algorithm 1~2
% sigma = 0.1; % This significantly influence the convergence rate!
% tau = .9/(L*sigma); 

% AHMOD
tau = 20;
sigma = 4 / (L * tau);

theta = 1;
Lamda = 0.16;% 0.08;
gamma = 0.7 * Lamda;
% initial valius
f = I_noisy; % zeros(M,N); % I_noisy;
g = grad(I_noisy) * 0;
f1 = f;
NIter = 500;

for k = 1:NIter
    fold = f;
%    g = ProxFS(g + sigma * grad(f1), sigma);  % An alternative way
    Temp = g + sigma * grad(f1);  % This seems to converge faster
    g = Temp./max(ones(size(Temp)),abs(Temp));
    f = ProxG(f + tau * div(g(:,:,1),g(:,:,2)), tau, Lamda);
    
%     % Algorithm 1~2
%     theta = 1;% 1/sqrt(1+2*gamma*tau);
%     tau = tau * theta;
%     sigma = sigma / theta;
%     tau_n(k) = tau;
%     f1 = f + theta * (f-fold);
    
    % AHMOD
    theta =  1;% 1/sqrt(1+2*gamma*tau);
    tau = tau * theta;
    sigma = sigma / theta;
    f1 = f;  
    
    Error(k) = norm(f1(:) - I_noisy(:))^2;
    Norm(k) = F(grad(f1));
    
    Obj(k) = 0.5 * Lamda * Error(k) + Norm(k); % denoising
end
figure;imshow(f1,[]);

figure;
loglog([1:k],Obj(1:k));
hold on;
loglog(Obj_appro);

% % loglog([1:k],40000000./[1:k].^2);

%%
F(grad(f1)),F(grad(X))


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
