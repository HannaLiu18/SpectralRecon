% Data: A, b, tao, NImage, gamma

function F = TV_denoising_Obj_Circle(X,Data)
    A = Data.A  * opDiag(Data.PreCon);
    b = Data.b;
    tao = Data.tao;
    NImage = Data.NImage;
    gamma = Data.gamma;
    Circle_Label = Data.Circle_Label;
    
    Diff = A * X - b;
    F = Diff' * Diff * tao;
    
    X_Square = zeros(NImage,NImage);
    X_Square(Circle_Label) = Data.PreCon(:) .* X;
    
    X_Square = reshape(X_Square,NImage,NImage); % Precondition
    X_extend = zeros(NImage + 2, NImage + 2);
    X_extend(2:end-1, 2:end-1) = X_Square;
    X_extend(1,:) = X_extend(2,:);
    X_extend(end,:) = X_extend(end-1,:);
    X_extend(:,1) = X_extend(:,2);
    X_extend(:,end) = X_extend(:,end-1);  
    
    G = zeros(NImage, NImage);
    
    % ANISOTROPIC
    Temp = X_extend(3:end,2:end-1) - X_Square;% X_extend(2:end-1,2:end-1);
    G = G + sqrt(Temp.^2 + gamma); % abs(Temp); % sqrt(Temp.^2 + gamma);
    Temp = X_extend(2:end-1,3:end) - X_Square;% X_extend(2:end-1,2:end-1);
    G = G + sqrt(Temp.^2 + gamma); % abs(Temp); % sqrt(Temp.^2 + gamma);
    F = sum(G(:)) + F;

%     % ISOTROPIC
%     Temp = X_extend(3:end,2:end-1) - X_extend(2:end-1,2:end-1);
%     G = G + Temp.^2;
%     Temp = X_extend(2:end-1,3:end) - X_extend(2:end-1,2:end-1);
%     G = G + Temp.^2 + gamma;
%     F = sum(sqrt(G(:))) + F;
end