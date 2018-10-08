% Data: A, b, tao, NImage, gamma

function F = TV_denoising_Obj_OS(X,Data)
    A = Data.A  * opDiag(Data.PreCon);
    b = Data.b;
    NImage = Data.NImage;
    Nb = Data.Nb;
    gamma = Data.gamma;
    
    Diff = A * X - b;
    F = Diff' * Diff;
    
    % 
    for k = 1:Nb
%         X_Square = reshape(Data.PreCon(:) .* X,Nb, NImage,NImage); % Precondition
%         X_Square = squeeze(X_Square(k,:,:));
        
        X_Square = reshape(Data.PreCon(:) .* X, NImage,NImage,Nb); % Precondition
        X_Square = squeeze(X_Square(:,:,k));
        
        
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
        F = Data.tao(k) * sum(G(:)) + F;
    end  
%     % ISOTROPIC
%     Temp = X_extend(3:end,2:end-1) - X_extend(2:end-1,2:end-1);
%     G = G + Temp.^2;
%     Temp = X_extend(2:end-1,3:end) - X_extend(2:end-1,2:end-1);
%     G = G + Temp.^2 + gamma;
%     F = sum(sqrt(G(:))) + F;
end

