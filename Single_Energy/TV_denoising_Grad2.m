% Data: NImage,gamma, A, b, tao
function G = TV_denoising_Grad2(X,Data)
    NImage = Data.NImage;
    gamma = Data.gamma;
    A = Data.A;
    b = Data.b;
    tao = Data.tao;
    
    X_Square = reshape(Data.PreCon(:) .* X,NImage,NImage); % Precondition
    X_extend = zeros(NImage + 2, NImage + 2);
    X_extend(2:end-1, 2:end-1) = X_Square;
    X_extend(1,:) = X_extend(2,:);
    X_extend(end,:) = X_extend(end-1,:);
    X_extend(:,1) = X_extend(:,2);
    X_extend(:,end) = X_extend(:,end-1);  
    
    G = zeros(NImage, NImage);
    
    % ISOTROPIC regu is smaller than ANISOTROPIC, therefore requires
    % smaller lamda for data fidity term.
    
%     % ISOTROPIC
%     Temp_One = X_extend(2:end-1,1:end-2) - X_extend(3:end,1:end-2);
%     Temp_Two = X_extend(2:end-1,2:end-1) - X_extend(2:end-1,1:end-2);    
% 
%     Temp_Three = X_extend(2:end-1,2:end-1) - X_extend(1:end-2,2:end-1);    
%     Temp_Four = X_extend(1:end-2,2:end-1) - X_extend(1:end-2,3:end);   
%     
%     Temp_Five = X_extend(2:end-1,2:end-1) - X_extend(3:end,2:end-1);  
%     Temp_Six = X_extend(2:end-1,2:end-1) - X_extend(2:end-1,3:end);  
%     
%     G = G + Temp_Two./sqrt(Temp_One.^2 + Temp_Two.^2 + gamma);
%     G = G + Temp_Three./sqrt(Temp_Three.^2 + Temp_Four.^2 + gamma);
%     G = G + (Temp_Five + Temp_Six)./sqrt(Temp_Five.^2 + Temp_Six.^2 + gamma);

%     % ANISOTROPIC (should be just 0????--before 20180903)
%     Temp = X_Square - X_extend(1:end-2,2:end-1);
%     G = G + gamma./(Temp.^2 + gamma).^(3/2);
%     Temp = X_Square - X_extend(3:end,2:end-1);
%     G = G + gamma./(Temp.^2 + gamma).^(3/2);
%     Temp = X_Square - X_extend(2:end-1,1:end-2);
%     G = G + gamma./(Temp.^2 + gamma).^(3/2);
%     Temp = X_Square - X_extend(2:end-1,3:end);
%     G = G + gamma./(Temp.^2 + gamma).^(3/2);
% 
%     % ANISOTROPIC ()
%     Temp = X_Square - X_extend(1:end-2,2:end-1);
%     G = G + Temp.^2./(Temp.^2 + gamma).^(3/2);
%     Temp = X_Square - X_extend(3:end,2:end-1);
%     G = G + Temp.^2./(Temp.^2 + gamma).^(3/2);
%     Temp = X_Square - X_extend(2:end-1,1:end-2);
%     G = G + Temp.^2./(Temp.^2 + gamma).^(3/2);
%     Temp = X_Square - X_extend(2:end-1,3:end);
%     G = G + Temp.^2./(Temp.^2 + gamma).^(3/2);
    
    G = G(:);
    
    % Add error term
    G_Diff = 2 * A' * (A * ones(NImage * NImage,1));
    G = tao * G_Diff + G * 4; %(before 20180903)
    % G = tao * G_Diff + G * 2;
end
