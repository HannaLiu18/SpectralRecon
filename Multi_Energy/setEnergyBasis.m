%% Generate MAC curves for decomposition
%  PhotonElectric and Compton for materials without K-edges
%  structEBas:
%  Output:
%       strName,pfEnergy,pfBase
%       pnBaseLoc: matrix of size structDet.nEnergyBin * length(pfEnergy)

function structEBas = setEnergyBasis()
    % Absorption effect
    structEBas(1).strName = 'PhotonElectric';
    structEBas(1).pfEnergy = [10:1:150]; %(keV)
    structEBas(1).pfBase = 1./(structEBas(1).pfEnergy).^3; 
    % Incoherent effect
    structEBas(2).strName = 'Compton';
    structEBas(2).pfEnergy = [10:1:150];
    alpha = structEBas(2).pfEnergy /510.975;
        structEBas(2).pfBase = (1+alpha)./alpha.^2 .*...
        (2 * (1+alpha)./(1+2*alpha) - log(1+2*alpha)./alpha) + ...
        1/2./alpha .* log(1+2*alpha) - (1+3*alpha)./(1+2*alpha).^2;
    % Contrast (Gadolinium)
    load('Gadolinium.mat');
    structEBas(3).strName = 'Gadolinium';
    structEBas(3).pfEnergy = x(:,1)'; 
    structEBas(3).pfBase = x(:,2)'; 


    % Scale base functions, otherwise fminsearch
    % may have problem (not gradient-based method)?
%     structEBas(1).pfBase = mean(structEBas(3).pfBase) ...
%         /mean(structEBas(1).pfBase) *  structEBas(1).pfBase;
%     structEBas(2).pfBase = mean(structEBas(3).pfBase) ...
%         /mean(structEBas(2).pfBase) *  structEBas(2).pfBase;
     structEBas(1).pfBase = structEBas(1).pfBase * 10^6;
     mean(structEBas(1).pfBase)
     mean(structEBas(2).pfBase)
     mean(structEBas(3).pfBase)
end
