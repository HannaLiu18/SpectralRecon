%% Load MAS curves from .mat files
%  structMas:
%  Input:
%       pfEnergy:   vector of energies for Mas data points 
%       nBase:      number of different materials
%       cfFilename: files saving Mas data,note the last two bases should be
%                   contrast and water(used in function setTrueImg())
%  Output:
%       pfMas:      mas curves, size of length(pfEnergy) * nBase
function structMas = setMasCurve()
    % Input
    structMas.pfEnergy = [15:1:105];
    structMas.nBase = 6;
    structMas.cfFilename = {'SoftTissue.mat','Heart.mat','CorticalBone.mat',...
       'YellowMarrow.mat','GdContrast.mat','Water.mat'};
    % Output
    for k = 1:structMas.nBase
        load(structMas.cfFilename{k});
        structMas.pfMas(:,k) = interp1(log(x(:,1)),log(x(:,2)),log(structMas.pfEnergy));
        structMas.pfMas(:,k) = exp(structMas.pfMas(:,k));
    end
end