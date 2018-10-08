%% Set up energy discreminating detector properties
%  structDet: 
%  Input:
%       fEnergyLow:   lower bound of energy detected
%       fEnergyHigh:  upper bound of energy detected   
%       fEnerghWidth: width of energy bins
%  Output:
%       pfEnergyLow   lower bound of energy detected for all bins
%       pfEnergyHigh  upper bound of energy detected for all bins
%       nEnergyBin    number of energy bins
function structDet = setDetector()
    % Input
    structDet.fEnergyLow = 15; %%% Worst choice...
    structDet.fEnergyHigh = 105; % 105;
    structDet.fEnerghWidth = 10;
    % Output
    structDet.pfEnergyLow = [structDet.fEnergyLow:structDet.fEnerghWidth:...
        structDet.fEnergyHigh - structDet.fEnerghWidth];
    structDet.pfEnergyHigh = [structDet.fEnergyLow + structDet.fEnerghWidth...
        :structDet.fEnerghWidth:structDet.fEnergyHigh];
    structDet.nEnergyBin = length(structDet.pfEnergyLow);
end


