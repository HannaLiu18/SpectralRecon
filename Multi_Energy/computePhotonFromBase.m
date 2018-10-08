%% Compute photons from base functions and coefficients
%  Output:
%       pfPhoton: vector length of structDet.nEnergyBin

function pfPhoton = computePhotonFromBase(A,structEBas,structDet,...
                                          structSrc,structProj) 
    pfPhoton = zeros(structDet.nEnergyBin,1);
    % Consider all energy bins
    for k = 1:structDet.nEnergyBin % for every energy bin 
        pfLAC = 0;
        for n = 1:size(structEBas,2) % for every energy basis
            nBaseLoc = (structEBas(n).pfEnergy >=structDet.pfEnergyLow(k)&...
                        structEBas(n).pfEnergy < structDet.pfEnergyHigh(k));
            pfLAC = pfLAC + structEBas(n).pfBase(nBaseLoc) * A(n);
        end % end energy basis
        pfProjminusExp = exp(-pfLAC);
        % find source data points inside the kth bin
        nSrcLoc = (structSrc.pfEnergy >= structDet.pfEnergyLow(k)&...
                   structSrc.pfEnergy < structDet.pfEnergyHigh(k));
        pfPhoton(k) = structSrc.fCurxSize / structProj.nChannel * ...
                structSrc.pfPhoton(nSrcLoc)' * pfProjminusExp';
    end
    pfPhoton = pfPhoton(:);
end