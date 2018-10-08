%% Compute simulated multi-energy projection data
%  structProj:
%  Output:
%       pfLac: matrix of size nChannel * nAngle * length(structMas.pfEnergy) 
%       pfPhoton: matrix of size nChannel * nAngle * structDet.nEnergyBin,
%                 simulation of detected photons with no noise
%       pfPhotonrnd: pfPhoton plus poission random noise 
function structProj = computeProjME(structImg,structProj,structSrc,...
                                    structGeo,structMas,structDet)
    % compute LAC projections at every energy
    for m = 1:length(structMas.pfEnergy) 
        [sinogram_id, sinogram] = astra_create_sino(...
            structImg.pfImgLac(:,:,m), structGeo.proj_id);
        structProj.pfLac(:,:,m) = sinogram';
    end % end energy
    
    % compute photon number and integration over energy
    structProj.pfPhoton = zeros(structProj.nChannel,structProj.nAngle,...
                                structDet.nEnergyBin);
    for k = 1:structDet.nEnergyBin % for every energy bin 
        % find MAC data points inside the kth bin
        nProjLoc = (structMas.pfEnergy >=structDet.pfEnergyLow(k)&...
                    structMas.pfEnergy < structDet.pfEnergyHigh(k));
        pfProjminusExp = exp(-structProj.pfLac(:,:,nProjLoc));
        % find source data points inside the kth bin
        nSrcLoc = find(structSrc.pfEnergy >= structDet.pfEnergyLow(k)&...
                       structSrc.pfEnergy < structDet.pfEnergyHigh(k));
    
        for n = 1: length(nSrcLoc)
            structProj.pfPhoton(:,:,k) = structProj.pfPhoton(:,:,k) +  ...
                structSrc.fCurxSize / structProj.nChannel * ...
                structSrc.pfPhoton(nSrcLoc(n)) * pfProjminusExp(:,:,n);
        end  
        PhotonIn = sum(structSrc.pfPhoton(nSrcLoc)) * structSrc.fCurxSize / structProj.nChannel;
        
        structProj.pfPhotonrnd(:,:,k) = poissrnd(structProj.pfPhoton(:,:,k));
        structProj.pfLACrnd(:,:,k) =  log(PhotonIn./max(structProj.pfPhotonrnd(:,:,k),10));
    end
    structProj.pfPhotonrnd = poissrnd(structProj.pfPhoton);

end

