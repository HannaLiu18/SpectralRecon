%% Generate true material image
%  structImg:
%  Input:
%       fFOV:     field of view(mm)
%       nPixel:   image size
%       nRegion:  number of different regions in true image
%       pfE:      matrix of size nRegion * 6, for definition see phantom
%                 function in matlab, elements in first column are additive 
%                 indices for different materials
%       nNoContr: number of regions without contrast
%       pfMassPer:contrast element concentrations in mass percent at each
%                 contrast regions
%  Output:
%       pfDensity:  densities for different regions, for contrast-water
%                   mixture regions, density need to be computed 
%       pfMassPerMatrix: matrix of size nRegion*structMas.nBase, each row
%                        records the mass percentage of different base 
%                        materials in one certain region
%       pfImgLabel: true image of size nPixel * nPixel,pixel values are
%                   indices for different materials
%       pfImgLac:   true image of size nPixel * nPixel *
%                   length(structMas.pfEnergy),pixel values are linear
%                   atten. coeff. at different energies
%  structContrast:
%  Input:
%       fMMCompound: molar mass of the contrast molecular 
%       fMMAtom:     molar mass of the contrast element
%       fRatio:      fMMCompound / fMMAtom
%       fMolarOS:    contrast molarity in original solution (mol/ml)
%       fDensityOS:  density of the original solution (g/ml)
function [structImg,structContrast] = setTrueImg(structMas)
    % Input
    structImg.fFOV = 300.0;    
    structImg.nPixel = 128;
    structImg.nRegion = 7;
    structImg.nNoContr = 4;
    % Regions are labeled from 1 to 7 actually
    structImg.pfE=[1,   0.9,    0.9,    0,      0,      0;   % Soft tissue
                   1,   0.2,    0.2,    0,      0.5,    0;   % Heart
                   2,   0.18,   0.18,   0,      -0.06,  0;   % Cortical bone
                   1,   0.14,   0.14,   0,      -0.06,  0;   % Yellow marrow
                   4,   0.1,    0.1,    -0.3,   -0.5,   0;   % Gadolinium1
                   5,   0.02,   0.025,	0.4,	-0.2,  18;   % Gadolinium2
                   6,   0.05,   0.06,   -0.45,   0.3,   40]; % Gadolinium2
    
    structImg.pfDensity(1:structImg.nNoContr) = [1.02,1.06,1.92,0.98];
    structImg.pfMassPer(1:structImg.nRegion - structImg.nNoContr) = ...
        [0.03,0.01,0.01];
    % Properties of the contrast
    structContrast.fMMCompound = 583;
    structContrast.fMMAtom = 157;
    structContrast.fMolarOS = 0.5*10^-3;
    structContrast.fDensityOS = 1.15;
    structContrast.fRatio = structContrast.fMMCompound / structContrast.fMMAtom;

    % Output
    % Compute mass persentage of contrast in the solution
    structImg.pfMassPerMatrix = eye(structMas.nBase);
    for k = 1: (structImg.nRegion - structImg.nNoContr)
        structImg.pfMassPerMatrix(k + structImg.nNoContr,structMas.nBase-1) = ...
            structImg.pfMassPer(k) * structContrast.fRatio;
        structImg.pfMassPerMatrix(k + structImg.nNoContr,structMas.nBase) = ...
            1 - structImg.pfMassPerMatrix(k + structImg.nNoContr,structMas.nBase-1);
    end
    % Compute contrast solution densities
    % Reference: E Roessl and R Proksa(2007),equation(9) and appendix A4
    for k = 1:(structImg.nRegion - structImg.nNoContr)
        x = structContrast.fMolarOS * structContrast.fMMAtom / structImg.pfMassPer(k)...
            - structContrast.fDensityOS; % A4
        
        structImg.pfDensity(k + structImg.nNoContr) = ...
            (structContrast.fDensityOS + x)/(1 + x);   %A1
    end
    
    % Images showing region indices
    structImg.pfImgLabel = phantom(structImg.pfE,structImg.nPixel);
    % Images showing lac at different energies
    structImg.pfImgLac = zeros(structImg.nPixel,structImg.nPixel,...
        length(structMas.pfEnergy));
    for m = 1:length(structMas.pfEnergy) % for every energy 
        for n = 1:structImg.nRegion % for every region
            % find pixels in this region
            fTemp = (structImg.pfImgLabel(:,:) == n);
            % compute MAC
            fMAC = structImg.pfMassPerMatrix(n,:) * structMas.pfMas(m,:)';
            % compute LAC
            fTemp = fTemp * fMAC * structImg.pfDensity(n) * 0.1;
            structImg.pfImgLac(:,:,m) = structImg.pfImgLac(:,:,m) + fTemp;
        end % end material
    end % end energy
end