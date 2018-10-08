%% X-ray tube spectral Generated from matlab package "spektr 3.0" 
%  Downloaded from http://istar.jhu.edu/downloads/
%  structSrc:
%  Input:
%       pfEnergy:   vector of energies for x-ray tube spectral
%       fCurxSize:  mm^2*mAs
%  Output:
%       pfPhoton:   photons emitted at different energies
function structSrc = setSource()
    % Input
    structSrc.pfEnergy = [1:150];  % (KeV)
    structSrc.fCurxSize = 100;     % Reasonable value?
    % Output
    structSrc.pfPhoton = spektrSpectrum(120, [0 0], 'TASMIP', 0);
    % Plot Figure: X-ray tube spectra
    figure;plot(structSrc.pfEnergy,structSrc.pfPhoton);grid on;
    title('X-ray tube spectra');
    xlabel('E[KeV]');ylabel('Photons/mm^2/mAs at 100 cm from the source');
end