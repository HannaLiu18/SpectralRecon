%% Set forward projection geometry(using astra-toolbox)
%  structGeo: 
%  Output:
%       vol_geom: object geometry
%       proj_geom:projectio geometry
%       proj_id:  projector id
%       W:        system matrix
%  structProj:(Currently consider parallel geometries)
%  Input:
%       nChannel: number of detectors
%       nAngle:   number of projections
%       fFOV:     field of view
%  Output:
%       fDetWidth:detector width
function [structGeo,structProj] = setGeometryFwd(structImg)
    % Input
    structProj.nChannel = 672 / 2;
    structProj.nAngle = 180; % 90;
    structProj.fFOV = 500;
    % Output
    structProj.fDetWidth = structProj.fFOV / structProj.nChannel;
    % Set geometries for astra
    %%%% Step One:Set volume and projection geometries %%%%
    structGeo.vol_geom = astra_create_vol_geom(...
         structImg.nPixel, structImg.nPixel,...
        -structImg.fFOV/2, structImg.fFOV/2,...
        -structImg.fFOV/2, structImg.fFOV/2);
    structGeo.proj_geom = astra_create_proj_geom('parallel', ...
        structProj.fDetWidth, structProj.nChannel, ...
        linspace(0,2*pi,structProj.nAngle));
    %%%% Step Two: Create projector %%%%
    structGeo.proj_id = astra_create_projector('linear', ...
        structGeo.proj_geom, structGeo.vol_geom);

    %%%% Step Three: Generate system matrix %%%%
    structGeo.W = opTomo('linear', structGeo.proj_geom, structGeo.vol_geom); 
end