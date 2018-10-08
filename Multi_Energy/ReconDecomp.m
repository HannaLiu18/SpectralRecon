function Image = ReconDecomp(structGeo,x)
    recon_id = astra_mex_data2d('create', '-vol', structGeo.vol_geom, 0);
    sinogram_id = astra_mex_data2d('create', '-sino', structGeo.proj_geom, 0);
    astra_mex_data2d('store', sinogram_id, x);

    cfg = astra_struct('FBP');
    cfg.ProjectorId = structGeo.proj_id;
    cfg.ProjectionDataId = sinogram_id;
    cfg.ReconstructionDataId = recon_id;
    fbp_id = astra_mex_algorithm('create', cfg);
    astra_mex_algorithm('run', fbp_id);
    Image = astra_mex_data2d('get', recon_id);


    % garbage disposal
    astra_mex_data2d('delete', sinogram_id, recon_id);
    astra_mex_algorithm('delete', fbp_id);
end