%% Compute the log-likelihood func for material decomposition in raw data
%  Note that only one channel data is used, maximum likihood of
%  multi-channel data might be impossible to optimize for fminsearch
%  Output:
%       f: function handle of maximum likelihood function
function f = decomplikelihood(A,structEBas,structDet,structSrc,...
                              structProj,nChannel,nAngle)
    % Compute photon numbers from base functions structEBas
    pfPhoton = computePhotonFromBase(A,structEBas,structDet,...
                                     structSrc,structProj);
    pfdata = structProj.pfPhotonrnd(nChannel,nAngle,:);pfdata = pfdata(:);
    f = sum(pfPhoton - pfdata .* log(pfPhoton));
end