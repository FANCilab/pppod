function [fpR, fpG] = mixFP(varargin, waveL)

nFP = numel(varargin);

for iF = 1:nFP
    
    fp_name = varargin{iF};
    
    % get excitation spectrum
    excSpectrum(:, iF) = prism.excitationSpectrum(fp_name, waveL);
    
    %get emission spectrum, weighted by bscope red channel 607/70nm
    [emsR(:,iF), ~, pmt] = prism.emissionSpectrum(fp_name, 572:642);

    %get emission spectrum, weighted by bscope green channel 525/50nm
    emsG(:,iF) = prism.emissionSpectrum(fp_name, 500:550);
    
    %calculate mixing coefficient in R and G channel
    fpR(:, iF) =  excSpectrum(:, iF).*nansum(emsR(:,iF));
    fpG(:, iF) =  excSpectrum(:, iF).*nansum(emsG(:,iF));
    
end


end