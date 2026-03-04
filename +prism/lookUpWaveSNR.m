function wSnr = lookUpWaveSNR(waveL, z)

saveTo = 'C:\Users\Federico\Google Drive\CarandiniLab\CarandiniLab_MATLAB\Data\rvRetinotopy';
load(fullfile(saveTo, ['SNRvsZ_', 'FR141', '_', '1.mat']), 'ww', 'wFit');

for iw = 1:numel(waveL)
    
    wSnr(iw) = wFit{ww == waveL(iw)}(z)/ wFit{ww == 890}(z);
    
end

 
end