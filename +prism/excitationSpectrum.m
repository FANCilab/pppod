function [acs, cs, nw] = excitationSpectrum(fp, queryW, doPlot)

if nargin <3
    doPlot = 0;
end

if nargin <2 || isempty(queryW)
    queryW = 680:1400;
end


if ispc
    root = 'C:\Users\Federico\Google Drive\CarandiniLab\CarandiniLab_MATLAB\FedericoBox\2P\PresynapticNetwoks\+prism\2Pspectra';
else
    root = '/Users/Federico/Google Drive/CarandiniLab/CarandiniLab_MATLAB/FedericoBox/2P/PresynapticNetwoks/+prism/2Pspectra';
    
end

file = fullfile(root, [fp, '_w.txt']);

if exist(file)
    
    fileID = fopen(file, 'r');
    
    head = textscan(fileID, '%s', 3);
    
    data = textscan(fileID, '%f %f %f');
    
    fclose(fileID);
    
    w = data{1};
    cs = data{2};
    acs = data{3};
    
%      nw = 680:1400;
%    iw = 680:1400;
%    ics = interp1(w, cs, iw, 'linear');
%    iacs = interp1(w', acs, iw, 'linear');
%    ics = gaussFilt(ics, 10); 
%    iacs = gaussFilt(iacs, 10);


    nw = queryW;
    cs = interp1(w, cs, nw, 'linear');
    acs = interp1(w', acs, nw, 'linear');


    if doPlot
    figure; 
    plot(nw, cs, 'k'); hold on; plot(nw, acs, 'r');
    xlabel('wavelength (nm)');
    ylabel('cs (black) -- acs(red)')
    formatAxes
    end
else
    
    warning('Missing data from FP requested. Try TagBFP dsRed2 mCherry boundGCaMP6 freeGCaMP6')
end

end

