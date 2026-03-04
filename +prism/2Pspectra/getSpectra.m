function [acs, cs, nw] = getSpectra(fp, doPlot)

if nargin <2
    doPlot = 0;
end

root = 'C:\Users\Federico\Google Drive\CarandiniLab\CarandiniLab_MATLAB\FedericoBox\2P\PresynapticNetwoks\+prism\2Pspectra';

file = fullfile(root, [fp, '_w.txt']);

if exist(file)
    
    fileID = fopen(file, 'r');
    
    head = textscan(fileID, '%s', 3);
    
    data = textscan(fileID, '%f %f %f');
    
    fclose(fileID);
    
    w = data{1};
    cs = data{2};
    acs = data{3};
    
    nw = 680:1400;
    cs = interp1(w, cs, nw, 'linear');
    acs = interp1(w', acs, nw, 'linear');

    cs = gaussFilt(cs',1);
    acs = gaussFilt(acs', 1);
    
    if doPlot
    figure; 
    plot(nw, cs, 'k'); hold on; plot(nw, acs, 'r');
    xlabel('wavelength (nm)');
    ylabel('cs (black) -- acs(red)')
    formatAxes
    end
else
    
    warning('Missing data from FP requested')
end

end

