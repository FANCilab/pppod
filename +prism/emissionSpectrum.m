function [ems, nw, pmt] = emissionSpectrum(fp, queryW, doPlot)

if nargin <3
    doPlot = 0;
end

if nargin <2 || isempty(queryW)
    queryW = 300:800;
end

if ispc
    root = 'C:\Users\Federico\Google Drive\CarandiniLab\CarandiniLab_MATLAB\FedericoBox\2P\PresynapticNetwoks\+prism\2Pspectra';
else
    root = '/Users/Federico/Google Drive/CarandiniLab/CarandiniLab_MATLAB/FedericoBox/2P/PresynapticNetwoks/+prism/2Pspectra';
    
end

pmt_file = fullfile(root, ['GaAsP_PMT', '.txt']);

if exist(pmt_file)
    
    fileID = fopen(pmt_file, 'r');
    
    head = textscan(fileID, '%s', 2);
    
    data = textscan(fileID, '%f %f');
    
    fclose(fileID);
    
    w = data{1};
    cs = data{2};
    cs = cs./max(cs);
    nw = queryW;
    pmt = interp1(w, cs, nw, 'linear');
    
else
    pmt = ones(numel(queryW), 1);
end

file = fullfile(root, [fp, '_ems.txt']);

if exist(file)
    
    fileID = fopen(file, 'r');
    
    head = textscan(fileID, '%s', 2);
    
    data = textscan(fileID, '%f %f');
    
    fclose(fileID);
    
    w = data{1};
    cs = data{2};
    
    nw = queryW;
    ems = interp1(w, cs, nw, 'linear');

%     cs = gaussFilt(cs',1);
%     acs = gaussFilt(acs', 1);
    
%     for iQ = 1: numel(queryW)
%         
%     ncs(iQ) = cs(nw == queryW(iQ));
%     nacs(iQ) = acs(nw == queryW(iQ));
%     nnw(iQ) = nw(nw == queryW(iQ));
%     end

    if doPlot
    figure; 
    plot(nw, ems, 'k'); hold on;
    xlabel('wavelength (nm)');
    ylabel('ems(black)')
    formatAxes
    end
else
    
    warning('Missing data from FP requested')
end

end

