function save_figure_compat(fig, fileName)
% SAVE_FIGURE_COMPAT Save figure using exportgraphics when available.

    [folderPath, ~, ext] = fileparts(fileName);
    if ~isempty(folderPath) && ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end

    ext = lower(ext);
    if exist('exportgraphics', 'file') == 2
        exportgraphics(fig, fileName, 'Resolution', 150);
        return
    end

    switch ext
        case '.png'
            print(fig, fileName, '-dpng', '-r150');
        case '.pdf'
            print(fig, fileName, '-dpdf', '-r150');
        case {'.jpg', '.jpeg'}
            print(fig, fileName, '-djpeg', '-r150');
        case '.tif'
            print(fig, fileName, '-dtiff', '-r150');
        otherwise
            saveas(fig, fileName);
    end
end
