function choices = sortROI(R, C,M)

nRois = size(R,3);

rsat = prctile(R(:), [5,95]);
csat = prctile(C(:), [5,95]);

pfig = figure;
drw = figure;

im =1;
choices = [];
% acc = 0;
n =[];

while im <= nRois
    
    
    % plotting and feedback to keep or discard peaks
    figure(pfig); clf
    
    gim = subplot(1,2,1);
    
    imagesc(R(:,:, im)); axis image; colormap(gim, gray); hold on
%     caxis([ 0 prctile(makeVec(R(:,:,im)), 99)])
    caxis(rsat)
    contour(any(M(:,:,im),3),1,'r', 'Linewidth', 0.2)
    
    rim = subplot(1,2,2);
    
    imagesc(C(:,:, im)); axis image; colormap(rim, gray); hold on
%     caxis([ 0 prctile(makeVec(C(:,:,im)), 99)])
    caxis(csat)

    
    contour(any(M(:,:,im),3),1,'r', 'Linewidth', 0.2)
    title(sprintf('roi # %d of %d', im, nRois))


fprintf('ROI %i:\n Inspect and press: \n n to redraw the current ROI before continue \n i to accept an IN \n p to accept a ePC \n x to leave unclassified \n d to discard \n u to undo last choice\n e to end here - discard the rest \n q to quit\n', im);

waitforbuttonpress;
c = get(gcf, 'CurrentCharacter');


switch c
    case 'q'
        return;
        
    case 'n'
        fprintf('\n Draw new ROI and then press enter twice');
        
        figure (drw); clf;
        new = imagesc(R(:,:,im)); axis image; colormap(gray);
        newroi = imfreehand(gca);
        wait(newroi);
        M(:,:,im) = repmat(createMask(newroi,new), [1 1 1]);
        pause;
        im = im;
        
    case 'e'
        n = [n, ones(1, nRois- im +1)];
        choices = [choices, repmat('d', [1, nRois- im+1])];
        im = nRois +1;
        
    case 'i'
        %             acc = acc +1;
        choices = [choices c];
        %             if numel(multiplets{im}) >1
        %                 waitforbuttonpress;
        %
        %                 n(im) = str2num(get(gcf, 'CurrentCharacter'));
        %             else
        %                 n(im) = 1;
        %             end
        %             accepted(acc,:) = neurons(multiplets{im}(n(im)),:);
        
        im = im+1;
        
    case 'p'
        %             acc = acc +1;
        choices = [choices c];
        %             if numel(multiplets{im}) >1
        %                 waitforbuttonpress;
        %
        %                 n(im) = str2num(get(gcf, 'CurrentCharacter'));
        %             else
        %                 n(im) = 1;
        %             end
        %             accepted(acc,:) = neurons(multiplets{im}(n(im)),:);
        
        im = im+1;
        
    case 'x'
        %             acc = acc +1;
        choices = [choices c];
        %             if numel(multiplets{im}) >1
        %                 waitforbuttonpress;
        %
        %                 n(im) = str2num(get(gcf, 'CurrentCharacter'));
        %             else
        %                 n(im) = 1;
        %             end
        %             accepted(acc,:) = neurons(multiplets{im}(n(im)),:);
        
        im = im+1;
        
    case 'd'
        choices = [choices 'd'];
        %             if numel(multiplets{im}) >1
        %                 waitforbuttonpress;
        %
        %                 n(im) = str2num(get(gcf, 'CurrentCharacter'));
        %             else
        %                 n(im) = 1;
        %             end
        im = im+1;
        
    case 'u'
        
        if choices(end) == 'i' || choices(end) == 'p' || choices(end) == 'x'
            im = im-1;
            %                 acc = acc -1;
            
        elseif choices(end) == 'd'
            im = im-1;
        end
        choices = choices(1:(end-1));
        %             n = n(1:(end-1));
end

end
end