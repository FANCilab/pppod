function [n2keep, choices, keeper] =sortROIred(multiplets, R, C, G, M)

nRois = size(R,3);

rsat = prctile(R(:), [5,95]);
csat = prctile(C(:), [5,95]);
gsat = prctile(G(:), [5,95]);

pfig = figure;
drw = figure;

im =1;
choices = [];
% acc = 0;
keeper =[];
while im <= numel(multiplets)
    
    if im == numel(multiplets)
        fprintf('LAST ROI!! \n');
    end
    % plotting and feedback to keep or discard peaks
    figure(pfig); clf
    nDuplic = numel(multiplets{im});
    
   
    for iD = 1:nDuplic
        
        rim{iD} = subplot(3, nDuplic, iD);
        
        imagesc(R(:,:, multiplets{im}(iD))); axis image; colormap(rim{iD}, gray); hold on
        caxis(rsat)
        contour(any(M(:,:,multiplets{im}(iD)),3),1,'r', 'Linewidth', 0.2)
        title('R');
        cim{iD} = subplot(3, nDuplic,nDuplic + iD);
        
        imagesc(C(:,:, multiplets{im}(iD))); axis image; colormap(cim{iD}, gray); hold on
        caxis(csat)
        
        contour(any(M(:,:,multiplets{im}(iD)),3),1,'r', 'Linewidth', 0.2)
        title('C');

        gim{iD} = subplot(3, nDuplic,2*nDuplic + iD);
        
        imagesc(G(:,:, multiplets{im}(iD))); axis image; colormap(gim{iD}, gray); hold on
        caxis(gsat)
        
        contour(any(M(:,:,multiplets{im}(iD)),3),1,'r', 'Linewidth', 0.2)
         title(sprintf('G, roi # %d of %d', im, numel(multiplets)))
    end
    
    fprintf('ROI %i:\n Inspect and press: \n n to redraw the current ROI before continue \n i to accept an IN \n p to accept a ePC \n x to leave unclassified \n d to discard \n u to undo last choice\n e to end here - discard the rest \n q to quit\n', im);
    
    waitforbuttonpress;
    c = get(gcf, 'CurrentCharacter');
    
    
    switch c
        case 'q'
            return;
            
        case 'n'
            fprintf('\n Draw new ROI and then press enter twice');
            
            figure (drw); clf;
            new = imagesc(mean(R(:,:, multiplets{im}), 3)); axis image; colormap(gray);
            newroi = imfreehand(gca);
            wait(newroi);
            multiM{im} = repmat(createMask(newroi,new), [1 1 nDuplic]);
            pause;
            im = im;
            
        case 'e'
            keeper = [keeper, zeros(1, numel(multiplets)- im +1)];
            choices = [choices, repmat('d', [1, numel(multiplets)- im+1])];
            im = numel(multiplets) +1;
            
        case 'i'
            %             acc = acc +1;
            choices = [choices c];
            if numel(multiplets{im}) >1
                waitforbuttonpress;
                
                keeper(im) = str2num(get(gcf, 'CurrentCharacter'));
            else
                keeper(im) = 1;
            end
            %             accepted(acc,:) = neurons(multiplets{im}(n(im)),:);
            
            im = im+1;
            
        case 'p'
            %             acc = acc +1;
            choices = [choices c];
            if numel(multiplets{im}) >1
                waitforbuttonpress;
                
                keeper(im) = str2num(get(gcf, 'CurrentCharacter'));
            else
                keeper(im) = 1;
            end
            %             accepted(acc,:) = neurons(multiplets{im}(n(im)),:);
            
            im = im+1;
            
        case 'x'
            %             acc = acc +1;
            choices = [choices c];
            
                keeper(im) = 0;
            %             accepted(acc,:) = neurons(multiplets{im}(n(im)),:);
            
            im = im+1;
            
        case 'd'
            choices = [choices 'd'];

            keeper(im) = 0;
            
            im = im+1;
            
        case 'u'
            
            if im >1 
            if choices(end) == 'i' || choices(end) == 'p' || choices(end) == 'x'
                im = im-1;
                %                 acc = acc -1;

            elseif choices(end) == 'd' && im >1
                im = im-1;
            end
            choices = choices(1:(end-1));
            keeper = keeper(1:(end-1));
            else
            fprintf('You are looking at the first ROI, cannot go back \n');
            end
    end
    
end

n2keep = zeros(nRois, 1);
for im = 1:numel(multiplets)
    if keeper(im)>0
    n2keep(multiplets{im}(keeper(im))) = 1;
    end
end

end