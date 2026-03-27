
function [rCh1, rCh2, vX, vY, dx, dy] =registerStackColumn(ch1, ch2, regCh, info)

if nargin <3
    regCh = 2;
end

if nargin <2 || isempty(ch2)
    regCh = 1;
    ch2 = [];
end

ch1 = single(mat2gray(ch1));
ch2 = single(mat2gray(ch2));

switch regCh
    case 1
        refMov = ch1;
    case 2
        refMov = ch2;
end
%%
%create a Gaussian filter for filtering registration frames
hGauss = fspecial('gaussian', [5 5], 1);

middleFrame = floor(size(refMov,3)/2);
dx = zeros(size(refMov,3),1);
dy = zeros(size(refMov,3),1);
dx(middleFrame) =0; dy(middleFrame) =0;

overfprintf(0, 'Registering...');
fDone = 0;
nMsgChars =0;
for t = (middleFrame+1): size(refMov,3)
    target = imfilter(refMov(:,:,t-1), hGauss, 'same', 'replicate');
    %find the best registration translation
    fftFrame = fft2(imfilter(refMov(:,:,t), hGauss, 'same', 'replicate'));
    output = img.dftregistration(fft2(target), fftFrame, 5);
    dx(t) = output(4) + dx(t-1);
    dy(t) = output(3) + dy(t-1);
    fDone = fDone+1;
    nMsgChars =overfprintf(nMsgChars, sprintf('%d', fDone));
end


for t = 1: middleFrame-1
    target = imfilter(refMov(:,:,middleFrame-t+1), hGauss, 'same', 'replicate');
    %find the best registration translation
    fftFrame = fft2(imfilter(refMov(:,:,middleFrame-t), hGauss, 'same', 'replicate'));
    output = img.dftregistration(fft2(target), fftFrame, 5);
    dx(middleFrame-t) = output(4) + dx(middleFrame-t+1);
    dy(middleFrame-t) = output(3) + dy(middleFrame-t+1);
    fDone = fDone+1;
    nMsgChars = overfprintf(nMsgChars, sprintf('%d', fDone));

end
fprintf(' Complete.\n');

ddx = [0; diff(dx)];
ddy = [0; diff(dy)];
dx(ddx>5 | ddy>5) = 0;
dy(ddx>5| ddy>5) =0;

[rCh1, vX, vY]=img.translate(ch1,dx, dy);

nCh = info.nChannels;

if nCh ==1

    %     zStack{iCh} = zstack.registerStackColumn(zStack{iCh});

    rCh1= uint16(mat2gray(rCh1)*(2^16-1));

    saveastiff(zStack{iCh}, fullfile(info.folderZstack , [info.expRef, '_zStackMean_reg.tif']));
else

    [rCh2, vX, vY]=img.translate(ch2,dx, dy);

    for iCh = 1:nCh


        switch iCh

            case 1
                rCh1= uint16(mat2gray(rCh1)*(2^16-1));
                tiff.saveastiff(rCh1, fullfile(info.folderZstack , [info.expRef, '_zStackMean_G_reg.tif']));
            case 2
                rCh2= uint16(mat2gray(rCh2)*(2^16-1));

                tiff.saveastiff(rCh2, fullfile(info.folderZstack , [info.expRef, '_zStackMean_R_reg.tif']));

        end
    end
end


