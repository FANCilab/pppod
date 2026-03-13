function [regMov, vX, vY] = register(mov)

[nx,ny, nt] = size(mov);

filtMov = imgaussfilt(mov); 

meanMov = makeVec(mean(filtMov, 3));

filtMov = reshape(filtMov, nx*ny, nt);

similar = corr(meanMov, filtMov);

thrs = prctile(similar, 90);

target = reshape(mean(filtMov(:, similar> thrs), 2), nx, ny);


hGauss = fspecial('gaussian', [5 5], 1);

% register
dx = zeros(nt, 1);
dy = zeros(nt, 1);

for t = 1: nt;

        fftFrame = fft2(imfilter(mov(:,:,t), hGauss, 'same', 'replicate'));
        output = dftregistration(fft2(target), fftFrame, 20);
        dx(t) = output(4) + dx(t);
        dy(t) = output(3) + dy(t);
        display(t);
end

[regMov, vX, vY]=img.translate(single(mov),dx, dy,'clip');


end