clear;

%% select movie 890
[nameG, path] = uigetfile('*.tif');

cd(path);

nameR = [nameG(1:end-5), 'R.tif'];

movG890 = img.loadFrames(nameG);
movR890 = img.loadFrames(nameR);
%% select movie 780
[nameG, path] = uigetfile('*.tif');

cd(path);

nameR = [nameG(1:end-5), 'R.tif'];

movG780 = img.loadFrames(nameG);
movR780 = img.loadFrames(nameR);

%% subtract G channel from red

samples = 1:1000:numel(movG890);
samplesG= cat(2, movG890(samples), movG780(samples));
samplesR = cat(2, movR890(samples), movR780(samples));
bleedG = robustfit(single(samplesG'), single(samplesR'));

R890 = mat2gray(double(movR890) - double(movG890).*bleedG(2) -bleedG(1));
R780 = mat2gray(double(movR780) - double(movG780).*bleedG(2) -bleedG(1));

%% register the red channels

nPlanes = size(R890);
hGauss = fspecial('gaussian', [5 5], 1);

for iPlane = 1:nPlanes
    
    target = imfilter(double(movG890(:,:,iPlane)), hGauss, 'same', 'replicate');
    %find the best registration translation
    fftFrame = fft2(imfilter(double(movG780(:,:,iPlane)), hGauss, 'same', 'replicate'));
    output = dftregistration(fft2(target), fftFrame, 20);
    dx(iPlane) = output(4);
    dy(iPlane) = output(3);
    iPlane
    
end

[R780, vX, vY] = img.translate(R780, dx, dy);



%% get coefficient for mCherry and dsRed
addpath('/Users/Federico/Google Drive/CarandiniLab/CarandiniLab_MATLAB/2Pspectra');
[mCherry, ~,  w] = getSpectra('mCherry',1);
[dsRed, ~, w] = getSpectra('DsRed2',1);

farRed780 = mCherry(w ==780);
farRed890 = mCherry(w == 890);

red890 = dsRed(w==890);
red780 = dsRed(w==780);


figure; plot(R890(samples), R780(samples),'.');

bleedR = robustfit(R890(samples), R780(samples));

mCherryMov = R780 - bleedR(2)*R890 -bleedR(1);
dsRedMov = R890;
%% solve linear combination

trueSig = pinv([red890/(red890 + red), farRed890; red780, farRed780])*[R890(:)'; R780(:)'];

trueSign = [R890(:)'; R780(:)']\[red890, farRed890; red780, farRed780];

dsRedMov = reshape(trueSig(1,:), size(R890));
mCherryMov = reshape(trueSig(2,:), size(R890));

%%

saveastiff(int16(dsRedMov*10000), [path(1:end-3), 'dsRed.tif']);
saveastiff(int16(mCherryMov*10000), [path(1:end-3), 'mCherry.tif']);


