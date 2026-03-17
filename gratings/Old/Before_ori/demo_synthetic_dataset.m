% DEMO_SYNTHETIC_DATASET Create a small synthetic dataset and run the workflow.
%
% This demo builds a simple experiment with active direction and size,
% generates synthetic time courses and response amplitudes for a few neurons,
% flattens stimulus conditions into the preferred nStim-by-repeat input
% layout, and runs the grating analysis workflow.

rng(1);

% -------------------------------------------------------------------------
% Define stimulus axes
% -------------------------------------------------------------------------
directions = 0:45:315;                 % 8 directions
sizes = [5 10 20];                     % 3 sizes (deg)
tfs = 2;                               % inactive TF
sfs = 0.08;                            % inactive SF
nRep = 6;
t = linspace(-1, 2, 80);               % 80 time samples
nNeurons = 5;

nDir = numel(directions);
nSize = numel(sizes);
nTf = numel(tfs);
nSf = numel(sfs);
nTime = numel(t);
nStim = nDir * nSize * nTf * nSf;

% -------------------------------------------------------------------------
% Build synthetic amplitude responses in canonical form:
%   amp = [dir size tf sf rep neuron]
%   tc  = [dir size tf sf rep time neuron]
% -------------------------------------------------------------------------
ampCanon = zeros(nDir, nSize, nTf, nSf, nRep, nNeurons);
tcCanon = zeros(nDir, nSize, nTf, nSf, nRep, nTime, nNeurons);

preferredDirections = [90 180 270 0 135];
preferredSizes = [10 20 5 10 20];
responsiveScale = [1.2 0.9 0.0 1.0 0.7];

kernel = exp(-0.5 * ((t - 0.35) / 0.18).^2);

for iNeuron = 1:nNeurons
    for iDir = 1:nDir
        dirDiff = abs(mod(directions(iDir) - preferredDirections(iNeuron) + 180, 360) - 180);
        dirGain = exp(-0.5 * (dirDiff / 45).^2);

        for iSize = 1:nSize
            sizeDiff = abs(log2(sizes(iSize) / preferredSizes(iNeuron)));
            sizeGain = exp(-0.5 * (sizeDiff / 0.6).^2);

            meanResp = responsiveScale(iNeuron) * 2.5 * dirGain * sizeGain;

            for iRep = 1:nRep
                trialAmp = meanResp + 0.35 * randn();
                ampCanon(iDir, iSize, 1, 1, iRep, iNeuron) = trialAmp;

                baseline = 0.08 * randn(1, nTime);
                tcCanon(iDir, iSize, 1, 1, iRep, :, iNeuron) = baseline + trialAmp * kernel;
            end
        end
    end
end

% -------------------------------------------------------------------------
% Optional blanks (not used by the requested workflow)
% blankAmp shape: [nBlanks, neuron]
% blankTc shape : [nBlanks, time, neuron]
% Here nBlanks = nSize * nTf * nSf * nRep = 18
% -------------------------------------------------------------------------
nBlanks = nSize * nTf * nSf * nRep;
blankAmp = 0.15 * randn(nBlanks, nNeurons);
blankTc = 0.08 * randn(nBlanks, nTime, nNeurons);

% -------------------------------------------------------------------------
% Pack input struct in the preferred flattened-stimulus format:
%   data.tc  = [nStim, nRep, nTime, nNeuron]
%   data.amp = [nStim, nRep, nNeuron]
% -------------------------------------------------------------------------
data = struct();
data.tc = reshape(tcCanon, [nStim, nRep, nTime, nNeurons]);
data.amp = reshape(ampCanon, [nStim, nRep, nNeurons]);
data.blankTc = blankTc;
data.blankAmp = blankAmp;
data.t = t;
data.directions = directions;
data.sizes = sizes;
data.tfs = tfs;
data.sfs = sfs;
data.activeParams = {'direction', 'size'};

% -------------------------------------------------------------------------
% Run workflow
% -------------------------------------------------------------------------
targetFolder = fullfile(pwd, 'demo_grating_output');
opts = struct('alpha', 0.05, 'saveExt', 'png', 'visible', 'off');

results = analyze_grating_experiment(data, targetFolder, opts); %#ok<NASGU>

disp('Demo complete. Results saved to:');
disp(targetFolder);
