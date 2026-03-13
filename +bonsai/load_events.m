function event = load_events(ops)


% to do for Heather

% - Check saving of timestamps
% - sparse noise has values -128, -1, 0
% - fix folders structure of saved data
% - implement automatic saving in correct folder
% - check rotary encoder initial values
% - add timestamp to blink_box
% - save stimulus metadata - position, duration, size
%- save screen info data - pixels size, refresh rate, mouse distance,
%   geometric ocnfiguration
% - map stimulus to degree

ops.root_dir = fullfile(ops.root_storage, ops.subject, ops.date, num2str(ops.exp));

% Preallocate structs to hold data
event = struct();

if exist(fullfile(ops.root_dir, 'SparseNoise_Log.bin'), 'file')
    stim_type = 'sparse_noise';
elseif exist(fullfile(ops.root_dir, 'metadata.csv'), 'file')
    stim_type = 'gratings';
end
    fprintf('Loading %s experiment', stim_type);

    fs = 1000; %harp sampling frequency;
%% photodiode
filepath = fullfile(ops.root_dir, 'photodiode.csv');
if exist(filepath, 'file')
    phd= readmatrix(filepath);
    event.photodiode.phd =phd(:,2);
    event.photodiode.time =phd(:,1);
    event.exp_start_t = event.photodiode.time(1);
    event.photodiode.time = event.photodiode.time  - event.exp_start_t;
    fprintf('Loaded: %s\n', filepath);
else
    fprintf('File not found: %s\n', filepath);
end

%% stim times (times_grid)
switch stim_type
    case 'sparse noise'
        filepath = fullfile(ops.root_dir, 'times_grid.csv');

        if exist(filepath, 'file')
            event.stim_times= readmatrix(filepath);
            fprintf('Loaded: %s\n', filepath);
        else
            fprintf('times_grid.csv not found: %s\n', filepath);
        end

    case 'gratings'
        filepath = fullfile(ops.root_dir, 'stimon.csv');

        if exist(filepath, 'file')
            stimon = readmatrix(filepath);
            stimont = stimon(:,1)- event.exp_start_t;
            stimont(isnan(stimont))=[];
            flips = diff(stimont)>0.1;
            onsetflip = logical([1; flips(1:end)]);  
            offsetflip = logical([flips(1:end); 1]);  

            event.stim_times.onset = stimont(onsetflip); % This would fail if the recording starts in the middle of a stimulus
            event.stim_times.offset = stimont(offsetflip);

           if numel(event.stim_times.onset) ~= numel(event.stim_times.offset)
                warning('You have not recorded all the stimuli onset and offset. Times might be wrong')
            end
            fprintf('Loaded: %s\n', filepath);
        else
            fprintf('stimon.csv not found: %s\n', filepath);
        end

end

%% flicker square
filepath = fullfile(ops.root_dir, 'blink_box.csv');
if exist(filepath, 'file')
    event.blink_command = readmatrix(filepath);

    fprintf('Loaded: %s\n', filepath);
else
    fprintf('File not found: %s\n', filepath);
end

%% rotary encoder
filepath = fullfile(ops.root_dir, 'encoder_log.csv');
if exist(filepath, 'file')
    rencoder= readmatrix(filepath);
    event.rencoder.pos = cumsum(rencoder(:,2));
    event.rencoder.time = rencoder(:,1)- event.exp_start_t;
    
    fprintf('Loaded: %s\n', filepath);
else
    fprintf('File not found: %s\n', filepath);
end


%% sparse noise stimulus
filepath = fullfile(ops.root_dir, 'SparseNoise_Log.bin');
if exist(filepath, 'file')

    fid = fopen(filepath, 'rb');
    stimulus = fread(fid, 'int8');  % Assuming int8 for SparseNoise_Log1
    fclose(fid);

    % TODO: save and load metadata about the stimulus. Currently hard
    % coded.
    try
            stimulus = reshape(stimulus, 10, 30, []);

    catch
    stimulus = reshape(stimulus, 7, 27, []);
    end
    stimulus(stimulus ==0) = +1;
    stimulus(stimulus ==-128) = 0;
    event.sparse_noise.frames = stimulus;

    fprintf('Loaded: %s\n', filepath);
else
    fprintf('File not found: %s\n', filepath);
        stim_type = 'gratings';

end
%% digital timestamps (mic frames)
filepath = fullfile(ops.root_dir, 'timestamped_digital.csv');

if exist(filepath, 'file')
    switch stim_type
        case 'sparse_noise'
            try
            dig_timestamps= readtable(filepath);
            frame_on = strcmp(dig_timestamps{:,3}, 'DI3');
            event.frame.on = table2array(dig_timestamps(frame_on,1))- event.exp_start_t;
            event.frame.off = table2array(dig_timestamps(~frame_on,1))- event.exp_start_t;
            event.frame.number = cumsum(double(frame_on));
            event.frame.times = table2array(dig_timestamps(:,1))- event.exp_start_t;

            fprintf('Loaded: %s\n', filepath);
            catch
            dig_timestamps= readtable(filepath);
            frame_on = strcmp(dig_timestamps{:,2}, 'DI3');
            event.frame.on = table2array(dig_timestamps(frame_on,1))- event.exp_start_t;
            event.frame.off = table2array(dig_timestamps(~frame_on,1))- event.exp_start_t;
            event.frame.number = cumsum(double(frame_on));
            event.frame.times = table2array(dig_timestamps(:,1))- event.exp_start_t;

            fprintf('Loaded: %s\n', filepath);
            end

        case 'gratings'
            dig_timestamps= readtable(filepath);
            frame_on = strcmp(dig_timestamps{:,2}, 'DI3');
            event.frame.on = table2array(dig_timestamps(frame_on,1))- event.exp_start_t;
            event.frame.off = table2array(dig_timestamps(~frame_on,1))- event.exp_start_t;
            event.frame.number = cumsum(double(frame_on));
            event.frame.times = table2array(dig_timestamps(:,1))- event.exp_start_t;
            fprintf('Loaded: %s\n', filepath);

    end
else
    fprintf('File not found: %s\n', filepath);
end



%% gratings stimuli
filepath = fullfile(ops.root_dir, 'metadata.csv');
if exist(filepath, 'file')

    grating = bonsai.parse_gratings_protocol(filepath);
    % grating.parameters= readmatrix(filepath);
    % grating.ang= grating.parameters(:,1);
    % grating.size = grating.parameters(:,4);
    % grating.tf = grating.parameters(:,3);
    % grating.sf = grating.parameters(:,3);
    % grating.contrast = grating.parameters(:,5);
    % grating.seq = 1+ grating.ang/30;
    % grating.labels = unique(grating.ang);
    % grating.labels = num2cell(grating.labels);
    % grating.is_stim = unique(grating.ang)<360;
    % grating.is_blanks = unique(grating.ang)>=360;
    % grating.labels(grating.is_blanks) = {'blank'};
   
    if isfield(event.stim_times, 'onset')
        grating.stimTimes.onset = event.stim_times.onset;
        grating.stimTimes.offset = event.stim_times.offset;
    else
        grating.stimTimes = phdStimTimes_LFR(event.photodiode);

    end

    grating.stimMatrix= bonsai.buildStimMatrix(grating.stimSequence, grating.stimTimes, event.frame.on);

    event.grating = grating;
    fprintf('Loaded: %s\n', filepath);
else
    fprintf('File not found: %s\n', filepath);
end



end


function times = phdStimTimes_LFR(photodiode)

phd = photodiode.phd;

phd=smoothdata(phd,'gaussian', 50);
phd = abs(phd-median(phd));
phd=(phd-min(phd))/(max(phd)-min(phd));
% deltas = smoothdata(abs([0; diff(phd)]),'rlowess',50);
dp = medfilt1(abs([0; diff(phd)]), 50, 'truncate');
dp = medfilt1(abs([0; diff(phd)]), 50, 'truncate');

thr=0.005; % using one threshold here

isHigh=dp>thr;
% deltas=[0; diff(above)];
risingEdges = find(diff(isHigh) == 1) + 1;
fallingEdges = find(diff(isHigh) == -1) + 1;

minDur = 800;

dur = risingEdges(2:end) - fallingEdges(1:end-1);
valid = find(dur>minDur);
risingEdges = [risingEdges(1); risingEdges(valid+1)];
fallingEdges = [fallingEdges(valid); fallingEdges(end)];

times.onset = photodiode.time(risingEdges);
times.offset = photodiode.time(fallingEdges);

figure; plot(photodiode.time, dp);
hold on;
plot(times.onset, 0.005*ones(size(times.onset)), 'or')
plot(times.offset, 0.005*ones(size(times.offset)), 'ob')

end

