
%close all;
clear;
clc;
%% Set path to relevant code

cd(code_repo);
addpath(genpath(code_repo));
addpath('C:\Users\User\Documents\Code\pppod');

%% database of the recordings, one entry per zstack
i = 0;

i = i+1;
db(i).mouse_name    = 'HS034'; % animal name
db(i).date          = '20260206'; % date of the recording
db(i).exp_n         = [2]; % all the experiments in the recording
db(i).expID         = 1; % the experiment you want to compute pixel map of
db(i).regto = 1;

i = i+1;
db(i).mouse_name    = 'HS034'; % animal name
db(i).date          = '20260210'; % date of the recording
db(i).exp_n         = [1]; % all the experiments in the recording
db(i).expID         = 1; % the experiment you want to compute pixel map of
db(i).regto = 1;

i = i+1;
db(i).mouse_name    = 'HS034'; % animal name
db(i).date          = '20260209'; % date of the recording
db(i).exp_n         = [1]; % all the experiments in the recording
db(i).expID         = 1; % the experiment you want to compute pixel map of
db(i).regto = 1;

%% Initialise options for registration

gcp;
tic
clear options;

options.nFramesPerChunk = 1024;
options.nFrames4TargetSelection = 100;
% options.channels = 1;  %% if just Red Ch was acquired
options.registrationChannel = 1;
options.doClipping = false;  %!!!!SUPER IMPORTANT!!!!
options.targetFrame= 'average';
options.doBidi = 0;
%% average and register the stack. Need to add code to save the relative 'info file' in the same folder
% tiffs with the G and R zstack are going to be saved in
% info.folderProcessed. check that it points to a directory on your PC!

for iDb = 3%:numel(db)
    
    %which channel to register to
    options.registrationChannel = db(iDb).regto;

    %extract info about the recording from tiff headers
    info= getExpInfo(db(iDb).mouse_name , db(iDb).date , db(iDb).exp_n(db(iDb).expID), 1);

    %registers each plane of the zstack individually
    [info(iDb), zStack{iDb}] = zstack.regZstack(info, options);

    %registers planes to each other
    [zStack{iDb}{1}, zStack{iDb}{2} ] = zstack.registerStackColumn(zStack{iDb}{1}, zStack{iDb}{2}, 1, info);

end

%%
%%