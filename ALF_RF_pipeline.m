%% ALF_RF_pipeline

clear;
clc;

addpath(fileparts(mfilename('fullpath')));

ALFroot = 'D:\Pipeline\Data';
resultsRoot = 'D:\Pipeline\Results';

s2pPaths = [
    "\\10.233.25.135\FANCiNAS1\Data\2P\HS035\Processed\20260512\1"
    "\\10.233.25.135\FANCiNAS1\Data\2P\NM023\Processed\20260209\1_2"
    "\\10.233.25.135\FANCiNAS1\Data\2P\NM023\Processed\20260428\1"
    "\\10.233.25.135\FANCiNAS1\Data\2P\NM034\Processed\20260506\1"
    "\\10.233.25.135\FANCiNAS1\Data\2P\NM035\Processed\20260513\1"
    "\\10.233.25.135\FANCiNAS1\Data\2P\NM035\Processed\20260515\2\suite2p"
    "\\10.233.25.135\FANCiNAS1\Data\2P\NM035\Processed\20260519\2\suite2p"
    "\\10.233.25.135\FANCiNAS1\Data\2P\NM035\Processed\20260520\1\suite2p"
    "\\10.233.25.135\FANCiNAS1\Data\2P\NM036\20260507\Processed\1"
    "\\10.233.25.135\FANCiNAS1\Data\2P\NM037\20260609\Processed"
    "\\10.233.25.135\FANCiNAS1\Data\2P\NM039\20260610\Processed\1"
    "\\10.233.25.135\FANCiNAS1\Data\2P\NM040\20260616\Processed\1"
    "\\10.233.25.135\FANCiNAS1\Data\2P\NM044\20260610\Processed\1"
];

modalities = [
    "deconv"
    "deconv"
    "deconv"
    "deconv"
    "traces"
    "traces"
    "traces"
    "traces"
    "deconv"
    "traces"
    "deconv"
    "deconv"
    "traces"
];

report = ALFexporter(s2pPaths, ALFroot);

for i = 1:height(report)
    if ~report.ok(i)
        continue
    end

    sessions = table( ...
        string(report.subject(i)), ...
        string(report.date(i)), ...
        string(report.session(i)), ...
        'VariableNames', {'subject', 'date', 'session'});

    loadedSession = ALFloader(ALFroot, sessions);
    ALF_main_analyse_RFs(loadedSession, modalities(report.index(i)), resultsRoot, false, true);
end
