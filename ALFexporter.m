function [sessions, report] = ALFexporter(s2pPaths, ALFroot)
%ALFEXPORTER Export selected Suite2p runs to ALF/ONE from MATLAB.
%
%   sessions = ALFexporter(s2pPaths, ALFroot)
%   [sessions, report] = ALFexporter(s2pPaths, ALFroot)
%
%   Inputs
%   ------
%   s2pPaths : char, string, string array, or cellstr
%       One or more Suite2p run folders. Each run folder should contain
%       db.npy, run.log, and plane*/ folders.
%
%   ALFroot : char or string scalar
%       Output root. Files are written below subject/date/session/alf.
%
%   Output
%   ------
%   sessions : struct array
%       Loader-ready struct array with string scalar fields subject, date, and
%       session. Pass this directly to ALFloader(ALFroot, sessions). Failed
%       exports are omitted.
%
%   report : table
%       One row per Suite2p path with ok/status/message/nWritten/writtenFiles.
%
%   Notes
%   -----
%   This function does NOT call main_export.py or rf_master_export.py.
%   It calls the Python exporter modules directly:
%       calcium_export.single_session_calcium_export(...)
%       wheel_export.single_session_wheel_export(...)
%       sparse_noise_export.single_session_sparse_noise_export(...)
%
%   For Suite2p runs concatenating multiple raw sessions, wheel and sparse-noise
%   export are restricted to the raw session containing SparseNoise_Log.bin,
%   using the same low-level helpers from calcium_export.py, wheel_export.py,
%   and sparse_noise_export.py.
%
%   Put ALFexporter.m in the same folder as these Python files, or make sure
%   that folder is reachable from the MATLAB current folder:
%       calcium_export.py
%       sparse_noise_export.py
%       wheel_export.py
%
%   Example
%   -------
%   s2pPaths = [
%       "D:\Data\2P\NM023\Processed\20260209\1"
%       "D:\Data\2P\NM034\Processed\20260301\1"
%   ];
%   [sessions, report] = ALFexporter(s2pPaths, "D:\Pipeline\Data");
%   ALF = ALFloader("D:\Pipeline\Data", sessions);

    if nargin ~= 2
        error('ALFexporter requires exactly two inputs: ALFexporter(s2pPaths, ALFroot).');
    end

    s2pPaths = normalizePathList(s2pPaths, 's2pPaths');
    ALFroot = normalizeScalarPath(ALFroot, 'ALFroot');

    if isempty(s2pPaths)
        error('s2pPaths is empty. Provide at least one Suite2p run folder.');
    end

    if ~isfolder(ALFroot)
        mkdir(ALFroot);
    end

    pythonFolder = findPythonExporterFolder();
    addPythonPath(pythonFolder);
    bridge = ensureDirectBridgeModule();

    n = numel(s2pPaths);
    ok = false(n, 1);
    status = strings(n, 1);
    subject = strings(n, 1);
    date = strings(n, 1);
    session = strings(n, 1);
    message = strings(n, 1);
    nCalciumWritten = zeros(n, 1);
    nWheelWritten = zeros(n, 1);
    nSparseNoiseWritten = zeros(n, 1);
    nWritten = zeros(n, 1);
    elapsedSeconds = zeros(n, 1);
    writtenFiles = cell(n, 1);

    fprintf('ALFexporter: output root = %s\n', ALFroot);
    fprintf('ALFexporter: Python module folder = %s\n', pythonFolder);
    fprintf('ALFexporter: exporting %d Suite2p run(s).\n', n);

    pathlib = py.importlib.import_module('pathlib');
    edgesPy = py.tuple(num2cell([-135.0, 135.0, -40.0, 40.0]));

    for i = 1:n
        runPath = char(s2pPaths(i));
        fprintf('\n[%d/%d] EXPORT %s\n', i, n, runPath);
        tStart = tic;

        try
            if ~isfolder(runPath)
                error('Suite2p run folder does not exist: %s', runPath);
            end

            result = bridge.export_sparse_noise_run_direct( ...
                pathlib.Path(runPath), ...
                pathlib.Path(ALFroot), ...
                pyargs( ...
                    'encoder_counts_per_revolution', 65536.0, ...
                    'sparse_noise_edges_deg', edgesPy));

            ok(i) = true;
            status(i) = pyDictString(result, 'status');
            subject(i) = pyDictString(result, 'subject');
            date(i) = pyDictString(result, 'date');
            session(i) = pyDictString(result, 'session');
            nCalciumWritten(i) = pyDictNumber(result, 'n_calcium');
            nWheelWritten(i) = pyDictNumber(result, 'n_wheel');
            nSparseNoiseWritten(i) = pyDictNumber(result, 'n_sparse_noise');
            nWritten(i) = pyDictNumber(result, 'n_written');
            writtenFiles{i} = pySequenceToStrings(result.get('written'));
            if status(i) == "skipped"
                message(i) = pyDictString(result, 'message');
            else
                message(i) = sprintf('Wrote %d file(s).', nWritten(i));
            end
            fprintf('[OK] %s/%s/%s | %s\n', subject(i), date(i), session(i), message(i));

        catch ME
            ok(i) = false;
            status(i) = "failed";
            nWritten(i) = 0;
            writtenFiles{i} = strings(0, 1);
            message(i) = string(ME.message);
            fprintf(2, '[FAILED] %s\n', ME.message);
        end

        elapsedSeconds(i) = toc(tStart);
    end

    report = table( ...
        (1:n)', ...
        s2pPaths(:), ...
        repmat(string(ALFroot), n, 1), ...
        ok, ...
        status, ...
        subject, ...
        date, ...
        session, ...
        nCalciumWritten, ...
        nWheelWritten, ...
        nSparseNoiseWritten, ...
        nWritten, ...
        elapsedSeconds, ...
        message, ...
        writtenFiles, ...
        'VariableNames', { ...
            'index', ...
            's2pPath', ...
            'ALFroot', ...
            'ok', ...
            'status', ...
            'subject', ...
            'date', ...
            'session', ...
            'nCalciumWritten', ...
            'nWheelWritten', ...
            'nSparseNoiseWritten', ...
            'nWritten', ...
            'elapsedSeconds', ...
            'message', ...
            'writtenFiles'});

    goodRows = ok & subject ~= "" & date ~= "" & session ~= "";
    sessions = localSessionStructArray(subject(goodRows), date(goodRows), session(goodRows));

    nSkipped = nnz(status == "skipped");
    fprintf('\nALFexporter summary: %d exported, %d skipped, %d failed.\n', ...
        nnz(ok) - nSkipped, nSkipped, nnz(~ok));
    fprintf('ALFexporter sessions for ALFloader: %d\n', numel(sessions));

    if any(~ok)
        warning('ALFexporter:SomeRunsFailed', ...
            '%d of %d Suite2p run(s) failed. Check report.message.', nnz(~ok), n);
    end
end

function sessions = localSessionStructArray(subject, date, session)
    if isempty(subject)
        sessions = localEmptySessionStructArray();
        return
    end

    sessionTable = table( ...
        string(subject(:)), ...
        string(date(:)), ...
        string(session(:)), ...
        'VariableNames', {'subject', 'date', 'session'});
    sessionTable = unique(sessionTable, 'rows', 'stable');

    sessions = table2struct(sessionTable);
end

function sessions = localEmptySessionStructArray()
    sessions = repmat(struct( ...
        'subject', string.empty(0, 0), ...
        'date', string.empty(0, 0), ...
        'session', string.empty(0, 0)), 0, 1);
end

function paths = normalizePathList(value, name)
    if ischar(value)
        paths = string(cellstr(value));
    elseif isstring(value)
        paths = value(:);
    elseif iscell(value)
        try
            paths = string(value(:));
        catch ME
            error('%s must be a char, string array, or cell array of path strings. %s', name, ME.message);
        end
    else
        error('%s must be a char, string array, or cell array of path strings.', name);
    end

    paths = strip(paths);
    paths = paths(paths ~= "");
end

function pathText = normalizeScalarPath(value, name)
    if ischar(value)
        pathText = string(value);
    elseif isstring(value) && isscalar(value)
        pathText = value;
    else
        error('%s must be a char or string scalar.', name);
    end

    pathText = strip(pathText);
    if pathText == ""
        error('%s is empty.', name);
    end

    pathText = char(pathText);
end

function pythonFolder = findPythonExporterFolder()
    thisFile = mfilename('fullpath');
    thisFolder = fileparts(thisFile);

    candidates = strings(0, 1);
    candidates(end + 1, 1) = string(thisFolder);
    candidates(end + 1, 1) = string(pwd);
    candidates(end + 1, 1) = string(fullfile(thisFolder, 'ONEpy'));
    candidates(end + 1, 1) = string(fullfile(pwd, 'ONEpy'));

    parent = string(thisFolder);
    for k = 1:5
        nextParent = string(fileparts(char(parent)));
        if nextParent == parent || nextParent == ""
            break;
        end
        parent = nextParent;
        candidates(end + 1, 1) = parent;
        candidates(end + 1, 1) = string(fullfile(char(parent), 'ONEpy'));
    end

    required = [
        "calcium_export.py"
        "sparse_noise_export.py"
        "wheel_export.py"
    ];

    candidates = unique(candidates, 'stable');
    for i = 1:numel(candidates)
        folder = char(candidates(i));
        if isempty(folder) || ~isfolder(folder)
            continue;
        end

        allPresent = true;
        for j = 1:numel(required)
            if ~isfile(fullfile(folder, char(required(j))))
                allPresent = false;
                break;
            end
        end

        if allPresent
            pythonFolder = folder;
            return;
        end
    end

    error(['Could not find calcium_export.py, sparse_noise_export.py, and wheel_export.py. ', ...
        'Put ALFexporter.m in the same folder as those files, or run MATLAB from that folder.']);
end

function addPythonPath(folder)
    folder = char(folder);
    sys = py.importlib.import_module('sys');
    currentPaths = cell(sys.path);

    alreadyPresent = false;
    for i = 1:numel(currentPaths)
        if strcmp(char(py.str(currentPaths{i})), folder)
            alreadyPresent = true;
            break;
        end
    end

    if ~alreadyPresent
        sys.path.insert(int32(0), folder);
    end
end

function bridge = ensureDirectBridgeModule()
    moduleName = 'alfexporter_direct_bridge';
    bridgeFolder = fullfile(tempdir, 'ALFexporter_pybridge');
    if ~isfolder(bridgeFolder)
        mkdir(bridgeFolder);
    end

    bridgePath = fullfile(bridgeFolder, [moduleName '.py']);
    code = directBridgePythonCode();

    fid = fopen(bridgePath, 'w');
    if fid < 0
        error('Could not write Python bridge module: %s', bridgePath);
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, '%s\n', char(code));
    clear cleaner

    addPythonPath(bridgeFolder);
    importlib = py.importlib.import_module('importlib');
    importlib.invalidate_caches();
    bridge = importlib.import_module(moduleName);
    bridge = importlib.reload(bridge);
end

function value = pyDictString(dictObj, key)
    value = string(char(py.str(dictObj.get(key))));
end

function value = pyDictNumber(dictObj, key)
    value = str2double(char(py.str(dictObj.get(key))));
end

function values = pySequenceToStrings(pySeq)
    try
        items = cell(pySeq);
    catch
        n = double(py.len(pySeq));
        items = cell(n, 1);
        for i = 1:n
            items{i} = pySeq{i};
        end
    end

    values = strings(numel(items), 1);
    for i = 1:numel(items)
        values(i) = string(char(py.str(items{i})));
    end
end

function code = directBridgePythonCode()
    lines = [
"from pathlib import Path"
"import json"
"import numpy as np"
"import calcium_export as cal"
"import wheel_export as wheel"
"import sparse_noise_export as sparse_noise"
""
"def _corrected_data_paths(run_path):"
"    run_db = cal._load_dict_npy(Path(run_path) / 'db.npy', 'run-level Suite2p db.npy')"
"    data_paths = cal._as_list(run_db.get('data_path'))"
"    if hasattr(cal, '_data_paths_with_known_edge_cases'):"
"        data_paths = cal._data_paths_with_known_edge_cases(data_paths, Path(run_path))"
"    return [str(value) for value in data_paths]"
""
"def _is_multi_session_data_path_list(data_paths):"
"    sessions = []"
"    for value in data_paths:"
"        parts = [part for part in str(value).replace('\\', '/').split('/') if part]"
"        if parts:"
"            sessions.append(parts[-1])"
"    return len(set(sessions)) > 1"
""
"def _raw_session_has_file(raw_session, filename):"
"    try:"
"        wheel._find_raw_file(raw_session.local_raw_path, filename)"
"    except FileNotFoundError:"
"        return False"
"    return True"
""
"def _sparse_noise_raw_session(run_path):"
"    run_path = Path(run_path)"
"    run_db = cal._load_dict_npy(run_path / 'db.npy', 'run-level Suite2p db.npy')"
"    raw_sessions = cal._raw_sessions_from_db(run_db, run_path)"
"    matches = [raw_session for raw_session in raw_sessions if _raw_session_has_file(raw_session, 'SparseNoise_Log.bin')]"
"    if not matches:"
"        raise FileNotFoundError('none of the raw sessions contains SparseNoise_Log.bin')"
"    if len(matches) > 1:"
"        raise ValueError('multiple raw sessions contain SparseNoise_Log.bin: ' + str([raw_session.session for raw_session in matches]))"
"    return matches[0]"
""
"def _export_wheel_raw_session(raw_session, out, *, encoder_counts_per_revolution):"
"    encoder_path = wheel._find_raw_file(raw_session.local_raw_path, 'encoder_log.csv')"
"    timestamps, position_counts = wheel._load_unwrapped_encoder_counts(encoder_path)"
"    timestamps = cal._to_session_relative_times("
"        timestamps,"
"        cal._session_time_zero(raw_session.local_raw_path),"
"        description='wheel timestamps from ' + str(encoder_path),"
"    )"
"    if timestamps.size == 0:"
"        raise ValueError('No wheel samples found in ' + str(encoder_path))"
"    if not np.all(np.isfinite(position_counts)):"
"        raise ValueError('Wheel positions contain non-finite values: ' + str(encoder_path))"
"    if timestamps.size != position_counts.size:"
"        raise ValueError('Wheel timestamps and positions differ in length for ' + str(encoder_path) + ': ' + str(timestamps.size) + ' vs ' + str(position_counts.size))"
"    if np.any(np.diff(timestamps) <= 0):"
"        raise ValueError('Wheel timestamps are not strictly increasing; refusing to drop duplicates: ' + str(encoder_path))"
"    scale = 2 * np.pi / float(encoder_counts_per_revolution)"
"    position = position_counts.astype(np.float64) * scale"
"    velocity = wheel._velocity(position, timestamps)"
"    alf_folder = Path(out) / raw_session.subject / raw_session.date / raw_session.session / 'alf'"
"    alf_folder.mkdir(parents=True, exist_ok=True)"
"    return ["
"        wheel._save_npy(alf_folder / '_ibl_wheel.position.npy', position),"
"        wheel._save_npy(alf_folder / '_ibl_wheel.timestamps.npy', timestamps.astype(np.float64)),"
"        wheel._save_npy(alf_folder / '_ibl_wheel.velocity.npy', velocity),"
"    ]"
""
"def _export_sparse_noise_raw_session(raw_session, out, *, stimulus_edges_deg):"
"    sparse_noise_path = wheel._find_raw_file(raw_session.local_raw_path, 'SparseNoise_Log.bin')"
"    times_grid_path = wheel._find_raw_file(raw_session.local_raw_path, 'times_grid.csv')"
"    stim_map = sparse_noise._load_sparse_noise_movie(sparse_noise_path)"
"    frame_times = sparse_noise._load_frame_times(times_grid_path, raw_session.local_raw_path)"
"    stim_map, frame_times = sparse_noise._align_movie_and_times(stim_map, frame_times)"
"    sparse_times, sparse_positions = sparse_noise._sparse_noise_events(stim_map, frame_times, stimulus_edges_deg)"
"    alf_folder = Path(out) / raw_session.subject / raw_session.date / raw_session.session / 'alf'"
"    raw_passive_folder = Path(out) / raw_session.subject / raw_session.date / raw_session.session / 'raw_passive_data'"
"    alf_folder.mkdir(parents=True, exist_ok=True)"
"    raw_passive_folder.mkdir(parents=True, exist_ok=True)"
"    return ["
"        sparse_noise._save_npy(alf_folder / '_ibl_passiveRFM.times.npy', frame_times.astype(np.float64)),"
"        sparse_noise._save_sparse_noise_raw_bin("
"            raw_passive_folder / '_iblrig_RFMapStim.raw.bin',"
"            sparse_noise_path,"
"            expected_shape=stim_map.shape,"
"        ),"
"        sparse_noise._save_npy(alf_folder / '_ibl_sparseNoise.times.npy', sparse_times.astype(np.float64)),"
"        sparse_noise._save_npy(alf_folder / '_ibl_sparseNoise.xy.npy', sparse_positions.astype(np.float64)),"
"    ]"
""
"def _alf_folders_from_written(written):"
"    folders = set()"
"    for path in written:"
"        parts = Path(path).parts"
"        for idx, part in enumerate(parts):"
"            if part == 'alf':"
"                folders.add(Path(*parts[: idx + 1]))"
"                break"
"    return sorted(folders)"
""
"def _load_time_vector(path):"
"    values = np.asarray(np.load(path, allow_pickle=False), dtype=np.float64).reshape(-1)"
"    if values.size == 0:"
"        raise ValueError('Time vector is empty: ' + str(path))"
"    if not np.all(np.isfinite(values)):"
"        raise ValueError('Time vector contains non-finite values: ' + str(path))"
"    if values.size > 1 and np.any(np.diff(values) <= 0):"
"        raise ValueError('Time vector is not strictly increasing: ' + str(path))"
"    return values"
""
"def _time_ranges_overlap(a, b):"
"    return max(float(a[0]), float(b[0])) <= min(float(a[-1]), float(b[-1]))"
""
"def _validate_alf_timebase(alf_folder):"
"    alf_folder = Path(alf_folder)"
"    stim_path = alf_folder / '_ibl_passiveRFM.times.npy'"
"    wheel_path = alf_folder / '_ibl_wheel.timestamps.npy'"
"    fov_time_paths = sorted(alf_folder.glob('FOV_*/mpci.times.npy'))"
"    stim_times = _load_time_vector(stim_path) if stim_path.is_file() else None"
"    wheel_times = _load_time_vector(wheel_path) if wheel_path.is_file() else None"
"    fov_times = [(path.parent.name, _load_time_vector(path)) for path in fov_time_paths]"
"    if stim_times is not None:"
"        for fov_name, mpci_times in fov_times:"
"            if not _time_ranges_overlap(stim_times, mpci_times):"
"                raise ValueError('No overlap between sparse-noise times and ' + fov_name + ' calcium times in ' + str(alf_folder) + ': stim=[' + format(stim_times[0], '.6g') + ', ' + format(stim_times[-1], '.6g') + '], mpci=[' + format(mpci_times[0], '.6g') + ', ' + format(mpci_times[-1], '.6g') + ']')"
"        if wheel_times is not None and not _time_ranges_overlap(stim_times, wheel_times):"
"            raise ValueError('No overlap between sparse-noise times and wheel times in ' + str(alf_folder) + ': stim=[' + format(stim_times[0], '.6g') + ', ' + format(stim_times[-1], '.6g') + '], wheel=[' + format(wheel_times[0], '.6g') + ', ' + format(wheel_times[-1], '.6g') + ']')"
""
"def _save_json(path, value):"
"    path = Path(path)"
"    path.parent.mkdir(parents=True, exist_ok=True)"
"    with path.open('w', encoding='utf-8') as file_obj:"
"        json.dump(value, file_obj, indent=2)"
"        file_obj.write(chr(10))"
"    return path"
""
"def _export_fanci_source_metadata(raw_session, run_path, alf_root):"
"    alf_folder = Path(alf_root) / raw_session.subject / raw_session.date / raw_session.session / 'alf'"
"    return ["
"        _save_json("
"            alf_folder / '_FANCi_source.rawDataPath.json',"
"            {'rawDataPath': str(raw_session.local_raw_path)},"
"        ),"
"        _save_json("
"            alf_folder / '_FANCi_source.s2pSessionPath.json',"
"            {'s2pSessionPath': str(Path(run_path))},"
"        ),"
"    ]"
""
"def _alf_export_present(raw_session, alf_root):"
"    session_folder = Path(alf_root) / raw_session.subject / raw_session.date / raw_session.session"
"    alf_folder = session_folder / 'alf'"
"    raw_passive_folder = session_folder / 'raw_passive_data'"
"    required_alf_files = ["
"        alf_folder / '_ibl_wheel.position.npy',"
"        alf_folder / '_ibl_wheel.timestamps.npy',"
"        alf_folder / '_ibl_wheel.velocity.npy',"
"        alf_folder / '_ibl_passiveRFM.times.npy',"
"        alf_folder / '_ibl_sparseNoise.times.npy',"
"        alf_folder / '_ibl_sparseNoise.xy.npy',"
"    ]"
"    required_raw_files = ["
"        raw_passive_folder / '_iblrig_RFMapStim.raw.bin',"
"    ]"
"    fov_folders = sorted(alf_folder.glob('FOV_*')) if alf_folder.is_dir() else []"
"    has_calcium = any("
"        fov_folder.is_dir()"
"        and (fov_folder / 'mpci.times.npy').is_file()"
"        and (fov_folder / 'mpci.ROIActivityF.npy').is_file()"
"        for fov_folder in fov_folders"
"    )"
"    return has_calcium and all(path.is_file() for path in required_alf_files + required_raw_files)"
""
"def export_sparse_noise_run_direct(run_path, alf_root, *, encoder_counts_per_revolution=65536, sparse_noise_edges_deg=(-135.0, 135.0, -40.0, 40.0)):"
"    run_path = Path(run_path)"
"    alf_root = Path(alf_root)"
"    sparse_noise_edges_deg = tuple(float(x) for x in sparse_noise_edges_deg)"
"    data_paths = _corrected_data_paths(run_path)"
"    raw_session = _sparse_noise_raw_session(run_path)"
"    if _alf_export_present(raw_session, alf_root):"
"        return {"
"            'status': 'skipped',"
"            'subject': raw_session.subject,"
"            'date': raw_session.date,"
"            'session': raw_session.session,"
"            'n_calcium': 0,"
"            'n_wheel': 0,"
"            'n_sparse_noise': 0,"
"            'n_written': 0,"
"            'written': [],"
"            'message': 'Export already present; skipped existing ALF session.',"
"        }"
""
"    calcium_written = list(cal.single_session_calcium_export(run_path, alf_root, stim_type='sparse_noise'))"
""
"    if _is_multi_session_data_path_list(data_paths):"
"        wheel_written = _export_wheel_raw_session("
"            raw_session,"
"            alf_root,"
"            encoder_counts_per_revolution=encoder_counts_per_revolution,"
"        )"
"        sparse_written = _export_sparse_noise_raw_session("
"            raw_session,"
"            alf_root,"
"            stimulus_edges_deg=sparse_noise_edges_deg,"
"        )"
"    else:"
"        wheel_written = list(wheel.single_session_wheel_export("
"            run_path,"
"            alf_root,"
"            encoder_counts_per_revolution=encoder_counts_per_revolution,"
"        ))"
"        sparse_written = list(sparse_noise.single_session_sparse_noise_export("
"            run_path,"
"            alf_root,"
"            stimulus_edges_deg=sparse_noise_edges_deg,"
"        ))"
""
"    source_written = _export_fanci_source_metadata(raw_session, run_path, alf_root)"
"    written = calcium_written + wheel_written + sparse_written + source_written"
"    for alf_folder in _alf_folders_from_written(written):"
"        _validate_alf_timebase(alf_folder)"
""
"    return {"
"        'status': 'ok',"
"        'subject': raw_session.subject,"
"        'date': raw_session.date,"
"        'session': raw_session.session,"
"        'n_calcium': len(calcium_written),"
"        'n_wheel': len(wheel_written),"
"        'n_sparse_noise': len(sparse_written),"
"        'n_written': len(written),"
"        'written': [str(path) for path in written],"
"    }"
    ];

    code = strjoin(lines, sprintf('\n'));
end
