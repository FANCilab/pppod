# masterPipeline_single

Runs export -> deconv RF mapping -> deconv plots -> traces RF mapping -> traces plots
for one subject/date/session (all FOVs in it). Open `masterPipeline_single.m`, edit the
configuration blocks at the top, and run it from MATLAB.

## Setup on a new machine
1. MATLAB with Parallel Computing Toolbox and a supported GPU (the mappers use `gpuArray` and `parfor`),
   plus Optimization Toolbox (`lsqcurvefit`) and Image Processing/Statistics as used by the existing code.
2. Python 3.x with `pip install -r requirements.txt`.
3. Edit these in `masterPipeline_single.m`:
   - `pythonExecutable` - full path to your python.exe (the default is one user's install).
   - `outputRoot` - where ALF exports, results, plots and logs are written (never under the raw/Suite2p roots).
   - `raw2pRoot`, `suite2pRoot` - network paths (read only).
   - `subject`, `date`, `session`.
   - `npyMatlabPath` can be left alone: npy-matlab is bundled in `third_party/npy-matlab`
     and put on the MATLAB path automatically. If the configured folder is missing you get a harmless warning.
4. If MATLAB already loaded a different Python in this session, restart MATLAB.

## Files used
- `masterPipeline_single.m`
- `main_mapRF_ALF_deconv_single.m`, `alfLoaders_deconv.m`
- `main_mapRF_ALF_traces_single.m`, `alfLoaders_traces.m`
- `RFPlot_single.py`
- `ONEpy/` (`rf_master_export_single.py`, `calcium_export.py`, `main_export.py`, `s2p_path_utils.py`,
  `sparse_noise_export.py`, `wheel_export.py`)
- MATLAB packages `+whiteNoise`, `+krnl`, `+colmaps`, `+s2pUtils`
- `third_party/npy-matlab` (BSD licence, from https://github.com/kwikteam/npy-matlab)
