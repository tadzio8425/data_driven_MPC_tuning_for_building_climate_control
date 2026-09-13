%% startup.m: Project path configuration
%
%  Run once after opening MATLAB in this folder. It adds the project
%  subdirectories to the path, so that RequestHandler, the MPC step functions
%  and the identification script resolve from anywhere, and sets the working
%  directory to the project root, so that the relative models/ and data/ paths
%  resolve.
%
%  Order of work:
%    1. startup
%    2. identification/ekf_sysid_2C     → models/rc_params_2C6z.mat
%    3. mpc/tracking_mpc_params_build   → models/tracking_mpc_params.mat
%    4. run_baselines                   → the menu of controllers and sweeps

root = fileparts(mfilename('fullpath'));
if isempty(root), root = pwd; end

addpath(fullfile(root, 'boptest'));
addpath(fullfile(root, 'identification'));
addpath(fullfile(root, 'mpc'));

osqp_path = fullfile(root, 'osqp', 'osqp-0.6.2-matlab-windows64');   % prebuilt MEX
if exist(osqp_path, 'dir')
    addpath(osqp_path);
else
    warning('startup:osqp', 'OSQP folder not found: %s', osqp_path);
end

cd(root);

fprintf('─────────────────────────────────────────────\n');
fprintf('  Project paths configured.\n');
fprintf('  Root: %s\n', root);
fprintf('  Added: boptest/  identification/  mpc/  osqp/\n');
fprintf('─────────────────────────────────────────────\n');
