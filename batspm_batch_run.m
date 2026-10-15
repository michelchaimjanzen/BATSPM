function batspm_batch_run(input_folder, output_root, run_stats, batch_mat, tpm_name)
% BATSPM_BATCH_RUN
% Automates the pipeline for an entire folder of .nii files.
%
% USAGE:
%   batspm_batch_run('Input', 'Results')
%   batspm_batch_run('Input', 'Results', true)
%   batspm_batch_run('Input', 'Results', false, M)
%   batspm_batch_run('Input', 'Results', false, M, 'Bat_tpm_6ch_intensity.nii')
%
% batch_mat: 3x3 rotation matrix found once manually for this protocol.
%   Applied after reslice via spm_get_space (header only, no resampling).
%   Same M applies to all subjects in the batch.
%
% tpm_name: optional TPM filename. Default: 'Bat_tpm_6ch.nii'
%
% WORKFLOW PER BATCH:
%   1. Open one resliced subject in SPM Display
%   2. Reorient until it matches the TPM
%   3. Read the Dir Cos — that's your M
%   4. Run: batspm_batch_run(InputFolder, NewResults, false, M)

if nargin < 3, run_stats = false;              end
if nargin < 4, batch_mat = [];                 end
if nargin < 5, tpm_name  = 'Bat_tpm_6ch.nii'; end

if ~exist(input_folder, 'dir'), error('Input folder not found: %s', input_folder); end
if ~exist(output_root,  'dir'), mkdir(output_root); end

file_list = dir(fullfile(input_folder, '*.nii'));
num_files = length(file_list);

fprintf('---------------------------------------------------\n');
fprintf('BATSPM BATCH MODE (%d subjects)\n', num_files);
fprintf('Atlas Statistics: %s\n', string(run_stats));
fprintf('TPM: %s\n', tpm_name);
if ~isempty(batch_mat)
    fprintf('Orientation matrix: PROVIDED (applied after reslice)\n');
else
    fprintf('Orientation matrix: none\n');
end
fprintf('---------------------------------------------------\n');

for i = 1:num_files
    filename  = file_list(i).name;
    full_path = fullfile(input_folder, filename);
    fprintf('\n>>> Processing %d/%d: %s\n', i, num_files, filename);
    try
        tic;
        batspm_run(full_path, output_root, run_stats, batch_mat, tpm_name);
        fprintf('   [OK] %.1fs\n', toc);
    catch ME
        fprintf('   [ERROR] %s failed: %s\n', filename, ME.message);
    end
end

fprintf('\nBatch Complete.\n');
end