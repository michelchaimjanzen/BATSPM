function batspm_run(T2_path, out_dir, do_stats, batch_mat, tpm_name)
% BATSPM_RUN - Main pipeline wrapper for Bat MRI Analysis
% =========================================================
% USAGE:
%   batspm_run(image, output)                        -> Standard run
%   batspm_run(image, output, true)                  -> Run + Atlas Statistics
%   batspm_run(image, output, false, M)              -> With orientation matrix
%   batspm_run(image, output, false, M, 'tpm.nii')  -> With specific TPM file
%
% batch_mat: optional 3x3 rotation matrix (unit vectors).
%   Found once manually for a protocol by checking one subject in SPM Display.
%   Applied AFTER reslice via spm_get_space (header only, no resampling).
%   Reslice always outputs clean identity Dir Cos first, then M is stamped
%   onto the header so segmentation sees the correct orientation.
%
% tpm_name: optional TPM filename (must be in BATSPM folder).
%   Default: 'Bat_tpm_6ch.nii'

    if nargin < 3, do_stats  = false;              end
    if nargin < 4, batch_mat = [];                 end
    if nargin < 5, tpm_name  = 'Bat_tpm_6ch.nii'; end

    [~, fname, fext] = fileparts(char(T2_path));
    if strcmpi(fext, '.gz') || strcmpi(fext, '.zip')
        error('STOP! Please unzip this file manually. This script only accepts .nii files.');
    end
    if ~strcmpi(fext, '.nii')
        error('STOP! The input file must be a .nii file. You provided: %s', fext);
    end

    fname_clean = regexprep(fname, '[^a-zA-Z0-9]', '');
    bat_dir     = fullfile(out_dir, fname_clean);
    subject_dir = bat_dir;
    if ~exist(bat_dir, 'dir'), mkdir(bat_dir); end
    new_T2_path = fullfile(bat_dir, [fname_clean '.nii']);
    
    if ~exist(new_T2_path, 'file')
        copyfile(T2_path, new_T2_path);
        fprintf('Moved data to workspace: %s\n', bat_dir);
    end
    
    batspm_dir = fileparts(mfilename('fullpath'));
    tpm_path   = fullfile(batspm_dir, tpm_name);
    if ~exist(tpm_path, 'file')
        error('TPM not found: %s', tpm_path);
    end
    fprintf('Using TPM: %s\n', tpm_name);
    V_tpm       = spm_vol(tpm_path);
    num_tissues = length(V_tpm);
    fprintf('TPM has %d tissue channels.\n', num_tissues);

    % --- STEP 4: Reslice ---
    % No batch_mat — always outputs clean identity-like Dir Cos.
    if exist('batspm_reslice_iso', 'file')
        batspm_reslice_iso(new_T2_path, 0.4);
    else
        warning('batspm_reslice_iso not found! Skipping geometry fix.');
    end

    % --- STEP 4.5: Apply orientation (header only, via spm_get_space) ---
    % Runs AFTER reslice. Does NOT resample — just stamps the correct
    % orientation matrix onto the header so SPM segmentation sees the
    % right axis directions. The voxel data stays exactly as resliced.
    if ~isempty(batch_mat)
        fprintf('   -> Applying batch orientation (header only)...\n');
        V_r  = spm_vol(new_T2_path); V_r = V_r(1);
        dims = V_r.dim(1:3);
        
        % Extract the TRUE voxel sizes (preserving negative signs like -0.4)
        original_zooms = sqrt(sum(V_r.mat(1:3,1:3).^2, 1)) .* sign(diag(V_r.mat(1:3,1:3)))';
        
        % Normalise columns just in case caller didn't
        col_norms        = sqrt(sum(batch_mat.^2, 1));
        R                = batch_mat ./ col_norms;
        
        % Build 4x4 affine: multiply rotation by the actual signed voxel sizes
        new_mat          = eye(4);
        new_mat(1:3,1:3) = R * diag(original_zooms);
        new_mat(1:3,4)   = -new_mat(1:3,1:3) * (dims(:) / 2);
        
        spm_get_space(new_T2_path, new_mat);
        
        fprintf('   -> Dir Cos set to:\n');
        fprintf('         [%.3f  %.3f  %.3f]\n', R(1,:));
        fprintf('         [%.3f  %.3f  %.3f]\n', R(2,:));
        fprintf('         [%.3f  %.3f  %.3f]\n', R(3,:));
    end
    % --- STEP 4.8: QC Orientation Check ---
    % Correlates resliced+oriented subject with the atlas.
    % If correlation < threshold, subject is flagged as outlier,
    % logged to outliers.txt in the output root, and pipeline skips
    % segmentation to avoid wasting time on a misoriented scan.
    % ---> EVERYTHINGS COMMENTED OUT FROM HERE
 %   atlas_path   = fullfile(batspm_dir, 'BatlasWithSkull.nii'); 
 %   qc_threshold = 0.15;
 %   if exist('batspm_check_orientation', 'file') && exist(atlas_path, 'file')
 %       [qc_score, is_outlier] = batspm_check_orientation(new_T2_path, atlas_path, qc_threshold);
 %       fprintf('   -> QC orientation score: %.3f (threshold: %.2f)\n', qc_score, qc_threshold);
 %       if is_outlier
 %           fprintf('   -> WARNING: Outlier detected! Skipping segmentation.\n');
 %           % Log to outliers.txt in output root
 %           outlier_log = fullfile(out_dir, 'outliers.txt');
 %           fid = fopen(outlier_log, 'a');
 %           fprintf(fid, '%s\t score=%.3f\t %s\n', fname_clean, qc_score, datestr(now));
 %           fclose(fid);
 %           fprintf('   -> Logged to: %s\n', outlier_log);
 %           return;  % Stop pipeline for this subject
 %       end
 %   else
 %       fprintf('   -> QC check skipped (batspm_check_orientation or atlas not found).\n');
 %   end

    % --- STEP 5: Segmentation ---
    fprintf('Running Bat Segmentation...\n');
    seg = batspm_segment(new_T2_path, tpm_path, num_tissues);
    
    % --- STEP 6: Normalization ---
    norm_res  = [];
    batspm_dir = fileparts(mfilename('fullpath'));
    if exist('batspm_normalize', 'file')
        norm_res = batspm_normalize(seg, subject_dir, batspm_dir);
    end
        
    % --- STEP 7: Jacobian ---
    jac = [];
    if exist('batspm_jacobian', 'file') && ~isempty(norm_res)
        jac = batspm_jacobian(seg, subject_dir);
        if do_stats
            atlas_path = fullfile(batspm_dir, 'Batlas_labels.nii');
            if exist('batspm_atlas_stats', 'file') && exist(atlas_path, 'file')
                batspm_atlas_stats(jac.det, atlas_path, subject_dir, fname_clean);
            else
                warning('Skipping Stats: Missing batspm_atlas_stats.m or Batlas.nii');
            end
        end
    end
    
    % --- STEP 8: Regional Volumes ---
    if exist('batspm_extract_volumes', 'file')
        seg_3class = seg;
        if isfield(seg_3class, 'c4'), seg_3class = rmfield(seg_3class, 'c4'); end
        if isfield(seg_3class, 'c5'), seg_3class = rmfield(seg_3class, 'c5'); end
        if isfield(seg_3class, 'c6'), seg_3class = rmfield(seg_3class, 'c6'); end
        volumes = batspm_extract_volumes(seg_3class, subject_dir);
        if exist('batspm_save_results', 'file')
            batspm_save_results(volumes, fname_clean, subject_dir);
        end
    end
    
    fprintf('Pipeline complete for %s.\n', fname_clean);
end

% =========================================================
% Sub-function: Segmentation
% =========================================================
function seg = batspm_segment(T2_path, tpm_path, num_tissues)
    spm('Defaults','fMRI');
    spm_jobman('initcfg');
    if nargin < 3, num_tissues = 3; end
    fprintf('DEBUG: Looking for file at: %s\n', T2_path);
    if exist(T2_path, 'file') ~= 2
        error('STOP! The file does not exist on the disk here: %s', T2_path);
    else
        fprintf('DEBUG: File found! Proceeding to SPM...\n');
    end
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.preproc.channel.vols     = {T2_path};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg  = 0.00000001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 15;
    matlabbatch{1}.spm.spatial.preproc.channel.write    = [0 1];
    ngaus_per_class = [3, 2, 2, 2, 2, 2];
    for k = 1:num_tissues
        matlabbatch{1}.spm.spatial.preproc.tissue(k).tpm   = {[tpm_path ',' num2str(k)]};
        matlabbatch{1}.spm.spatial.preproc.tissue(k).ngaus = ngaus_per_class(min(k, numel(ngaus_per_class)));
        if k <= 3
            matlabbatch{1}.spm.spatial.preproc.tissue(k).native = [1 0];
        else
            matlabbatch{1}.spm.spatial.preproc.tissue(k).native = [0 0];
        end
        matlabbatch{1}.spm.spatial.preproc.tissue(k).warped = [0 0];
    end
    matlabbatch{1}.spm.spatial.preproc.warp.mrf     = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg     = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg  = '';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm    = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp    = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.write   = [1 1];
    spm_jobman('run', matlabbatch);
    [pth, nam, ext] = fileparts(T2_path);
    seg.T2         = T2_path;
    seg.biascorr   = fullfile(pth, ['m'  nam ext]);
    seg.gm         = fullfile(pth, ['c1' nam ext]);
    seg.wm         = fullfile(pth, ['c2' nam ext]);
    seg.csf        = fullfile(pth, ['c3' nam ext]);
    seg.defforward = fullfile(pth, ['y_'  nam ext]);
    seg.definverse = fullfile(pth, ['iy_' nam ext]);
end

% =========================================================
% Sub-function: Normalization
% =========================================================
function norm_res = batspm_normalize(seg, out_dir, toolbox_dir)
    fprintf('Running Normalization...\n');
    spm('Defaults','fMRI');
    spm_jobman('initcfg');
    def_field      = seg.defforward;
    images_to_warp = {seg.T2; seg.gm; seg.wm; seg.csf};
    tpm_bb = [-9.900, -10.575, -15.150; 10.050, 10.425, 15.000];
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.normalise.write.subj.def      = {def_field};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = images_to_warp;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox  = [0.5 0.5 0.5];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb   = tpm_bb;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    spm_jobman('run', matlabbatch);
    [p, n, e] = fileparts(seg.T2);  norm_res.wT2  = fullfile(p, ['w' n e]);
    [p, n, e] = fileparts(seg.gm);  norm_res.wgm  = fullfile(p, ['w' n e]);
    [p, n, e] = fileparts(seg.wm);  norm_res.wwm  = fullfile(p, ['w' n e]);
    [p, n, e] = fileparts(seg.csf); norm_res.wcsf = fullfile(p, ['w' n e]);
    if nargin > 2 && exist(toolbox_dir, 'dir')
        fprintf('   -> Cleanup: Applying specific tissue masks...\n');
        cleanup_jobs = {
            norm_res.wT2,  'Bat_mask_Brain.nii';
            norm_res.wgm,  'Bat_mask_1.nii';
            norm_res.wwm,  'Bat_mask_2.nii';
            norm_res.wcsf, 'Bat_mask_3.nii';
        };
        for k = 1:size(cleanup_jobs, 1)
            img_path  = cleanup_jobs{k, 1};
            mask_path = fullfile(toolbox_dir, cleanup_jobs{k, 2});
            if exist(img_path, 'file') && exist(mask_path, 'file')
                apply_mask_safe(img_path, mask_path);
            end
        end
        fprintf('      Cleanup complete.\n');
    end
end

function apply_mask_safe(img_path, mask_path)
    V_img  = spm_vol(img_path);  V_img  = V_img(1);
    V_mask = spm_vol(mask_path); V_mask = V_mask(1);
    V_img.dim  = V_img.dim(1:3);
    V_mask.dim = V_mask.dim(1:3);
    if ~isequal(V_img.dim, V_mask.dim) || ~isequal(V_img.mat, V_mask.mat)
        flags.mean = false; flags.which = 1; flags.interp = 0; flags.prefix = 'tmp_';
        spm_reslice({img_path, mask_path}, flags);
        [p,n,e] = fileparts(mask_path);
        mask_path  = fullfile(p, ['tmp_' n e]);
        V_mask     = spm_vol(mask_path); V_mask = V_mask(1);
        V_mask.dim = V_mask.dim(1:3);
        temp_created = true;
    else
        temp_created = false;
    end
    Y_clean = spm_read_vols(V_img) .* (spm_read_vols(V_mask) > 0.1);
    spm_write_vol(V_img, Y_clean);
    if temp_created, delete(mask_path); end
end

% =========================================================
% Sub-function: Jacobian
% =========================================================
function jac = batspm_jacobian(seg, out_dir)
    fprintf('Calculating Jacobian...\n');
    def_field = seg.defforward;
    [def_path, def_name, ~] = fileparts(def_field);
    matlabbatch = [];
    matlabbatch{1}.spm.util.defs.comp{1}.def        = {def_field};
    matlabbatch{1}.spm.util.defs.out{1}.savejac.savedir.saveusr = {def_path};
    try
        spm('Defaults', 'fMRI'); spm_jobman('initcfg'); spm_jobman('run', matlabbatch);
    catch ME
        warning('batspm_jacobian: SPM jobman failed: %s', ME.message);
        jac = []; return;
    end
    final_name = ['j_' def_name '.nii'];
    final_path = fullfile(out_dir, final_name);
    for i = 1:2
        candidates = {fullfile(def_path, final_name), fullfile(pwd, final_name)};
        if exist(candidates{i}, 'file')
            if ~strcmp(candidates{i}, final_path), movefile(candidates{i}, final_path); end
            fprintf('   -> Jacobian saved: %s\n', final_path);
            jac.det = final_path; return;
        end
    end
    warning('batspm_jacobian: Could not locate Jacobian output.'); jac = [];
end

% =========================================================
% Sub-function: Save Results
% =========================================================
function batspm_save_results(volumes, bat_name, out_dir)
    csv_file  = fullfile(out_dir, [bat_name '_volumes.csv']);
    flag_str  = ''; if isfield(volumes, 'flag_summary'), flag_str = volumes.flag_summary; end
    skull_vol = NaN; if isfield(volumes, 'Skull'), skull_vol = volumes.Skull; end
    T = table({bat_name}, volumes.GM, volumes.WM, volumes.CSF, volumes.TIV, skull_vol, {flag_str}, ...
        'VariableNames', {'ID','GM_mm3','WM_mm3','CSF_mm3','TIV_mm3','Skull_mm3','QC_flags'});
    writetable(T, csv_file);
    fprintf('Results saved to: %s\n', csv_file);
end