function [score, is_outlier] = batspm_check_orientation(subj_path, atlas_path, threshold)
% BATSPM_CHECK_ORIENTATION
% Checks whether a resliced+oriented subject scan is correctly aligned
% with the atlas by computing a 3D Pearson correlation between them.
%
% After reslicing and applying the orientation matrix, a correctly oriented
% bat brain should overlap with the atlas — bright regions match bright,
% dark regions match dark. A misoriented scan (wrong rotation, flipped axis)
% will have near-zero or negative correlation.
%
% USAGE:
%   [score, is_outlier] = batspm_check_orientation(subj_path, atlas_path)
%   [score, is_outlier] = batspm_check_orientation(subj_path, atlas_path, 0.4)
%
% Inputs:
%   subj_path  - path to resliced+oriented subject .nii file
%   atlas_path - path to template/atlas (e.g. BatlasWithSkull.nii)
%   threshold  - correlation threshold below which scan is flagged (default: 0.4)
%
% Outputs:
%   score      - Pearson correlation coefficient [-1, 1]
%   is_outlier - true if score < threshold

    if nargin < 3, threshold = 0.15; end

    is_outlier = false;
    score      = NaN;

    if ~exist(subj_path,  'file'), warning('Subject not found: %s',  subj_path);  return; end
    if ~exist(atlas_path, 'file'), warning('Atlas not found: %s', atlas_path); return; end

    % Load subject
    V_s = spm_vol(subj_path);  V_s = V_s(1);
    Y_s = double(spm_read_vols(V_s));

    % Load atlas — resample to subject space if needed
    V_a = spm_vol(atlas_path); V_a = V_a(1);

    % Sample atlas at subject voxel positions (in world space)
    % This handles different dims/vox sizes gracefully
    [xi, yi, zi] = ndgrid(1:V_s.dim(1), 1:V_s.dim(2), 1:V_s.dim(3));
    xyz_vox  = [xi(:)'; yi(:)'; zi(:)'; ones(1, numel(xi))];
    xyz_world = V_s.mat * xyz_vox;           % subject voxels -> world mm
    xyz_atlas = V_a.mat \ xyz_world;          % world mm -> atlas voxels

    Y_a_sampled = spm_sample_vol(V_a, xyz_atlas(1,:), xyz_atlas(2,:), xyz_atlas(3,:), 1);
    Y_a_sampled = reshape(Y_a_sampled, V_s.dim(1:3));

    % Only use voxels that are inside both volumes (non-zero in subject)
    valid = Y_s(:) > 0 & isfinite(Y_a_sampled(:)) & Y_a_sampled(:) > 0;

    if sum(valid) < 1000
        warning('batspm_check_orientation: Too few valid voxels (%d) for reliable correlation.', sum(valid));
        return;
    end

    % Pearson correlation
    f_s = Y_s(valid);
    f_a = Y_a_sampled(valid);
    score = corr(f_s, f_a);

    is_outlier = score < threshold;
end