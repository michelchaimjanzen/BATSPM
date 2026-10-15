function batspm_reslice_iso(nii_path, target_res, force_z_thick, batch_mat)
% BATSPM_RESLICE_ISO (TPM-ORIGIN VERSION v7)
% 1. Checks for "Swapped Header" (Pencil Effect) and fixes it.
% 2. Finds the Brain boundary for cropping.
% 3. Clamps huge ghosts (>50mm).
% 4. Reslices to isotropic.
% 5. Centers subject at world (0,0,0) to match TPM center.
% 6. If batch_mat provided: output grid and header use that orientation,
%    so voxel data and header always agree. Otherwise uses TPM diagonal signs.
%
% batch_mat: optional 3x3 rotation matrix (unit vectors, no scale).
%   The output image will have this orientation baked into both the
%   voxel data AND the header — they will always match.

    if nargin < 2, target_res = 0.4;  end
    if nargin < 3, force_z_thick = []; end
    if nargin < 4, batch_mat = [];     end

    % --- LOCATE TPM (to read its diagonal signs as fallback) ---
    batspm_dir = fileparts(mfilename('fullpath'));
    tpm_path = fullfile(batspm_dir, 'Bat_tpm_6ch.nii');
    if ~exist(tpm_path, 'file')
        tpm_path = fullfile(batspm_dir, 'Bat_tpm.nii');
    end
    if ~exist(tpm_path, 'file')
        tpm_path = fullfile(batspm_dir, 'TPM.nii');
    end

    if exist(tpm_path, 'file')
        V_tpm     = spm_vol(tpm_path);
        tpm_signs = sign(diag(V_tpm(1).mat(1:3, 1:3)));
        fprintf('   -> TPM axis orientations: [%+d %+d %+d]\n', ...
            tpm_signs(1), tpm_signs(2), tpm_signs(3));
    else
        warning('TPM not found. Using default positive orientation.');
        tpm_signs = [1; 1; 1];
    end

    % --- Determine output rotation matrix ---
    % If batch_mat provided: use it (normalise rows to unit length just in case)
    % Otherwise: fall back to diagonal TPM signs (identity-like)
    if ~isempty(batch_mat)
        % Normalise each column to unit length (in case caller didn't)
        col_norms = sqrt(sum(batch_mat.^2, 1));
        R_out = batch_mat ./ col_norms;
        fprintf('   -> Using provided batch_mat as output orientation.\n');
    else
        R_out = diag(tpm_signs);
        fprintf('   -> Using TPM diagonal signs as output orientation.\n');
    end

    % 1. Load Data
    V = spm_vol(nii_path);
    V = V(1);
    Y = double(spm_read_vols(V));

    fprintf('Processing %s...\n', spm_file(nii_path, 'filename'));

    % --- STEP 1.5: HEADER AUTO-REPAIR (The "Pencil" Fix) ---
    dims      = V.dim(1:3);
    vox_sizes = sqrt(sum(V.mat(1:3,1:3).^2));

    [~, min_dim_idx] = min(dims);
    [~, max_vox_idx] = max(vox_sizes);

    if min_dim_idx ~= max_vox_idx
        fprintf('   -> DETECTED HEADER MISMATCH (Pencil Effect).\n');
        fprintf('      Swapping voxel sizes: Axis %d (was %.2f) <-> Axis %d (was %.2f)\n', ...
            min_dim_idx, vox_sizes(min_dim_idx), max_vox_idx, vox_sizes(max_vox_idx));

        R       = V.mat(1:3, 1:3);
        norm_R  = R ./ vox_sizes;

        new_vox_sizes             = vox_sizes;
        new_vox_sizes(min_dim_idx) = vox_sizes(max_vox_idx);
        new_vox_sizes(max_vox_idx) = vox_sizes(min_dim_idx);

        V.mat(1:3, 1:3) = norm_R .* new_vox_sizes;
        fprintf('      New Physical FOV: %.1f x %.1f x %.1f mm\n', dims .* new_vox_sizes);
    end

    % --- OVERRIDE (Manual force if needed) ---
    if ~isempty(force_z_thick)
        V.mat(3,3) = force_z_thick * sign(V.mat(3,3));
    end

    % 2. FIND THE BRAIN (for cropping boundary only)
    Y_smooth = smooth3(Y, 'box', 3);
    [max_val, max_idx] = max(Y_smooth(:));
    
    % FIX 1: Lower threshold from 0.15 to 0.05 so it doesn't ignore the dimmer frontal lobe
    mask = Y_smooth > (max_val * 0.10);

    CC = bwconncomp(mask);
    found_brain = false;
    for i = 1:CC.NumObjects
        if ismember(max_idx, CC.PixelIdxList{i})
            clean_mask = false(size(mask));
            clean_mask(CC.PixelIdxList{i}) = true;
            found_brain = true;
            break;
        end
    end
    if ~found_brain, clean_mask = mask; end

    % 3. CALCULATE BOUNDING BOX IN VOXEL SPACE
    [ind_x, ind_y, ind_z] = ind2sub(size(clean_mask), find(clean_mask));
    
    % FIX 2: Increase padding from 5 to 15 voxels for a massive safety buffer
    pad     = 10;
    min_ijk = [min(ind_x)-pad, min(ind_y)-pad, min(ind_z)-pad];
    max_ijk = [max(ind_x)+pad, max(ind_y)+pad, max(ind_z)+pad];

    % --- SAFETY CLAMP (Expanded Limit) ---
    current_vox  = sqrt(sum(V.mat(1:3,1:3).^2));
    box_size_mm  = (max_ijk - min_ijk) .* current_vox;
    
    % FIX 3: Increase maximum bounding box size from 50mm to 80mm
    max_allowed_mm = 70;
    
    for d = 1:3
        if box_size_mm(d) > max_allowed_mm
            center_idx   = round((min_ijk(d) + max_ijk(d)) / 2);
            half_width   = round((max_allowed_mm / current_vox(d)) / 2);
            min_ijk(d)   = center_idx - half_width;
            max_ijk(d)   = center_idx + half_width;
        end
    end
    min_ijk = max(min_ijk, 1);
    max_ijk = min(max_ijk, V.dim(1:3));

    % 4. CONVERT BOUNDING BOX CORNERS TO WORLD SPACE
    [crop_x, crop_y, crop_z] = ndgrid(min_ijk(1):max_ijk(1), ...
                                      min_ijk(2):max_ijk(2), ...
                                      min_ijk(3):max_ijk(3));

    corners_x = V.mat(1,1)*crop_x + V.mat(1,2)*crop_y + V.mat(1,3)*crop_z + V.mat(1,4);
    corners_y = V.mat(2,1)*crop_x + V.mat(2,2)*crop_y + V.mat(2,3)*crop_z + V.mat(2,4);
    corners_z = V.mat(3,1)*crop_x + V.mat(3,2)*crop_y + V.mat(3,3)*crop_z + V.mat(3,4);

    % 5. BUILD OUTPUT GRID ALONG R_out AXES
    % Project all corner world coords onto each R_out axis to get extent
    % This ensures the output grid spans the brain along the target axes
    world_corners = [corners_x(:)'; corners_y(:)'; corners_z(:)'];  % 3 x N

    % Project corners onto each output axis
    proj = R_out' * world_corners;   % 3 x N  (coords along each output axis)

    min_proj = min(proj, [], 2);     % 3x1
    max_proj = max(proj, [], 2);     % 3x1

    % Build 1D grids along each output axis
    ax1 = min_proj(1) : target_res : max_proj(1);
    ax2 = min_proj(2) : target_res : max_proj(2);
    ax3 = min_proj(3) : target_res : max_proj(3);

    [g1, g2, g3] = ndgrid(ax1, ax2, ax3);

    % Convert output grid coords back to world space
    % world = R_out * [g1; g2; g3]
    xq_mm = R_out(1,1)*g1 + R_out(1,2)*g2 + R_out(1,3)*g3;
    yq_mm = R_out(2,1)*g1 + R_out(2,2)*g2 + R_out(2,3)*g3;
    zq_mm = R_out(3,1)*g1 + R_out(3,2)*g2 + R_out(3,3)*g3;

    % 6. RESAMPLE: world coords -> source voxel coords -> interpolate
    F     = griddedInterpolant(Y, 'linear', 'none');
    M_inv = inv(V.mat);

    src_x = M_inv(1,1)*xq_mm + M_inv(1,2)*yq_mm + M_inv(1,3)*zq_mm + M_inv(1,4);
    src_y = M_inv(2,1)*xq_mm + M_inv(2,2)*yq_mm + M_inv(2,3)*zq_mm + M_inv(2,4);
    src_z = M_inv(3,1)*xq_mm + M_inv(3,2)*yq_mm + M_inv(3,3)*zq_mm + M_inv(3,4);

    Y_new = F(src_x, src_y, src_z);
    Y_new(isnan(Y_new)) = 0;

    % 7. SAVE: header uses R_out scaled to target_res, centered at origin
    V_new       = V;
    V_new.dim   = size(Y_new);
    V_new.dim   = V_new.dim(1:3);
    V_new.fname = nii_path;

    % Output affine: R_out directions scaled to target_res
    V_new.mat           = eye(4);
    V_new.mat(1:3, 1:3) = R_out * target_res;

    % Center at world (0,0,0)
    subj_center_vox    = V_new.dim(:) / 2;
    V_new.mat(1:3, 4)  = -V_new.mat(1:3,1:3) * subj_center_vox;

    spm_write_vol(V_new, Y_new);

    fprintf('   -> Output orientation: R_out axes.\n');
    fprintf('   -> Subject centered at world (0,0,0). Origin voxel: [%.1f %.1f %.1f]\n', ...
        subj_center_vox);
    fprintf('   -> Success! Resliced to %.1f mm isotropic. (Dims: %d %d %d)\n', ...
        target_res, V_new.dim(1), V_new.dim(2), V_new.dim(3));
end