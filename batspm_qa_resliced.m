function batspm_qa_resliced(results_dir)
% BATSPM_QA_RESLICED
    if nargin < 1
        results_dir = uigetdir(pwd, 'Select your BATSPM Results Folder');
        if results_dir == 0
            fprintf('QA canceled.\n'); return;
        end
    end

    fprintf('\n=== BATSPM QA: Resliced Viewer ===\n');
    search_pattern = fullfile(results_dir, '**', 'T2*.nii');
    file_list = dir(search_pattern);
    
    if isempty(file_list)
        warning('No resliced T2 files found in %s', results_dir);
        return;
    end
    
    all_files = cell(length(file_list), 1);
    for i = 1:length(file_list)
        all_files{i} = fullfile(file_list(i).folder, file_list(i).name);
    end
    
    num_files = length(all_files);
    batch_size = 6; 
    fprintf('Found %d Resliced T2 files to review.\n', num_files);
    
    for i = 1:batch_size:num_files
        end_idx = min(i + batch_size - 1, num_files);
        current_batch = all_files(i:end_idx);
        
        fprintf('Displaying images %d to %d...\n', i, end_idx);
        spm_check_registration(char(current_batch));
        
        % FIX: Force MATLAB to draw the graphics and bring window to front
        drawnow;
        try figure(spm_figure('GetWin', 'Graphics')); catch; end
        
        if end_idx < num_files
            input('Press [ENTER] to view the next 6 bats (or Ctrl+C to quit)...', 's');
        else
            fprintf('Finished reviewing all files.\n');
        end
    end
end