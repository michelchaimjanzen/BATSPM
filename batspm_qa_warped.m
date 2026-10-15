function batspm_qa_warped(results_dir)
% BATSPM_QA_WARPED
    if nargin < 1
        results_dir = uigetdir(pwd, 'Select your BATSPM Results Folder');
        if results_dir == 0
            fprintf('QA canceled.\n'); return;
        end
    end

    fprintf('\n=== BATSPM QA: Warped Viewer ===\n');
    search_pattern = fullfile(results_dir, '**', 'w*.nii');
    file_list = dir(search_pattern);
    
    if isempty(file_list)
        warning('No normalized wT2 files found in %s', results_dir);
        return;
    end
    
    all_files = cell(length(file_list), 1);
    for i = 1:length(file_list)
        all_files{i} = fullfile(file_list(i).folder, file_list(i).name);
    end
    
    num_files = length(all_files);
    batch_size = 6; 
    fprintf('Found %d Normalized wT2 files to review.\n', num_files);
    
    for i = 1:batch_size:num_files
        end_idx = min(i + batch_size - 1, num_files);
        current_batch = all_files(i:end_idx);
        
        fprintf('Displaying images %d to %d...\n', i, end_idx);
        
        % 1. Find or create the SPM Graphics window
        fg = spm_figure('GetWin', 'Graphics');
        
        % 2. Load the 6 images into Check Reg
        spm_check_registration(char(current_batch));
        
        % 3. Force the SPM window to the absolute front of your screen
        figure(fg);
        drawnow;
        
        if end_idx < num_files
            input('Press [ENTER] in this Command Window to view the next 6 bats...', 's');
        else
            fprintf('Finished reviewing all files.\n');
        end
    end
    ends