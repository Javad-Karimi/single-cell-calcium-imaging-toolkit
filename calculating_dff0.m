function calculating_dff0(animal_ID, root_directory, dates, sessions)

baseline_percentile_win = 32; % 32 seconds
baseline_percentile = 10;
file_name = 'Fall.mat';

% Inputs:
%       animal_ID               a string containing the animal ID
%       root_directory          a string containing the root directory of
%                               the data
%       dates                   a cell containig the dates of all the
%                               recording days
%       sessions                a cell containing the name of the session
%                               recorded in each day
%       reference_date          a string containg the date of the registeration-reference recording

for day_counter = 1:length(dates)
    load_path = fullfile(root_directory,animal_ID, dates{day_counter},sessions{day_counter}, ...
        'suite2p', 'plane0', 'tau_1point5_SparseMode_True_ThresholdScaling_0point9');
    load(fullfile(load_path,file_name), 'F','Fneu', 'ops','spks');
    camera_srate = ops.fs;
    
    save_path = fullfile(root_directory_save,animal_ID, dates{day_counter});
    
    F_corrected = bsxfun(@minus, F, mean(F,2)) - 0.7*bsxfun(@minus, Fneu, mean(Fneu,2)); % nuropil correction
    F_corrected = bsxfun(@plus, F_corrected, mean(F,2));
    clearvars F Fneu
    
    figure;
    neurons = randi(size(F_corrected,1),1,4);
    for neuron_counter = 1:length(neurons)
        subplot(4,1,neuron_counter);
        plot(F_corrected(neurons(neuron_counter),:));
        title(['neuron ', num2str(neurons(neuron_counter))]);
    end
    
    baseline = zeros(size(F_corrected));
    parfor frame_counter = 1:size(F_corrected,2)
        if frame_counter - round(baseline_percentile_win*camera_srate) > 0
            baseline(:,frame_counter) = ...
                prctile(F_corrected(:, frame_counter + (-round(baseline_percentile_win*camera_srate):0)),baseline_percentile,2);
        else
            baseline(:,frame_counter) = prctile(F_corrected(:,1:frame_counter),baseline_percentile,2);
        end
    end
    dff0 = (F_corrected - baseline)./baseline;
    
    dff0_denoised = zeros(size(dff0));
    parfor cell_counter = 1:size(dff0,1)
        dff0_denoised(cell_counter,:) = single(wden(dff0(cell_counter,:),'sqtwolog','s','mln',5,'sym4'));
    end
    dff0_denoised = single(dff0_denoised);
    
    % visualizing the signals calculated for a random set of neurons
    figure;
    neurons = randi(size(dff0_denoised,1),1,12);
    for neuron_counter = 1:length(neurons)
        subplot(3,4,neuron_counter);
        yyaxis left;
        plot(F_corrected(neurons(neuron_counter),:));
        hold on;
        plot(baseline(neurons(neuron_counter),:),'g');
        yyaxis right;
        plot(dff0_denoised(neurons(neuron_counter),:));
        title(['neuron ', num2str(neurons(neuron_counter))]);
    end
    
    figure;
    neurons = randi(size(dff0_denoised,1),1,12);
    for neuron_counter = 1:length(neurons)
        subplot(3,4,neuron_counter);
        plot(dff0(neurons(neuron_counter),:),'b');
        hold on;
        plot(dff0_denoised(neurons(neuron_counter),:),'r');
        title(['neuron ', num2str(neurons(neuron_counter))]);
    end
    
    save(fullfile(save_path, 'dff0.mat'),'dff0', 'camera_srate', '-v7.3');
    save(fullfile(save_path, 'dff0_denoised.mat'),'dff0_denoised', 'camera_srate', '-v7.3');
    save(fullfile(save_path, 'spks.mat'),'spks', 'camera_srate', '-v7.3');
end

end