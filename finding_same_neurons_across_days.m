function finding_same_neurons_across_days(animal_ID, root_directory, dates, sessions, reference_date)
% It is an implementaiot nof the method described in "Reactivation predicts
% the consolidation of unbiased long-term cognitive maps"

microns_per_pixel = 835.76/800;
maximal_distance = 12; % cell-pairs that are more than 12 micrometers apart are assumed to be different cells
file_name = 'Fall.mat';
center_outliers_threshold = 10;

% Inputs:
%       animal_ID               a string containing the animal ID
%       root_directory          a string containing the root directory of
%                               the data
%       dates                   a cell containig the dates of all the
%                               recording days
%       sessions                a cell containing the name of the session
%                               recorded in each day
%       reference_date          a string containg the date of the registeration-reference recording
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
days_data = struct;
for day_counter = 1:length(dates)
    days_data.(['day_', num2str(day_counter)]).date = dates{day_counter};
    load_path = fullfile(root_directory,animal_ID, dates{day_counter},sessions{day_counter}, ...
        'suite2p', 'plane0', 'tau_1point5_SparseMode_True_ThresholdScaling_0point9');
    data = load(fullfile(load_path,file_name), 'ops','stat');
    
    save_path = fullfile(root_directory,animal_ID, dates{day_counter});
    
    MeanImg = data.ops.meanImg;
    days_data.(['day_', num2str(day_counter)]).MeanImg = MeanImg;
    if ~strcmp(dates{day_counter},reference_date)
        mytransform = load(fullfile(save_path, 'mytform.mat'));
        mytransform = mytransform.mytform;
        days_data.(['day_', num2str(day_counter)]).MeanImg_registered = ...
            imwarp(MeanImg,mytransform,'OutputView',imref2d(size(MeanImg)),'Interp','cubic');
    else
        days_data.(['day_', num2str(day_counter)]).MeanImg_registered = [];
    end
    
    neuron_pixles_indx = cellfun(@(z) z.ipix, data.stat,'UniformOutput', false);
    clearvars data
    xxx = cell(1,length(neuron_pixles_indx));
    yyy = cell(1,length(neuron_pixles_indx));
    date_temp = dates{day_counter};
    parfor cell_counter = 1:length(neuron_pixles_indx)
        img = zeros(size(MeanImg));
        img = img(:);
        img(neuron_pixles_indx{cell_counter}) = 1;
        img = reshape(img, size(MeanImg));
        img = img';
        if ~strcmp(date_temp,reference_date)
            img_regisered = imwarp(img,mytransform,'OutputView',imref2d(size(img)),'Interp','nearest');
        else
            img_regisered = img;
        end
        [yyy{cell_counter}, xxx{cell_counter}] = find(img_regisered == 1);
    end
    clearvars img img_regisered date_temp
    days_data.(['day_', num2str(day_counter)]).centers_x = cellfun(@(z) round(trimmean(z, center_outliers_threshold)), xxx)';
    days_data.(['day_', num2str(day_counter)]).centers_y = cellfun(@(z) round(trimmean(z, center_outliers_threshold)), yyy)';
    days_data.(['day_', num2str(day_counter)]).neurons_x = xxx;
    days_data.(['day_', num2str(day_counter)]).neurons_y = yyy;
    
end
clearvars day_counter xxx yyy centers_x centers_y mytransform load_path save_path

D_centers = cell(length(dates), length(dates));
D_Jaccard = cell(length(dates), length(dates));
for day_counter_1 = 1:length(dates)-1
    centers_x_day_1 = days_data.(['day_', num2str(day_counter_1)]).centers_x;
    centers_y_day_1 = days_data.(['day_', num2str(day_counter_1)]).centers_y;
    
    for day_counter_2 = day_counter_1 + 1 : length(dates)
        centers_x_day_2 = days_data.(['day_', num2str(day_counter_2)]).centers_x;
        centers_y_day_2 = days_data.(['day_', num2str(day_counter_2)]).centers_y;
        
        % plotting the registeration results
        figure;
        ax1 = subplot(1,2,1);
        if ~strcmp(dates{day_counter_1},reference_date)
            aaa = days_data.(['day_', num2str(day_counter_1)]).MeanImg_registered;
        else
            aaa = days_data.(['day_', num2str(day_counter_1)]).MeanImg;
        end
        
        for ii = 1:length(centers_x_day_1)
            if ~isnan(centers_y_day_1(ii)) && ~isnan(centers_x_day_1(ii))
                aaa(centers_y_day_1(ii), centers_x_day_1(ii)) = max(aaa(:));
            else
                continue
            end
        end
        imagesc(aaa); colormap('gray');
        title(['day ', num2str(day_counter_1)]);
        ax2 = subplot(1,2,2);
        if ~strcmp(dates{day_counter_2},reference_date)
            aaa = days_data.(['day_', num2str(day_counter_2)]).MeanImg_registered;
        else
            aaa = days_data.(['day_', num2str(day_counter_2)]).MeanImg;
        end
        for ii = 1:length(centers_x_day_2)
            if ~isnan(centers_y_day_2(ii)) && ~isnan(centers_x_day_2(ii))
                aaa(centers_y_day_2(ii), centers_x_day_2(ii)) = max(aaa(:));
            else
                continue
            end
        end
        imagesc(aaa); colormap('gray');
        title(['day ', num2str(day_counter_2)]);
        linkaxes([ax1,ax2],'xy');
        clearvars aaa
        savefig(gcf,fullfile(root_directory, animal_ID, ['registeration_days_', num2str(day_counter_1), '_and_', ...
            num2str(day_counter_2), '_', animal_ID, '.fig']),'compact');
        
        % calculating euclidean and Jaccard distances
        D_centers{day_counter_1, day_counter_2} = single(pdist2([centers_x_day_1, centers_y_day_1], ...
            [centers_x_day_2, centers_y_day_2],'euclidean'));
        D_Jaccard_temp = zeros(size(D_centers{day_counter_1, day_counter_2}),'single');
        D_centers_temp = D_centers{day_counter_1, day_counter_2};
        xxx_1 = days_data.(['day_', num2str(day_counter_1)]).neurons_x;
        yyy_1 = days_data.(['day_', num2str(day_counter_1)]).neurons_y;
        xxx_2 = days_data.(['day_', num2str(day_counter_2)]).neurons_x;
        yyy_2 = days_data.(['day_', num2str(day_counter_2)]).neurons_y;
        for cell_counter_1 = 1:size(D_centers_temp,1)
            parfor cell_counter_2 = 1:size(D_centers_temp,2)
                overlap_size = size(intersect([yyy_1{cell_counter_1}, xxx_1{cell_counter_1}], ...
                    [yyy_2{cell_counter_2}, xxx_2{cell_counter_2}], 'rows'),1);
                if ~isempty(overlap_size)
                    Jaccard_coef = overlap_size / (length(yyy_1{cell_counter_1}) + ...
                        length(yyy_2{cell_counter_2}) - overlap_size);
                else
                    Jaccard_coef = 0;
                end
                D_Jaccard_temp(cell_counter_1,cell_counter_2) = 1 - Jaccard_coef;
            end
        end
        D_Jaccard{day_counter_1, day_counter_2} = D_Jaccard_temp;
        
    end
end
clearvars D_Jaccard_temp D_centers_temp xxx_1 yyy_1 xxx_2 yyy_2 centers_x_day_1 centers_y_day_1 centers_x_day_2 ...
    centers_y_day_2 Jaccard_coef
save(fullfile(root_directory, animal_ID, ['finding_same_neurons_across_days_', animal_ID, '_results.mat']), 'days_data', ...
    'D_Jaccard', 'D_centers', 'dates', '-v7.3');

D_Jaccard_temp = [];
for day_counter_1 = 1:length(dates)-1
    for day_counter_2 = day_counter_1 + 1 : length(dates)
        D_centers_temp = D_centers{day_counter_1, day_counter_2};
        D_Jaccard_temp = [D_Jaccard_temp; ...
            reshape(D_Jaccard{day_counter_1, day_counter_2}(D_centers_temp <= round(maximal_distance/microns_per_pixel)),[],1)];
    end
end
D_Jaccard_temp = D_Jaccard_temp(D_Jaccard_temp < 1);
clearvars D_centers_temp

rng('default'); % for reproducibility

phat_beta_init = mle(double(D_Jaccard_temp(D_Jaccard_temp > 0.5)),'distribution','beta');
figure;
histogram(D_Jaccard_temp(D_Jaccard_temp > 0.5), 'BinWidth',0.02, 'Normalization', 'pdf');
hold on;
aaa = pdf('beta', D_Jaccard_temp(D_Jaccard_temp > 0.5),phat_beta_init(1),phat_beta_init(2));
plot(D_Jaccard_temp(D_Jaccard_temp > 0.5),aaa,'.r');


phat_gamma_init = mle(double(D_Jaccard_temp(D_Jaccard_temp < 0.5)),'distribution','gamma');
figure;
histogram(D_Jaccard_temp(D_Jaccard_temp < 0.5), 'BinWidth',0.02, 'Normalization', 'pdf');
hold on;
aaa =  pdf('gamma', D_Jaccard_temp(D_Jaccard_temp < 0.5),phat_gamma_init(1),phat_gamma_init(2));
plot(D_Jaccard_temp(D_Jaccard_temp < 0.5),aaa,'.m');

mixed_model = @(data,k,phat_beta_a,phat_beta_b,phat_gamma_a,phat_gamma_b) k*pdf('beta', data,phat_beta_a,phat_beta_b) + ...
    (1-k)*pdf('gamma', data,phat_gamma_a,phat_gamma_b);
phat_mixed_model = mle(D_Jaccard_temp,'pdf',mixed_model,'start',[0.5, phat_beta_init(1), phat_beta_init(2), ...
    phat_gamma_init(1),phat_gamma_init(2)]);

figure;
histogram(D_Jaccard_temp, 'BinWidth',0.01, 'Normalization', 'pdf');
hold on;
plot(D_Jaccard_temp, mixed_model(D_Jaccard_temp, phat_mixed_model(1),phat_mixed_model(2),phat_mixed_model(3),phat_mixed_model(4),phat_mixed_model(5)),'r.');

prior_beta = phat_mixed_model(1);
prior_gamma = 1 - prior_beta;
phat_beta = [phat_mixed_model(2),phat_mixed_model(3)]; % beta is for different cells
phat_gamma = [phat_mixed_model(4),phat_mixed_model(5)]; % beta is for same cells
posterior_prob = cell(size(D_Jaccard));
posterior_prob_temp = [];
for day_counter_1 = 1:length(dates)-1
    for day_counter_2 = day_counter_1 + 1 : length(dates)
        Jaccard_given_same = pdf('gamma',D_Jaccard{day_counter_1,day_counter_2}(:),phat_gamma(1),phat_gamma(2));
        Jaccard_given_different = pdf('beta',D_Jaccard{day_counter_1,day_counter_2}(:),phat_beta(1),phat_beta(2));
        
        posterior_prob{day_counter_1,day_counter_2} = Jaccard_given_same*prior_gamma ./ (Jaccard_given_same*prior_gamma + ...
            Jaccard_given_different*prior_beta);
        posterior_prob{day_counter_1,day_counter_2} = reshape(posterior_prob{day_counter_1,day_counter_2}, ...
            size(D_Jaccard{day_counter_1,day_counter_2},1),[]);
        posterior_prob{day_counter_1,day_counter_2}(D_centers{day_counter_1,day_counter_2} > round(maximal_distance/microns_per_pixel) | ...
            D_Jaccard{day_counter_1,day_counter_2} >= 1) = 0;
        posterior_prob_temp = [posterior_prob_temp; reshape(posterior_prob{day_counter_1,day_counter_2}...
            (D_centers{day_counter_1,day_counter_2} <= round(maximal_distance/microns_per_pixel) & ...
            D_Jaccard{day_counter_1,day_counter_2} < 1),[],1)];
    end
end

figure;
histogram(posterior_prob_temp, 'BinWidth',0.01, 'Normalization', 'probability');

% constructing the augmented similarity matrix for hirerachical clustering
for day_counter_1 = 1:length(dates)-1
    for day_counter_2 = day_counter_1 + 1 : length(dates)
        posterior_prob{day_counter_2, day_counter_1} = posterior_prob{day_counter_1,day_counter_2}';
    end
end

num_cells_per_day = [];
for day_counter = 1:length(dates)
    if day_counter ~= length(dates)
        posterior_prob{day_counter, day_counter} = eye(size(posterior_prob{day_counter, day_counter+1},1),'single');
    else
        posterior_prob{day_counter, day_counter} = eye(size(posterior_prob{day_counter, day_counter-1},1),'single');
    end
    
    num_cells_per_day = [num_cells_per_day; size(posterior_prob{day_counter, day_counter},1)];
end

posterior_prob_augmented = cell2mat(posterior_prob);
figure; imagesc(posterior_prob_augmented); colormap('gray'); caxis([0.99,1]);

dissimilarity_matrix = 10 - posterior_prob_augmented;
dissimilarity_matrix = dissimilarity_matrix(tril(dissimilarity_matrix,-1) ~= 0);
dissimilarity_matrix = dissimilarity_matrix - 9;
tree = linkage(dissimilarity_matrix');

% figure;
% dendrogram(tree,0);
T = cluster(tree,'cutoff',0.01,'criterion','distance');
figure; histogram(T, 'BinWidth', 1);

[N,~] = histcounts(T, 'BinWidth', 1, 'BinLimits', [min(T), max(T)]);
sum(N > length(dates))
figure; histogram(N, 'BinWidth', 1);

cluster_members = cell(max(T),1);
cluster_members_days = cell(max(T),1);
num_cells_per_day_cumulative = cumsum(num_cells_per_day);
for cluster_counter = 1:max(T)
    cluster_members{cluster_counter} = [];
    cluster_members_days{cluster_counter} = [];
    aaa = find(T == cluster_counter);
    for cell_counter = 1:length(aaa)
        bbb = find(aaa(cell_counter) <=  num_cells_per_day_cumulative, 1);
        cluster_members_days{cluster_counter} = [cluster_members_days{cluster_counter}, bbb];
        if bbb > 1
            cluster_members{cluster_counter} = [cluster_members{cluster_counter}, aaa(cell_counter) - num_cells_per_day_cumulative(bbb-1)];
        else
            cluster_members{cluster_counter} = [cluster_members{cluster_counter}, aaa(cell_counter)];
        end
    end
end
clearvars aaa bbb

save(fullfile(root_directory_save, animal_ID, ['finding_same_neurons_across_days_', animal_ID, '_results.mat']), ...
    'T', 'cluster_members_days', 'cluster_members', '-append');

% visulaizing a set of randomly chosen neurons that are common across all
% the recording days

cluster_num = randperm(length(cluster_members), 10);

for cluster_counter = 1:length(cluster_num)
    figure;
    ax = zeros(1,length(dates));
    for day_counter = 1:length(dates)
        if ~strcmp(dates{day_counter},reference_date)
            aaa = days_data.(['day_', num2str(day_counter)]).MeanImg_registered;
        else
            aaa = days_data.(['day_', num2str(day_counter)]).MeanImg;
        end
        
        ax(day_counter) = subplot(2,3,day_counter);
        imagesc(aaa); colormap('gray'); axis square; hold on; zoom on;
        caxis([100,0.6*max(aaa(:))]);
        title(['day ', num2str(day_counter)]);
    end
    linkaxes(ax,'xy');
    
    cell_counter = 1;
    while cell_counter <= length(cluster_members{cluster_num(cluster_counter)})
        if ~strcmp(days_data.(['day_', num2str(cluster_members_days{cluster_num(cluster_counter)}(cell_counter))]).date,reference_date)
            aaa = days_data.(['day_', num2str(cluster_members_days{cluster_num(cluster_counter)}(cell_counter))]).MeanImg_registered;
        else
            aaa = days_data.(['day_', num2str(cluster_members_days{cluster_num(cluster_counter)}(cell_counter))]).MeanImg;
        end
        centers_x = days_data.(['day_', num2str(cluster_members_days{cluster_num(cluster_counter)}(cell_counter))]).centers_x;
        centers_y = days_data.(['day_', num2str(cluster_members_days{cluster_num(cluster_counter)}(cell_counter))]).centers_y;
        
        subplot(2,3,cluster_members_days{cluster_num(cluster_counter)}(cell_counter));
        imagesc(aaa); colormap('gray'); axis square; hold on;
        title({['day ', num2str(cluster_members_days{cluster_num(cluster_counter)}(cell_counter))], ['cluster ', num2str(cluster_num(cluster_counter))]});
        
        cell_counter_temp = cell_counter;
        while cell_counter_temp <= length(cluster_members{cluster_num(cluster_counter)}) && ...
                cluster_members_days{cluster_num(cluster_counter)}(cell_counter_temp) == cluster_members_days{cluster_num(cluster_counter)}(cell_counter)
            if ~isnan(centers_y(cluster_members{cluster_num(cluster_counter)}(cell_counter_temp))) && ...
                    ~isnan(centers_x(cluster_members{cluster_num(cluster_counter)}(cell_counter_temp)))
                text(centers_x(cluster_members{cluster_num(cluster_counter)}(cell_counter_temp)),centers_y(cluster_members{cluster_num(cluster_counter)}(cell_counter_temp)), '*', ...
                    'Color', [1 0 0], 'FontSize',15);
            end
            cell_counter_temp = cell_counter_temp + 1;
        end
        cell_counter = cell_counter_temp;
    end
    
    
    figure;
    ax = zeros(1,length(dates));
    for day_counter = 1:length(dates)
        if ~strcmp(dates{day_counter},reference_date)
            aaa = days_data.(['day_', num2str(day_counter)]).MeanImg_registered;
        else
            aaa = days_data.(['day_', num2str(day_counter)]).MeanImg;
        end
        
        ax(day_counter) = subplot(2,3,day_counter);
        imagesc(aaa); colormap('gray'); axis square; hold on; zoom on;
        caxis([100,0.6*max(aaa(:))]);
        title(['day ', num2str(day_counter)]);
    end
    linkaxes(ax,'xy');
    
    cell_counter = 1;
    while cell_counter <= length(cluster_members{cluster_counter})
        if ~strcmp(days_data.(['day_', num2str(cluster_members_days{cluster_counter}(cell_counter))]).date,reference_date)
            aaa = days_data.(['day_', num2str(cluster_members_days{cluster_counter}(cell_counter))]).MeanImg_registered;
        else
            aaa = days_data.(['day_', num2str(cluster_members_days{cluster_counter}(cell_counter))]).MeanImg;
        end
        neurons_x = days_data.(['day_', num2str(cluster_members_days{cluster_counter}(cell_counter))]).neurons_x;
        neurons_y = days_data.(['day_', num2str(cluster_members_days{cluster_counter}(cell_counter))]).neurons_y;
        
        
        cell_counter_temp = cell_counter;
        while cell_counter_temp <= length(cluster_members{cluster_counter}) && ...
                cluster_members_days{cluster_counter}(cell_counter_temp) == cluster_members_days{cluster_counter}(cell_counter)
            for pixel_counter = 1:length(neurons_x{cluster_members{cluster_counter}(cell_counter_temp)})
                if ~isnan(neurons_x{cluster_members{cluster_counter}(cell_counter_temp)}(pixel_counter)) && ...
                        ~isnan(neurons_y{cluster_members{cluster_counter}(cell_counter_temp)}(pixel_counter))
                    aaa(neurons_y{cluster_members{cluster_counter}(cell_counter_temp)}(pixel_counter),...
                        neurons_x{cluster_members{cluster_counter}(cell_counter_temp)}(pixel_counter)) = ...
                        max(aaa(:));
                end
            end
            cell_counter_temp = cell_counter_temp + 1;
        end
        
        subplot(2,3,cluster_members_days{cluster_counter}(cell_counter));
        imagesc(aaa); colormap('gray'); axis square; hold on;
        title({['day ', num2str(cluster_members_days{cluster_counter}(cell_counter))], ['cluster ', num2str(cluster_counter)]});
        
        cell_counter = cell_counter_temp;
    end
    
    clearvars neurons_x neurons_y centers_x centers_y aaa
end

end

