function ridge_regression(img_stk, regressors,srate,Mask )
% Inputs:
%   img_stk: a stack of frames with dimentions pixels(128) by pixels(128) by frames
%   regressors: a matrix whose rows are regressors time-series (e.g. animal facial movement)
%   srate: imaging sampling rate
%   Mask: a binary matrix determining the boundary of the imaging window

% orthogonalizing the regressors
[regressors,~] = qr(regressors',0);
regressors = regressors';

% constructing the desing matrix
tic
kernel_length = round(srate*2);
design_mat = [];
parfor ii = 1:size(regressors,1)
    design_mat = [design_mat, design_matrix(regressors(ii,:),kernel_length,1)];
end
toc
figure; plot(any(design_mat,2));

img_mat = reshape(img_stk,128^2,size(img_stk,3));
clearvars img_stk
NonZeroPixelsIndex = find(reshape(Mask,128*128,1) == 1);
img_mat = img_mat(NonZeroPixelsIndex,:);

[U,S,V] = svd(img_mat, 'econ');
M = 1;
EV = 0;
while EV < 0.95
    EV = sum(diag(S(1:M,1:M)).^2) / sum(diag(S).^2);
    M = M + 1;
end

% analytical solution for ridge regression (the ridgeMML function is written by Matt Kaufman, 2018. mattkaufman@uchicago.edu)
[L, B, convergenceFailures] = ridgeMML((S(1:M,1:M)*V(motion_frames,1:M)')', design_mat, 0);

% cross-validation
PC_count = 300; % or = M
fold_count = 10;
onservation_indx = randperm(size(design_mat,1));
fold_size = floor(length(onservation_indx)/fold_count);
% L = cell(1,fold_count);
beta = cell(1,fold_count);
% convergenceFailures = cell(1,fold_count);
test_onservation_indx = cell(1,fold_count);
V_recon = zeros(size(V,1), PC_count, 'single');

for fold_counter = 1:fold_count
    if fold_counter ~= fold_count
        test_onservation_indx{fold_counter} = onservation_indx((fold_counter-1)*fold_size + (1:fold_size));
    elseif fold_counter == fold_count
        test_onservation_indx{fold_counter} = onservation_indx((fold_counter-1)*fold_size+1:end);
    end
    training_onservation_indx = setdiff(onservation_indx, test_onservation_indx{fold_counter});
    if fold_counter == 1
        [L, beta{fold_counter}, convergenceFailures] = ridgeMML((S(1:PC_count,1:PC_count)*V(training_onservation_indx,1:PC_count)')', ...
                                               design_mat(training_onservation_indx,:), 0);
    else
        [~, beta{fold_counter}] = ridgeMML((S(1:PC_count,1:PC_count)*V(training_onservation_indx,1:PC_count)')', ...
                                               design_mat(training_onservation_indx,:), 0, L);        
    end
    V_recon(test_onservation_indx{fold_counter},:) = design_mat(test_onservation_indx{fold_counter},:)*beta{fold_counter}(2:end,:) + ...
                                                   repmat(beta{fold_counter}(1,:),length(test_onservation_indx{fold_counter}),1);
end

% explained variance
covV = cov((S(1:PC_count,1:PC_count)*V(:,1:PC_count)')');
covV_recon = cov(V_recon);
cCovV = bsxfun(@minus, V_recon, mean(V_recon,1))' * bsxfun(@minus,(S(1:PC_count,1:PC_count)*V(:,1:PC_count)')',...
          mean(S(1:PC_count,1:PC_count)*V(:,1:PC_count)',2)') / (size(V, 1) - 1);
figure; plot(diag(cCovV ./ sqrt(abs(covV .* covV_recon))));
covP = sum((U(:,1:PC_count) * cCovV) .* U(:,1:PC_count), 2);
varP1 = sum((U(:,1:PC_count) * covV) .* U(:,1:PC_count), 2);
varP2 = sum((U(:,1:PC_count) * covV_recon) .* U(:,1:PC_count), 2);
stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5;
corrMat = covP ./ stdPxPy;
mean(corrMat)
std(corrMat)
EV = corrMat .^ 2;
% cross-validated explained_variance and correlation coeffcient maps
EV_map = zeros(128^2,1,'single');
EV_map(NonZeroPixelsIndex,:) = EV;
EV_map = reshape(EV_map,128,128);
figure; imagesc(EV_map);
colormap('jet');

cv_corr_map = zeros(128^2,1,'single');
cv_corr_map(NonZeroPixelsIndex,:) = corrMat;
cv_corr_map = reshape(cv_corr_map,128,128);
figure; imagesc(cv_corr_map);
colormap('jet');

% kernel maps
beta_pixel = U(:,1:PC_count) * Beta(2:end,:)';

for ROI_counter = 1:8
    beta_pixel_reshaped = zeros(128^2,1);
    beta_pixel_reshaped(NonZeroPixelsIndex) = mean(beta_pixel(:,(ROI_counter-1)*kernel_length + (1:kernel_length)),2);
    beta_pixel_reshaped = reshape(beta_pixel_reshaped,128,128);
    figure;
    imagesc(beta_pixel_reshaped);
    title(['ROI ', num2str(ROI_counter)]);
end
% kernel movies
figure;
ROI_counter = 2;
beta_pixel_reshaped = zeros(128^2,1);
for time_counter = 1:kernel_length
    beta_pixel_reshaped(NonZeroPixelsIndex) = beta_pixel(:,(ROI_counter-1)*kernel_length + time_counter);
    beta_pixel_reshaped = reshape(beta_pixel_reshaped,128,128);
    imagesc(beta_pixel_reshaped);
    caxis([-0.005,0.005]);
    colormap('jet');
    title({['ROI ', num2str(ROI_counter)]; ['frame = ', num2str(time_counter)]});
    pause(0.1);    
end
% reconstructing imaging stacks
session_counter = 1;
first_frame = 1 + single(session_counter>1)*sum(nFrames(1:session_counter-1));
last_frame = first_frame + nFrames(session_counter)-1;
V_recon = design_mat(first_frame:last_frame,:)*Beta(2:end,:) + repmat(Beta(1,:),nFrames(session_counter),1);
img_mat_recon = U(:,1:PC_count)*V_recon';
%
figure;
plot(mean(img_mat_recon,1));
%
img_stk_recon = zeros(128^2,nFrames(session_counter),'single');
img_stk_recon(NonZeroPixelsIndex,:) = img_mat_recon;
clearvars img_mat_recon
img_stk_recon = reshape(img_stk_recon,128,128,nFrames(session_counter));

end