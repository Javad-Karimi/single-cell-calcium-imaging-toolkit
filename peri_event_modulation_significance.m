function peri_event_modulation_significance(peri_event_signals, before_event, significance_threshold, srate)
% finding modulated neurons around specific events based on Jadhav et al. 2016
% Inputs:
%   peri_event_signals: neuron by time by trials (events)
%   before_event: the period before and after events for determination of
%                   significant modulation
%   significance_threshold: percentile value 
%   srate: sampling rate

rng('default');
after_event = before_event;
permutations_quant = 1000;

events_quant = size(peri_event_signals,3);
cell_quant = size(peri_event_signals,1);
midpoint = round(size(peri_event_signals,2)/2);
sig_shuffled = zeros(size(peri_event_signals,1),size(peri_event_signals,2),permutations_quant, 'single');

parfor permutation_counter = 1:permutations_quant
    sig_temp = zeros(size(peri_event_signals),'single');
    shift = 2.5*before_event*(rand(events_quant,1,'single')-0.5); % up to 500 ms shift
    shift = round(shift*srate);
    for ripple_counter = 1:events_quant
        sig_temp(:,:,ripple_counter) = circshift(peri_event_signals(:,:,ripple_counter),shift(ripple_counter),2);
    end
    sig_shuffled(:,:,permutation_counter) = squeeze(mean(sig_temp,3));
end
clearvars temp_var sig_temp
sig_shuffled_after_events = sig_shuffled(:,midpoint+(0:round(after_event*srate)),:);
sig_shuffled_before_events = sig_shuffled(:,midpoint+(-round(before_event*srate):0),:);

peri_event_signals = squeeze(mean(peri_event_signals,3));
baseline = median(peri_event_signals(:,1:midpoint-round(2*before_event*camera_srate)),2);
sig_after_events = peri_event_signals(:,midpoint+(0:round(after_event*camera_srate)));
sig_before_events = peri_event_signals(:,midpoint+(-round(before_event*camera_srate):0));

modulation_metric_after_events = mean((sig_after_events - squeeze(mean(sig_shuffled_after_events,3))).^2,2);
modulation_metric_before_events = mean((sig_before_events - squeeze(mean(sig_shuffled_before_events,3))).^2,2);

modulation_metric_shuffled_after_events = bsxfun(@minus, sig_shuffled_after_events, squeeze(mean(sig_shuffled_after_events,3)));
modulation_metric_shuffled_after_events = squeeze(mean(modulation_metric_shuffled_after_events.^2, 2));

modulation_metric_shuffled_before_events = bsxfun(@minus, sig_shuffled_before_events, squeeze(mean(sig_shuffled_before_events,3)));
modulation_metric_shuffled_before_events = squeeze(mean(modulation_metric_shuffled_before_events.^2, 2));

modulation_flag_after_events = modulation_metric_after_events > prctile(modulation_metric_shuffled_after_events,significance_threshold,2);
positive_modulation_flag_after_events = (mean(sig_after_events,2) > baseline) & modulation_flag_after_events;
negative_modulation_flag_after_events = (mean(sig_after_events,2) < baseline) & modulation_flag_after_events;

modulation_flag_before_ripples = modulation_metric_before_events > prctile(modulation_metric_shuffled_before_events,significance_threshold,2);
positive_modulation_flag_before_events = (mean(sig_before_events,2) > baseline) & modulation_flag_before_ripples;
negative_modulation_flag_before_events = (mean(sig_before_events,2) < baseline) & modulation_flag_before_ripples;

end