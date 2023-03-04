function principal_angles_cos = principal_angles_analysis(peri_event_signals_ii, peri_event_signals_jj, before_event, srate, latent_dim)
% principal angles analysis for two sets of peri-event signals implemeted
% based on the article "Cortical population activity within a preserved
% neural manifold underlies multiple motor behaviors", Nature
% communications, 2018

% Inputs:
%   peri_event_signals_ii: neuron by time by trials (events)
%   peri_event_signals_jj: neuron by time by trials (events)
%   before_event: the period before and after events for calculatio of
%                   pricipial angles
%   srate: sampling rate

after_event = before_event;

midpoint = round(size(peri_event_signals_ii,2)/2);
peri_event_signals_ii = peri_event_signals_ii(:,midpoint+(-round(before_event*srate):round(after_event*srate)),:);
peri_event_signals_ii = bsxfun(@minus, peri_event_signals_ii, mean(peri_event_signals_ii,2));
peri_event_signals_ii = reshape(peri_event_signals_ii, size(peri_event_signals_ii,1), []);

midpoint = round(size(peri_event_signals_jj,2)/2);
peri_event_signals_jj = peri_event_signals_jj(:,midpoint+(-round(before_event*srate):round(after_event*srate)),:);
peri_event_signals_jj = bsxfun(@minus, peri_event_signals_jj, mean(peri_event_signals_jj,2));
peri_event_signals_jj = reshape(peri_event_signals_jj, size(peri_event_signals_jj,1), []);

[U_ii,~,~] = svd(peri_event_signals_ii,'econ');
[U_jj,~,~] = svd(peri_event_signals_jj,'econ');

[~, principal_angles_cos ,~] = svd(U_ii(:,1:latent_dim)' * U_jj(:,1:latent_dim), 'econ');

figure;
bar(acosd(diag(principal_angles_cos)));
xlabel('component number');
ylabel('principal angle (degrees)');

end