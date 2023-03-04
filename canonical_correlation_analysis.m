function canonical_correlation_analysis(peri_event_signals_ii, peri_event_signals_jj, before_event,latent_dim)
% canonical correlation analysis for two sets of peri-event signals implemented based
% on the article "Long-term stability of cortical population dynamics underlying
% consistent behavior", NN, 2020

% Inputs:
%   peri_event_signals_ii: neuron by time by trials (events)
%   peri_event_signals_jj: neuron by time by trials (events)
%   before_event: the period before and after events for calculation of
%                   canonical components
%   srate: sampling rate
%   latent_dim: number of canonical components to be considered

after_event = before_event;

midpoint = round(size(peri_event_signals_ii,2)/2);
peri_event_signals_ii = peri_event_signals_ii(:,midpoint+(-round(before_event*camera_srate):round(after_event*camera_srate)),:);
peri_event_signals_ii = bsxfun(@minus, peri_event_signals_ii, mean(peri_event_signals_ii,2));
peri_event_signals_ii = reshape(peri_event_signals_ii, size(peri_event_signals_ii,1), []);

[U_ii,S_ii,V_ii] = svd(peri_event_signals_ii,'econ');

midpoint = round(size(peri_event_signals_jj,2)/2);
peri_event_signals_jj = peri_event_signals_jj(:,midpoint+(-round(before_event*camera_srate):round(after_event*camera_srate)),:);
peri_event_signals_jj = bsxfun(@minus, peri_event_signals_jj, mean(peri_event_signals_jj,2));
peri_event_signals_jj = reshape(peri_event_signals_jj, size(peri_event_signals_jj,1), []);

[U_jj,S_jj,V_jj] = svd(peri_event_signals_jj,'econ');

L_ii = S_ii(1:latent_dim,1:latent_dim)*V_ii(:,1:latent_dim)';
L_jj = S_jj(1:latent_dim,1:latent_dim)*V_jj(:,1:latent_dim)';

[Q_ii,R_ii] = qr(L_ii',0);
Q_ii = single(Q_ii);
R_ii = single(R_ii);

[Q_jj,R_jj] = qr(L_jj',0);
Q_jj = single(Q_jj);
R_jj = single(R_jj);

[U,~,V] = svd(Q_ii(:,1:m)' * Q_jj(:,1:m), 'econ');

M_ii = inv(R_ii(1:m,1:m)) * U;
M_jj = inv(R_jj(1:m,1:m)) * V;

L_ii_aligned = M_ii' * L_ii;
L_jj_aligned = M_jj' * L_jj;

L_ii_to_jj = inv(M_jj)' * M_ii' * L_ii;
L_jj_to_ii = inv(M_ii)' * M_jj' * L_jj;

PC_num = 1;
figure;
yyaxis left
plot(L_ii(PC_num,:));
yyaxis right
plot(L_jj(PC_num,:));

corrcoef(L_ii(PC_num,:), L_jj(PC_num,:))
%%
PC_num = 1;
figure;
yyaxis left
plot(L_ii_aligned(PC_num,:));
yyaxis right
plot(L_jj_aligned(PC_num,:));

corrcoef(L_ii_aligned(PC_num,:), L_jj_aligned(PC_num,:))
%%
PC_num = 1;
figure;
yyaxis left
plot(L_ii(PC_num,:));
yyaxis right
plot(L_jj_to_ii(PC_num,:));

corrcoef(L_ii(PC_num,:), L_jj_to_ii(PC_num,:))
%%
PC_num = 1;
figure;
yyaxis left
plot(L_jj(PC_num,:));
yyaxis right
plot(L_ii_to_jj(PC_num,:));

corrcoef(L_jj(PC_num,:), L_ii_to_jj(PC_num,:))
%%
figure;
plot(mean(U_ii(:,1:m)*L_ii,1));
hold on;
plot(mean(U_ii(:,1:m)*L_jj_to_ii,1),'r');

corrcoef(mean(U_ii(:,1:m)*L_ii,1), mean(U_ii(:,1:m)*L_jj_to_ii,1))
%%
figure;
plot(mean(U_jj(:,1:m)*L_jj,1));
hold on;
plot(mean(U_jj(:,1:m)*L_ii_to_jj,1),'r');

corrcoef(mean(U_jj(:,1:m)*L_jj,1), mean(U_jj(:,1:m)*L_ii_to_jj,1))
%%
aaa = U_ii(:,1:m)*L_ii;
bbb = U_ii(:,1:m)*L_jj_to_ii;
neuron_num = 300;
figure;
plot(aaa(neuron_num,:));
hold on;
plot(bbb(neuron_num,:));

corrcoef(aaa(neuron_num,:), bbb(neuron_num,:))
%%
corr_mat = zscore(aaa,1,2) * zscore(bbb,1,2)' / size(aaa,2);
figure; plot(diag(corr_mat));
figure; histogram(diag(corr_mat));
%%
aaa = U_jj(:,1:m)*L_jj;
bbb = U_jj(:,1:m)*L_ii_to_jj;
neuron_num = 300;
figure;
plot(aaa(neuron_num,:));
hold on;
plot(bbb(neuron_num,:));

corrcoef(aaa(neuron_num,:), bbb(neuron_num,:))
%%
corr_mat = zscore(aaa,1,2) * zscore(bbb,1,2)' / size(aaa,2);
figure; plot(diag(corr_mat));
figure; histogram(diag(corr_mat));
%%
aaa = U_jj(:,1:m)*L_jj;
bbb = U_ii(:,1:m)*L_ii;
neuron_num = 300;
figure;
plot(aaa(neuron_num,:));
hold on;
plot(bbb(neuron_num,:));

corrcoef(aaa(neuron_num,:), bbb(neuron_num,:))
%%
corr_mat = zscore(aaa,1,2) * zscore(bbb,1,2)' / size(aaa,2);
figure; plot(diag(corr_mat));
figure; histogram(diag(corr_mat),'BinWidth',0.1);
%%
corr_mat = zscore(peri_event_signals_ii,1,2) * zscore(peri_event_signals_jj,1,2)' / size(peri_event_signals_ii,2);
figure; plot(diag(corr_mat));
figure; histogram(diag(corr_mat),'BinWidth',0.1);

end
