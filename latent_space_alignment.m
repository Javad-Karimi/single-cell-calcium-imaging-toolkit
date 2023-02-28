function [s, O] = latent_space_alignment(A_ref, A_nonref, s_ref, s_nonref, B)
% It is an implementation of the algorithm described in the article "Stabilization of a
% brain–computer interface via the alignment of low-dimensional spaces of neural activity"

% Inputs:
% A_ref: factor loading matrix for the reference day
% A_nonref: factor loading matrix for the non-reference day
% s_ref: reference-day-identity of common neurons between the reference and
%        non-reference days
% s_nonref: non-reference-day-identity of common neurons between the reference and
%        non-reference days
% B: the number of neurons used for calculation of the alignment matrix

% Outputs:
% s: a vector containing the cell numbers used in calculation of alignment
%    matrix
% O: alignment matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_ref = A_ref(s_ref,:);
A_nonref = A_nonref(s_nonref,:);

% threshold_1 = prctile(sum(A_ref.^2, 2), 5);
% s_1 = find(sum(A_ref.^2, 2) < threshold_1);
%
% threshold_2 = prctile(sum(A_nonref.^2, 2), 5);
% s_2 = find(sum(A_nonref.^2, 2) < threshold_2);
%
% s = setdiff(1:size(A_ref,1), union(s_1,s_2));
s = 1:size(A_ref,1);


if length(s) < B
    [~,~,transform] = procrustes(A_ref(s,:),A_nonref(s,:),'scaling', false, 'reflection', 'best');
    O = transform.T;
    O = O';
    
    Delta = A_ref(s,:) - A_nonref(s,:) * O';
    disp([num2str(length(s)), ';', 'error = ', num2str(norm(Delta,2))]);
    
end

jj = [];
while length(s) > B
    s = setdiff(s, s(jj));
    % [~, T] = rotatefactors(A_nonref(s,:),'Method','procrustes','Target',A_ref(s,:));
    [~,~,transform] = procrustes(A_ref(s,:),A_nonref(s,:),'scaling', false, 'reflection', 'best');
    O = transform.T;
    % O = T;
    O = O';
    Delta = A_ref(s,:) - A_nonref(s,:) * O';
    [~, jj] = max(sum(Delta.^2,2));
    disp([num2str(s(jj)), ' ; ', num2str(length(s)), ';', 'error = ', num2str(norm(Delta,2))]);
end

end