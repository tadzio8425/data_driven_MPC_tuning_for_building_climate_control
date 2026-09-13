function kpi = load_kpi(filepath)
%LOAD_KPI  Read the BOPTEST KPI struct out of a *_simdata.mat file.
%
%  Returns [] when the file is absent, so the aggregating scripts can report a
%  missing cell instead of failing. The PI and MPC runners nest their logs in a
%  sim struct while the RL and BO runners save theirs at the top level, but
%  every one of them also saves kpi at the top level, which is what is read here.

    kpi = [];
    if ~exist(filepath, 'file'), return; end
    S = load(filepath, 'kpi');
    if isfield(S, 'kpi') && isstruct(S.kpi), kpi = S.kpi; end
end
