function run_all_figures(mode)
% RUN_ALL_FIGURES Draw Figures 1--7 of Satouri_Persister_Final.pdf.
% RUN_ALL_FIGURES('recompute') Recalculate the included numerical results.
% Figure 8 uses the author's separate code and is not included here.
if nargin<1, mode='paper'; end
root=fileparts(mfilename('fullpath'));
switch lower(mode)
    case 'paper'
        loaded=load(fullfile(root,'data','paper_data.mat'),'d'); d=loaded.d;
    case 'recompute'
        d=recompute_manuscript_data();
        out=fullfile(root,'results'); if ~exist(out,'dir'), mkdir(out); end
        save(fullfile(out,'all_recomputed.mat'),'d');
    otherwise
        error('Use run_all_figures(''paper'') or run_all_figures(''recompute'').');
end
outdir=fullfile(root,['figures_' lower(mode)]);
for number=1:7, make_figure(number,d,outdir); end
fprintf('Figures 1--7 saved as PDF, PNG and editable FIG files in:\n%s\n',outdir);
fprintf('Figure 8 is handled separately by the author.\n');
end
