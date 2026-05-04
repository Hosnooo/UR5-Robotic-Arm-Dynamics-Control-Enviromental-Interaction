%% export_unified_parameterization_latex_v2.m
% Exports LaTeX from results_unified produced by
% unified_4dof_parameterization_pipeline.m
%
% Outputs:
%   - Full parameterization: Pi_full, Y_full
%   - Reduced parameterization: Pi_r, Y_r
%   - Base parameterization: Theta, Y_b
%   - Optional report-ready summary paragraph
%   - Optional dependency relations
%
% Notes:
% - Uses results_unified (new robust pipeline)
% - Handles wide matrices (24+ columns)
% - Uses safer beautification for dq, ddq, trig shorthands
% - Does not save files unless enabled

clear; clc;

%% ============================================================
% 0) User settings
% ============================================================
use_workspace_first      = true;
mat_path                 = fullfile('parametric_results', 'unified_4dof_parameterization_results.mat');

display_latex_blocks     = true;
display_report_summary   = true;
display_dependency_tex   = true;
display_row_equations    = false;

save_tex_files           = true;
outdir                   = 'parametric_results';
snippet_suffix           = '_snippet';

max_matrix_cols_floor    = 80;

scale_full_Y_to_width    = true;
scale_reduced_Y_to_width = true;
scale_base_Y_to_width    = true;
scale_vector_to_width    = false;
scale_roweq_to_width     = true;

% vectors as row-transpose in report
full_pi_as_rowT          = true;
reduced_pi_as_rowT       = true;
theta_as_rowT            = true;

full_pi_name             = '\boldsymbol{\Pi}';
full_Y_name              = '\mathbf{Y}';

reduced_pi_name          = '\boldsymbol{\Pi}_{r}';
reduced_Y_name           = '\mathbf{Y}_{r}';

theta_name               = '\boldsymbol{\Theta}';
base_Y_name              = '\mathbf{Y}_{b}';

%% ============================================================
% 1) Load results_unified
% ============================================================
R = [];

if use_workspace_first
    try
        hasVar = evalin('base', 'exist(''results_unified'',''var'')');
        if hasVar
            R = evalin('base', 'results_unified');
        end
    catch
    end
end

if isempty(R)
    if ~isfile(mat_path)
        error('Could not find results_unified in workspace or MAT-file: %s', mat_path);
    end

    S = load(mat_path);
    if ~isfield(S, 'results_unified')
        error('MAT-file does not contain variable "results_unified".');
    end
    R = S.results_unified;
end

%% ============================================================
% 2) Validate fields
% ============================================================
mustHaveTop = {'full','reduced','base','check'};
for k = 1:numel(mustHaveTop)
    if ~isfield(R, mustHaveTop{k})
        error('results_unified is missing field "%s".', mustHaveTop{k});
    end
end

mustHaveFull = {'Pi','Y'};
for k = 1:numel(mustHaveFull)
    if ~isfield(R.full, mustHaveFull{k})
        error('results_unified.full is missing field "%s".', mustHaveFull{k});
    end
end

mustHaveReduced = {'Pi','Y','b','check'};
for k = 1:numel(mustHaveReduced)
    if ~isfield(R.reduced, mustHaveReduced{k})
        error('results_unified.reduced is missing field "%s".', mustHaveReduced{k});
    end
end

mustHaveBase = {'Theta','Y_base','check'};
for k = 1:numel(mustHaveBase)
    if ~isfield(R.base, mustHaveBase{k})
        error('results_unified.base is missing field "%s".', mustHaveBase{k});
    end
end

%% ============================================================
% 3) Extract
% ============================================================
Pi_full = R.full.Pi(:);
Y_full  = R.full.Y;

Pi_red  = R.reduced.Pi(:);
Y_red   = R.reduced.Y;
b_red   = R.reduced.b(:);
chk_red = R.reduced.check(:);

Theta   = R.base.Theta(:);
Y_base  = R.base.Y_base;
chk_base= R.base.check(:);

if isfield(R.check,'ind_cols')
    ind_cols = R.check.ind_cols;
else
    ind_cols = [];
end

if isfield(R.check,'dep_cols')
    dep_cols = R.check.dep_cols;
else
    dep_cols = [];
end

if isfield(R.base,'W_sym')
    W_sym = R.base.W_sym;
else
    W_sym = sym([]);
end

%% ============================================================
% 4) Build LaTeX blocks
% ============================================================
full_pi_expr   = latexNamedVectorExpr(Pi_full, full_pi_name, full_pi_as_rowT);
full_Y_expr    = latexNamedMatrixExpr(Y_full, full_Y_name);

red_pi_expr    = latexNamedVectorExpr(Pi_red, reduced_pi_name, reduced_pi_as_rowT);
red_Y_expr     = latexNamedMatrixExpr(Y_red, reduced_Y_name);

theta_expr      = latexNamedVectorExpr(Theta, theta_name, theta_as_rowT);
base_Y_expr     = latexNamedMatrixExpr(Y_base, base_Y_name);

b_expr          = latexNamedVectorExpr(b_red, 'b_r', false);
chk_red_expr    = latexNamedVectorExpr(chk_red, '\varepsilon_r', false);
chk_base_expr   = latexNamedVectorExpr(chk_base, '\varepsilon_b', false);

full_pi_tex   = latexDisplayBlock(full_pi_expr,   max(numel(Pi_full),1), max_matrix_cols_floor, scale_vector_to_width);
full_Y_tex    = latexDisplayBlock(full_Y_expr,    size(Y_full,2),        max_matrix_cols_floor, scale_full_Y_to_width);

red_pi_tex    = latexDisplayBlock(red_pi_expr,    max(numel(Pi_red),1),  max_matrix_cols_floor, scale_vector_to_width);
red_Y_tex     = latexDisplayBlock(red_Y_expr,     size(Y_red,2),         max_matrix_cols_floor, scale_reduced_Y_to_width);

theta_tex     = latexDisplayBlock(theta_expr,     max(numel(Theta),1),   max_matrix_cols_floor, scale_vector_to_width);
base_Y_tex    = latexDisplayBlock(base_Y_expr,    size(Y_base,2),        max_matrix_cols_floor, scale_base_Y_to_width);

b_tex         = latexDisplayBlock(b_expr,         1,                     max_matrix_cols_floor, true);
chk_red_tex   = latexDisplayBlock(chk_red_expr,   1,                     max_matrix_cols_floor, true);
chk_base_tex  = latexDisplayBlock(chk_base_expr,  1,                     max_matrix_cols_floor, true);

%% ============================================================
% 5) Row equations (optional)
% ============================================================
row_eq_red = cell(size(Y_red,1),1);
for i = 1:size(Y_red,1)
    row_eq_red{i} = latexDisplayBlock( ...
        latexRowEquationExpr(i, Y_red(i,:), Pi_red, b_red(i), '\tau'), ...
        size(Y_red,2), max_matrix_cols_floor, scale_roweq_to_width);
end

row_eq_base = cell(size(Y_base,1),1);
for i = 1:size(Y_base,1)
    row_eq_base{i} = latexDisplayBlock( ...
        latexRowEquationExpr(i, Y_base(i,:), Theta, sym(0), '\tau'), ...
        size(Y_base,2), max_matrix_cols_floor, scale_roweq_to_width);
end

%% ============================================================
% 6) Report summary and dependency relations
% ============================================================
summary_tex    = buildReportSummaryTex(R);
dependency_tex = buildDependencyRelationsTex(dep_cols, ind_cols, W_sym);

%% ============================================================
% 7) Display
% ============================================================
if display_latex_blocks
    fprintf('\n================ FULL Pi (LaTeX) ================\n%s\n', full_pi_tex);
    fprintf('\n================ FULL Y (LaTeX) ================\n%s\n', full_Y_tex);

    fprintf('\n================ REDUCED Pi_r (LaTeX) ================\n%s\n', red_pi_tex);
    fprintf('\n================ REDUCED Y_r (LaTeX) ================\n%s\n', red_Y_tex);

    fprintf('\n================ BASE Theta (LaTeX) ================\n%s\n', theta_tex);
    fprintf('\n================ BASE Y_b (LaTeX) ================\n%s\n', base_Y_tex);

    fprintf('\n================ b_r (LaTeX) ================\n%s\n', b_tex);
    fprintf('\n================ reduced check (LaTeX) ================\n%s\n', chk_red_tex);
    fprintf('\n================ base check (LaTeX) ================\n%s\n', chk_base_tex);
end

if display_report_summary
    fprintf('\n================ REPORT SUMMARY (LaTeX) ================\n%s\n', summary_tex);
end

if display_dependency_tex
    fprintf('\n================ DEPENDENCY RELATIONS (LaTeX) ================\n%s\n', dependency_tex);
end

if display_row_equations
    for i = 1:numel(row_eq_red)
        fprintf('\n================ Reduced row equation %d (LaTeX) ================\n%s\n', i, row_eq_red{i});
    end
    for i = 1:numel(row_eq_base)
        fprintf('\n================ Base row equation %d (LaTeX) ================\n%s\n', i, row_eq_base{i});
    end
end

%% ============================================================
% 8) Optional save
% ============================================================
if save_tex_files
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    writeTextFile(fullfile(outdir, ['Pi_full', snippet_suffix, '.tex']),       full_pi_tex);
    writeTextFile(fullfile(outdir, ['Y_full', snippet_suffix, '.tex']),        full_Y_tex);
    writeTextFile(fullfile(outdir, ['Pi_reduced', snippet_suffix, '.tex']),    red_pi_tex);
    writeTextFile(fullfile(outdir, ['Y_reduced', snippet_suffix, '.tex']),     red_Y_tex);
    writeTextFile(fullfile(outdir, ['Theta_base', snippet_suffix, '.tex']),    theta_tex);
    writeTextFile(fullfile(outdir, ['Y_base', snippet_suffix, '.tex']),        base_Y_tex);
    writeTextFile(fullfile(outdir, ['b_reduced', snippet_suffix, '.tex']),     b_tex);
    writeTextFile(fullfile(outdir, ['check_reduced', snippet_suffix, '.tex']), chk_red_tex);
    writeTextFile(fullfile(outdir, ['check_base', snippet_suffix, '.tex']),    chk_base_tex);
    writeTextFile(fullfile(outdir, ['report_summary', snippet_suffix, '.tex']), summary_tex);
    writeTextFile(fullfile(outdir, ['dependency_relations', snippet_suffix, '.tex']), dependency_tex);

    fprintf('\nSaved LaTeX snippets in: %s\n', outdir);
end

%% ============================================================
% 9) Put outputs in workspace
% ============================================================
latex_export = struct();

latex_export.full.Pi      = full_pi_tex;
latex_export.full.Y       = full_Y_tex;

latex_export.reduced.Pi   = red_pi_tex;
latex_export.reduced.Y    = red_Y_tex;
latex_export.reduced.b    = b_tex;
latex_export.reduced.chk  = chk_red_tex;

latex_export.base.Theta   = theta_tex;
latex_export.base.Y       = base_Y_tex;
latex_export.base.chk     = chk_base_tex;

latex_export.summary      = summary_tex;
latex_export.dependencies = dependency_tex;
latex_export.row_eq_red   = row_eq_red;
latex_export.row_eq_base  = row_eq_base;

assignin('base', 'latex_export', latex_export);

%% ============================================================
% Local function: named matrix
% ============================================================
function s = latexNamedMatrixExpr(M, nameStr)
[nr, nc] = size(M);
rows = cell(nr,1);

for i = 1:nr
    entries = cell(1,nc);
    for j = 1:nc
        entries{j} = beautifyDynLatex(latex(M(i,j)));
    end
    rows{i} = strjoin(entries, ' & ');
end

body = strjoin(rows, latexRowBreak());
s = sprintf('%s = \\begin{bmatrix}\n%s\n\\end{bmatrix}', nameStr, body);
end

%% ============================================================
% Local function: named vector
% ============================================================
function s = latexNamedVectorExpr(v, nameStr, asRowTranspose)
v = v(:);

if asRowTranspose
    entries = cell(1,numel(v));
    for k = 1:numel(v)
        entries{k} = beautifyDynLatex(latex(v(k)));
    end
    body = strjoin(entries, ' & ');
    s = sprintf('%s = \\begin{bmatrix}\n%s\n\\end{bmatrix}^{\\!\\top}', nameStr, body);
else
    entries = cell(numel(v),1);
    for k = 1:numel(v)
        entries{k} = beautifyDynLatex(latex(v(k)));
    end
    body = strjoin(entries, latexRowBreak());
    s = sprintf('%s = \\begin{bmatrix}\n%s\n\\end{bmatrix}', nameStr, body);
end
end

%% ============================================================
% Local function: unnamed column vector
% ============================================================
function s = latexUnnamedColumnVectorExpr(v)
v = v(:);
entries = cell(numel(v),1);
for k = 1:numel(v)
    entries{k} = beautifyDynLatex(latex(v(k)));
end
body = strjoin(entries, latexRowBreak());
s = sprintf('\\begin{bmatrix}\n%s\n\\end{bmatrix}', body);
end

%% ============================================================
% Local function: unnamed matrix
% ============================================================
function s = latexUnnamedMatrixExpr(M)
[nr, nc] = size(M);
rows = cell(nr,1);

for i = 1:nr
    entries = cell(1,nc);
    for j = 1:nc
        entries{j} = beautifyDynLatex(latex(M(i,j)));
    end
    rows{i} = strjoin(entries, ' & ');
end

body = strjoin(rows, latexRowBreak());
s = sprintf('\\begin{bmatrix}\n%s\n\\end{bmatrix}', body);
end

%% ============================================================
% Local function: row equation
% ============================================================
function s = latexRowEquationExpr(idx, Yrow, Pvec, bscalar, tauName)
lhs  = sprintf('%s_{%d}', tauName, idx);
Yltx = latexUnnamedMatrixExpr(Yrow);
Pltx = latexUnnamedColumnVectorExpr(Pvec);

if isAlways(bscalar == 0, 'Unknown', 'false')
    rhs = sprintf('%s\\,%s', Yltx, Pltx);
else
    rhs = sprintf('%s\\,%s - %s', Yltx, Pltx, beautifyDynLatex(latex(bscalar)));
end

s = sprintf('%s = %s', lhs, rhs);
end

%% ============================================================
% Local function: display block
% ============================================================
function s = latexDisplayBlock(expr, numCols, floorCols, scaleToWidth)
expr = fixLatexRowBreaks(expr);
preamble = sprintf('\\setcounter{MaxMatrixCols}{%d}\n', max(numCols, floorCols));

if scaleToWidth
    body = sprintf('\\[\n\\resizebox{\\linewidth}{!}{$\\displaystyle %s$}\n\\]', expr);
else
    body = sprintf('\\[\n%s\n\\]', expr);
end

s = [preamble, body];
end

%% ============================================================
% Local function: report summary
% ============================================================
function s = buildReportSummaryTex(R)

nfull_r = size(R.full.Y,1);
nfull_c = size(R.full.Y,2);

nred_r  = size(R.reduced.Y,1);
nred_c  = size(R.reduced.Y,2);

nbase_r = size(R.base.Y_base,1);
nbase_c = size(R.base.Y_base,2);

rankY   = R.check.rank_Y_stack;
rankYb  = R.check.rank_Y_base_stack;

s = sprintf([ ...
'\\paragraph{Parameterization summary.}\n' ...
'The reduced 4-DOF symbolic model was parameterized after midpoint substitution. ' ...
'The full regressor has dimensions $%d \\times %d$. ' ...
'After removing zero columns, the reduced regressor has dimensions $%d \\times %d$. ' ...
'The stacked sampled reduced regressor has numerical rank $%d$, so the reduced regressor is not full column rank. ' ...
'After eliminating linear dependence, the base regressor has dimensions $%d \\times %d$ and stacked numerical rank $%d$, which confirms full column rank of the final base model.\n\n' ...
'The exact checks produced: reconstruction = %s, $b=0$ = %s, dependence verification $Y_{\\mathrm{dep}}=Y_{\\mathrm{ind}}W$ = %s, and base-model reconstruction = %s.\n'], ...
nfull_r, nfull_c, ...
nred_r, nred_c, ...
rankY, ...
nbase_r, nbase_c, rankYb, ...
passfail(R.check.reconstruction_ok), ...
passfail(R.check.b_zero_exact), ...
passfail(R.check.dep_relation_ok), ...
passfail(R.check.base_ok));
end

%% ============================================================
% Local function: dependency relations LaTeX
% ============================================================
function s = buildDependencyRelationsTex(dep_cols, ind_cols, W)
if isempty(dep_cols)
    s = 'No dependent columns were detected.';
    return;
end

lines = cell(numel(dep_cols),1);

for j = 1:numel(dep_cols)
    pieces = {};
    for i = 1:numel(ind_cols)
        cij = W(i,j);
        if isAlways(cij == 0, 'Unknown', 'false')
            continue;
        end

        cij_tex = beautifyDynLatex(latex(cij));

        if isAlways(cij == 1, 'Unknown', 'false')
            pieces{end+1} = sprintf('Y_{:, %d}', ind_cols(i)); %#ok<AGROW>
        elseif isAlways(cij == -1, 'Unknown', 'false')
            pieces{end+1} = sprintf('-Y_{:, %d}', ind_cols(i)); %#ok<AGROW>
        else
            pieces{end+1} = sprintf('%s\\,Y_{:, %d}', cij_tex, ind_cols(i)); %#ok<AGROW>
        end
    end

    rhs = strjoin(pieces, ' + ');
    rhs = strrep(rhs, '+ -', '- ');

    lines{j} = sprintf('\\[%s = %s\\]', sprintf('Y_{:, %d}', dep_cols(j)), rhs);
end

s = strjoin(lines, sprintf('\n'));
end

%% ============================================================
% Local function: row break
% ============================================================
function rb = latexRowBreak()
rb = [char(92), char(92), sprintf('\n')];
end

%% ============================================================
% Local function: beautify LaTeX
% ============================================================
function s = beautifyDynLatex(s)
if isstring(s)
    s = char(s);
end

s = regexprep(s, '\\left', '');
s = regexprep(s, '\\right', '');

% ddq and dq
s = regexprep(s, 'ddq([0-9]+)', '\\ddot{q}_{$1}');
s = regexprep(s, 'dq([0-9]+)',  '\\dot{q}_{$1}');
s = regexprep(s, 'q([0-9]+)',   'q_{$1}');

% nicer trig shorthand
s = strrep(s, '\cos(q_{2})', 'c_{2}');
s = strrep(s, '\sin(q_{2})', 's_{2}');
s = strrep(s, '\cos(q_{3})', 'c_{3}');
s = strrep(s, '\sin(q_{3})', 's_{3}');
s = strrep(s, '\cos(q_{2}+q_{3})', 'c_{23}');
s = strrep(s, '\sin(q_{2}+q_{3})', 's_{23}');
s = strrep(s, '\cos(q_{2}+q_{3}+q_{4})', 'c_{234}');
s = strrep(s, '\sin(q_{2}+q_{3}+q_{4})', 's_{234}');

% doubled-angle shorthand
s = strrep(s, '\cos(2 q_{2})', 'c_{2(2)}');
s = strrep(s, '\sin(2 q_{2})', 's_{2(2)}');
s = strrep(s, '\cos(2 q_{2}+q_{3})', 'c_{2q_{2}+q_{3}}');
s = strrep(s, '\sin(2 q_{2}+q_{3})', 's_{2q_{2}+q_{3}}');
s = strrep(s, '\cos(2 q_{2}+2 q_{3})', 'c_{2(23)}');
s = strrep(s, '\sin(2 q_{2}+2 q_{3})', 's_{2(23)}');
s = strrep(s, '\cos(2 q_{2}+2 q_{3}+2 q_{4})', 'c_{2(234)}');
s = strrep(s, '\sin(2 q_{2}+2 q_{3}+2 q_{4})', 's_{2(234)}');

% remove ugly double spaces
s = regexprep(s, ' +', ' ');
end

%% ============================================================
% Local function: fix row breaks
% ============================================================
function s = fixLatexRowBreaks(s)
if isstring(s)
    s = char(s);
end

nl   = sprintf('\n');
good = [char(92), char(92), nl];
bad  = [char(92), nl];

token = '__ROWBREAK_TOKEN__';
while contains(s, token)
    token = [token, '_X'];
end

s = strrep(s, good, token);
s = strrep(s, bad, good);
s = strrep(s, token, good);
end

%% ============================================================
% Local function: pass/fail text
% ============================================================
function s = passfail(tf)
if tf
    s = 'PASS';
else
    s = 'FAIL';
end
end

%% ============================================================
% Local function: write text file
% ============================================================
function writeTextFile(filename, txt)
fid = fopen(filename, 'w');
if fid == -1
    error('Could not open file for writing: %s', filename);
end
fprintf(fid, '%s', txt);
fclose(fid);
end