%% unified_4dof_parameterization_pipeline.m
% One robust script for:
%   1) symbolic 4-DOF dynamics derivation (Craig MDH)
%   2) midpoint substitution BEFORE parameterization
%   3) exact monomial linear parameterization tau = Y*Pi - b
%   4) zero-column pruning -> reduced Y, Pi
%   5) sampled-rank/base-parameter check
%
% Prints only:
%   - full Y and Pi
%   - reduced Y and Pi
%   - checker results
%
% Optional:
%   - base regressor Y_base and base parameter vector Theta
%
% Notes:
% - This fixes the upstream issue by rebuilding the regressor from tau_model
%   AFTER midpoint substitution, instead of substituting into an old regressor.
% - Default behavior is display-oriented and does not save files unless enabled.

clear; clc;

%% ============================================================
% 0) User settings
% ============================================================
use_midpoint_for_parameterization = true;   % critical fix
clean_steps                      = 40;

display_full_param               = true;
display_reduced_param            = true;
display_checker                  = true;
display_base_model               = true;

save_results_to_mat              = true;
outdir                           = 'parametric_results';
mat_filename                     = 'unified_4dof_parameterization_results.mat';

rank_tol                         = 1e-9;
num_samples                      = 160;
t0                               = 0;
tf                               = 12;
g_value                          = 9.81;

rational_tol                     = 1e-10;
numeric_zero_tol                 = 1e-9;

if save_results_to_mat && ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% ============================================================
% 1) Symbolic variables
% ============================================================
syms q1 q2 q3 q4 dq1 dq2 dq3 dq4 ddq1 ddq2 ddq3 ddq4 real
syms d1 a2 a3 d4 g real
syms m1 m2 m3 m4 real

syms Ixx1 Iyy1 Izz1 real
syms Ixx2 Iyy2 Izz2 real
syms Ixx3 Iyy3 Izz3 real
syms Ixx4 Iyy4 Izz4 real

syms lc1 lc2 lc3 lc4 real

q   = [q1; q2; q3; q4];
dq  = [dq1; dq2; dq3; dq4];
ddq = [ddq1; ddq2; ddq3; ddq4];
n   = 4;

%% ============================================================
% 2) Craig modified DH data
% ============================================================
a_prev = [ sym(0), sym(0), a2, a3 ];

alpha_prev = sym(zeros(1,4));
alpha_prev(2) = sym(pi)/2;

d_i     = [ d1, sym(0), sym(0), d4 ];
theta_i = [ q1, q2, q3, q4 ];

%% ============================================================
% 3) COM locations and midpoint substitutions
% ============================================================
c = cell(1,n);
c{1} = [0;   0;   lc1];
c{2} = [lc2; 0;   0  ];
c{3} = [lc3; 0;   0  ];
c{4} = [0;   0;   lc4];

mid_subs = [lc1,   lc2,   lc3,   lc4];
mid_vals = [d1/2,  a2/2,  a3/2,  d4/2];

%% ============================================================
% 4) Diagonal inertia tensors
% ============================================================
Ibody = cell(1,n);
Ibody{1} = diag([Ixx1, Iyy1, Izz1]);
Ibody{2} = diag([Ixx2, Iyy2, Izz2]);
Ibody{3} = diag([Ixx3, Iyy3, Izz3]);
Ibody{4} = diag([Ixx4, Iyy4, Izz4]);

m = [m1; m2; m3; m4];

%% ============================================================
% 5) Forward kinematics
% ============================================================
A  = cell(1,n);
T0 = cell(1,n);

T = sym(eye(4));
for i = 1:n
    A{i}  = craigMDH(a_prev(i), alpha_prev(i), d_i(i), theta_i(i));
    T     = T * A{i};
    T0{i} = T;
end

%% ============================================================
% 6) Origins and z-axes
% ============================================================
o = cell(1,n+1);
z = cell(1,n+1);

o{1} = sym([0;0;0]);
z{1} = sym([0;0;1]);

for i = 1:n
    o{i+1} = T0{i}(1:3,4);
    z{i+1} = T0{i}(1:3,3);
end

%% ============================================================
% 7) COM positions and Jacobians
% ============================================================
pC = cell(1,n);
R0 = cell(1,n);
Jv = cell(1,n);
Jw = cell(1,n);

for i = 1:n
    R0{i} = T0{i}(1:3,1:3);
    pC{i} = R0{i} * c{i} + T0{i}(1:3,4);

    Jv{i} = sym(zeros(3,n));
    Jw{i} = sym(zeros(3,n));

    for j = 1:i
        Jw{i}(:,j) = z{j};
        Jv{i}(:,j) = cross(z{j}, pC{i} - o{j});
    end
end

%% ============================================================
% 8) Dynamics: M, C, G, tau
% ============================================================
M = sym(zeros(n,n));
V_total = sym(0);

for i = 1:n
    I0i = R0{i} * Ibody{i} * R0{i}.';
    M = M + m(i) * (Jv{i}.' * Jv{i}) + (Jw{i}.' * I0i * Jw{i});
    V_total = V_total + m(i) * g * pC{i}(3);
end

G = jacobian(V_total, q).';

dM = cell(1,n);
for i = 1:n
    dM{i} = diff(M, q(i));
end

C = sym(zeros(n,n));
for k = 1:n
    for j = 1:n
        ckj = sym(0);
        for i = 1:n
            c_ijk = sym(1)/2 * ( dM{i}(k,j) + dM{j}(k,i) - dM{k}(i,j) );
            ckj   = ckj + c_ijk * dq(i);
        end
        C(k,j) = ckj;
    end
end

tau = M*ddq + C*dq + G;

%% ============================================================
% 9) Choose model used for parameterization
%     IMPORTANT: midpoint is applied BEFORE parameterization
% ============================================================
M_model   = M;
C_model   = C;
G_model   = G;
tau_model = tau;

if use_midpoint_for_parameterization
    M_model   = subs(M_model,   mid_subs, mid_vals);
    C_model   = subs(C_model,   mid_subs, mid_vals);
    G_model   = subs(G_model,   mid_subs, mid_vals);
    tau_model = subs(tau_model, mid_subs, mid_vals);

    phys_atoms = [ ...
        d1; a2; a3; d4; ...
        m1; m2; m3; m4; ...
        Ixx1; Iyy1; Izz1; ...
        Ixx2; Iyy2; Izz2; ...
        Ixx3; Iyy3; Izz3; ...
        Ixx4; Iyy4; Izz4];
else
    phys_atoms = [ ...
        d1; a2; a3; d4; ...
        lc1; lc2; lc3; lc4; ...
        m1; m2; m3; m4; ...
        Ixx1; Iyy1; Izz1; ...
        Ixx2; Iyy2; Izz2; ...
        Ixx3; Iyy3; Izz3; ...
        Ixx4; Iyy4; Izz4];
end

M_model   = cleanSym(M_model,   clean_steps);
C_model   = cleanSym(C_model,   clean_steps);
G_model   = cleanSym(G_model,   clean_steps);
tau_model = cleanSym(tau_model, clean_steps);

%% ============================================================
% 10) Exact monomial parameterization of tau_model
% ============================================================
[Pi_full, Y_full, b_full, check_full] = buildMonomialRegressor(tau_model, phys_atoms, clean_steps);

Pi_full    = cleanSym(Pi_full,    clean_steps);
Y_full     = cleanSym(Y_full,     clean_steps);
b_full     = cleanSym(b_full,     clean_steps);
check_full = cleanSym(check_full, clean_steps);

%% ============================================================
% 11) Zero-column pruning
% ============================================================
keep_col = false(1, size(Y_full,2));
for j = 1:size(Y_full,2)
    keep_col(j) = ~isZeroSymArrayStrong(Y_full(:,j), clean_steps);
end

Y_red   = cleanSym(Y_full(:, keep_col), clean_steps);
Pi_red  = cleanSym(Pi_full(keep_col),   clean_steps);
b_red   = cleanSym(b_full,              clean_steps);
tau_red = tau_model;

check_red = cleanSym(tau_red - (Y_red*Pi_red - b_red), clean_steps);

%% ============================================================
% 12) Checker on reduced parameterization
% ============================================================
dyn_vars = [q; dq; ddq];

reconstruction_ok = isZeroSymArrayStrong(check_red, clean_steps);

b_zero_exact = isZeroSymArrayStrong(b_red, clean_steps);

Pi_vars    = stableUniqueSym(symvar(Pi_red));
Pi_has_dyn = symIntersect(Pi_vars, dyn_vars);
Pi_dyn_ok  = isempty(Pi_has_dyn);

Y_vars            = stableUniqueSym(symvar(Y_red));
known_constants   = sym('g');
Y_has_raw_params  = symIntersect(Y_vars, phys_atoms);
Y_only_dyn_plus_g = isempty(symSetdiff(Y_vars, [dyn_vars; known_constants]));

zero_cols_after_prune = findZeroColumnsStrong(Y_red, clean_steps);
zero_cols_ok          = isempty(zero_cols_after_prune);

%% ============================================================
% 13) Sampled stacked-rank test on reduced Y
% ============================================================
t_samples = linspace(t0, tf, num_samples);

q0   = [ 0.30; -0.80;  0.60;  0.20];
A1   = [ 0.70;  0.60;  0.50;  0.40];
A2   = [ 0.20;  0.15;  0.18;  0.12];
w1   = [ 0.80;  1.10;  0.95;  1.30];
w2   = [ 1.70;  1.40;  1.90;  1.60];
phi1 = [ 0.10;  0.70;  1.10;  1.40];
phi2 = [ 0.30;  0.90;  1.30;  1.80];

extra_Y_syms = symSetdiff(Y_vars, dyn_vars);

extra_vals = nan(numel(extra_Y_syms),1);
for k = 1:numel(extra_Y_syms)
    namek = char(extra_Y_syms(k));
    switch namek
        case 'g'
            extra_vals(k) = g_value;
        otherwise
            error(['Y still contains the symbol "', namek, ...
                '". Add a numeric value for it before rank testing.']);
    end
end

nY = size(Y_red,1);
pY = size(Y_red,2);
Y_stack = zeros(nY*num_samples, pY);

for k = 1:num_samples
    t = t_samples(k);

    qk   = q0 + A1.*sin(w1.*t + phi1) + A2.*cos(w2.*t + phi2);
    dqk  =      A1.*w1.*cos(w1.*t + phi1) - A2.*w2.*sin(w2.*t + phi2);
    ddqk =    - A1.*(w1.^2).*sin(w1.*t + phi1) - A2.*(w2.^2).*cos(w2.*t + phi2);

    subs_syms = [q; dq; ddq; extra_Y_syms];
    subs_vals = [qk; dqk; ddqk; extra_vals];

    Yk = double(subs(Y_red, subs_syms, subs_vals));
    Y_stack((k-1)*nY+1:k*nY, :) = Yk;
end

rank_Y_stack  = rank(Y_stack, rank_tol);
full_col_rank = (rank_Y_stack == pY);
sv            = svd(Y_stack);

%% ============================================================
% 14) Base-parameter reduction
% ============================================================
[~, R, E] = qr(Y_stack, 0);
diagR = abs(diag(R));
if isempty(diagR)
    rank_qr = 0;
else
    rank_qr = sum(diagR > rank_tol);
end

ind_cols = sort(E(1:rank_qr));
dep_cols = setdiff(1:pY, ind_cols);

Y_ind  = Y_red(:, ind_cols);
Pi_ind = Pi_red(ind_cols);

Y_dep  = Y_red(:, dep_cols);
Pi_dep = Pi_red(dep_cols);

W_num                = [];
W_sym                = sym([]);
dep_relation_resid   = sym([]);
dep_relation_ok      = false;

Theta       = sym([]);
Y_base      = Y_ind;
base_check  = sym([]);
base_ok     = false;

if ~isempty(dep_cols)
    Yind_stack = Y_stack(:, ind_cols);
    Ydep_stack = Y_stack(:, dep_cols);

    W_num = Yind_stack \ Ydep_stack;
    W_sym = numericMatrixToSymbolicSimple(W_num, rational_tol);

    dep_relation_resid = cleanSym(Y_dep - Y_ind*W_sym, clean_steps);
    dep_relation_ok    = isZeroSymArrayStrong(dep_relation_resid, clean_steps);

    if dep_relation_ok
        Theta      = cleanSym(Pi_ind + W_sym*Pi_dep, clean_steps);
        base_check = cleanSym(tau_red - (Y_base*Theta - b_red), clean_steps);
        base_ok    = isZeroSymArrayStrong(base_check, clean_steps);
    end
else
    dep_relation_ok = true;
    Theta      = Pi_ind;
    base_check = cleanSym(tau_red - (Y_base*Theta - b_red), clean_steps);
    base_ok    = isZeroSymArrayStrong(base_check, clean_steps);
end

%% ============================================================
% 14.1) Checks on the NEW base parameterization itself
% ============================================================
Theta_vars     = stableUniqueSym(symvar(Theta));
Theta_has_dyn  = symIntersect(Theta_vars, dyn_vars);
Theta_dyn_ok   = isempty(Theta_has_dyn);

Ybase_vars            = stableUniqueSym(symvar(Y_base));
Ybase_has_raw_params  = symIntersect(Ybase_vars, phys_atoms);
Ybase_only_dyn_plus_g = isempty(symSetdiff(Ybase_vars, [dyn_vars; known_constants]));

zero_cols_base = findZeroColumnsStrong(Y_base, clean_steps);
zero_cols_base_ok = isempty(zero_cols_base);

% Sampled stacked regressor for the base model
if isempty(ind_cols)
    Y_base_stack = Y_stack;
else
    Y_base_stack = Y_stack(:, ind_cols);
end

rank_Y_base_stack  = rank(Y_base_stack, rank_tol);
full_col_rank_base = (rank_Y_base_stack == size(Y_base,2));
sv_base            = svd(Y_base_stack);

% Optional conditioning indicators
if isempty(sv_base)
    min_sv_base   = NaN;
    cond_Y_base   = NaN;
else
    min_sv_base = sv_base(end);
    if min_sv_base > 0
        cond_Y_base = sv_base(1) / min_sv_base;
    else
        cond_Y_base = Inf;
    end
end

%% ============================================================
% 15) Print only requested outputs
% ============================================================
if display_full_param
    disp(' ');
    disp('================ FULL PARAMETERIZATION ================');
    disp('================ Pi_full ================');
    disp(Pi_full);
    disp('================ Y_full ================');
    disp(Y_full);
end

if display_reduced_param
    disp(' ');
    disp('================ REDUCED PARAMETERIZATION (ZERO COLUMNS REMOVED) ================');
    disp('================ Pi_reduced ================');
    disp(Pi_red);
    disp('================ Y_reduced ================');
    disp(Y_red);
end

if display_checker
    fprintf('\n============================================================\n');
    fprintf('CHECKER RESULTS\n');
    fprintf('============================================================\n');
    fprintf('Midpoint used before parameterization : %s\n', yesno(use_midpoint_for_parameterization));
    fprintf('Full Y size                           : %d x %d\n', size(Y_full,1), size(Y_full,2));
    fprintf('Reduced Y size                        : %d x %d\n', size(Y_red,1), size(Y_red,2));
    fprintf('------------------------------------------------------------\n');
    fprintf('Exact reconstruction                  : %s\n', passfail(reconstruction_ok));
    fprintf('b == 0 (exact)                        : %s\n', passfail(b_zero_exact));
    fprintf('Pi contains no q,dq,ddq               : %s\n', passfail(Pi_dyn_ok));
    fprintf('Y contains no physical atoms          : %s\n', passfail(isempty(Y_has_raw_params)));
    fprintf('Y only has q,dq,ddq,(known const)     : %s\n', passfail(Y_only_dyn_plus_g));
    fprintf('No zero columns remain                : %s\n', passfail(zero_cols_ok));
    fprintf('------------------------------------------------------------\n');
    fprintf('Stacked Y size                        : %d x %d\n', size(Y_stack,1), size(Y_stack,2));
    fprintf('rank(Y_stack)                         : %d\n', rank_Y_stack);
    fprintf('Number of reduced columns             : %d\n', pY);
    fprintf('Full column rank?                     : %s\n', passfail(full_col_rank));
    fprintf('Independent columns                   : ');
    disp(ind_cols);
    fprintf('Dependent columns                     : ');
    disp(dep_cols);
    fprintf('Exact symbolic Y_dep = Y_ind*W ?      : %s\n', passfail(dep_relation_ok));
    fprintf('Base model exact?                     : %s\n', passfail(base_ok));
    fprintf('============================================================\n');

    fprintf('Theta contains no q,dq,ddq           : %s\n', passfail(Theta_dyn_ok));
    fprintf('Y_base contains no physical atoms    : %s\n', passfail(isempty(Ybase_has_raw_params)));
    fprintf('Y_base only has q,dq,ddq,(const)     : %s\n', passfail(Ybase_only_dyn_plus_g));
    fprintf('No zero columns remain in Y_base     : %s\n', passfail(zero_cols_base_ok));
    fprintf('rank(Y_base_stack)                   : %d\n', rank_Y_base_stack);
    fprintf('Y_base full column rank?             : %s\n', passfail(full_col_rank_base));
    fprintf('min singular value of Y_base_stack   : %.6e\n', min_sv_base);
    fprintf('cond(Y_base_stack)                   : %.6e\n', cond_Y_base);

    fprintf('\nSingular values of stacked Y:\n');
    disp(sv.');

    if display_base_model && dep_relation_ok
        disp(' ');
        disp('================ BASE PARAMETERIZATION ================');
        disp('================ Y_base ================');
        disp(Y_base);
        disp('================ Theta ================');
        disp(Theta);
        disp('================ base check ================');
        disp(base_check);

        fprintf('\n================ COLUMN DEPENDENCY RELATIONS ================\n');
        printColumnRelations(dep_cols, ind_cols, W_sym);
    end
end

%% ============================================================
% 16) Save / workspace
% ============================================================
results_unified = struct();

results_unified.meta.use_midpoint_for_parameterization = use_midpoint_for_parameterization;
results_unified.meta.clean_steps                       = clean_steps;
results_unified.meta.mid_subs                          = mid_subs;
results_unified.meta.mid_vals                          = mid_vals;

results_unified.model.M    = M_model;
results_unified.model.C    = C_model;
results_unified.model.G    = G_model;
results_unified.model.tau  = tau_model;

results_unified.full.Pi    = Pi_full;
results_unified.full.Y     = Y_full;
results_unified.full.b     = b_full;
results_unified.full.check = check_full;

results_unified.reduced.Pi    = Pi_red;
results_unified.reduced.Y     = Y_red;
results_unified.reduced.b     = b_red;
results_unified.reduced.check = check_red;
results_unified.reduced.keep_col = keep_col;

results_unified.check.reconstruction_ok     = reconstruction_ok;
results_unified.check.b_zero_exact          = b_zero_exact;
results_unified.check.Pi_dyn_ok             = Pi_dyn_ok;
results_unified.check.Y_has_raw_params      = Y_has_raw_params;
results_unified.check.Y_only_dyn_plus_g     = Y_only_dyn_plus_g;
results_unified.check.zero_cols_ok          = zero_cols_ok;
results_unified.check.rank_Y_stack          = rank_Y_stack;
results_unified.check.full_col_rank         = full_col_rank;
results_unified.check.singular_values       = sv;
results_unified.check.ind_cols              = ind_cols;
results_unified.check.dep_cols              = dep_cols;
results_unified.check.dep_relation_ok       = dep_relation_ok;
results_unified.check.base_ok               = base_ok;

results_unified.base.Y_base   = Y_base;
results_unified.base.Theta    = Theta;
results_unified.base.W_num    = W_num;
results_unified.base.W_sym    = W_sym;
results_unified.base.check    = base_check;

results_unified.check.Theta_dyn_ok          = Theta_dyn_ok;
results_unified.check.Theta_has_dyn         = Theta_has_dyn;
results_unified.check.Ybase_has_raw_params  = Ybase_has_raw_params;
results_unified.check.Ybase_only_dyn_plus_g = Ybase_only_dyn_plus_g;
results_unified.check.zero_cols_base_ok     = zero_cols_base_ok;
results_unified.check.zero_cols_base        = zero_cols_base;
results_unified.check.rank_Y_base_stack     = rank_Y_base_stack;
results_unified.check.full_col_rank_base    = full_col_rank_base;
results_unified.check.singular_values_base  = sv_base;
results_unified.check.min_sv_base           = min_sv_base;
results_unified.check.cond_Y_base           = cond_Y_base;

assignin('base', 'results_unified', results_unified);
assignin('base', 'Pi_full_new',     Pi_full);
assignin('base', 'Y_full_new',      Y_full);
assignin('base', 'Pi_reduced_new',  Pi_red);
assignin('base', 'Y_reduced_new',   Y_red);
assignin('base', 'Y_base_new',      Y_base);
assignin('base', 'Theta_base_new',  Theta);

if save_results_to_mat
    save(fullfile(outdir, mat_filename), 'results_unified', '-v7.3');
    fprintf('\nSaved unified results to: %s\n', fullfile(outdir, mat_filename));
end

%% ============================================================
% Local function: Craig modified DH transform
% ============================================================
function A = craigMDH(a_prev, alpha_prev, d_i, theta_i)
ct = cos(theta_i);
st = sin(theta_i);
ca = cos(alpha_prev);
sa = sin(alpha_prev);

A = [ ct,      -st,       0,      a_prev;
    st*ca,   ct*ca,   -sa,   -d_i*sa;
    st*sa,   ct*sa,    ca,    d_i*ca;
    0,         0,       0,        1   ];
end

%% ============================================================
% Local function: cleanup
% ============================================================
function expr_out = cleanSym(expr_in, steps)
expr_out = simplify(expr_in, 'Steps', steps);
try
    expr_out = simplifyFraction(expr_out);
catch
end
expr_out = simplify(expr_out, 'Steps', steps);
end

%% ============================================================
% Local function: strong symbolic zero check
% ============================================================
function tf = isZeroSymArrayStrong(expr, steps)
expr = expr(:);
tf = true;
for k = 1:numel(expr)
    ek = simplify(expr(k), 'Steps', steps);
    try
        ek = simplifyFraction(ek);
    catch
    end
    ek = simplify(ek, 'Steps', steps);

    if isAlways(ek == 0, 'Unknown', 'false')
        continue;
    end
    if strcmp(strtrim(char(ek)), '0')
        continue;
    end
    tf = false;
    return;
end
end

%% ============================================================
% Local function: build exact monomial regressor tau = Y*Pi - b
% ============================================================
function [Pi_lin, Y_lin, b_lin, check_lin] = buildMonomialRegressor(tau_vec, phys_vars, clean_steps)

tau_vec = tau_vec(:);

mono_all = sym(zeros(0,1));
for k = 1:numel(tau_vec)
    expr_k = expand(tau_vec(k));
    [~, terms_k] = coeffs(expr_k, phys_vars.');
    terms_k = terms_k(:);

    for t = 1:numel(terms_k)
        mono_all(end+1,1) = terms_k(t); %#ok<AGROW>
    end
end

Pi_lin = stableUniqueSym(mono_all);

keep = true(size(Pi_lin));
for k = 1:numel(Pi_lin)
    if isAlways(Pi_lin(k) == 0, 'Unknown', 'false') || ...
            isAlways(Pi_lin(k) == 1, 'Unknown', 'false')
        keep(k) = false;
    end
end
Pi_lin = Pi_lin(keep);

Pi_lin = sortMonomialsForSub(Pi_lin);

if isempty(Pi_lin)
    Y_lin = sym(zeros(numel(tau_vec),0));
    b_lin = -tau_vec;
    check_lin = simplify(tau_vec - (Y_lin*Pi_lin - b_lin), 'Steps', clean_steps);
    return;
end

p_tmp   = sym('p_tmp_', [numel(Pi_lin), 1], 'real');
tau_sub = tau_vec;

for k = 1:numel(Pi_lin)
    tau_sub = subs(tau_sub, Pi_lin(k), p_tmp(k));
end

[Y_lin, b_lin] = equationsToMatrix(tau_sub, p_tmp);
check_lin = simplify(tau_vec - (Y_lin*Pi_lin - b_lin), 'Steps', clean_steps);
end

%% ============================================================
% Local function: stable unique symbolic vector
% ============================================================
function u = stableUniqueSym(v)
v = v(:);
u = sym(zeros(0,1));
seen = {};

for k = 1:numel(v)
    key = char(v(k));
    if ~any(strcmp(seen, key))
        seen{end+1} = key; %#ok<AGROW>
        u(end+1,1) = v(k); %#ok<AGROW>
    end
end
end

%% ============================================================
% Local function: sort monomials for safe substitution
% ============================================================
function v_sorted = sortMonomialsForSub(v)
if isempty(v)
    v_sorted = v;
    return;
end

score = zeros(numel(v),1);
for k = 1:numel(v)
    s = char(v(k));
    score(k) = 100*numel(symvar(v(k))) + length(s);
end

[~, idx] = sort(score, 'descend');
v_sorted = v(idx);
end

%% ============================================================
% Local function: find zero columns strongly
% ============================================================
function idx = findZeroColumnsStrong(M, steps)
idx = [];
for j = 1:size(M,2)
    if isZeroSymArrayStrong(M(:,j), steps)
        idx(end+1) = j; %#ok<AGROW>
    end
end
end

%% ============================================================
% Local function: symbolic intersection
% ============================================================
function out = symIntersect(a, b)
a = stableUniqueSym(a);
b = stableUniqueSym(b);

out = sym([]);
b_names = cell(numel(b),1);
for i = 1:numel(b)
    b_names{i} = char(b(i));
end

for i = 1:numel(a)
    if any(strcmp(char(a(i)), b_names))
        out(end+1,1) = a(i); %#ok<AGROW>
    end
end
end

%% ============================================================
% Local function: symbolic set difference a \ b
% ============================================================
function out = symSetdiff(a, b)
a = stableUniqueSym(a);
b = stableUniqueSym(b);

out = sym([]);
b_names = cell(numel(b),1);
for i = 1:numel(b)
    b_names{i} = char(b(i));
end

for i = 1:numel(a)
    if ~any(strcmp(char(a(i)), b_names))
        out(end+1,1) = a(i); %#ok<AGROW>
    end
end
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
% Local function: yes/no text
% ============================================================
function s = yesno(tf)
if tf
    s = 'YES';
else
    s = 'NO';
end
end

%% ============================================================
% Local function: convert numeric matrix to simple symbolic constants
% ============================================================
function Wsym = numericMatrixToSymbolicSimple(Wnum, tol)
Wsym = sym(zeros(size(Wnum)));

for i = 1:size(Wnum,1)
    for j = 1:size(Wnum,2)
        x = Wnum(i,j);

        if abs(x) < tol
            Wsym(i,j) = sym(0);
            continue;
        end

        if abs(x - round(x)) < tol
            Wsym(i,j) = sym(round(x));
            continue;
        end

        [n,d] = rat(x, tol);
        if abs(x - n/d) < tol
            Wsym(i,j) = sym(n)/sym(d);
        else
            Wsym(i,j) = sym(x);
        end
    end
end
end

%% ============================================================
% Local function: print dependency relations
% ============================================================
function printColumnRelations(dep_cols, ind_cols, W)
if isempty(dep_cols)
    fprintf('No dependent columns.\n');
    return;
end

for j = 1:numel(dep_cols)
    pieces = {};
    for i = 1:numel(ind_cols)
        cij = W(i,j);
        if isAlways(cij == 0, 'Unknown', 'false')
            continue;
        end

        if isAlways(cij == 1, 'Unknown', 'false')
            pieces{end+1} = sprintf('col%d', ind_cols(i)); %#ok<AGROW>
        elseif isAlways(cij == -1, 'Unknown', 'false')
            pieces{end+1} = sprintf('-col%d', ind_cols(i)); %#ok<AGROW>
        else
            pieces{end+1} = sprintf('(%s)*col%d', char(cij), ind_cols(i)); %#ok<AGROW>
        end
    end

    if isempty(pieces)
        rhs = '0';
    else
        rhs = strjoin(pieces, ' + ');
        rhs = strrep(rhs, '+ -', '- ');
    end

    fprintf('col%d = %s\n', dep_cols(j), rhs);
end
end