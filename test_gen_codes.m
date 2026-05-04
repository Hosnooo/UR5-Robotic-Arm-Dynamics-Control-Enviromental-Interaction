clc; clear; close all;
addpath('maple_gen')

%% Test settings
% Use a fixed seed so the same random states can be reproduced.
rng(1);

% Run this many random states through the checks.
NSAMPLES = 20;

% Run the timing section at the end.
DO_TIMING = true;

% Repeat each timed call and report the median.
NREPS = 3;

% Run the skew-symmetry test for Mdot - 2C.
CHECK_SKEW = true;

% Run h checks if UR5_h exists.
CHECK_H = (exist('UR5_h','file') == 2);

% Run forward/inverse dynamics round-trip if UR5_fdyn exists.
CHECK_FDYN = (exist('UR5_fdyn','file') == 2);

% Sample states inside these limits.
q_lim   = pi   * ones(6,1);
dq_lim  = 1.5  * ones(6,1);
ddq_lim = 3.0  * ones(6,1);

%% Tolerances
% These tolerances are tight enough to catch real issues but loose enough
% to avoid false failures from floating-point noise.
tol.fk_lastrow = 1e-12;
tol.R_orth     = 1e-10;
tol.detR       = 1e-10;
tol.M_sym      = 1e-9;
tol.M_pd       = 1e-10;
tol.C_zero     = 1e-12;
tol.h_zero     = 1e-12;
tol.tau        = 1e-8;
tol.g          = 1e-10;
tol.h          = 1e-8;
tol.roundtrip  = 1e-8;
tol.skew       = 1e-5;
tol.power      = 1e-7;

% Use this finite-difference step for Mdot.
epsMdot = 1e-6;

%% Load robot parameters
[p, Pi] = UR5_params();
n = 6;

fprintf('============================================================\n');
fprintf('UR5 symbolic dynamics validation\n');
fprintf('Samples: %d\n', NSAMPLES);
fprintf('CHECK_H: %d | CHECK_FDYN: %d | CHECK_SKEW: %d\n', CHECK_H, CHECK_FDYN, CHECK_SKEW);
fprintf('============================================================\n\n');

%% Zero-state smoke test
% This section confirms that every main function runs, returns the right
% sizes, and produces finite numbers at q = 0, dq = 0, ddq = 0.
q   = zeros(n,1);
dq  = zeros(n,1);
ddq = zeros(n,1);

T0e = UR5_fkine(q, Pi);
[Jv, Jw, Jg] = UR5_jacobian_geometric(q, Pi);
[M, C, G] = UR5_MCG(q, dq, Pi);
tau = UR5_idyn(q, dq, ddq, Pi);

assert(all(size(T0e) == [4 4]), 'T0e must be 4x4');
assert(all(size(Jv)  == [3 6]), 'Jv must be 3x6');
assert(all(size(Jw)  == [3 6]), 'Jw must be 3x6');
assert(all(size(Jg)  == [6 6]), 'Jg must be 6x6');
assert(all(size(M)   == [6 6]), 'M must be 6x6');
assert(all(size(C)   == [6 6]), 'C must be 6x6');
assert(all(size(G)   == [6 1]), 'G must be 6x1');
assert(all(size(tau) == [6 1]), 'tau must be 6x1');

assert(all(isfinite(T0e(:))), 'T0e contains non-finite values');
assert(all(isfinite(Jv(:))),  'Jv contains non-finite values');
assert(all(isfinite(Jw(:))),  'Jw contains non-finite values');
assert(all(isfinite(Jg(:))),  'Jg contains non-finite values');
assert(all(isfinite(M(:))),   'M contains non-finite values');
assert(all(isfinite(C(:))),   'C contains non-finite values');
assert(all(isfinite(G(:))),   'G contains non-finite values');
assert(all(isfinite(tau(:))), 'tau contains non-finite values');

%% Track the worst error seen across all random tests
maxerr.fk_lastrow = 0;
maxerr.R_orth     = 0;
maxerr.detR       = 0;
maxerr.M_sym      = 0;
maxerr.C_zero     = 0;
maxerr.h_zero     = 0;
maxerr.tau        = 0;
maxerr.g          = 0;
maxerr.h          = 0;
maxerr.roundtrip  = 0;
maxerr.skew       = 0;
maxerr.power      = 0;

% Track the smallest eigenvalue of M and the state that caused it.
minEigM = inf;
q_minEig  = nan(n,1);
dq_minEig = nan(n,1);
ddq_minEig = nan(n,1);

%% Randomized validation loop
for k = 1:NSAMPLES
    % Draw a random joint position, velocity, and acceleration.
    q   = rand_range_vec(-q_lim,   q_lim);
    dq  = rand_range_vec(-dq_lim,  dq_lim);
    ddq = rand_range_vec(-ddq_lim, ddq_lim);

    %% Kinematics checks
    % Check that the forward kinematics transform is a valid homogeneous transform.
    T0e = UR5_fkine(q, Pi);
    [Jv, Jw, Jg] = UR5_jacobian_geometric(q, Pi);

    assert(all(isfinite(T0e(:))), 'T0e has non-finite values');
    assert(all(isfinite(Jv(:))),  'Jv has non-finite values');
    assert(all(isfinite(Jw(:))),  'Jw has non-finite values');
    assert(all(isfinite(Jg(:))),  'Jg has non-finite values');

    R = T0e(1:3,1:3);

    % The last row of a homogeneous transform must be [0 0 0 1].
    maxerr.fk_lastrow = max(maxerr.fk_lastrow, norm(T0e(4,:) - [0 0 0 1], inf));

    % A valid rotation matrix must satisfy R'R = I.
    maxerr.R_orth = max(maxerr.R_orth, norm(R.'*R - eye(3), 'fro'));

    % A valid rotation matrix must satisfy det(R) = +1.
    maxerr.detR = max(maxerr.detR, abs(det(R) - 1));

    %% Dynamics calls
    M = UR5_M(q, Pi);
    C = UR5_C(q, dq, Pi);
    G = UR5_G(q, Pi);
    tau_id = UR5_idyn(q, dq, ddq, Pi);

    assert(all(isfinite(M(:))),      'M has non-finite values');
    assert(all(isfinite(C(:))),      'C has non-finite values');
    assert(all(isfinite(G(:))),      'G has non-finite values');
    assert(all(isfinite(tau_id(:))), 'tau_id has non-finite values');

    %% Inertia matrix checks
    % M must be symmetric.
    maxerr.M_sym = max(maxerr.M_sym, norm(M - M.', 'fro'));

    % M must be positive definite.
    % The smallest eigenvalue must stay positive at every tested state.
    Msym = 0.5*(M + M.');
    ev = eig(Msym);
    thisMinEig = min(ev);

    if thisMinEig < minEigM
        minEigM = thisMinEig;
        q_minEig = q;
        dq_minEig = dq;
        ddq_minEig = ddq;
    end

    %% Coriolis matrix zero-velocity check
    % C(q,0) must be zero for the standard Coriolis matrix definition used in
    % tau = M*ddq + C*dq + G.
    C0 = UR5_C(q, zeros(n,1), Pi);
    maxerr.C_zero = max(maxerr.C_zero, norm(C0, 'fro'));

    %% Inverse dynamics consistency check
    % UR5_idyn must match the direct formula tau = M*ddq + C*dq + G.
    tau_ref = M*ddq + C*dq + G;
    maxerr.tau = max(maxerr.tau, norm(tau_id - tau_ref, inf));

    %% Gravity consistency check
    % At zero velocity and zero acceleration, inverse dynamics must return gravity only.
    tau_g = UR5_idyn(q, zeros(n,1), zeros(n,1), Pi);
    maxerr.g = max(maxerr.g, norm(tau_g - G, inf));

    %% h consistency check
    % h must match C*dq.
    % h must be zero when dq = 0.
    if CHECK_H
        h = UR5_h(q, dq, Pi);
        h0 = UR5_h(q, zeros(n,1), Pi);

        assert(all(isfinite(h(:))), 'h has non-finite values');

        maxerr.h = max(maxerr.h, norm(h - C*dq, inf));
        maxerr.h_zero = max(maxerr.h_zero, norm(h0, inf));
    end

    %% Forward/inverse dynamics round-trip check
    % If tau comes from inverse dynamics, forward dynamics must recover ddq.
    if CHECK_FDYN
        qdd_fd = UR5_fdyn(q, dq, tau_id, Pi);
        assert(all(isfinite(qdd_fd(:))), 'qdd_fd has non-finite values');

        maxerr.roundtrip = max(maxerr.roundtrip, norm(qdd_fd - ddq, inf));
    end

    %% Skew-symmetry check
    % For a correct rigid-body model, Mdot - 2C must be skew-symmetric.
    % The scalar dq'*(Mdot - 2C)*dq must also be zero.
    if CHECK_SKEW
        Mdot = approx_Mdot(q, dq, Pi, epsMdot);
        N = Mdot - 2*C;

        maxerr.skew = max(maxerr.skew, norm(N + N.', 'fro'));
        maxerr.power = max(maxerr.power, abs(dq.' * N * dq));
    end
end

%% Summary
fprintf('-------------------- Summary --------------------\n');
print_result('FK last row [0 0 0 1]',      maxerr.fk_lastrow, tol.fk_lastrow);
print_result('Rotation orthogonality',     maxerr.R_orth,     tol.R_orth);
print_result('det(R) = 1',                 maxerr.detR,       tol.detR);
print_result('M symmetry',                 maxerr.M_sym,      tol.M_sym);
print_result('min eig(M) > 0',             minEigM,           tol.M_pd, 'greater');
print_result('C(q,0) = 0',                 maxerr.C_zero,     tol.C_zero);
print_result('tau = M*ddq + C*dq + G',     maxerr.tau,        tol.tau);
print_result('tau(q,0,0) = G(q)',          maxerr.g,          tol.g);

if CHECK_H
    print_result('h = C*dq',               maxerr.h,          tol.h);
    print_result('h(q,0) = 0',             maxerr.h_zero,     tol.h_zero);
else
    fprintf('[SKIP] h checks (UR5_h not found)\n');
end

if CHECK_FDYN
    print_result('fdyn/idyn round-trip',   maxerr.roundtrip,  tol.roundtrip);
else
    fprintf('[SKIP] fdyn/idyn round-trip (UR5_fdyn not found)\n');
end

if CHECK_SKEW
    print_result('Mdot - 2C skew symmetry', maxerr.skew,      tol.skew);
    print_result('dq''*(Mdot-2C)*dq = 0',   maxerr.power,     tol.power);
else
    fprintf('[SKIP] skew-symmetry check\n');
end
fprintf('-------------------------------------------------\n');

%% Final pass/fail decision
all_ok = true;
all_ok = all_ok && (maxerr.fk_lastrow <= tol.fk_lastrow);
all_ok = all_ok && (maxerr.R_orth     <= tol.R_orth);
all_ok = all_ok && (maxerr.detR       <= tol.detR);
all_ok = all_ok && (maxerr.M_sym      <= tol.M_sym);
all_ok = all_ok && (minEigM           >  tol.M_pd);
all_ok = all_ok && (maxerr.C_zero     <= tol.C_zero);
all_ok = all_ok && (maxerr.tau        <= tol.tau);
all_ok = all_ok && (maxerr.g          <= tol.g);

if CHECK_H
    all_ok = all_ok && (maxerr.h      <= tol.h);
    all_ok = all_ok && (maxerr.h_zero <= tol.h_zero);
end

if CHECK_FDYN
    all_ok = all_ok && (maxerr.roundtrip <= tol.roundtrip);
end

if CHECK_SKEW
    all_ok = all_ok && (maxerr.skew  <= tol.skew);
    all_ok = all_ok && (maxerr.power <= tol.power);
end

if all_ok
    fprintf('\nPASS: all kinematics and dynamics checks passed.\n');
else
    fprintf('\nFAIL: one or more checks failed.\n');

    if minEigM <= tol.M_pd
        fprintf('\nMass matrix positive-definiteness failed.\n');
        fprintf('Smallest eigenvalue of M: %.6e\n', minEigM);
        fprintf('State that produced the smallest eigenvalue:\n');
        fprintf('q   = [% .6f  % .6f  % .6f  % .6f  % .6f  % .6f]^T\n', q_minEig);
        fprintf('dq  = [% .6f  % .6f  % .6f  % .6f  % .6f  % .6f]^T\n', dq_minEig);
        fprintf('ddq = [% .6f  % .6f  % .6f  % .6f  % .6f  % .6f]^T\n', ddq_minEig);
    end
end

%% Timing benchmark
if DO_TIMING
    fprintf('\n-------------------- Timing ---------------------\n');

    % Use one random state for timing.
    q   = rand_range_vec(-q_lim,   q_lim);
    dq  = rand_range_vec(-dq_lim,  dq_lim);
    ddq = rand_range_vec(-ddq_lim, ddq_lim);
    tau = UR5_idyn(q, dq, ddq, Pi);

    % Time each function separately.
    t_fk   = bench_ms(@() UR5_fkine(q, Pi),                    NREPS);
    t_jac  = bench_ms(@() UR5_jacobian_geometric(q, Pi),       NREPS);
    t_M    = bench_ms(@() UR5_M(q, Pi),                        NREPS);
    t_C    = bench_ms(@() UR5_C(q, dq, Pi),                    NREPS);
    t_G    = bench_ms(@() UR5_G(q, Pi),                        NREPS);
    t_MCG  = bench_ms(@() UR5_MCG(q, dq, Pi),                  NREPS);
    t_idyn = bench_ms(@() UR5_idyn(q, dq, ddq, Pi),            NREPS);

    fprintf('UR5_fkine               : %9.3f ms\n', t_fk);
    fprintf('UR5_jacobian_geometric  : %9.3f ms\n', t_jac);
    fprintf('UR5_M                   : %9.3f ms\n', t_M);
    fprintf('UR5_C                   : %9.3f ms\n', t_C);
    fprintf('UR5_G                   : %9.3f ms\n', t_G);
    fprintf('UR5_MCG                 : %9.3f ms\n', t_MCG);
    fprintf('UR5_idyn                : %9.3f ms\n', t_idyn);

    if CHECK_H
        t_h = bench_ms(@() UR5_h(q, dq, Pi), NREPS);
        fprintf('UR5_h                   : %9.3f ms\n', t_h);
    end

    if CHECK_FDYN
        t_fdyn = bench_ms(@() UR5_fdyn(q, dq, tau, Pi), NREPS);
        fprintf('UR5_fdyn                : %9.3f ms\n', t_fdyn);
    end

    fprintf('-------------------------------------------------\n');
end

%% Local functions
function x = rand_range_vec(a, b)
    x = a + (b-a).*rand(size(a));
end

function Mdot = approx_Mdot(q, dq, Pi, epsMdot)
    % Approximate Mdot along the direction dq using central difference.
    if norm(dq) < 1e-14
        Mdot = zeros(6,6);
        return;
    end

    Mp = UR5_M(q + epsMdot*dq, Pi);
    Mm = UR5_M(q - epsMdot*dq, Pi);
    Mdot = (Mp - Mm) / (2*epsMdot);
end

function t_ms = bench_ms(fh, nreps)
    % Warm up once, then report the median runtime in milliseconds.
    fh();
    t = zeros(nreps,1);

    for i = 1:nreps
        tic;
        fh();
        t(i) = toc;
    end

    t_ms = 1e3 * median(t);
end

function print_result(name, value, tol, mode)
    % mode = 'less'    means value <= tol must hold
    % mode = 'greater' means value > tol must hold
    if nargin < 4
        mode = 'less';
    end

    switch mode
        case 'less'
            ok = (value <= tol);
            limit_text = sprintf('(tol %.1e)', tol);
        case 'greater'
            ok = (value > tol);
            limit_text = sprintf('(min %.1e)', tol);
        otherwise
            error('Unknown mode');
    end

    if ok
        tag = 'PASS';
    else
        tag = 'FAIL';
    end

    fprintf('[%s] %-26s : %.3e   %s\n', tag, name, value, limit_text);
end