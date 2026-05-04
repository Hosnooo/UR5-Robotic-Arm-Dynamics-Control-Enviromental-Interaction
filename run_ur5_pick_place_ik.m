%% run_ur5_pick_place_ik.m
% UR5 numerical IK pick-and-place animation with video export.
%
% This script:
%   1. Loads the generated UR5 kinematic model.
%   2. Builds a pick-and-place Cartesian task.
%   3. Solves the task using utils/inverse_kinematics.m.
%   4. Animates the UR5 motion.
%   5. Saves the result as an MP4 video.

clear; clc; close all;

%% Paths
addpath('maple_gen');
addpath('utils');

USE_FIG_HELPERS = (exist('figureoptscall','file') == 2) && ...
                  (exist('saveFigureAsPDF','file') == 2);

if USE_FIG_HELPERS
    figureoptscall;
end

%% Load UR5 parameters
[~, Pi] = UR5_params();

%% Video settings
videoFile = fullfile('ur5_pick_place_ik_animation.mp4');
fps = 30;

%% Numerical IK settings
ikOpts = struct();
ikOpts.maxIter = 120;
ikOpts.positionTolerance = 2e-4;
ikOpts.lambda = 0.035;
ikOpts.alpha = 0.9;
ikOpts.maxStep = deg2rad(4);
ikOpts.useOrientation = false;  % Position-only IK is robust for this demo.

%% Define reachable task locations from known UR5 configurations
% These joint configurations are used only to define reachable Cartesian
% pick/place points. The trajectory itself is solved using numerical IK.
qHome = deg2rad([  30;  -80;  105; -115;  -90;    0]);
qPickSeed = deg2rad([  45;  -95;  115; -110;  -90;   20]);
qPlaceSeed = deg2rad([ -35;  -80;  105; -115;  -90;  -20]);

THome = UR5_fkine(qHome, Pi);
TPick = UR5_fkine(qPickSeed, Pi);
TPlace = UR5_fkine(qPlaceSeed, Pi);

pHome = THome(1:3,4);
pPick = TPick(1:3,4);
pPlace = TPlace(1:3,4);

liftHeight = 0.14;
pPickAbove = pPick + [0; 0; liftHeight];
pPlaceAbove = pPlace + [0; 0; liftHeight];

%% Build smooth Cartesian task trajectory
% Object modes:
%   0 = object waiting at pick location
%   1 = object attached to the end-effector
%   2 = object placed at final location

segments = {
    pHome,       pPickAbove,   50, 0, 'Move above object';
    pPickAbove,  pPick,        30, 0, 'Lower to grasp';
    pPick,       pPickAbove,   35, 1, 'Lift object';
    pPickAbove,  pPlaceAbove,  65, 1, 'Carry object';
    pPlaceAbove, pPlace,       30, 1, 'Lower to place';
    pPlace,      pPlaceAbove,  35, 2, 'Release object';
    pPlaceAbove, pHome,        50, 2, 'Return home'
};

pDes = [];
objectMode = [];
phaseIndex = [];

phaseNames = strings(size(segments,1),1);

for i = 1:size(segments,1)
    pA = segments{i,1};
    pB = segments{i,2};
    nFrames = segments{i,3};
    mode = segments{i,4};
    phaseNames(i) = string(segments{i,5});

    pts = smoothCartesianSegment(pA, pB, nFrames);

    if i > 1
        pts = pts(:,2:end); % Avoid duplicate frames at segment boundaries.
    end

    pDes = [pDes, pts]; %#ok<AGROW>
    objectMode = [objectMode, mode*ones(1,size(pts,2))];
    phaseIndex = [phaseIndex, i*ones(1,size(pts,2))];
end

nTotal = size(pDes,2);

%% Solve IK along the Cartesian path
Q = zeros(6,nTotal);
ikPositionError = zeros(1,nTotal);
ikConverged = false(1,nTotal);

qPrev = qHome;

fprintf('Solving numerical IK for %d trajectory frames...\n', nTotal);

for k = 1:nTotal
    [qSol, info] = inverse_kinematics(pDes(:,k), qPrev, Pi, ikOpts);

    Q(:,k) = qSol;
    qPrev = qSol;

    ikPositionError(k) = info.positionErrorNorm;
    ikConverged(k) = info.converged;
end

fprintf('IK solved.\n');
fprintf('Maximum position error: %.6f m\n', max(ikPositionError));
fprintf('Mean position error   : %.6f m\n', mean(ikPositionError));
fprintf('Converged frames      : %d / %d\n', nnz(ikConverged), nTotal);

%% Pre-compute link positions for animation limits
linkPositions = zeros(3,7,nTotal);
TeeAll = zeros(4,4,nTotal);

for k = 1:nTotal
    [Tlist, Tee] = UR5_fkall(Q(:,k), Pi);
    TeeAll(:,:,k) = Tee;

    for j = 1:7
        linkPositions(:,j,k) = Tlist{j}(1:3,4);
    end
end

allPts = reshape(linkPositions, 3, []);
allPts = [allPts, pDes, pPick, pPlace, pPickAbove, pPlaceAbove];

padding = 0.18;
xyzMin = min(allPts, [], 2) - padding;
xyzMax = max(allPts, [], 2) + padding;

xyzMin(3) = min(xyzMin(3), -0.02);
xyzMax(3) = max(xyzMax(3), 0.65);

%% Prepare video writer
try
    writerObj = VideoWriter(videoFile, 'MPEG-4');
catch
    videoFile = fullfile(scriptDir, 'ur5_pick_place_ik_animation.avi');
    writerObj = VideoWriter(videoFile, 'Motion JPEG AVI');
end

writerObj.FrameRate = fps;
open(writerObj);

%% Animate and save video
fig = figure('Color', 'w', 'Position', [100 100 1200 800]);

cubeSide = 0.06;
cubePickCenter = pPick - [0; 0; cubeSide/2];
cubePlaceCenter = pPlace - [0; 0; cubeSide/2];

fprintf('Writing video: %s\n', videoFile);

for k = 1:nTotal
    clf(fig);
    hold on; grid on; axis equal;

    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');

    xlim([xyzMin(1), xyzMax(1)]);
    ylim([xyzMin(2), xyzMax(2)]);
    zlim([xyzMin(3), xyzMax(3)]);

    view(170, 25);
    set(gca, 'FontSize', 11);

    drawFloor(xyzMin, xyzMax);

    % Full desired path and completed path
    plot3(pDes(1,:), pDes(2,:), pDes(3,:), '--', ...
        'Color', [0.65 0.65 0.65], 'LineWidth', 1.0);
    plot3(pDes(1,1:k), pDes(2,1:k), pDes(3,1:k), ...
        'Color', [0.0 0.45 0.74], 'LineWidth', 2.0);

    % Pick and place markers
    plot3(pPick(1), pPick(2), pPick(3), 'o', ...
        'MarkerSize', 8, 'MarkerFaceColor', [0.47 0.67 0.19], ...
        'MarkerEdgeColor', 'k');
    plot3(pPlace(1), pPlace(2), pPlace(3), 's', ...
        'MarkerSize', 8, 'MarkerFaceColor', [0.85 0.33 0.10], ...
        'MarkerEdgeColor', 'k');

    % Object state
    if objectMode(k) == 0
        cubeCenter = cubePickCenter;
    elseif objectMode(k) == 1
        pe = TeeAll(1:3,4,k);
        cubeCenter = pe - [0; 0; cubeSide/2];
    else
        cubeCenter = cubePlaceCenter;
    end

    drawCube(cubeCenter, cubeSide, [1.0 0.8 0.2], 0.95);

    % Robot and gripper
    currentLinks = linkPositions(:,:,k);
    drawRobot(currentLinks);
    drawGripper(TeeAll(:,:,k), objectMode(k) == 1);

    phaseText = char(phaseNames(phaseIndex(k)));
    title(sprintf('UR5 Numerical IK Pick-and-Place | %s | Frame %d/%d', ...
        phaseText, k, nTotal), 'Interpreter', 'none');

    drawnow;
    writeVideo(writerObj, getframe(fig));
end

close(writerObj);

fprintf('Done. Video saved to:\n%s\n', videoFile);

%% Local helper functions

function pts = smoothCartesianSegment(pA, pB, nFrames)
%SMOOTHCARTESIANSEGMENT Quintic time-scaling between two 3D points.

    tau = linspace(0, 1, nFrames);
    s = 10*tau.^3 - 15*tau.^4 + 6*tau.^5;
    pts = pA(:) + (pB(:) - pA(:)) .* s;
end

function drawRobot(P)
%DRAWROBOT Draw UR5 link chain from precomputed frame origins.

    plot3(P(1,:), P(2,:), P(3,:), '-o', ...
        'LineWidth', 4.0, ...
        'MarkerSize', 7, ...
        'MarkerFaceColor', [0.1 0.1 0.1], ...
        'MarkerEdgeColor', 'k', ...
        'Color', [0.0 0.25 0.70]);

    plot3(P(1,1), P(2,1), P(3,1), 's', ...
        'MarkerSize', 12, ...
        'MarkerFaceColor', [0.2 0.2 0.2], ...
        'MarkerEdgeColor', 'k');
end

function drawGripper(Tee, isClosed)
%DRAWGRIPPER Draw a simple two-finger gripper at the end-effector frame.

    R = Tee(1:3,1:3);
    p = Tee(1:3,4);

    if isClosed
        jawWidth = 0.028;
    else
        jawWidth = 0.075;
    end

    fingerLength = 0.07;

    sideAxis = R(:,2);
    approachAxis = -R(:,3);

    base1 = p + 0.5*jawWidth*sideAxis;
    base2 = p - 0.5*jawWidth*sideAxis;

    tip1 = base1 + fingerLength*approachAxis;
    tip2 = base2 + fingerLength*approachAxis;

    plot3([base1(1) tip1(1)], [base1(2) tip1(2)], [base1(3) tip1(3)], ...
        'k-', 'LineWidth', 3.0);
    plot3([base2(1) tip2(1)], [base2(2) tip2(2)], [base2(3) tip2(3)], ...
        'k-', 'LineWidth', 3.0);

    plot3([base1(1) base2(1)], [base1(2) base2(2)], [base1(3) base2(3)], ...
        'k-', 'LineWidth', 2.0);

    plot3(p(1), p(2), p(3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
end

function drawCube(center, sideLength, faceColor, faceAlpha)
%DRAWCUBE Draw a cube centered at center.

    c = center(:).';
    s = sideLength/2;

    vertices = [
        -s -s -s;
         s -s -s;
         s  s -s;
        -s  s -s;
        -s -s  s;
         s -s  s;
         s  s  s;
        -s  s  s] + c;

    faces = [
        1 2 3 4;
        5 6 7 8;
        1 2 6 5;
        2 3 7 6;
        3 4 8 7;
        4 1 5 8];

    patch('Vertices', vertices, ...
          'Faces', faces, ...
          'FaceColor', faceColor, ...
          'FaceAlpha', faceAlpha, ...
          'EdgeColor', 'k', ...
          'LineWidth', 1.0);
end

function drawFloor(xyzMin, xyzMax)
%DRAWFLOOR Draw a simple ground plane at z = 0.

    x1 = xyzMin(1);
    x2 = xyzMax(1);
    y1 = xyzMin(2);
    y2 = xyzMax(2);

    floorVertices = [
        x1 y1 0;
        x2 y1 0;
        x2 y2 0;
        x1 y2 0];

    floorFaces = [1 2 3 4];

    patch('Vertices', floorVertices, ...
          'Faces', floorFaces, ...
          'FaceColor', [0.93 0.93 0.93], ...
          'FaceAlpha', 0.35, ...
          'EdgeColor', [0.75 0.75 0.75]);
end
