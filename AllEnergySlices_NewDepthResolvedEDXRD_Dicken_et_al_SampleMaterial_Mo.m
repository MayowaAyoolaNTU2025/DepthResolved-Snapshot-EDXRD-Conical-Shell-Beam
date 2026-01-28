clear; clc; close all;

dataFolder = 'C:\Users\padra\Desktop\MATLABMcXtrace\SampleMaterial\Mo_15mm_750nE';

Nx = 80;
Ny = 80;

pixelSize = 0.25e-3;    % meters
f = 0.1;               % meters
A = 0.4;
phi_deg = 3.92;
phi = deg2rad(phi_deg);

% Finding all .epsd files

files = dir(fullfile(dataFolder, 'HEXITEC.epsd_*'));
numFiles = numel(files);

if numFiles == 0
    error('No EPSD files found.');
end

disp(['Found ', num2str(numFiles), ' EPSD files']);

% Detector Geometry

x = ((1:Nx) - (Nx+1)/2) * pixelSize;
y = ((1:Ny) - (Ny+1)/2) * pixelSize;

[X, Y] = meshgrid(x, y);
R = sqrt(X.^2 + Y.^2);

% Sample distance map
zMap = (A .* R) ./ (R + f * tan(phi));

% Scattering angle map
twothetaMap = atan(R ./ f) + phi;
twothetaMap(twothetaMap <= 0) = NaN;
thetaMap = twothetaMap/2;

% Binning Axes

% d-spacing axis (Å)
dMin = 1.0;
dMax = 2.5;
numDBins = 300;
dAxis = linspace(dMin, dMax, numDBins);

% z axis (m)
zMin = min(zMap(:));
zMax = max(zMap(:));
numZBins = 200;
zAxis = linspace(zMin, zMax, numZBins);
zAxis_mm = zAxis * 1e3;

% Accumulator
intensity_dz = zeros(numDBins, numZBins);

% For total intensity image
totalImage = zeros(Nx, Ny);

% Loop over .epsd files
disp('Processing EPSD files...');

for k = 1:numFiles

    filename = fullfile(dataFolder, files(k).name);

    % ---- Extract Energy from Header ----
    fid = fopen(filename, 'r');
    if fid == -1
        warning('Cannot open %s', files(k).name);
        continue;
    end

    E = NaN;

    while ~feof(fid)
        line = fgetl(fid);
        if contains(line, 'title: PSD monitor Energy slice')
            tokens = regexp(line, '([\d.]+)\s*keV', 'tokens');
            if ~isempty(tokens)
                E = str2double(tokens{1}{1});   % keV
                break;
            end
        end
    end
    fclose(fid);

    if isnan(E)
        warning('Energy not found in %s', files(k).name);
        continue;
    end

    % ---- Load Intensity Data ----
    data = readmatrix(filename, ...
        'FileType', 'text', ...
        'CommentStyle', '#');

    [rows, ~] = size(data);

    if rows == 3*Nx
        intensity = data(1:Nx, :);
    elseif rows == Nx
        intensity = data;
    else
        warning('Unexpected dimensions in %s', files(k).name);
        continue;
    end

    % ---- Accumulate Total Intensity Image ----
    totalImage = totalImage + intensity;

    % ---- Compute d-spacing Map for This Energy ----
    dMap = 6.2 ./ (E .* sin(thetaMap));   % Å

    % ---- Flatten and Accumulate into (d, z) ----
    dVals = dMap(:);
    zVals = zMap(:);
    iVals = intensity(:);

    valid = isfinite(dVals) & isfinite(zVals) & iVals > 0;

    dVals = dVals(valid);
    zVals = zVals(valid);
    iVals = iVals(valid);

    dIdx = discretize(dVals, dAxis);
    zIdx = discretize(zVals, zAxis);

    for n = 1:length(iVals)
        if ~isnan(dIdx(n)) && ~isnan(zIdx(n))
            intensity_dz(dIdx(n), zIdx(n)) = ...
                intensity_dz(dIdx(n), zIdx(n)) + iVals(n);
        end
    end

end

disp('All files processed.');

% Total Intensity Image

figure;
imagesc(totalImage);
axis image;
colorbar;
title('Intensity Image (Mo.lau)');
xlim([20 60]);
ylim([20 60]);
xticks(20:20:60);
yticks(20:20:60);
xlabel('Pixel X');
ylabel('Pixel Y');

% d-Spacing vs Sample Distance Heatmap
figure;
imagesc(zAxis_mm, dAxis, intensity_dz);
set(gca, 'YDir', 'reverse');
colorbar;

xlabel('Sample distance z (mm)');
ylabel('d-spacing (Å)');
title('d-spacing vs Sample Distance (Mo.lau)');
xlim([80 180]);

% Material-Specific Diffractograms (z-sliced)

exposureTime = 300; % 300 seconds

Counts_per_second = intensity_dz / 300; 

% d-spacing range
d_low  = 1.0;   % Å
d_high = 2.5;   % Å
dMask = (dAxis >= d_low) & (dAxis <= d_high);
d_selected = dAxis(dMask);

% z ranges (mm)
zlow_mm = 100;
zhigh_mm = 200;

zMask = (zAxis_mm >= zlow_mm) & (zAxis_mm <= zhigh_mm);

% Integrate intensity
I_d = sum(Counts_per_second(dMask, zMask), 2);

% Plotting

figure;
plot(d_selected, I_d, 'LineWidth', 2);
grid on;

xlabel('d-spacing (Å)');
ylabel('Counts per second');

legend( ...
    'Mo.lau', ...
    'Location', 'best' );

title('1D Diffractogram (Mo.lau)');
