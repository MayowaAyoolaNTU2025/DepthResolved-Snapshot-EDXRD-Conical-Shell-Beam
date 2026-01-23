
clear; clc; close all;

dataFolder = 'C:\Users\padra\Desktop\MATLABMcXtrace';

Nx = 80;
Ny = 80;
numBins = 750;

pixelSize = 0.25e-3;    
f = 0.1;                
A = 0.4;
phi_deg = 3.92;
phi = deg2rad(phi_deg);

exposureTime = 300;       

Emin = 1;
Emax = 130;
energy = linspace(Emin, Emax, numBins);

% Preallocating 3D Data Cube
detectorCube = zeros(Nx, Ny, numBins);

% Loading .epsd files
disp('Loading EPSD files into 3D cube...');

for k = 1:numBins

    fileIndex = k - 1;
    filename = fullfile(dataFolder, sprintf('HEXITEC.epsd_%03d', fileIndex));

    if exist(filename, 'file')

        data = readmatrix(filename, 'FileType', 'text', 'CommentStyle', '#');

        [rows, cols] = size(data);

        if rows == 3*Nx
            intensityData = data(1:Nx, :);
        elseif rows == Nx
            intensityData = data;
        else
            error('Unexpected EPSD dimensions in %s: %d x %d', filename, rows, cols);
        end

        detectorCube(:,:,k) = intensityData;

    else
        warning('Missing file: %s', filename);
    end
end

disp('All files loaded.');
disp(['Cube size: ', mat2str(size(detectorCube))]);

% Generation of Total Intensity Image

totalImage = sum(detectorCube, 3);

figure;
imagesc(totalImage);
axis image;
colorbar;
title('Integrated Intensity Image');
xlabel('Pixel X');
ylabel('Pixel Y');

% Creation of Detector Coordinates

x = ((1:Nx) - (Nx+1)/2) * pixelSize;
y = ((1:Ny) - (Ny+1)/2) * pixelSize;

[X, Y] = meshgrid(x, y);

R = sqrt(X.^2 + Y.^2);

% Sample Distance, z = A * r / (r + f * tan(phi))

zMap = (A .* R) ./ (R + f * tan(phi));

%Scattering Angle 

twothetaMap = (atan(R ./ f) + phi);   % radians
twothetaMap(twothetaMap <= 0) = NaN;

% d-spacing AND z bin axes

% d-spacing axis (Å)
dMin = 0.5;
dMax = 6.0;
numDBins = 150;
dAxis = linspace(dMin, dMax, numDBins);

% z axis (m)
zMin = min(zMap(:));
zMax = max(zMap(:));
numZBins = 100;
zAxis = linspace(zMin, zMax, numZBins);
zAxis_mm = zAxis * 1e3;

% Accumulating Intensity into (d, z) space

intensity_dz = zeros(numDBins, numZBins);

disp('Generating d-spacing vs sample distance (z) map...');

for k = 1:numBins

    E = energy(k);                     % keV
    slice = detectorCube(:,:,k);

    dMap = 6.2 ./ (E .* sin(twothetaMap));

    % Flatten
    dVals = dMap(:);
    zVals = zMap(:);
    iVals = slice(:);

    % Valid values only
    valid = isfinite(dVals) & isfinite(zVals) & iVals > 0;

    dVals = dVals(valid);
    zVals = zVals(valid);
    iVals = iVals(valid);

    % Bin indices
    dIdx = discretize(dVals, dAxis);
    zIdx = discretize(zVals, zAxis);

    % Accumulate intensity
    for n = 1:length(iVals)
        if ~isnan(dIdx(n)) && ~isnan(zIdx(n))
            intensity_dz(dIdx(n), zIdx(n)) = ...
                intensity_dz(dIdx(n), zIdx(n)) + iVals(n);
        end
    end
end

disp('d–z intensity map completed.');

% d-Spacing vs Sample Distance HeatMap

figure;
imagesc(zAxis_mm, dAxis, intensity_dz);
set(gca, 'YDir', 'normal');
colorbar;

xlabel('Sample distance z (mm)');
xlim([50 200]);

ylabel('d-spacing (Å)');
title('d-spacing vs Sample Distance (z) Intensity Map');

% Selection of Ring Radii

disp('Please select TWO points on the diffraction rings in the image.');

figure; 
imagesc(totalImage); 
axis image; 
colorbar;
title('Click TWO points on the diffraction rings');

% ginput returns [x, y] coordinates of mouse clicks in axes units
[clickX, clickY] = ginput(2);

% Convert clicked pixel coordinates to distances from image center
centerX = (Nx+1)/2;
centerY = (Ny+1)/2;

ringRadius1 = sqrt((clickX(1)-centerX)^2 + (clickY(1)-centerY)^2) * pixelSize;
ringRadius2 = sqrt((clickX(2)-centerX)^2 + (clickY(2)-centerY)^2) * pixelSize;

disp(['Selected ring radii: ', num2str(ringRadius1), ' m and ', num2str(ringRadius2), ' m']);

tolerance = 0.001;

ringMask1 = (R > ringRadius1 - tolerance) & (R < ringRadius1 + tolerance);
ringMask2 = (R > ringRadius2 - tolerance) & (R < ringRadius2 + tolerance);

% Integrate Intensity vs Energy for Each Ring

intensity1 = zeros(1, numBins);
intensity2 = zeros(1, numBins);

for k = 1:numBins
    imageSlice = detectorCube(:,:,k);
    intensity1(k) = sum(imageSlice(ringMask1));
    intensity2(k) = sum(imageSlice(ringMask2));
end

% Counts per second

cps1 = intensity1 / exposureTime;
cps2 = intensity2 / exposureTime;

% Calculation of d-spacing

twoTheta1 = atan(ringRadius1 / f) + phi;
twoTheta2 = atan(ringRadius2 / f) + phi;

theta1 = twoTheta1 / 2;
theta2 = twoTheta2 / 2;

d1 = 6.199 ./ (energy .* sin(theta1));
d2 = 6.199 ./ (energy .* sin(theta2));

% Smoothening Data

cps1s = smooth(cps1, 5);
cps2s = smooth(cps2, 5);

% 1D Diffractograms

figure;
plot(d1, cps1s, 'LineWidth', 2); hold on;
plot(d2, cps2s, 'LineWidth', 2);
xlabel('d-spacing (Å)');
ylabel('Counts per second');
title('Diffractograms for C_diamond.laz and Si.laz');
legend('C_diamond.laz', 'Si.laz');
grid on;
xlim([0 15]);
%xticks(0:1:3); 
