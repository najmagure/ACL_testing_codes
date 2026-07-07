%% 1. Load Specimen
[file, path] = uigetfile({'*.stl;*.STL','STL files (*.stl)'; '*.*','All files'}, 'Select STL file');

fname = fullfile(path,file);

TR = stlread(fname);
V = TR.Points;
F = TR.ConnectivityList;

%% 2. PCA Main Axis 
[coeff, ~, ~, ~, ~, mu] = pca(V);

n = coeff(:,1)' / norm(coeff(:,1));

%% 3. Cut Positions
perc = [20, 50, 80];

proj = V * n'; 

lengthAlongN = max(proj) - min(proj);
minp = min(proj);
p_coords = minp + (perc(:)/100) * lengthAlongN;

P0_proj = mu * n';

%% 4. Cross-Section Calculation
num = numel(p_coords);
areas = zeros(num,1);
cutPoints = zeros(num,3);
loopsByCut = cell(num,1);

for k = 1:num
    Ppos = mu + (p_coords(k) - mu*n') * n; 
    cutPoints(k,:) = Ppos;

    [areas(k), loopsByCut{k}, ~] = crossSectionAreaFromMesh(V, F, Ppos, n);
end

%% 5. Display Results
T = table(perc(:), p_coords, areas, ...
    'VariableNames', {'Percent', 'ProjCoord', 'Area'});
disp(T);

%% 6. Plot Setup
figure;
trisurf(F, V(:,1), V(:,2), V(:,3), ...
    'FaceAlpha', 0.6, ...
    'EdgeColor', 'none', ...
    'DisplayName', 'Specimen');
hold on;
axis equal;
colors = lines(num);

%% 7. In-Plane Directions
if abs(n(3)) < 0.9
    tmp = [0 0 1];
else
    tmp = [0 1 0];
end

u = cross(n,tmp); u = u / norm(u);
v = cross(n,u); v = v / norm(v);

scale = 1.2 * max(range(V));

%% 8. PCA Axis Visualization
nLine = [mu - 0.55*lengthAlongN*n; mu + 0.55*lengthAlongN*n];

plot3(nLine(:,1), nLine(:,2), nLine(:,3), 'k--', ...
    'LineWidth', 2, ...
    'DisplayName', 'n axis');

%% 9. Draw Cuts
for k = 1:num
    Ppos = cutPoints(k,:);
    loops3D = loopsByCut{k};

    for L = 1:numel(loops3D)
        pts = loops3D{L};

        if L == 1
            plot3(pts(:,1), pts(:,2), pts(:,3), '-', ...
                'Color', colors(k,:), ...
                'LineWidth', 1.5, ...
                'DisplayName', sprintf('%d%% cut', perc(k)));
        else
            plot3(pts(:,1), pts(:,2), pts(:,3), '-', ...
                'Color', colors(k,:), ...
                'LineWidth', 1.5, ...
                'HandleVisibility', 'off');
        end
    end

    [su,sv] = meshgrid([-1 1], [-1 1]);
    planeCorners = Ppos + (su(:)*scale).*u + (sv(:)*scale).*v;

    patch('Vertices', planeCorners, ...
        'Faces', [1 2 4 3], ...
        'FaceColor', colors(k,:), ...
        'FaceAlpha', 0.15, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end

%% 10. Final Format
title('Specimen Cross Sections from PCA Axis');
xlabel('X'); ylabel('Y'); zlabel('Z');
legend('show', 'Location', 'best');
grid on;
view(3);
camlight;
lighting gouraud;