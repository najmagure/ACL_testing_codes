%% 1. Load Specimen
[file, path] = uigetfile({'*.stl;*.STL','STL files (*.stl)'; '*.*','All files'}, 'Select STL file');
fname = fullfile(path,file);

TR = stlread(fname);
V = TR.Points;
F = TR.ConnectivityList;

%% 2. PCA Main Axis
% First principal component = direction of greatest elongation
[coeff, ~, ~, ~, ~, mu] = pca(V);
n = coeff(:,1)' / norm(coeff(:,1));

%% 3. Cut Positions
% Convert percentages along the main axis to positions in projected space
perc = [20, 50, 80];
proj = V * n';
lengthAlongN = max(proj) - min(proj);
minp = min(proj);
p_coords = minp + (perc(:)/100) * lengthAlongN;

%% 4. Cross-Section Calculation
num = numel(p_coords);
areas = zeros(num,1);
cutPoints = zeros(num,3);
loopsByCut = cell(num,1);

for k = 1:num
    Ppos = mu + (p_coords(k) - mu*n') * n;   % point on cutting plane along n
    cutPoints(k,:) = Ppos;
    [areas(k), loopsByCut{k}, ~] = crossSectionAreaFromMesh(V, F, Ppos, n);
end

%% 5. Display Results
T = table(perc(:), p_coords, areas, 'VariableNames', {'Percent', 'ProjCoord', 'Area'});
disp(T);

%% 6. Plot Setup
figure;
trisurf(F, V(:,1), V(:,2), V(:,3), 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Specimen');
hold on; axis equal;
colors = lines(num);

%% 7. In-Plane Directions
% Build an orthonormal basis (u,v) spanning each cutting plane
tmp = [0 0 1];
if abs(n(3)) >= 0.9
    tmp = [0 1 0];
end
u = cross(n,tmp); u = u / norm(u);
v = cross(n,u);   v = v / norm(v);

scale = 1.2 * max(range(V));

%% 8. PCA Axis Visualization
nLine = [mu - 0.55*lengthAlongN*n; mu + 0.55*lengthAlongN*n];
plot3(nLine(:,1), nLine(:,2), nLine(:,3), 'k--', 'LineWidth', 2, 'DisplayName', 'n axis');

%% 9. Draw Cuts
for k = 1:num
    Ppos = cutPoints(k,:);
    loops3D = loopsByCut{k};

    for L = 1:numel(loops3D)
        pts = loops3D{L};
        visFlag = 'off';
        if L == 1, visFlag = 'on'; end
        plot3(pts(:,1), pts(:,2), pts(:,3), '-', 'Color', colors(k,:), ...
            'LineWidth', 1.5, 'DisplayName', sprintf('%d%% cut', perc(k)), 'HandleVisibility', visFlag);
    end

    % translucent square patch to visualize the cutting plane
    [su,sv] = meshgrid([-1 1], [-1 1]);
    planeCorners = Ppos + (su(:)*scale).*u + (sv(:)*scale).*v;
    patch('Vertices', planeCorners, 'Faces', [1 2 4 3], 'FaceColor', colors(k,:), ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

%% 10. Final Format
title('Specimen Cross Sections from PCA Axis');
xlabel('X'); ylabel('Y'); zlabel('Z');
legend('show', 'Location', 'best');
grid on; view(3); camlight; lighting gouraud;

%% 11. Flat 2D Projections of Each Cut
figure('Name','Flat 2D Cross-Section Projections');

for k = 1:num
    subplot(1, num, k);
    hold on; axis equal; grid on; box on;

    Ppos = cutPoints(k,:);
    loops3D = loopsByCut{k};

    allXY = zeros(0,2);
    for L = 1:numel(loops3D)
        pts3 = loops3D{L};
        XY = [(pts3 - Ppos) * u', (pts3 - Ppos) * v'];
        plot(XY(:,1), XY(:,2), '-', 'Color', colors(k,:), 'LineWidth', 1.5); % gap left as-is, not force-closed
        allXY = [allXY; XY];
    end

    if size(allXY,1) >= 3
        [ellX, ellY] = fitEquivalentEllipse2D(allXY);
        plot(ellX, ellY, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
    end

    title(sprintf('%d%% cut (Area = %.3f)', perc(k), areas(k)));
    xlabel('u'); ylabel('v');

    % add 10% padding around the shape
    x_bounds = xlim; y_bounds = ylim;
    x_pad = range(x_bounds) * 0.1;
    y_pad = range(y_bounds) * 0.1;
    xlim(x_bounds + [-x_pad, x_pad]);
    ylim(y_bounds + [-y_pad, y_pad]);
end

sgtitle('Flat 2D Projections of Cross Sections');

%% Reference Ellipse Fit (point-based second moments, no shoelace)
function [ellX, ellY] = fitEquivalentEllipse2D(XY)
c = mean(XY,1);
Xc = XY - c;
Cov = (Xc' * Xc) / size(Xc,1);

[evec, eval] = eig(Cov);
[lam, order] = sort(diag(eval), 'descend');
evec = evec(:, order);

a_axis = sqrt(2*max(lam(1),0));
b_axis = sqrt(2*max(lam(2),0));
phi = atan2(evec(2,1), evec(1,1));   % major-axis angle

t = linspace(0, 2*pi, 100)';
ellX = c(1) + a_axis*cos(t)*cos(phi) - b_axis*sin(t)*sin(phi);
ellY = c(2) + a_axis*cos(t)*sin(phi) + b_axis*sin(t)*cos(phi);
end

%% Cross Sectional Area from Mesh
function [areaTotal, loops3D, areas] = crossSectionAreaFromMesh(V, F, P0, n)
n = n(:)' / norm(n);
tol = 1e-12;
edges = [1 2; 2 3; 3 1];

% Step 1: find where the cutting plane crosses each triangle's edges
segments = zeros(size(F,1)*2,6);
segcount = 0;
for i = 1:size(F,1)
    tri = V(F(i,:),:);
    d = (tri - P0) * n';         % signed distance of each vertex to the plane
    pts = zeros(0,3);
    for e = 1:3
        a = edges(e,1); b = edges(e,2);
        da = d(a); db = d(b);
        if abs(da) < tol && abs(db) < tol
            pts = [pts; tri(a,:); tri(b,:)];
        elseif abs(da) < tol
            pts = [pts; tri(a,:)];
        elseif abs(db) < tol
            pts = [pts; tri(b,:)];
        elseif da * db < 0
            t = da / (da - db);
            pts = [pts; tri(a,:) + t*(tri(b,:) - tri(a,:))]; %#ok<AGROW>
        end
    end
    if size(pts,1) >= 2
        pts = unique(round(pts,12),'rows');
        if size(pts,1) == 2
            segcount = segcount + 1;
            segments(segcount,:) = [pts(1,:) pts(2,:)];
        elseif size(pts,1) > 2
            for k = 1:size(pts,1)-1
                segcount = segcount + 1;
                segments(segcount,:) = [pts(k,:) pts(k+1,:)];
            end
        end
    end
end
segments = segments(1:segcount,:);

if isempty(segments)
    areaTotal = 0; loops3D = {}; areas = [];
    return;
end

% Step 2: build unique point list + edge index list for stitching
P_list = unique(round([segments(:,1:3); segments(:,4:6)],12),'rows','stable');
toIndex = @(pnt) find(all(abs(P_list - pnt) < 1e-8, 2), 1);
E = zeros(size(segments,1),2);
for k = 1:size(segments,1)
    E(k,1) = toIndex(segments(k,1:3));
    E(k,2) = toIndex(segments(k,4:6));
end

% Step 3: stitch edge segments into closed loops (walk forward, then backward)
used = false(size(E,1),1);
loopsIdx = {};
while any(~used)
    e = find(~used,1); used(e) = true;
    chain = E(e,:);
    startIdx = chain(1); endIdx = chain(2);
    extended = true;
    while extended
        extended = false;
        for ee = find(~used)'
            if E(ee,1) == endIdx
                used(ee) = true; endIdx = E(ee,2); chain = [chain, endIdx]; extended = true; break;
            elseif E(ee,2) == endIdx
                used(ee) = true; endIdx = E(ee,1); chain = [chain, endIdx]; extended = true; break;
            end
        end
        if ~extended
            for ee = find(~used)'
                if E(ee,2) == startIdx
                    used(ee) = true; startIdx = E(ee,1); chain = [startIdx, chain]; extended = true; break;
                elseif E(ee,1) == startIdx
                    used(ee) = true; startIdx = E(ee,2); chain = [startIdx, chain]; extended = true; break;
                end
            end
        end
    end
    loopsIdx{end+1} = chain;
end

% Step 4: project each loop into the plane's (u,v) basis and compute area
u_guess = P_list(1,:) - P0;
if norm(u_guess) < 1e-8
    u_guess = P_list(min(2,size(P_list,1)),:) - P0;
end
u = u_guess - n*(n*u_guess'); u = u / norm(u);
v = cross(n,u); v = v / norm(v);

numLoops = numel(loopsIdx);
loops3D = cell(numLoops,1);
areas = zeros(numLoops,1);
for k = 1:numLoops
    pts3 = P_list(loopsIdx{k},:);
    if isequal(pts3(1,:), pts3(end,:))
        pts3(end,:) = [];   % drop duplicate closing point
    end
    loops3D{k} = pts3;

    XY = [(pts3 - P0) * u', (pts3 - P0) * v'];
    c = mean(XY,1);
    ang = atan2(XY(:,2)-c(2), XY(:,1)-c(1));
    [~, order] = sort(ang);
    XYs = XY(order,:);
    areas(k) = polyarea(XYs(:,1), XYs(:,2));
end
areaTotal = sum(areas);
end