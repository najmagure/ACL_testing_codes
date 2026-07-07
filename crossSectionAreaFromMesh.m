function [areaTotal, loops3D, areas] = crossSectionAreaFromMesh(V, F, P0, n)
% crossSectionAreaFromMesh  Cross-sectional area of triangular mesh by plane
% [areaTotal, loops3D, areas] = crossSectionAreaFromMesh(V,F,P0,n)
% V: Nx3 vertex list
% F: Mx3 face indices
% P0: 1x3 point on plane
% n:  1x3 plane normal (not necessarily unit)
% areaTotal: scalar total cross-sectional area
% loops3D: cell array of Kx3 loop point arrays (3D coordinates)
% areas: vector of area per loop (same order as loops3D)

n = n(:)' / norm(n);
tol = 1e-12;
edges = [1 2; 2 3; 3 1];

% collect intersection segments
segments = zeros(size(F,1)*2,6);
segcount = 0;
for i = 1:size(F,1)
    tri = V(F(i,:),:);                 % 3x3
    d = (tri - P0) * n';               % signed distances
    pts = zeros(0,3);
    for e = 1:3
        a = edges(e,1); b = edges(e,2);
        da = d(a); db = d(b);
        if abs(da) < tol && abs(db) < tol
            pts = [pts; tri(a,:); tri(b,:)]; %#ok<AGROW>
        elseif abs(da) < tol
            pts = [pts; tri(a,:)]; %#ok<AGROW>
        elseif abs(db) < tol
            pts = [pts; tri(b,:)]; %#ok<AGROW>
        elseif da * db < 0
            t = da / (da - db);
            pnt = tri(a,:) + t*(tri(b,:) - tri(a,:));
            pts = [pts; pnt]; %#ok<AGROW>
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

% no intersection
if isempty(segments)
    areaTotal = 0;
    loops3D = {};
    areas = [];
    return;
end

% unique points and adjacency
P_list = [segments(:,1:3); segments(:,4:6)];
P_list = unique(round(P_list,12),'rows','stable');

toIndex = @(pnt) find(all(abs(P_list - pnt) < 1e-8, 2), 1);
E = zeros(size(segments,1),2);
for k = 1:size(segments,1)
    p1 = segments(k,1:3); p2 = segments(k,4:6);
    E(k,1) = toIndex(p1); E(k,2) = toIndex(p2);
end

% stitch segments into loops
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
    if chain(1) ~= chain(end)
        chain = [chain, chain(1)];
    end
    loopsIdx{end+1} = chain;
end

% compute areas by projecting into plane basis
% construct in-plane orthonormal basis u,v
u_guess = P_list(1,:) - P0;
if norm(u_guess) < 1e-8
    u_guess = P_list(min(2,size(P_list,1)),:) - P0;
end
u = u_guess - n*(n*u_guess');
u = u / norm(u);
v = cross(n,u); v = v / norm(v);

numLoops = numel(loopsIdx);
loops3D = cell(numLoops,1);
areas = zeros(numLoops,1);
for k = 1:numLoops
    idxs = loopsIdx{k};
    pts3 = P_list(idxs,:);
    if isequal(pts3(1,:), pts3(end,:)), pts3(end,:) = []; end
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