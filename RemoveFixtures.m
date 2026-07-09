%selecting the file to use
[file, path] = uigetfile({'*.stl;*.STL','STL files (*.stl)'; '*.*','All files'}, 'Select STL file');

fname = fullfile(path,file);

%% Read STL

TR = stlread(fname);

V = TR.Points;
F = TR.ConnectivityList;

% PCA axis

[coeff,~,~,~,~] = pca(V);

axisVec = coeff(:,1);
axisVec = axisVec ./ norm(axisVec);

%% Coordinate along PCA axis


center = mean(V,1);

t = (V-center)*axisVec;

%% Plot mesh

figure

trisurf(F,...
        V(:,1),V(:,2),V(:,3),...
        'FaceAlpha',0.2,...
        'EdgeColor','none');

hold on
axis equal
camlight
lighting gouraud

%% PCA axis line

L = linspace(min(t),max(t),200)';

axisPts = center + L*axisVec';

plot3(axisPts(:,1),...
      axisPts(:,2),...
      axisPts(:,3),...
      'r','LineWidth',3)
title('Inspect PCA Axis')

figure('Name','Select Ligament Region')

histogram(t,200)

xlabel('Position Along PCA Axis')
ylabel('Vertex Count')

title({'Drag the red lines'; ...
       'Double-click each line when finished'})

hold on

xmin = min(t);
xmax = max(t);

x1 = xmin + 0.3*(xmax-xmin);
x2 = xmin + 0.7*(xmax-xmin);

h1 = drawline('Position',[x1 0; x1 max(ylim)], ...
              'Color','r','LineWidth',2);

h2 = drawline('Position',[x2 0; x2 max(ylim)], ...
              'Color','r','LineWidth',2);
%Read Locations ======
p1 = h1.Position;
p2 = h2.Position;

t1 = p1(1,1);
t2 = p2(1,1);

if t1 > t2
    tmp = t1;
    t1 = t2;
    t2 = tmp;
end

fprintf('Selected region:\n')
fprintf('t1 = %.3f\n',t1)
fprintf('t2 = %.3f\n',t2)

%%Preview retained
keepVertex = (t >= t1) & (t <= t2);

figure

scatter3(V(:,1),V(:,2),V(:,3), ...
    2,[0.8 0.8 0.8]);

hold on

scatter3(V(keepVertex,1), ...
    V(keepVertex,2), ...
    V(keepVertex,3), ...
    8,'r','filled')

axis equal
grid on

title('Red = Retained Vertices')

insideCount = sum(keepVertex(F),2);

keepFace = insideCount >= 2;

Fkeep = F(keepFace,:);

figure

trisurf(F,...
    V(:,1),V(:,2),V(:,3), ...
    'FaceAlpha',0.1,...
    'EdgeColor','none')

hold on

trisurf(Fkeep,...
    V(:,1),V(:,2),V(:,3), ...
    'FaceColor','cyan',...
    'EdgeColor','none')

axis equal
camlight
lighting gouraud

title('Retained Region')
%added edits
