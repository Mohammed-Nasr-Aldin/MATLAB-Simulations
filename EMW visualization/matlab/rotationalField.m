%% F = (-y, x, 0): 3D field, its 2D planes, and curl   (dark theme)
clear; clc; close all;

%% ---------- field + numerical div / curl ----------
[x, y, z] = meshgrid(-2:0.5:2, -2:0.5:2, -2:1:2);
U = -y;  V = x;  W = zeros(size(x));

divF = divergence(x, y, z, U, V, W);
[cx, cy, cz] = curl(x, y, z, U, V, W);

divMax = max(abs(divF(:)));
cxM = mean(cx(:));  cyM = mean(cy(:));  czM = mean(cz(:));

% Stokes check on a circle of radius 1 in the XY plane
th   = linspace(0, 2*pi, 400);
xc   = cos(th);  yc = sin(th);
circ = trapz(th, (-yc).*(-yc) + (xc).*(xc));      % integral of F . dl

fprintf('max|div F| = %.3f\n', divMax);
fprintf('curl F = (%.2f, %.2f, %.2f)\n', cxM, cyM, czM);
fprintf('Circulation on r=1 circle = %.4f ; curl_z*area = %.4f\n', circ, czM*pi);

% colours
cField = [0.00 0.85 1.00];   % cyan   : F
cCurl  = [1.00 0.65 0.00];   % orange : curl F
cXY    = [0.40 1.00 0.40];   % green  : plane z = 0
cXZ    = [1.00 0.90 0.30];   % yellow : plane y = y0
cYZ    = [1.00 0.45 0.45];   % red    : plane x = x0

y0 = 1;  x0 = 1;             % offset planes (so the arrows are visible)

%% =====================================================
%% FIGURE 1 : 3D field + one 2D plot per plane
%% =====================================================
fig1 = figure('Color','k','Position',[20 20 1500 1000], 'WindowState','maximized', ...
              'InvertHardcopy','off','Name','Field and its planes');

% ---- (1) 3D field with the three planes drawn in it ----
ax1 = cellax(fig1, 1);
quiver3(ax1, x, y, z, U, V, W, 'Color', cField, 'LineWidth', 1); hold(ax1,'on');
patch(ax1,'XData',[-2 2 2 -2],'YData',[-2 -2 2 2],'ZData',[0 0 0 0], ...
      'FaceColor',cXY,'FaceAlpha',0.15,'EdgeColor',cXY);
patch(ax1,'XData',[-2 2 2 -2],'YData',y0*ones(1,4),'ZData',[-2 -2 2 2], ...
      'FaceColor',cXZ,'FaceAlpha',0.15,'EdgeColor',cXZ);
patch(ax1,'XData',x0*ones(1,4),'YData',[-2 2 2 -2],'ZData',[-2 -2 2 2], ...
      'FaceColor',cYZ,'FaceAlpha',0.15,'EdgeColor',cYZ);
axis(ax1,'equal'); grid(ax1,'on'); view(ax1,35,25); darkax(ax1);
xlabel(ax1,'x'); ylabel(ax1,'y'); zlabel(ax1,'z');
title(ax1,'3D field  F = (-y, x, 0)','Color','w');
cellnote(fig1, 1, {'Field rotates counter-clockwise around the z-axis.', ...
    sprintf('Slicing planes:  z = 0 (green),  y = %g (yellow),  x = %g (red).', y0, x0), ...
    'Each plane is shown in 2D in the other three plots.'});

% ---- 2D grid used for all planes ----
[a, b] = meshgrid(-2:0.25:2, -2:0.25:2);
one = ones(size(a));  zer = zeros(size(a));

% ---- (2) plane XY (z = 0) ----
ax2 = cellax(fig1, 2);
quiver(ax2, a, b, -b, a, 1.2, 'Color', cXY, 'LineWidth', 1); hold(ax2,'on');
plot(ax2, xc, yc, '--', 'Color', 'w', 'LineWidth', 1.2);       % Stokes loop
axis(ax2,'equal'); grid(ax2,'on'); darkax(ax2); xlim(ax2,[-2.5 2.5]); ylim(ax2,[-2.5 2.5]);
xlabel(ax2,'x'); ylabel(ax2,'y'); title(ax2,'Plane XY  (z = 0)','Color','w');
cellnote(fig1, 2, { ...
    'In-plane F = (-y, x): circulates counter-clockwise.', ...
    sprintf('div F = %.2f   (no source / sink)', divMax), ...
    sprintf('curl_z = %.2f   (strong circulation in this plane)', czM), ...
    sprintf('Stokes check (dashed circle, r = 1): loop integral of F.dl = %.3f ,  curl_z x area = %.3f', circ, czM*pi)});

% ---- (3) plane XZ (y = y0) ----
ax3 = cellax(fig1, 3);
quiver(ax3, a, b, -y0*one, zer, 1.2, 'Color', cXZ, 'LineWidth', 1);
axis(ax3,'equal'); grid(ax3,'on'); darkax(ax3); xlim(ax3,[-2.5 2.5]); ylim(ax3,[-2.5 2.5]);
xlabel(ax3,'x'); ylabel(ax3,'z'); title(ax3,sprintf('Plane XZ  (y = %g)', y0),'Color','w');
cellnote(fig1, 3, { ...
    sprintf('In-plane F = (%g, 0): uniform, points toward -x.', -y0), ...
    'Out-of-plane component: Fy = x (not visible in this view).', ...
    sprintf('div F = %.2f', divMax), ...
    sprintf('curl_y = %.2f   (no circulation in this plane)', cyM)});

% ---- (4) plane YZ (x = x0) ----
ax4 = cellax(fig1, 4);
quiver(ax4, a, b, x0*one, zer, 1.2, 'Color', cYZ, 'LineWidth', 1);
axis(ax4,'equal'); grid(ax4,'on'); darkax(ax4); xlim(ax4,[-2.5 2.5]); ylim(ax4,[-2.5 2.5]);
xlabel(ax4,'y'); ylabel(ax4,'z'); title(ax4,sprintf('Plane YZ  (x = %g)', x0),'Color','w');
cellnote(fig1, 4, { ...
    sprintf('In-plane F = (%g, 0): uniform, points toward +y.', x0), ...
    'Out-of-plane component: Fx = -y (not visible in this view).', ...
    sprintf('div F = %.2f', divMax), ...
    sprintf('curl_x = %.2f   (no circulation in this plane)', cxM)});

%% =====================================================
%% FIGURE 2 : field + curl (big)  and the results panel
%% =====================================================
fig2 = figure('Color','k','Position',[20 20 1500 1000], 'WindowState','maximized', ...
              'InvertHardcopy','off','Name','F and curl F');

% ---- (1) field ----
bx1 = axes('Parent', fig2, 'Position', [0.05 0.42 0.42 0.50]);
quiver3(bx1, x, y, z, U, V, W, 'Color', cField, 'LineWidth', 1);
axis(bx1,'equal'); grid(bx1,'on'); view(bx1,35,25); darkax(bx1);
xlabel(bx1,'x'); ylabel(bx1,'y'); zlabel(bx1,'z');
title(bx1,'F = (-y, x, 0)','Color','w');
note(fig2, [0.05 0.28 0.42 0.09], {'Cyan arrows: the field F = (-y, x, 0).', ...
    'The arrows circulate around the z-axis: the field rotates.', ...
    sprintf('div F = %.2f everywhere (flat plane at zero): no source, no sink.', divMax)});

% ---- (2) curl ----
bx2 = axes('Parent', fig2, 'Position', [0.54 0.42 0.42 0.50]);
quiver3(bx2, x, y, z, cx, cy, cz, 'Color', cCurl, 'LineWidth', 1);
axis(bx2,'equal'); grid(bx2,'on'); view(bx2,35,25); darkax(bx2);
xlabel(bx2,'x'); ylabel(bx2,'y'); zlabel(bx2,'z');
title(bx2,'curl F','Color','w');
note(fig2, [0.54 0.28 0.42 0.09], {sprintf('curl F = (%.2f, %.2f, %.2f)', cxM, cyM, czM), ...
    'Constant and pointing along +z (right-hand rule of the rotation).', ...
    'Same value at every point: uniform rotation.'});

% ---- (3) results panel (two columns) ----
bx3 = axes('Parent', fig2, 'Position', [0.05 0.02 0.91 0.22]);
axis(bx3, 'off'); set(bx3, 'Color', 'k');
left = { ...
    'RESULTS', ...
    '', ...
    'F = (-y, x, 0)', ...
    '', ...
    'div F = dFx/dx + dFy/dy + dFz/dz = 0 + 0 + 0 = 0', ...
    sprintf('   numerical: max|div F| = %.2g   (zero everywhere, flat plane z = 0)', divMax), ...
    '', ...
    'curl F = (dFz/dy - dFy/dz,  dFx/dz - dFz/dx,  dFy/dx - dFx/dy)', ...
    '       = (0, 0, 1 - (-1)) = (0, 0, 2)', ...
    sprintf('   numerical: curl F = (%.2f, %.2f, %.2f)', cxM, cyM, czM)};
right = { ...
    'MEANING', ...
    '', ...
    '- div = 0: field lines neither start nor end anywhere.', ...
    '', ...
    '- curl = 2 z: the field circulates; a small paddle wheel', ...
    '  would spin about +z with angular velocity curl/2 = 1 rad/s.', ...
    '', ...
    sprintf('- Stokes: loop integral of F.dl = %.3f = curl_z x pi r^2  (r = 1).', circ)};
text(bx3, 0.00, 1.0, left,  'Units','normalized', 'VerticalAlignment','top', ...
     'Color','w', 'FontSize',11, 'Interpreter','none');
text(bx3, 0.52, 1.0, right, 'Units','normalized', 'VerticalAlignment','top', ...
     'Color','w', 'FontSize',11, 'Interpreter','none');

%% ---------- optional: save images with the black background ----------
% exportgraphics(fig1, 'plot1_field_planes.png', 'BackgroundColor','k', 'Resolution',300);
% exportgraphics(fig2, 'plot2_F_curl.png',       'BackgroundColor','k', 'Resolution',300);

%% ---------- local functions ----------
function darkax(ax)
set(ax, 'Color','k', 'XColor','w', 'YColor','w', 'ZColor','w', ...
        'GridColor','w', 'GridAlpha',0.3, 'Box','on');
end

function ax = cellax(fig, k)          % big plot area of cell k (1..4); text goes below it
ll = [0 0.5; 0.5 0.5; 0 0; 0.5 0];
p  = ll(k,:);
ax = axes('Parent', fig, 'Position', [p(1)+0.07, p(2)+0.15, 0.40, 0.31]);
end

function cellnote(fig, k, str)        % text box under the plot of cell k
ll = [0 0.5; 0.5 0.5; 0 0; 0.5 0];
p  = ll(k,:);
note(fig, [p(1)+0.05, p(2)+0.005, 0.43, 0.095], str);
end

function note(fig, pos, str)          % generic text box (outside the axes)
annotation(fig, 'textbox', pos, 'String', str, 'Color','w', 'BackgroundColor','k', ...
    'EdgeColor',[0.5 0.5 0.5], 'Interpreter','none', 'FontSize',10, ...
    'VerticalAlignment','middle', 'FitBoxToText','off');
end