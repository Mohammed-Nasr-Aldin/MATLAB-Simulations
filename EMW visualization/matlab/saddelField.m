%% F = (x, -y, 0): div F = 0 AND curl F = 0   (white theme + HQ export)
%  F = grad(phi), phi = (x^2 - y^2)/2  -> curl = 0 and div = 0
clear; clc; close all;

dpi = 600;                   % image resolution (use 300 if you run out of memory)

%% ---------- field + numerical div / curl ----------
[x, y, z] = meshgrid(-2:0.5:2, -2:0.5:2, -2:1:2);
U = x;  V = -y;  W = zeros(size(x));

divF = divergence(x, y, z, U, V, W);
[cx, cy, cz] = curl(x, y, z, U, V, W);

clean = @(v) round(v, 10) + 0;                    % avoids printing -0.00
divF = clean(divF);  cx = clean(cx);  cy = clean(cy);  cz = clean(cz);
divMax  = max(abs(divF(:)));
cxM = mean(cx(:));  cyM = mean(cy(:));  czM = mean(cz(:));
cmag    = sqrt(cx.^2 + cy.^2 + cz.^2);
curlMax = max(cmag(:));

% Checks on a circle of radius 1 in the XY plane
th   = linspace(0, 2*pi, 400);
xc   = cos(th);  yc = sin(th);
circ = clean(trapz(th, xc.*(-yc) + (-yc).*xc));   % circulation: integral of F . dl
flux = clean(trapz(th, xc.*xc + (-yc).*yc));      % flux:        integral of F . n dl

fprintf('max|div F| = %.3f\n', divMax);
fprintf('curl F = (%.2f, %.2f, %.2f)\n', cxM, cyM, czM);
fprintf('Circulation on r=1 circle = %.4f ; curl_z*area = %.4f\n', circ, czM*pi);
fprintf('Flux on r=1 circle        = %.4f ; div*area    = %.4f\n', flux, 0);

% colours (dark enough for a white background)
cField = [0.00 0.45 0.74];   % blue   : F
cCurl  = [0.85 0.33 0.10];   % orange : curl F
cXY    = [0.10 0.60 0.20];   % green  : plane z = 0
cXZ    = [0.85 0.65 0.00];   % gold   : plane y = y0
cYZ    = [0.85 0.15 0.15];   % red    : plane x = x0

% colormaps: zero is light gray so dots stay visible on white
n = 128;  t = linspace(0,1,n)';
gry = [0.85 0.85 0.85];
cmDiv  = [(1-t)*[0.15 0.35 0.85] + t*gry ; (1-t)*gry + t*[0.85 0.15 0.15]];   % blue-gray-red
cmCurl = (1-t)*gry + t*cCurl;                                                  % gray-orange

y0 = 1;  x0 = 1;             % offset planes (so the arrows are visible)

outDir = fullfile(pwd, 'output_images');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% =====================================================
%% FIGURE 1 : 3D field + one 2D plot per plane
%% =====================================================
fig1 = figure('Color','w','Position',[20 20 1500 1000],'Visible','off', ...
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
axis(ax1,'equal'); grid(ax1,'on'); view(ax1,35,25); lightax(ax1);
xlabel(ax1,'x'); ylabel(ax1,'y'); zlabel(ax1,'z');
title(ax1,'3D field  F = (x, -y, 0)','Color','k');
cellnote(fig1, 1, {'Field flows outward along x and inward along y (saddle flow).', ...
    sprintf('Slicing planes:  z = 0 (green),  y = %g (gold),  x = %g (red).', y0, x0), ...
    'Each plane is shown in 2D in the other three plots.'});

% ---- 2D grid used for all planes ----
[a, b] = meshgrid(-2:0.25:2, -2:0.25:2);
zer = zeros(size(a));

% ---- (2) plane XY (z = 0) ----
ax2 = cellax(fig1, 2);
quiver(ax2, a, b, a, -b, 1.2, 'Color', cXY, 'LineWidth', 1); hold(ax2,'on');
plot(ax2, xc, yc, '--', 'Color', 'k', 'LineWidth', 1.2);
axis(ax2,'equal'); grid(ax2,'on'); lightax(ax2); xlim(ax2,[-2.5 2.5]); ylim(ax2,[-2.5 2.5]);
xlabel(ax2,'x'); ylabel(ax2,'y'); title(ax2,'Plane XY  (z = 0)','Color','k');
cellnote(fig1, 2, { ...
    'In-plane F = (x, -y): saddle flow, stagnation point at the origin.', ...
    sprintf('div F = %.2f   (dFx/dx = +1 cancels dFy/dy = -1)', divMax), ...
    sprintf('curl_z = %.2f   (no circulation in this plane)', czM), ...
    sprintf('Dashed circle r = 1:  circulation = %.3f ,  net flux = %.3f', circ, flux)});

% ---- (3) plane XZ (y = y0) ----
ax3 = cellax(fig1, 3);
quiver(ax3, a, b, a, zer, 1.2, 'Color', cXZ, 'LineWidth', 1);
axis(ax3,'equal'); grid(ax3,'on'); lightax(ax3); xlim(ax3,[-2.5 2.5]); ylim(ax3,[-2.5 2.5]);
xlabel(ax3,'x'); ylabel(ax3,'z'); title(ax3,sprintf('Plane XZ  (y = %g)', y0),'Color','k');
cellnote(fig1, 3, { ...
    'In-plane F = (x, 0): horizontal, pointing away from x = 0.', ...
    sprintf('Out-of-plane component: Fy = -y = %g (constant, not visible here).', -y0), ...
    sprintf('div F = %.2f', divMax), ...
    sprintf('curl_y = %.2f   (no circulation in this plane)', cyM)});

% ---- (4) plane YZ (x = x0) ----
ax4 = cellax(fig1, 4);
quiver(ax4, a, b, -a, zer, 1.2, 'Color', cYZ, 'LineWidth', 1);
axis(ax4,'equal'); grid(ax4,'on'); lightax(ax4); xlim(ax4,[-2.5 2.5]); ylim(ax4,[-2.5 2.5]);
xlabel(ax4,'y'); ylabel(ax4,'z'); title(ax4,sprintf('Plane YZ  (x = %g)', x0),'Color','k');
cellnote(fig1, 4, { ...
    'In-plane F = (-y, 0): horizontal, pointing toward y = 0.', ...
    sprintf('Out-of-plane component: Fx = x = %g (constant, not visible here).', x0), ...
    sprintf('div F = %.2f', divMax), ...
    sprintf('curl_x = %.2f   (no circulation in this plane)', cxM)});

savehq(fig1, outDir, 'plot1_field_planes', dpi);

%% =====================================================
%% FIGURE 2 : field + div F (colour map) + results panel
%% =====================================================
fig2 = figure('Color','w','Position',[20 20 1500 1000],'Visible','off', ...
              'InvertHardcopy','off','Name','F and div F');

% ---- (1) field ----
bx1 = axes('Parent', fig2, 'Position', [0.05 0.42 0.40 0.50]);
quiver3(bx1, x, y, z, U, V, W, 'Color', cField, 'LineWidth', 1);
axis(bx1,'equal'); grid(bx1,'on'); view(bx1,35,25); lightax(bx1);
xlabel(bx1,'x'); ylabel(bx1,'y'); zlabel(bx1,'z');
title(bx1,'F = (x, -y, 0)','Color','k');
note(fig2, [0.05 0.28 0.40 0.09], {'Blue arrows: the field F = (x, -y, 0).', ...
    'Arrows do not circulate: they stretch along x and squeeze along y.', ...
    'What spreads out in x comes back in through y.'});

% ---- (2) div F : every grid point coloured by its divergence ----
bx2 = axes('Parent', fig2, 'Position', [0.54 0.42 0.36 0.50]);
scatter3(bx2, x(:), y(:), z(:), 45, divF(:), 'filled', ...
         'MarkerEdgeColor', [0.2 0.2 0.2]);
colormap(bx2, cmDiv);  clim(bx2, [-1 1]);
axis(bx2,'equal'); grid(bx2,'on'); view(bx2,35,25); lightax(bx2);
xlabel(bx2,'x'); ylabel(bx2,'y'); zlabel(bx2,'z');
title(bx2,'div F  (colour = value of div F at each point)','Color','k');
cb = colorbar(bx2, 'Position', [0.93 0.45 0.015 0.44]);
cb.Color = 'k';  cb.Label.String = 'div F';  cb.Label.Color = 'k';
set(bx2, 'Position', [0.54 0.42 0.36 0.50]);
note(fig2, [0.54 0.28 0.42 0.09], {sprintf('div F = %.2f at every point (max|div F| = %.2g)', 0, divMax), ...
    'All dots have the same neutral colour = 0 on the colour scale (-1 ... +1).', ...
    'No red (source) and no blue (sink) anywhere: incompressible flow.'});

% ---- (3) results panel (two columns) ----
bx3 = axes('Parent', fig2, 'Position', [0.05 0.02 0.91 0.22]);
axis(bx3, 'off'); set(bx3, 'Color', 'w');
left = { ...
    'RESULTS', ...
    '', ...
    'F = (x, -y, 0)', ...
    '', ...
    'div F = dFx/dx + dFy/dy + dFz/dz = 1 + (-1) + 0 = 0', ...
    sprintf('   numerical: max|div F| = %.2g   (zero everywhere)', divMax), ...
    '', ...
    'curl F = (dFz/dy - dFy/dz,  dFx/dz - dFz/dx,  dFy/dx - dFx/dy)', ...
    '       = (0 - 0, 0 - 0, 0 - 0) = (0, 0, 0)', ...
    sprintf('   numerical: curl F = (%.2f, %.2f, %.2f)', cxM, cyM, czM)};
right = { ...
    'MEANING', ...
    '', ...
    '- div = 0: no sources or sinks, the flow is incompressible.', ...
    '', ...
    '- curl = 0: no rotation, a small paddle wheel would not spin.', ...
    '', ...
    '- F = grad(phi) with phi = (x^2 - y^2)/2 (a potential exists), and', ...
    '  laplacian(phi) = 1 - 1 = 0 (phi is harmonic).', ...
    sprintf('- Check on r = 1 circle: circulation = %.3f , net flux = %.3f.', circ, flux)};
text(bx3, 0.00, 1.0, left,  'Units','normalized', 'VerticalAlignment','top', ...
     'Color','k', 'FontSize',11, 'Interpreter','none');
text(bx3, 0.52, 1.0, right, 'Units','normalized', 'VerticalAlignment','top', ...
     'Color','k', 'FontSize',11, 'Interpreter','none');

savehq(fig2, outDir, 'plot2_F_div', dpi);

%% =====================================================
%% FIGURE 3 : the 3D field alone
%% =====================================================
fig3 = figure('Color','w','Position',[20 20 1200 900],'Visible','off', ...
              'InvertHardcopy','off','Name','3D field only');
ax5 = axes('Parent', fig3, 'Position', [0.08 0.10 0.84 0.82]);
quiver3(ax5, x, y, z, U, V, W, 'Color', cField, 'LineWidth', 1.2);
axis(ax5,'equal'); grid(ax5,'on'); view(ax5,35,25); lightax(ax5);
xlabel(ax5,'x'); ylabel(ax5,'y'); zlabel(ax5,'z');
title(ax5,'3D vector field  F = (x, -y, 0)','Color','k','FontSize',14);

savehq(fig3, outDir, 'plot3_field_3D', dpi);

%% =====================================================
%% FIGURE 4 : curl F alone (colour = |curl F|)
%% =====================================================
fig4 = figure('Color','w','Position',[20 20 1200 900],'Visible','off', ...
              'InvertHardcopy','off','Name','curl F only');
ax6 = axes('Parent', fig4, 'Position', [0.08 0.14 0.76 0.78]);
scatter3(ax6, x(:), y(:), z(:), 55, cmag(:), 'filled', 'MarkerEdgeColor', [0.2 0.2 0.2]);
colormap(ax6, cmCurl);  clim(ax6, [0 1]);
axis(ax6,'equal'); grid(ax6,'on'); view(ax6,35,25); lightax(ax6);
xlabel(ax6,'x'); ylabel(ax6,'y'); zlabel(ax6,'z');
title(ax6, sprintf('curl F = (%.2f, %.2f, %.2f)  at every point', cxM, cyM, czM), ...
      'Color','k','FontSize',14);
cb2 = colorbar(ax6, 'Position', [0.89 0.25 0.02 0.55]);
cb2.Color = 'k';  cb2.Label.String = '|curl F|';  cb2.Label.Color = 'k';
set(ax6, 'Position', [0.08 0.14 0.76 0.78]);
note(fig4, [0.08 0.02 0.76 0.07], {sprintf('max |curl F| over the grid = %.2g : irrotational field.', curlMax), ...
    'The curl is the zero vector everywhere, so no arrows can be drawn; dots mark the grid.'});

savehq(fig4, outDir, 'plot4_curl', dpi);

%% ---------- done: show all figures ----------
set([fig1 fig2 fig3 fig4], 'Visible', 'on');
fprintf('\nImages saved in: %s\n', outDir);

%% ---------- local functions ----------
function lightax(ax)
set(ax, 'Color','w', 'XColor','k', 'YColor','k', 'ZColor','k', ...
        'GridColor','k', 'GridAlpha',0.15, 'Box','on');
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
annotation(fig, 'textbox', pos, 'String', str, 'Color','k', 'BackgroundColor','w', ...
    'EdgeColor',[0.5 0.5 0.5], 'Interpreter','none', 'FontSize',10, ...
    'VerticalAlignment','middle', 'FitBoxToText','off');
end

function savehq(fig, folder, name, dpi)   % PNG at high resolution, white background
pngFile = fullfile(folder, [name '.png']);
drawnow;
try
    exportgraphics(fig, pngFile, 'Resolution', dpi, 'BackgroundColor', 'w');
    fprintf('Saved: %s\n', pngFile);
catch ME
    warning('exportgraphics failed for "%s" (%s). Trying print...', name, ME.message);
    try
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, pngFile, '-dpng', sprintf('-r%d', dpi));
        fprintf('Saved (print): %s\n', pngFile);
    catch ME2
        warning('Could not save "%s": %s', name, ME2.message);
    end
end
end