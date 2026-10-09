%% EM plane wave in free space: E, H, curl, div, grad|F| and the wave in 3D
%  E(r,t) = E0 cos(k n.r - w t) p            p : unit polarisation, p . n = 0
%  H(r,t) = (E0/eta) cos(k n.r - w t) h      h = n x p
%  Lengths are in wavelengths (lambda = 1 -> k = 2*pi), so derivatives are "per wavelength".
clear; clc; close all;

%% ---------------- parameters (change here) ----------------
theta = 50;       % polar angle of n from +z   [deg]
phi   = 35;       % azimuth angle of n from +x [deg]
psi   = 30;       % polarisation angle of E0 in the plane perpendicular to n [deg]
E0    = 1;        % E amplitude [V/m]
eta   = 120*pi;   % intrinsic impedance [ohm] (free space)
wt    = 0;        % omega*t : phase of the snapshot [rad]
L     = 0.75;     % half size of the plotted cube [wavelengths]
dpi   = 400;      % image resolution (600 = max quality, bigger files)
s0    = 0;        % position of the transverse plane along n [wavelengths]

%% ---------------- propagation direction n and polarisation ----------------
n  = [sind(theta)*cosd(phi), sind(theta)*sind(phi), cosd(theta)];   % unit vector
e1 = [cosd(theta)*cosd(phi), cosd(theta)*sind(phi), -sind(theta)];  % orthonormal to n
e2 = [-sind(phi), cosd(phi), 0];                                    % orthonormal to n and e1
p  = cosd(psi)*e1 + sind(psi)*e2;                                   % direction of E
h  = cross(n, p);                                                   % direction of H
k  = 2*pi;
H0 = E0/eta;

fprintf('n = (%.3f, %.3f, %.3f)\n', n);
fprintf('p = (%.3f, %.3f, %.3f)   (E direction)\n', p);
fprintf('h = (%.3f, %.3f, %.3f)   (H direction)\n', h);
fprintf('|n.p| = %.1e ,  |p x h - n| = %.1e\n', abs(dot(n,p)), norm(cross(p,h)-n));

%% ---------------- field = amplitude vector * g(phase) ----------------
AE = E0*p;   AH = H0*h;
gC = @(ph) cos(ph);
gS = @(ph) sin(ph);
gD = @(ph) -sin(ph);
gG = @(ph) sin(ph).*sign(cos(ph));          % d|cos|/d(phase)  (minus sign is in the amplitude)
clean = @(v) round(v, 12) + 0;

%% ---------------- numerical check (interior points of a fine grid) ----------------
[xf,yf,zf] = meshgrid(linspace(-L, L, 41));
PHf = k*(n(1)*xf + n(2)*yf + n(3)*zf) - wt;
[cxn,cyn,czn] = curl(xf,yf,zf, AE(1)*cos(PHf), AE(2)*cos(PHf), AE(3)*cos(PHf));
dvn = divergence(xf,yf,zf, AE(1)*cos(PHf), AE(2)*cos(PHf), AE(3)*cos(PHf));
ii  = 4:38;   sel = @(a) a(ii,ii,ii);
ACE = -k*E0*h;
dC  = sqrt(sel(cxn - ACE(1)*sin(PHf)).^2 + sel(cyn - ACE(2)*sin(PHf)).^2 + sel(czn - ACE(3)*sin(PHf)).^2);
fprintf('curl E: max|numerical - analytic| / (k E0) = %.2e  (finite-difference error only)\n', max(dC(:))/(k*E0));
fprintf('div  E: max|numerical| / (k E0)            = %.2e  (analytic value = 0)\n', max(abs(sel(dvn)),[],'all')/(k*E0));

%% ---------------- colours + shared settings ----------------
cE  = [0.00 0.45 0.74];   cCE = [0.85 0.33 0.10];   cGE = [0.10 0.60 0.20];
cH  = [0.80 0.10 0.10];   cCH = [0.55 0.20 0.70];   cGH = [0.50 0.35 0.10];

W = struct('n',n,'p',p,'h',h,'k',k,'wt',wt,'L',L,'dx',0.25, ...
           'cXY',[0.10 0.60 0.20],'cXZ',[0.85 0.65 0.00],'cYZ',[0.40 0.40 0.85], ...
           'cE',cE,'cH',cH);
W.s0 = s0;
W.sub = sprintf(['n = (%.3f, %.3f, %.3f)     p = (%.3f, %.3f, %.3f)     ' ...
                 'h = n x p = (%.3f, %.3f, %.3f)     \\omega t = %.2f rad     lengths in wavelengths'], ...
                 n, p, h, wt);

outDir = fullfile(pwd, 'EM_images');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ---------------- the 8 fields (each gets its own figure with 4 graphs) ----------------
mk = @(file,lab,ttl,A,g,kind,col,cb,cl) struct('file',file,'lab',lab,'title',ttl, ...
        'A',A,'g',g,'kind',kind,'col',col,'cb',cb,'clim',cl);
ph = 'k n\cdot r - \omega t';

F(1) = mk('plot1_E','E', ['Electric field    E = E_0 cos(' ph ')  p'], ...
          AE, gC, 'vec', cE, '|E|  [V/m]', [0 1]);
F(2) = mk('plot2_curlE','curl E', ['curl E = \nabla\times E = -k E_0 sin(' ph ')  h      (Faraday: -\mu \partial H/\partial t)'], ...
          -k*E0*h, gS, 'vec', cCE, '|curl E|  [V/m per wavelength]', [0 1]);
F(3) = mk('plot3_divE','div E', ['div E = \nabla\cdot E = -k (n\cdot E_0) sin(' ph ') = 0      (n \perp E)'], ...
          clean(k*dot(n,AE)), gD, 'scal', cE, 'div E  [V/m per wavelength]', [-k*E0 k*E0]);
F(4) = mk('plot4_gradabsE','grad|E|', ['grad|E| = \nabla|E| = -k E_0 sign(cos) sin(' ph ')  n      (|E| = E_0 |cos|)'], ...
          -k*E0*n, gG, 'vec', cGE, '|grad|E||  [V/m per wavelength]', [0 1]);

F(5) = mk('plot5_H','H', ['Magnetic field    H = (E_0/\eta) cos(' ph ')  h'], ...
          AH, gC, 'vec', cH, '|H|  [A/m]', [0 1]);
F(6) = mk('plot6_curlH','curl H', ['curl H = \nabla\times H = +k H_0 sin(' ph ')  p      (Ampere: \epsilon \partial E/\partial t)'], ...
          k*H0*p, gS, 'vec', cCH, '|curl H|  [A/m per wavelength]', [0 1]);
F(7) = mk('plot7_divH','div H', ['div H = \nabla\cdot H = -k (n\cdot H_0) sin(' ph ') = 0      (n \perp H)'], ...
          clean(k*dot(n,AH)), gD, 'scal', cH, 'div H  [A/m per wavelength]', [-k*H0 k*H0]);
F(8) = mk('plot8_gradabsH','grad|H|', ['grad|H| = \nabla|H| = -k H_0 sign(cos) sin(' ph ')  n      (|H| = H_0 |cos|)'], ...
          -k*H0*n, gG, 'vec', cGH, '|grad|H||  [A/m per wavelength]', [0 1]);

%% ---------------- draw + save ----------------
figs = gobjects(numel(F)+1, 1);
for i = 1:numel(F)
    figs(i) = plotField(W, F(i), outDir, dpi);
end
figs(end) = plotWave(W, outDir, dpi);
figs(end+1) = plotWaveTransverse(W, outDir, dpi);

set(figs, 'Visible', 'on');
fprintf('\nImages saved in: %s\n', outDir);

%% =====================================================================
%% local functions
%% =====================================================================
function fig = plotWaveTransverse(W, outDir, dpi)
fig = figure('Color','w','Position',[20 20 1500 1000],'Visible','off', ...
             'InvertHardcopy','off','Name','plot10_wave_transverse');
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');
title(tl, 'Plane wave and the transverse plane (perpendicular to n):   E only,   H only,   E and H together', ...
      'Color','k', 'FontWeight','bold', 'FontSize',13);
subtitle(tl, [W.sub sprintf('     plane at s_0 = %.2f \\lambda', W.s0)], ...
         'Color',[0.25 0.25 0.25], 'FontSize',10);

Ls = 1.5;  Ad = 0.45;  Wf = 0.8;                   % axis half length, drawn amplitude, plane half size
s  = linspace(-Ls, Ls, 301)';   sa = (-Ls:1/12:Ls)';
c  = cos(W.k*s  - W.wt);        ca = cos(W.k*sa - W.wt);
axisPts = s*W.n;
Etip = axisPts + Ad*c*W.p;      Htip = axisPts + Ad*c*W.h;
Pa = sa*W.n;   Es = Ad*ca*W.p;  Hs = Ad*ca*W.h;

c0 = cos(W.k*W.s0 - W.wt);                         % phase factor on the transverse plane
C0 = W.s0*W.n;                                     % centre of the plane
sg = sign(c0) + (c0 == 0);                         % +1 / -1 (never 0)

% ---------- (1) 3D wave + transverse plane ----------
ax = nexttile(tl, 1);  hold(ax, 'on');
h1 = plot3(ax, axisPts(:,1), axisPts(:,2), axisPts(:,3), 'k-', 'LineWidth',1.5);
h2 = plot3(ax, Etip(:,1), Etip(:,2), Etip(:,3), 'Color',W.cE, 'LineWidth',2);
h3 = plot3(ax, Htip(:,1), Htip(:,2), Htip(:,3), 'Color',W.cH, 'LineWidth',2);
quiver3(ax, Pa(:,1),Pa(:,2),Pa(:,3), Es(:,1),Es(:,2),Es(:,3), 0, 'Color',W.cE, 'LineWidth',0.8);
quiver3(ax, Pa(:,1),Pa(:,2),Pa(:,3), Hs(:,1),Hs(:,2),Hs(:,3), 0, 'Color',W.cH, 'LineWidth',0.8);
a = [-1;1;1;-1];  b = [-1;-1;1;1];
V = C0 + Wf*(a*W.p + b*W.h);                       % the transverse plane (spanned by p and h)
patch(ax,'XData',V(:,1),'YData',V(:,2),'ZData',V(:,3),'FaceColor',[0.55 0.55 0.55], ...
      'FaceAlpha',0.18,'EdgeColor',[0.3 0.3 0.3],'LineStyle','--','LineWidth',1.2);
quiver3(ax, C0(1),C0(2),C0(3), 0.5*c0*W.p(1),0.5*c0*W.p(2),0.5*c0*W.p(3), 0, ...
        'Color',W.cE, 'LineWidth',3, 'MaxHeadSize',0.8);
quiver3(ax, C0(1),C0(2),C0(3), 0.5*c0*W.h(1),0.5*c0*W.h(2),0.5*c0*W.h(3), 0, ...
        'Color',W.cH, 'LineWidth',3, 'MaxHeadSize',0.8);
e = Ls*W.n;
quiver3(ax, e(1),e(2),e(3), 0.6*W.n(1),0.6*W.n(2),0.6*W.n(3), 0, 'Color','k', 'LineWidth',3, 'MaxHeadSize',1.5);
text(ax, e(1)+0.7*W.n(1), e(2)+0.7*W.n(2), e(3)+0.7*W.n(3), 'n', 'Color','k', 'FontWeight','bold', 'FontSize',13);
axis(ax,'equal'); grid(ax,'on'); view(ax,35,25); lightax(ax);
xlabel(ax,'x/\lambda'); ylabel(ax,'y/\lambda'); zlabel(ax,'z/\lambda');
legend(ax, [h1 h2 h3], {'propagation axis n','E','H'}, 'Location','northeast');
title(ax, '3D wave  +  transverse plane (dashed)');

% ---------- common drawing of the transverse plane (u along p, v along h) ----------
[Ug, Vg] = meshgrid(-0.5:0.25:0.5);
Aa = 0.2*c0;                                       % arrow length on the grid
labU = 'u  (along p : E direction)';   labV = 'v  (along h : H direction)';

% ---------- (2) E only ----------
ax = nexttile(tl, 2);  hold(ax, 'on');
plot(ax, [-0.4 0.4], [0 0], 'k:', 'LineWidth',1);                       % oscillation range of E
quiver(ax, Ug-Aa/2, Vg, Aa*ones(size(Ug)), zeros(size(Ug)), 0, 'Color',W.cE, 'LineWidth',1.6, 'MaxHeadSize',0.6);
text(ax, 0, -0.72, sprintf('E = %.2f E_0 at this instant;  it oscillates along u between -E_0 and +E_0', c0), ...
     'HorizontalAlignment','center', 'Color','k', 'FontSize',10);
finishPlane(ax, Wf, labU, labV, 'Transverse plane:  E only');

% ---------- (3) H only ----------
ax = nexttile(tl, 3);  hold(ax, 'on');
plot(ax, [0 0], [-0.4 0.4], 'k:', 'LineWidth',1);                       % oscillation range of H
quiver(ax, Ug, Vg-Aa/2, zeros(size(Ug)), Aa*ones(size(Ug)), 0, 'Color',W.cH, 'LineWidth',1.6, 'MaxHeadSize',0.6);
text(ax, 0, -0.72, sprintf('H = %.2f H_0 at this instant;  it oscillates along v between -H_0 and +H_0', c0), ...
     'HorizontalAlignment','center', 'Color','k', 'FontSize',10);
finishPlane(ax, Wf, labU, labV, 'Transverse plane:  H only');

% ---------- (4) E and H together ----------
ax = nexttile(tl, 4);  hold(ax, 'on');
plot(ax, [-0.4 0.4], [0 0], 'k:', 'LineWidth',1);
plot(ax, [0 0], [-0.4 0.4], 'k:', 'LineWidth',1);
quiver(ax, 0, 0, 0.5*c0, 0, 0, 'Color',W.cE, 'LineWidth',3, 'MaxHeadSize',0.5);
quiver(ax, 0, 0, 0, 0.5*c0, 0, 'Color',W.cH, 'LineWidth',3, 'MaxHeadSize',0.5);
plot(ax, 0.08*sg*[1 1 0], 0.08*sg*[0 1 1], 'k-', 'LineWidth',1.3);      % right-angle mark
text(ax, 0.55*c0*sg + 0.05*sg, -0.05*sg, 'E', 'Color',W.cE, 'FontWeight','bold', 'FontSize',14);
text(ax, 0.05*sg, 0.57*c0*sg + 0.05*sg, 'H', 'Color',W.cH, 'FontWeight','bold', 'FontSize',14);
plot(ax, 0.55, 0.55, 'ko', 'MarkerSize',18, 'LineWidth',1.6);           % S out of the page
plot(ax, 0.55, 0.55, 'k.', 'MarkerSize',16);
text(ax, 0.55, 0.38, 'S = E x H', 'HorizontalAlignment','center', 'Color','k', 'FontWeight','bold');
text(ax, 0.55, 0.28, '(toward you, +n)', 'HorizontalAlignment','center', 'Color','k', 'FontSize',9);
text(ax, 0, -0.72, 'E \perp H \perp n  and  the fields are uniform over the whole plane', ...
     'HorizontalAlignment','center', 'Color','k', 'FontSize',10);
finishPlane(ax, Wf, labU, labV, 'Transverse plane:  E and H together');

savehq(fig, outDir, 'plot10_wave_transverse', dpi);
end

function finishPlane(ax, Wf, labU, labV, ttl)
axis(ax, 'equal');  xlim(ax, [-Wf Wf]);  ylim(ax, [-Wf Wf]);
grid(ax, 'on');  lightax(ax);
xlabel(ax, labU);  ylabel(ax, labV);  title(ax, ttl);
end

function fig = plotField(W, F, outDir, dpi)
fig = figure('Color','w','Position',[20 20 1500 1000],'Visible','off', ...
             'InvertHardcopy','off','Name',F.file);
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');
title(tl,    F.title, 'Color','k', 'FontWeight','bold', 'FontSize',13);
subtitle(tl, W.sub,   'Color',[0.25 0.25 0.25], 'FontSize',10);

isVec  = strcmp(F.kind, 'vec');
idx    = [1 2; 1 3; 2 3];
planes = {'XY  (z = 0)', 'XZ  (y = 0)', 'YZ  (x = 0)'};
lab    = {'x/\lambda','y/\lambda'; 'x/\lambda','z/\lambda'; 'y/\lambda','z/\lambda'};
L = W.L;

% ---------- (1) 3D ----------
ax = nexttile(tl, 1);
[xs,ys,zs] = meshgrid(-L:W.dx:L);
G = F.g(W.k*(W.n(1)*xs + W.n(2)*ys + W.n(3)*zs) - W.wt);
if isVec
    quiver3(ax, xs, ys, zs, F.A(1)*G, F.A(2)*G, F.A(3)*G, 1, 'Color',F.col, 'LineWidth',1);
else
    scatter3(ax, xs(:), ys(:), zs(:), 60, F.A*G(:), 'filled', 'MarkerEdgeColor',[0.3 0.3 0.3]);
    colormap(ax, divmap());  clim(ax, F.clim);
    cb = colorbar(ax);  cb.Label.String = F.cb;
end
hold(ax, 'on');
sq = [-L L L -L];  sq2 = [-L -L L L];  z4 = zeros(1,4);
patch(ax,'XData',sq,'YData',sq2,'ZData',z4,'FaceColor',W.cXY,'FaceAlpha',0.10,'EdgeColor',W.cXY);
patch(ax,'XData',sq,'YData',z4,'ZData',sq2,'FaceColor',W.cXZ,'FaceAlpha',0.10,'EdgeColor',W.cXZ);
patch(ax,'XData',z4,'YData',sq,'ZData',sq2,'FaceColor',W.cYZ,'FaceAlpha',0.10,'EdgeColor',W.cYZ);
q = 1.15*L*W.n;                                           % propagation direction n
quiver3(ax, 0,0,0, q(1),q(2),q(3), 0, 'Color','k', 'LineWidth',2.5, 'MaxHeadSize',0.6);
text(ax, 1.08*q(1), 1.08*q(2), 1.08*q(3), 'n', 'Color','k', 'FontWeight','bold', 'FontSize',13);
axis(ax,'equal'); grid(ax,'on'); view(ax,35,25); lightax(ax);
lim = 1.3*L;  xlim(ax,[-lim lim]); ylim(ax,[-lim lim]); zlim(ax,[-lim lim]);
xlabel(ax,'x/\lambda'); ylabel(ax,'y/\lambda'); zlabel(ax,'z/\lambda');
title(ax, ['3D  -  ' F.lab]);

% ---------- (2..4) the three faces ----------
gf = linspace(-L, L, 121);   ga = linspace(-L, L, 13);
[Af,Bf] = meshgrid(gf);      [Aa,Ba] = meshgrid(ga);
for pl = 1:3
    ax = nexttile(tl, pl+1);
    [Xf,Yf,Zf] = pl2xyz(pl, Af, Bf);
    Gf = F.g(W.k*(W.n(1)*Xf + W.n(2)*Yf + W.n(3)*Zf) - W.wt);
    if isVec
        nA = norm(F.A);
        imagesc(ax, gf, gf, nA*abs(Gf));  set(ax,'YDir','normal');
        colormap(ax, pastel(F.col));  clim(ax, [0 nA]);
        hold(ax, 'on');
        [Xa,Ya,Za] = pl2xyz(pl, Aa, Ba);
        Ga = F.g(W.k*(W.n(1)*Xa + W.n(2)*Ya + W.n(3)*Za) - W.wt);
        quiver(ax, Aa, Ba, F.A(idx(pl,1))*Ga, F.A(idx(pl,2))*Ga, 0.9, 'Color',F.col, 'LineWidth',1.1);
    else
        imagesc(ax, gf, gf, F.A*Gf);  set(ax,'YDir','normal');
        colormap(ax, divmap());  clim(ax, F.clim);
    end
    cb = colorbar(ax);  cb.Label.String = F.cb;
    axis(ax,'equal'); xlim(ax,[-L L]); ylim(ax,[-L L]); grid(ax,'on'); lightax(ax);
    xlabel(ax, lab{pl,1}); ylabel(ax, lab{pl,2});
    title(ax, [F.lab '  on plane  ' planes{pl}]);
end
savehq(fig, outDir, F.file, dpi);
end

function fig = plotWave(W, outDir, dpi)
fig = figure('Color','w','Position',[20 20 1500 1000],'Visible','off', ...
             'InvertHardcopy','off','Name','plot9_wave');
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');
title(tl, 'Plane wave:   E (blue)   and   H (red)   along the propagation axis n,    S = E \times H  (black arrow)', ...
      'Color','k', 'FontWeight','bold', 'FontSize',13);
subtitle(tl, [W.sub '     (E and H drawn with the same length; H is scaled by \eta)'], ...
         'Color',[0.25 0.25 0.25], 'FontSize',10);

Ls = 1.5;   Ad = 0.45;                         % half length of the axis, drawn amplitude
s  = linspace(-Ls, Ls, 301)';   sa = (-Ls:1/12:Ls)';
c  = cos(W.k*s  - W.wt);        ca = cos(W.k*sa - W.wt);
axisPts = s*W.n;
Etip = axisPts + Ad*c*W.p;      Htip = axisPts + Ad*c*W.h;
Pa = sa*W.n;   Es = Ad*ca*W.p;  Hs = Ad*ca*W.h;

% ---------- (1) 3D ----------
ax = nexttile(tl, 1);  hold(ax, 'on');
h1 = plot3(ax, axisPts(:,1), axisPts(:,2), axisPts(:,3), 'k-', 'LineWidth',1.5);
h2 = plot3(ax, Etip(:,1), Etip(:,2), Etip(:,3), 'Color',W.cE, 'LineWidth',2);
h3 = plot3(ax, Htip(:,1), Htip(:,2), Htip(:,3), 'Color',W.cH, 'LineWidth',2);
quiver3(ax, Pa(:,1),Pa(:,2),Pa(:,3), Es(:,1),Es(:,2),Es(:,3), 0, 'Color',W.cE, 'LineWidth',0.8);
quiver3(ax, Pa(:,1),Pa(:,2),Pa(:,3), Hs(:,1),Hs(:,2),Hs(:,3), 0, 'Color',W.cH, 'LineWidth',0.8);
a = [-1;1;1;-1];  b = [-1;-1;1;1];  Wf = 0.7;           % one wavefront (plane perpendicular to n)
V = 0*W.n + Wf*(a*W.p + b*W.h);
patch(ax,'XData',V(:,1),'YData',V(:,2),'ZData',V(:,3),'FaceColor',[0.5 0.5 0.5], ...
      'FaceAlpha',0.12,'EdgeColor',[0.4 0.4 0.4]);
e = Ls*W.n;
quiver3(ax, e(1),e(2),e(3), 0.6*W.n(1),0.6*W.n(2),0.6*W.n(3), 0, 'Color','k', 'LineWidth',3, 'MaxHeadSize',1.5);
text(ax, e(1)+0.7*W.n(1), e(2)+0.7*W.n(2), e(3)+0.7*W.n(3), 'S = E x H', 'Color','k', 'FontWeight','bold');
axis(ax,'equal'); grid(ax,'on'); view(ax,35,25); lightax(ax);
xlabel(ax,'x/\lambda'); ylabel(ax,'y/\lambda'); zlabel(ax,'z/\lambda');
legend(ax, [h1 h2 h3], {'propagation axis n','E','H'}, 'Location','northeast');
title(ax, '3D wave');

% ---------- (2..4) projections on the three faces ----------
idx    = [1 2; 1 3; 2 3];
planes = {'XY','XZ','YZ'};
lab    = {'x/\lambda','y/\lambda'; 'x/\lambda','z/\lambda'; 'y/\lambda','z/\lambda'};
for pl = 1:3
    ax = nexttile(tl, pl+1);  hold(ax, 'on');
    i1 = idx(pl,1);  i2 = idx(pl,2);
    plot(ax, axisPts(:,i1), axisPts(:,i2), 'k-', 'LineWidth',1.5);
    plot(ax, Etip(:,i1), Etip(:,i2), 'Color',W.cE, 'LineWidth',2);
    plot(ax, Htip(:,i1), Htip(:,i2), 'Color',W.cH, 'LineWidth',2);
    quiver(ax, Pa(:,i1), Pa(:,i2), Es(:,i1), Es(:,i2), 0, 'Color',W.cE, 'LineWidth',0.8);
    quiver(ax, Pa(:,i1), Pa(:,i2), Hs(:,i1), Hs(:,i2), 0, 'Color',W.cH, 'LineWidth',0.8);
    axis(ax,'equal'); grid(ax,'on'); lightax(ax);
    xlabel(ax, lab{pl,1}); ylabel(ax, lab{pl,2});
    title(ax, ['Projection on plane  ' planes{pl}]);
end
savehq(fig, outDir, 'plot9_wave', dpi);
end

function [X,Y,Z] = pl2xyz(pl, A, B)       % map 2D plane coordinates to x,y,z
O = zeros(size(A));
switch pl
    case 1, X = A; Y = B; Z = O;          % XY plane, z = 0
    case 2, X = A; Y = O; Z = B;          % XZ plane, y = 0
    case 3, X = O; Y = A; Z = B;          % YZ plane, x = 0
end
end

function cm = divmap()                    % blue - gray - red, zero = gray
t = linspace(0,1,128)';  gry = [0.85 0.85 0.85];
cm = [(1-t)*[0.15 0.35 0.85] + t*gry ; (1-t)*gry + t*[0.85 0.15 0.15]];
end

function cm = pastel(col)                 % white -> light version of col
t = linspace(0,1,128)';
cm = (1-t)*[1 1 1] + t*(0.45*[1 1 1] + 0.55*col);
end

function lightax(ax)
set(ax, 'Color','w', 'XColor','k', 'YColor','k', 'ZColor','k', ...
        'GridColor','k', 'GridAlpha',0.15, 'Box','on', 'Layer','top');
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