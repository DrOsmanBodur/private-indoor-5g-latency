function make_fig5_trial_latency()
% Figure 5 of the article: trial-level mean-latency distributions.
%   (a) empirical CDF of trial-level mean latency by target UDP rate
%   (b) proportions of 60-second trials whose mean latency exceeded 10, 20,
%       30 and 50 ms, with 95% Wilson confidence intervals
% Both panels are written side by side to one self-contained SVG file; the
% article shows them as Figure 5a and 5b.
%
% Input:   data/latency_210.csv
% Output:  fig5_trial_latency_ecdf_exceedance.svg in this folder

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));

Tmain = readtable(fullfile(root, 'data', 'latency_210.csv'));
bwList = sort(unique(Tmain.B_Mbps));
nBW = numel(bwList);
latByBW = cell(nBW,1);
for k = 1:nBW
    latByBW{k} = Tmain.latency_ms(Tmain.B_Mbps == bwList(k));
end

thresholds = [10 20 30 50];
exceedPropMat = nan(nBW, numel(thresholds));
exceedCiLoMat = nan(nBW, numel(thresholds));
exceedCiHiMat = nan(nBW, numel(thresholds));
for k = 1:nBW
    x = latByBW{k}; n = numel(x);
    for t = 1:numel(thresholds)
        cnt = sum(x > thresholds(t));
        exceedPropMat(k,t) = cnt/n;
        [lo,hi] = wilsonCI(cnt, n);
        exceedCiLoMat(k,t) = lo; exceedCiHiMat(k,t) = hi;
    end
end

bwLabels = arrayfun(@(v) sprintf('%g Mbps', v), bwList, 'UniformOutput', false);
bwColors = {'#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#084594'};
threshColors = {'#1b9e77','#d95f02','#7570b3','#e7298a'};

captionLines = {};   % the caption is set in the article

svg = buildCombinedFigure(latByBW, bwLabels, bwColors, thresholds, exceedPropMat, exceedCiLoMat, exceedCiHiMat, threshColors, captionLines);

outSvg = fullfile(here, 'fig5_trial_latency_ecdf_exceedance.svg');
fid = fopen(outSvg, 'w'); fprintf(fid, '%s', svg); fclose(fid);
fprintf('Wrote %s\n', outSvg);
end

%% ============ LOCAL FUNCTIONS ============

function p = empPercentile(x, q)
x = sort(x(:)); n = numel(x);
if n == 1, p = x(1); return; end
h = (q/100)*(n-1) + 1; lo = floor(h); hi = ceil(h);
if lo == hi, p = x(lo); else, p = x(lo) + (h-lo)*(x(hi)-x(lo)); end
end

function [lo,hi] = wilsonCI(x, n)
z = 1.959963984540054;
if n == 0, lo = NaN; hi = NaN; return; end
phat = x/n; denom = 1 + z^2/n;
center = phat + z^2/(2*n);
adj = z*sqrt((phat*(1-phat) + z^2/(4*n))/n);
lo = max(0, (center - adj)/denom); hi = min(1, (center + adj)/denom);
end

function s = svgRect(x,y,w,h,fill,stroke,strokeW,opacity)
if nargin < 8, opacity = 1; end
s = sprintf('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="%s" stroke="%s" stroke-width="%.2f" opacity="%.2f"/>\n', x,y,w,h,fill,stroke,strokeW,opacity);
end
function s = svgLine(x1,y1,x2,y2,stroke,strokeW)
s = sprintf('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="%.2f"/>\n', x1,y1,x2,y2,stroke,strokeW);
end
function s = svgText(x,y,txt,size,anchor,weight,color,rotate,style)
if nargin < 5, anchor='middle'; end
if nargin < 6, weight='normal'; end
if nargin < 7, color='#1a1a1a'; end
if nargin < 8, rotate=0; end
if nargin < 9, style='normal'; end
txt = strrep(txt,'&','&amp;');
rotAttr = ''; if rotate~=0, rotAttr = sprintf(' transform="rotate(%d %.2f %.2f)"', rotate,x,y); end
s = sprintf('<text x="%.2f" y="%.2f" font-family="Arial,Helvetica,sans-serif" font-size="%d" text-anchor="%s" font-weight="%s" font-style="%s" fill="%s"%s>%s</text>\n', ...
    x,y,size,anchor,weight,style,color,rotAttr,txt);
end
function s = svgCircle(cx,cy,r,fill,opacity)
if nargin<5, opacity=1; end
s = sprintf('<circle cx="%.2f" cy="%.2f" r="%.2f" fill="%s" opacity="%.2f"/>\n', cx,cy,r,fill,opacity);
end
function s = svgPolyline(xs,ys,stroke,strokeW)
pts=''; for i=1:numel(xs), pts=[pts sprintf('%.2f,%.2f ',xs(i),ys(i))]; end %#ok<AGROW>
s = sprintf('<polyline points="%s" fill="none" stroke="%s" stroke-width="%.2f"/>\n', pts, stroke, strokeW);
end
function step = niceStep(rough)
if rough <= 0, step = 1; return; end
mag = 10^floor(log10(rough)); resid = rough/mag;
if resid<1.5, step=1*mag; elseif resid<3, step=2*mag; elseif resid<7, step=5*mag; else, step=10*mag; end
end

function content = panelEcdf(groupData, labels, colors, panelW, panelH, xLabel)
% Returns inner SVG markup (no outer <svg> tag) for the ECDF panel, local coords (0,0)-(panelW,panelH).
plotLeft = 70; plotRight = panelW - 155; plotTop = 15; plotBottom = panelH - 45;
allVals = vertcat(groupData{:}); dataMin=min(allVals); dataMax=max(allVals);
rng_=dataMax-dataMin; pad=rng_*0.05; xlo=max(0,dataMin-pad); xhi=dataMax+pad;
xToPix = @(v) plotLeft + (v-xlo)/(xhi-xlo)*(plotRight-plotLeft);
yToPix = @(v) plotBottom - v*(plotBottom-plotTop);
content = '';
content = [content svgLine(plotLeft, plotTop, plotLeft, plotBottom, '#333333', 1.2)];
content = [content svgLine(plotLeft, plotBottom, plotRight, plotBottom, '#333333', 1.2)];
for f = 0:0.2:1
    yp = yToPix(f);
    content = [content svgLine(plotLeft-5, yp, plotRight, yp, '#e5e5e5', 0.8)]; %#ok<AGROW>
    content = [content svgText(plotLeft-10, yp+4, sprintf('%.1f',f), 11, 'end')]; %#ok<AGROW>
end
step = niceStep((xhi-xlo)/6);
for xv = ceil(xlo/step)*step : step : xhi
    xp = xToPix(xv);
    content = [content svgLine(xp, plotBottom, xp, plotBottom+6, '#333333', 1)]; %#ok<AGROW>
    content = [content svgText(xp, plotBottom+22, sprintf('%g',xv), 11, 'middle')]; %#ok<AGROW>
end
content = [content svgText(18, (plotTop+plotBottom)/2, 'Empirical CDF', 12, 'middle', 'bold', '#1a1a1a', -90)];
content = [content svgText((plotLeft+plotRight)/2, panelH-8, xLabel, 12, 'middle', 'bold')];
for g = 1:numel(groupData)
    x = sort(groupData{g}(:)); n = numel(x);
    xs=[]; ys=[];
    for i=1:n, xs=[xs,x(i),x(i)]; ys=[ys,(i-1)/n,i/n]; end %#ok<AGROW>
    col = colors{mod(g-1,numel(colors))+1};
    content = [content svgPolyline(arrayfun(xToPix,xs), arrayfun(yToPix,ys), col, 2)]; %#ok<AGROW>
    legY = plotTop + (g-1)*20 + 8;
    content = [content svgLine(plotRight+15, legY, plotRight+38, legY, col, 3)]; %#ok<AGROW>
    content = [content svgText(plotRight+44, legY+4, labels{g}, 11, 'start')]; %#ok<AGROW>
end
end

function content = panelExceedance(bwLabels, thresholds, propMat, ciLoMat, ciHiMat, colors, panelW, panelH)
[nBW, nT] = size(propMat);
plotLeft = 65; plotRight = panelW - 20; plotTop = 15; plotBottom = panelH - 95;
yToPix = @(v) plotBottom - v*(plotBottom-plotTop);
content = '';
content = [content svgLine(plotLeft, plotTop, plotLeft, plotBottom, '#333333', 1.2)];
content = [content svgLine(plotLeft, plotBottom, plotRight, plotBottom, '#333333', 1.2)];
maxY = min(1, max(ciHiMat(:))*1.15);
for f = 0:0.1:ceil(maxY*10)/10
    yp = yToPix(f);
    content = [content svgLine(plotLeft-5, yp, plotRight, yp, '#e5e5e5', 0.8)]; %#ok<AGROW>
    content = [content svgText(plotLeft-10, yp+4, sprintf('%.1f',f), 10, 'end')]; %#ok<AGROW>
end
content = [content svgText(16, (plotTop+plotBottom)/2, 'Proportion of trials exceeding threshold', 11, 'middle', 'bold', '#1a1a1a', -90)];
groupW = (plotRight-plotLeft)/nBW; barW = groupW/(nT+1.5);
for k = 1:nBW
    gx0 = plotLeft + groupW*(k-1);
    for t = 1:nT
        xc = gx0 + barW*(t+0.25);
        p = propMat(k,t); lo = ciLoMat(k,t); hi = ciHiMat(k,t);
        yb = yToPix(0); yp = yToPix(p);
        col = colors{mod(t-1,numel(colors))+1};
        content = [content svgRect(xc-barW*0.4, min(yp,yb), barW*0.8, abs(yb-yp), col, '#333333', 0.8)]; %#ok<AGROW>
        yLo = yToPix(lo); yHi = yToPix(hi);
        content = [content svgLine(xc, yLo, xc, yHi, '#222222', 1.2)]; %#ok<AGROW>
        content = [content svgLine(xc-4, yLo, xc+4, yLo, '#222222', 1.2)]; %#ok<AGROW>
        content = [content svgLine(xc-4, yHi, xc+4, yHi, '#222222', 1.2)]; %#ok<AGROW>
    end
    content = [content svgText(gx0+groupW/2, plotBottom+20, bwLabels{k}, 11, 'middle')]; %#ok<AGROW>
end
content = [content svgText((plotLeft+plotRight)/2, plotBottom+38, 'Target bandwidth', 12, 'middle', 'bold')];
legY = plotBottom + 62;
for t = 1:nT
    col = colors{mod(t-1,numel(colors))+1};
    lx = plotLeft + (t-1)*135;
    content = [content svgRect(lx, legY-11, 14, 14, col, '#333333', 0.8)]; %#ok<AGROW>
    content = [content svgText(lx+20, legY, sprintf('> %d ms', thresholds(t)), 11, 'start')]; %#ok<AGROW>
end
end

function svg = buildCombinedFigure(latByBW, bwLabels, bwColors, thresholds, propMat, ciLoMat, ciHiMat, threshColors, captionLines)
panelH = 470;
panelW_a = 660; % ECDF (needs room for legend)
panelW_b = 780; % exceedance (needs room for 6 bandwidth groups + legend)
gap = 30;
topPad = 30; labelRowH = 26; captionLineH = 20; bottomPad = 15;
W = panelW_a + gap + panelW_b;
H = topPad + labelRowH + panelH + bottomPad + numel(captionLines)*captionLineH + 10;

svg = sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" font-family="Arial,Helvetica,sans-serif">\n', W,H,W,H);
svg = [svg sprintf('<rect x="0" y="0" width="%d" height="%d" fill="white"/>\n', W,H)];

xA = 0; xB = panelW_a + gap; yPanels = topPad + labelRowH;
svg = [svg svgText(xA+10, topPad+18, '(a)', 15, 'start', 'bold')];
svg = [svg svgText(xB+10, topPad+18, '(b)', 15, 'start', 'bold')];

contentA = panelEcdf(latByBW, bwLabels, bwColors, panelW_a, panelH, 'Trial-level mean latency (ms)');
contentB = panelExceedance(bwLabels, thresholds, propMat, ciLoMat, ciHiMat, threshColors, panelW_b, panelH);

svg = [svg sprintf('<svg x="%d" y="%d" width="%d" height="%d" viewBox="0 0 %d %d">\n%s</svg>\n', xA, yPanels, panelW_a, panelH, panelW_a, panelH, contentA)];
svg = [svg sprintf('<svg x="%d" y="%d" width="%d" height="%d" viewBox="0 0 %d %d">\n%s</svg>\n', xB, yPanels, panelW_b, panelH, panelW_b, panelH, contentB)];

capY = yPanels + panelH + 30;
for i = 1:numel(captionLines)
    isFirst = (i==1);
    isCaveat = (i>=4);
    weight = 'normal'; color = '#1a1a1a'; sz = 12;
    if isFirst, weight = 'bold'; end
    if isCaveat, color = '#555555'; sz = 11; end
    svg = [svg svgText(10, capY + (i-1)*20, captionLines{i}, sz, 'start', weight, color)]; %#ok<AGROW>
end

svg = [svg '</svg>'];
end
