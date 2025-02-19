%% raincloud_plot - plots a combination of half-violin, boxplot,  and raw
% datapoints (1d scatter).
% Use as h = raincloud_plot(X), where X is a data vector is a cell array of handles for the various figure parts.
% See below for optional inputs.
% Based on https://micahallen.org/2018/03/15/introducing-raincloud-plots/
% Inspired by https://m.xkcd.com/1967/
% v1 - Written by Tom Marshall. www.tomrmarshall.com
% v2 - Updated inputs to be more flexible - Micah Allen 12/08/2018
%
% Thanks to Jacob Bellmund for some improvements


function [h, u] = raincloud_plot_smooth(X, varargin)

% ---------------------------- INPUT ----------------------------
%
% X - vector of data to be plotted, required.
%
% --------------------- OPTIONAL ARGUMENTS ----------------------
%
% color             - color vector for rainclouds (default gray, i.e. = [.5 .5 .5])
% band_width         - band_width of smoothing kernel (default = 1)
% density_type       - choice of density algo ('ks' or 'rath'). Default = 'ks'
% box_on             - logical to turn box plots on/off (default = 0)
% box_dodge          - logical to turn on/off box plot dodging (default = 0)
% box_dodge_amount    - mutiplicative value to increase dodge amount (default = 0)
% alpha             - scalar positive value to increase cloud alpha (defalut = 1)
% dot_dodge_amount    - scalar value to increase dot dodge amounts (defalut =0.6)
% box_col_match       - logical to set it so that boxes match the colour of clouds (default = 0)
% line_width         - scalar value to set global line width (default = 2)
% lwr_bnd        - mutiplicative value to increase spacing at bottom of plot(default = 1)
% bxcl           - color of box outline
% bxfacecl       - color of box face
%
% ---------------------------- OUTPUT ----------------------------
% h - figure handle to change more stuff
% u - parameter from kernel density estimate
%
% ------------------------ EXAMPLE USAGE -------------------------
%
% h = raincloud('X', myData, 'box_on', 1, 'color', [0.5 0.5 0.5])
%

%% check all the inputs and if they do not exist then revert to default settings
% input parsing settings
p = inputParser;
p.CaseSensitive = true;
p.Parameters;
p.Results;
p.KeepUnmatched = true;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

% set the desired and optional input arguments
addRequired(p, 'X', @isnumeric);
addOptional(p, 'color', [0.5 0.5 0.5], @isnumeric)
addOptional(p, 'band_width', [])
addOptional(p, 'density_type', 'ks', @ischar)
addOptional(p, 'box_on', 1, @isnumeric)
addOptional(p, 'box_dodge', 0, @isnumeric)
addOptional(p, 'box_dodge_amount', 1, @isnumeric)
addOptional(p, 'alpha', 1, validScalarPosNum)
addOptional(p, 'dot_dodge_amount', 0.4, @isnumeric)
addOptional(p, 'box_col_match', 0, @isnumeric)
addOptional(p, 'line_width', 2, validScalarPosNum)
addOptional(p, 'line_color', 0, @isnumeric)
addOptional(p, 'lwr_bnd', 1, @isnumeric)
addOptional(p, 'bxcl', [0 0 0], @isnumeric)
addOptional(p, 'bxfacecl', [1 1 1], @isnumeric)
addOptional(p, 'cloud_edge_col', [0 0 0], @isnumeric)


% parse the input
parse(p,X,varargin{:});
% then set/get all the inputs out of this structure
X                   = p.Results.X;
color               = p.Results.color;
density_type        = p.Results.density_type;
box_on              = p.Results.box_on;
box_dodge           = p.Results.box_dodge;
box_dodge_amount    = p.Results.box_dodge_amount;
alpha               = p.Results.alpha;
dot_dodge_amount    = p.Results.dot_dodge_amount;
box_col_match       = p.Results.box_col_match;
line_width          = p.Results.line_width;
line_color          = p.Results.line_color;
lwr_bnd             = p.Results.lwr_bnd;
bxcl                = p.Results.bxcl;
bxfacecl            = p.Results.bxfacecl;
cloud_edge_col      = p.Results.cloud_edge_col;
band_width          = p.Results.band_width;

% calculate kernel density
switch density_type
    case 'ks'

        counts = hist(X,300);
   B =    smoothdata([ 0 0 0 0 0 0 0 0 0 0 counts/length(X) 0 0 0 0 0 0 0 0 0 0 ],'gaussian',15);
    xx=linspace(min(X),max(X),300);
    steps=xx(2)-xx(1);
    C=linspace(min(X)-10*steps, max(X)+10*steps,320);
end

% density plot
h{1} = plot(B); hold on

% [a,b]=findpeaks(B);
%         c=[a;b];
%         cc=sortrows(c')
%  hold on
%         for s=1:size(cc,1)
%         plot(C(cc(s,2)),cc(s,1), 'o','Color', 'r' )
%         end
        
        
 % set(h{1}, 'FaceColor', color);
% set(h{1}, 'EdgeColor', cloud_edge_col);
set(h{1}, 'LineWidth', line_width);
set(h{1}, 'Color', line_color);
% set(h{1}, 'FaceAlpha', alpha);

ylim([0 max(B)+0.005])

% make some space under the density plot for the boxplot and raindrops
yl = get(gca, 'YLim');
set(gca, 'YLim', [-yl(2)*lwr_bnd yl(2)]);
yline(0)
% width of boxplot
wdth = yl(2) * 0.25;
 %set(gca, 'XScale', 'log')
 xlim([1000 10000]);

% jitter for raindrops
jit = (rand(size(X)) - 0.5) * wdth;

% info for making boxplot
quartiles   = quantile(X, [0.25 0.75 0.5]);
iqr         = quartiles(2) - quartiles(1);
Xs          = sort(X);
whiskers(1) = min(Xs(Xs > (quartiles(1) - (1.5 * iqr))));
whiskers(2) = max(Xs(Xs < (quartiles(2) + (1.5 * iqr))));
Y           = [quartiles whiskers];

% raindrops
if box_dodge
    drops_pos = (jit * 0.6) - yl(2) * dot_dodge_amount;
else
    drops_pos = jit - yl(2) / 2;
end

%h{2} = scatter(X, drops_pos);

h{2} = scatter_kde_horizontal(X, drops_pos);
     %colorbar
h{2}.SizeData = 10;
% h{2}.MarkerFaceColor = color;
% h{2}.MarkerEdgeColor = 'k';

if box_on
    
    if box_col_match
        
        bxcl = color;
        
    end
    
    if box_dodge
        box_pos = [Y(1) ((-yl(2) * box_dodge_amount) - (wdth * 0.3)) Y(2) - Y(1) (wdth * 0.6)];
        
        % median line
        h{4} = line([Y(3) Y(3)], [((-yl(2) * box_dodge_amount) - (wdth * 0.3)) ((-yl(2) * box_dodge_amount) + (wdth * 0.3))], 'col', bxcl, 'LineWidth', line_width);
        
        % whiskers
        h{5} = line([Y(2) Y(5)], [(-yl(2) * box_dodge_amount) (-yl(2) * box_dodge_amount)], 'col', bxcl, 'LineWidth', line_width);
        h{6} = line([Y(1) Y(4)], [(-yl(2) * box_dodge_amount) (-yl(2) * box_dodge_amount)], 'col', bxcl, 'LineWidth', line_width);
    else
        box_pos = [Y(1) -yl(2)/2-(wdth * 0.5) Y(2)-Y(1) wdth];
        % median line
        h{4} = line([Y(3) Y(3)], [-yl(2)/2-(wdth * 0.5) -yl(2) / 2 + (wdth * 0.5)], 'col', bxcl, 'LineWidth', line_width);
        
        % whiskers
        h{5} = line([Y(2) Y(5)], [-yl(2)/2 -yl(2)/2], 'col', bxcl, 'LineWidth', line_width);
        h{6} = line([Y(1) Y(4)], [-yl(2)/2 -yl(2)/2], 'col', bxcl, 'LineWidth', line_width);
    end
    % 'box' of 'boxplot'
    h{3} = rectangle('Position', box_pos);
    set(h{3}, 'EdgeColor', bxcl)
    set(h{3}, 'LineWidth', line_width);
    %set(h{3}, 'FaceColor', bxfacecl);
    % could also set 'FaceColor' here, etc
%    ylim('tickaligned') % set the ylim to capture all data
end