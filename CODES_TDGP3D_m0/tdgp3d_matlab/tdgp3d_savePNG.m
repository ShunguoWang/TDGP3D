function tdgp3d_savePNG(fig, file)
% TDGP3D_SAVEPNG  Save a figure of the tdgp3d_* plotting functions as a PNG (150 dpi), then show it.
%
%   fig = figure('Visible', 'off', 'Position', [...]);   % build the figure invisible
%   ...
%   tdgp3d_savePNG(fig, 'name.png')                       % saved at the size it was built with
%
% Some of these figures are larger than a laptop screen (tdgp3d_stationPDF with 100 sites asks
% for about 1670 x 1550 pixels). Shown on the screen, such a figure is shrunk to fit it, and in the
% R2025 MATLAB desktop print(fig, file, '-dpng') of it can stop with "Output file ... was not
% created" although the file is written, with an axes toolbar in it. So the figure is saved while
% invisible, at the size it was built with (remembered, so that a later save -- tdgp3d_plotNodes
% adds nodes to the figures of tdgp3d_plotModel and saves them again -- gets the same size), with
% the axes toolbars hidden; the image comes from print(..., '-RGBImage') and is written with
% imwrite. The figure is shown afterwards.

sz = getappdata(fig, 'tdgp3dSize');
if isempty(sz)
    sz = get(fig, 'Position');
    setappdata(fig, 'tdgp3dSize', sz);
end
set(fig, 'Visible', 'off');
pos = get(fig, 'Position');
set(fig, 'Position', [pos(1:2) sz(3:4)]);

ax = findall(fig, 'Type', 'axes');
tb = cell(numel(ax), 1);
for i = 1:numel(ax)
    if isprop(ax(i), 'Toolbar') && ~isempty(ax(i).Toolbar)
        tb{i} = ax(i).Toolbar.Visible;
        ax(i).Toolbar.Visible = 'off';
    end
end
drawnow;
img = print(fig, '-RGBImage', '-r150');
for i = 1:numel(ax)
    if ~isempty(tb{i})
        ax(i).Toolbar.Visible = tb{i};
    end
end

imwrite(img, file);
set(fig, 'Visible', 'on');

end
