%PLOTMULTIORBITS Convert the orbital elements of several satellites to
%the inertial p-q orbital plane
%
% Input:
%
%   ID : integer or array of integers
%       NORAD ID(s) of the desired satellite(s)
%       Examples: 25544, [25544, 07376]
%
%   drawMilestones : integer
%       turn milestone markers and axes on or off
%
% Output:
%   ax : map axes
%       axes containing a rotatable globe and satellite orbit(s)
%
% Authors:
%   Huanlin Dai
%   Meg Noah
%
% Last Modified September 1, 2021

%% IDs of Satellites

IDs = [25544, 07376];
drawMilestones = 1;

%% Create Figure and Globe Axes

grs80 = referenceEllipsoid('grs80', 'm');

figure('Renderer', 'opengl', 'Color', 'k', 'POSITION', [0 0 1000 1000])
ax = axesm('globe', 'Geoid', grs80);
ax.Position = [0 0 1 1];
axis equal
set(ax,'visible','off')
view(3)

load topo
h1 = geoshow(topo60c, topo60cR, 'DisplayType', 'texturemap');
demcmap(topo60c)
land = shaperead('landareas', 'UseGeoCoords', true);
h2 = plotm([land.Lat], [land.Lon], 'Color', 'black');
rivers = shaperead('worldrivers','UseGeoCoords',true);
h3 = plotm([rivers.Lat],[rivers.Lon],'Color','blue');

%% Calculate and Plot Orbits onto Globe Axes

for id = IDs
    addOrbit(ax, id, drawMilestones)
end

