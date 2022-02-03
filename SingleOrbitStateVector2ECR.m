%% *MultipleOrbitStateVectors2ECR.m*
%% *Purpose*
%  Converting the orbital elements of several satellites to the inertial
%  p-q orbital plane.
%% *History*
%  When       Who    What
%  ---------- ------ --------------------------------------------------
%  2021/07/06 hdai   editing code for personal use

%% Inputs
ID = 07376; % ISS Zarya: 25544; Molniya 3-50: 25847; Molniya 2-10: 07376
Hourf = 0; % [hrs] Desired time after TLE epoch
tstep = 0.1; % [hrs] time step for time evals
tend = 48; % [hrs] end of time evals
drawMileStones = 1;

%% *Standard Gravitational Parameter*
% The standard gravitational parameter $$\mu$$ of a celestial body is the
% product of the gravitational constant G and the mass M of the body. 
mu = 3.98618e14; % [m3/s2] Earth's geocentric gravitational constant

%% *Get the Satellite TLE*
[TLE] = getSatelliteTLE(ID);

%% *Convert to Orbital Elements*
[OE] = TLE2OrbitalElements(TLE);
fprintf(1,['Kepler Elements for satelliteID %d epoch %s:\n' ...
    '\ta [m] = %f\n\te = %f\n\ti [deg] = %f\n\tomega [deg] = %f\n' ...
    '\tOmega [deg] = %f\n\tM [deg] = %f\n'], floor(OE.satelliteID), ...
    datestr(OE.epoch),OE.a_km*1e3, OE.e, OE.i_deg, OE.omega_deg, ...
    OE.Omega_deg, OE.M_deg);
a_m = OE.a_km*1e3;
e = OE.e;
M_deg = OE.M_deg;

period = (2*pi*a_m^1.5)/sqrt(mu);
%% *Orbital Plane Coordinates*
%  p_m - [m] coordinate along axis through center and periapsis
%  q_m - [m] coordinate passing through focus and perpendicular to p-axis
%  E90 - Eccentric anomaly of 90 degrees of true anomaly
%  dpdt_m_per_s = [rad/s] p component velocity
%  dqdt_m_per_s = [rad/s] q component velocity

n_rad_per_s = sqrt(mu/a_m^3);  % [rad/s] mean motion
n_deg_per_s = rad2deg(n_rad_per_s); % [deg/s] mean motion
M_rad = deg2rad(M_deg);

E_rad = M_rad; 
dE = 99999;
eps = 1e-6; % [rad] control precision of Newton's method solution
while (abs(dE) > eps)
    dE = (E_rad - e * sin(E_rad) - M_rad)/(1 - e * cos(E_rad));
    E_rad = E_rad - dE;
end
p_m = a_m*(cos(E_rad) - e);
q_m = a_m*sqrt(1 - e^2)*sin(E_rad);
E90 = 2*atand(sqrt((1-e)/(1+e)));
dMdt_rad_per_s = n_rad_per_s;
dEdt_rad_per_s = dMdt_rad_per_s/(1 - e*cos(E_rad));
dpdt_m_per_s = -a_m*sin(E_rad)*dEdt_rad_per_s;
dqdt_m_per_s = a_m*cos(E_rad)*dEdt_rad_per_s*sqrt(1 - e^2);
E_deg_epoch = rad2deg(E_rad); 

%% *Rotate To ECI*
Rz_Omega = [ ...
    [cosd(OE.Omega_deg) sind(OE.Omega_deg) 0]; ...
    [-sind(OE.Omega_deg) cosd(OE.Omega_deg) 0]; ...
    [0 0 1]];
Rx_i = [ ...
    [1 0 0]; ...
    [0 cosd(OE.i_deg) sind(OE.i_deg)]; ...
    [0 -sind(OE.i_deg) cosd(OE.i_deg)]];
Rz_omega = [ ...
    [cosd(OE.omega_deg) sind(OE.omega_deg) 0]; ...
    [-sind(OE.omega_deg) cosd(OE.omega_deg) 0]; ...
    [0 0 1]];

% time of epoch
[Year,Month,Day,H,M,S] = datevec(OE.epoch);
HourUTC = H + M/60.0 + S/3600.0; 
jd = juliandate(Year,Month,Day,HourUTC,0,0);
jd0 = juliandate(Year,Month,Day,0,0,0);
% from time in Julian centuries from J2000
T = (jd - 2451545.0d0)./36525.0d0;
D0 = (jd0 - 2451545.0d0);
% % [deg] GMST = Sidereal Time at Greenwich
GMST = 6.697374558 + 0.06570982441908*D0 + 1.00273790935*HourUTC + ...
    0.000026*T^2;
% % [deg] Sidereal Time
%ST_deg = 100.46061837 + 36000.770053608*T + 0.000387933*T^2 - ...
    %T^3/38710000;
%GST_deg = ST_deg + 1.00273790935*HourUTC*15;
%LST_deg = GST_deg - ObsLongitude_degW;
% [rad] Earth rotation angle in relation to UT1
Theta_deg = 100.460618375 + 36000.770053608336*T + 0.0003879333*T^2 + ...
    15*H + M/4 + mod(S/240,360);

Rz_hour = [ ...
    [cosd(Theta_deg) sind(Theta_deg) 0]; ...
    [-sind(Theta_deg) cosd(Theta_deg) 0]; ...
    [0 0 1]];

% position of satellite at epoch in the orbit pq plane
r_pq = [p_m q_m 0]';

% position of satellite at epoch in ECI coordinates
omega_deg = OE.omega_deg;
Omega_deg = OE.Omega_deg;
i_deg = OE.i_deg;

%{
px = cosd(omega_deg)*cosd(Omega_deg) - sind(omega_deg)*cosd(i_deg)*sind(Omega_deg);
py = cosd(omega_deg)*sind(Omega_deg) + sind(omega_deg)*cosd(i_deg)*cosd(Omega_deg);
pz = sind(omega_deg)*sind(i_deg);

qx = -sind(omega_deg)*cosd(Omega_deg) - cosd(omega_deg)*cosd(i_deg)*sind(Omega_deg);
qy = -sind(omega_deg)*sind(Omega_deg) + cosd(omega_deg)*cosd(i_deg)*cosd(Omega_deg);
qz = cosd(omega_deg)*sind(i_deg);

wx = sind(i_deg)*cosd(omega_deg);
wy = -sind(i_deg)*sind(Omega_deg);
wz = cosd(i_deg);

r_ECI = [(p_m*px+q_m*qx) (p_m*py+q_m*qy) (p_m*pz+q_m*qz)];
%}

r_ECI = inv(Rz_Omega)*inv(Rx_i)*inv(Rz_omega)*r_pq;
% r_ECR = Rz_hour*r_ECI;
% r_LLA = ecef2lla(r_ECI');
r_LLA = eci2lla(r_ECI',datevec(datenum(Year, Month, Day, H, M, S)),'IAU-2000/2006');

% disp(num2str(r_LLA(1)));
% disp(num2str(wrapTo180(atan2d(r_ECI(2),r_ECI(1))-Theta_deg)));
% disp(datestr(OE.epoch));

Evals = 0:1:360.0; % [deg] values of the eccentric anomaly around orbit 
Orbit_p = a_m*(cosd(Evals)-e); % [m] orbit positions
Orbit_q = a_m*sqrt(1 - e^2)*sind(Evals); % [m] orbit positions

% [s] time since epoch along orbit
deltaT_s = ((Evals-E_deg_epoch) - e*sind(Evals-E_deg_epoch))/n_deg_per_s;
%  qmark - [m] q on orbit
qmarkp = a_m*(cosd(E90)-e);
qmarkq = a_m*sqrt(1-e^2)*sind(E90);
qmark = (inv(Rz_Omega)*inv(Rx_i)*inv(Rz_omega)*[qmarkp qmarkq 0]')';
nqmark = (inv(Rz_Omega)*inv(Rx_i)*inv(Rz_omega)*[qmarkp -qmarkq 0]')';
% r_ECR_orbit = Rz_hour*r_ECI_orbit';
% r_LLA_orbit = ecef2lla(r_ECR_orbit');
% only do one at a time... <sigh>
Orbit_ECI = zeros(numel(deltaT_s),3);
Orbit_LLA = zeros(numel(deltaT_s),3);
for ipt = 1:size(Orbit_ECI,1)
    r_pq = [Orbit_p(ipt) Orbit_q(ipt) 0]';
    Orbit_ECI(ipt,:) = (inv(Rz_Omega)*inv(Rx_i)*inv(Rz_omega)*r_pq)'; %[Rz_Omega*Rx_i*Rz_omega*r_pq]';
    lla = eci2lla(Orbit_ECI(ipt,:),datevec(datenum(Year, Month, Day, H, M, S+deltaT_s(ipt))),'IAU-2000/2006');
    Orbit_LLA(ipt,:) = lla;
end

%% *Plot Earth*
grs80 = referenceEllipsoid('grs80', 'm');

longs = mod(Theta_deg:30:Theta_deg+330, 360)-180; %rotated longitudes
figure('Color', 'k', 'POSITION', [0 0 1000 1000])
ax = axesm('globe', 'Geoid', grs80);
    %'Grid', 'on', 'GlineWidth', 1, ...
    %'GlineStyle', '-', 'Gcolor', [0.5 0.5 0.5], 'Galtitude', 100, ...
    %'MLineLocation', longs);
    % Meridian labeling doesn't seem compatible with rotated globes :(
    %'MeridianLabel', 'on', 'MLabelParallel', 'equator', ...
    %'MLabelLocation', [-30, 0, 30]);
ax.Position = [0 0 1 1];
axis equal
view(3)

load topo
h1 = geoshow(topo, topolegend, 'DisplayType', 'texturemap');
demcmap(topo)
land = shaperead('landareas', 'UseGeoCoords', true);
h2 = plotm([land.Lat], [land.Lon], 'Color', 'black');
rivers = shaperead('worldrivers','UseGeoCoords',true);
h3 = plotm([rivers.Lat],[rivers.Lon],'Color','blue');
rotate(h1, [0,0,1], Theta_deg) %rotating map
rotate(h2, [0,0,1], Theta_deg) %rotating map
rotate(h3, [0,0,1], Theta_deg) %rotating map
set(gca,'visible','off')

%% *Plot 3-D Cartesian Milestones*

if drawMileStones == 1 % Draw periapsis, apoapsis, milestone anomalies, etc.
    % Orbit periapsis
    plot3(ax, Orbit_ECI(1,1), Orbit_ECI(1,2), Orbit_ECI(1,3), ...
        's', 'Color', 'b', 'Markersize', 5, 'MarkerFaceColor', '#D9FFFF')

    % 90 degree eccentric anomaly / line from 90 degrees to -90 degrees EA
    plot3(ax, Orbit_ECI(90,1), Orbit_ECI(90,2), Orbit_ECI(90,3), ...
        'h', 'Color', 'b', 'Markersize', 8, 'MarkerFaceColor', '#D9FFFF')
    ninetonine = [Orbit_ECI(90,1), Orbit_ECI(90,2), Orbit_ECI(90,3); ...
                  Orbit_ECI(270,1), Orbit_ECI(270,2), Orbit_ECI(270,3)]';
    line(ax, ninetonine(1,:), ninetonine(2,:), ninetonine(3,:), ...
        'Color', 'cyan', 'LineWidth', 1)

    % 90 degrees true anomaly / line from 90 degrees to -90 degrees t.a.
    plot3(ax, qmark(1), qmark(2), qmark(3), ...
        'p', 'Color', 'b', 'Markersize', 8, 'MarkerFaceColor', '#D9FFFF')
    q2nq = [nqmark(1), nqmark(2), nqmark(3); ...
           qmark(1),qmark(2),qmark(3)]';
    line(ax, q2nq(1,:), q2nq(2,:), q2nq(3,:), ...
        'Color', 'cyan', 'LineWidth', 1)

    % Line from periapsis to apoapsis
    p2a = [Orbit_ECI(1,1),Orbit_ECI(1,2),Orbit_ECI(1,3); ...
           Orbit_ECI(181,1),Orbit_ECI(181,2),Orbit_ECI(181,3)]';
    line(ax, p2a(1,:), p2a(2,:), p2a(3,:), ...
        'Color', 'cyan', 'LineWidth', 1)
end

% 3-D Graph Labels and Settings
xlabel('ECI x [m]')
ylabel('ECI y [m]')
zlabel('ECI z [m]')
wish = 0.9*ones(1,3);
set(gca, 'Xcolor', wish, 'Ycolor', wish, 'Zcolor',wish);
xticks('auto')
yticks('auto')
zticks('auto')
maxr = sqrt((Orbit_ECI(181,1))^2 +(Orbit_ECI(181,2))^2 + (Orbit_ECI(181,3))^2);
xlim([-maxr,maxr])
ylim([-maxr,maxr])
zlim([-maxr,maxr])
title('Satellite Orbit in ECI Coordinates');
grid on
movegui(gca, 'center')

%% *Plot 2-D Ground Track*
%{
figure('color','k');
earth = imread('ear0xuu2.jpg');
lv= size(earth,1);
lh= size(earth,2);
lats =  (1:lv)*180/lv - 90;
lons =  (1:lh)*360/lh - 180;
image(lons, -lats, earth)
hold on;
set(gca,'ydir','normal');
grid on
plot(Orbit_LLA(:,2),Orbit_LLA(:,1),'.r');
plot(r_LLA(2), r_LLA(1), 'p', 'MarkerFaceColor', [1 0 0.7], ...
    'MarkerEdgeColor', 'none', 'Markersize', 14);
set(gca,'XTick',-180:30:180);
set(gca,'YTick',-90:30:90);
set(gca,'Xcolor',0.9*ones(1,3));
set(gca,'Ycolor',0.9*ones(1,3));
title(['Ground Track for Epoch ' datestr(OE.epoch)])
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
axis equal
axis tight
set(gca,'fontweight','bold');
movegui(gca, 'east')
%}

%% 3-D Animation calculations
tevals = 0:tstep:tend; % [hrs] times to evaluate animation at
Tperiod = (2*pi/sqrt(mu))*a_m^(1.5); % [s] period of orbit
Mevals = M_rad+(2*pi/Tperiod)*tevals*3600; % [rad] mean anomalies for all times
Eevals = zeros(numel(tevals),1);
for k = 1:numel(tevals)
    Eevals(k) = Mevals(k); 
    dEev = 99999;
    eps = 1e-6; % [rad] control precision of Newton's method solution
    while (abs(dEev) > eps)
        dEev = (Eevals(k) - e * sin(Eevals(k)) - Mevals(k))/(1 - e * cos(Eevals(k)));
        Eevals(k) = Eevals(k) - dEev;
    end
end
pevals = a_m*(cos(Eevals) - e);
qevals = a_m*sqrt(1 - e^2)*sin(Eevals);
Thetaevals = 100.460618375 + 36000.770053608336*T + 0.0003879333*T^2 + ...
    15*(H+tevals) + M/4 + mod(S/240,360);
r_pqevals = [pevals'; qevals'; zeros(1,numel(pevals))]';
Orbit_pevals = a_m*(cos(Eevals)-e); % [m] orbit positions
Orbit_qevals = a_m*sqrt(1 - e^2)*sin(Eevals); % [m] orbit positions
Orbit_ECIevals = zeros(numel(tevals),3);
for ipteval = 1:size(Orbit_ECIevals,1)
    r_pqeval = [Orbit_pevals(ipteval) Orbit_qevals(ipteval) 0]';
    Orbit_ECIevals(ipteval,:) = (inv(Rz_Omega)*inv(Rx_i)*inv(Rz_omega)*r_pqeval)';
end
orbtail = ceil(period/(3600*tstep));

%% Animate
tvallen = length(tevals);
orbplots = gobjects(numel(tevals),3);
for t = 1:tvallen
    satpos = plot3(ax, Orbit_ECIevals(t,1), Orbit_ECIevals(t,2), ...
            Orbit_ECIevals(t,3), 'o', 'Color', 'b', 'Markersize', 3, ...
            'MarkerFaceColor', '#D9FFFF');
    c2satline = line(ax, [Orbit_ECIevals(t,1),0], ...
        [Orbit_ECIevals(t,2),0], [Orbit_ECIevals(t,3),0], 'Color', ...
        'yellow', 'LineWidth', 1);
    orbplots(t) = plot3(ax, Orbit_ECIevals(t,1), ...
                        Orbit_ECIevals(t,2), ...
                        Orbit_ECIevals(t,3), ...
                        '.', 'Color', 'r', 'Markersize', 1);
    rotate(h1, [0,0,1], 15*tstep) %rotating map
    rotate(h2, [0,0,1], 15*tstep) %rotating map
    rotate(h3, [0,0,1], 15*tstep) %rotating map
    pause(0.04)
    if t ~= tvallen
        delete(satpos)
        delete(c2satline)  
    end
    if t > orbtail
        delete(orbplots(t-orbtail))
    end
end

% % Satellite orbit
% plot3(ax, Orbit_ECI(:,1), Orbit_ECI(:,2), Orbit_ECI(:,3), '-o', ...
%     'Color', 'r', 'Markersize', 2)


