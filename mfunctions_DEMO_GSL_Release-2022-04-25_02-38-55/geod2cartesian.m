function X = geod2cartesian(lat,lon,H,EARTH)
% geod2cartesian converts geodetic WGS84 coordinates
% latitude, lonitude [rad] and height [m] above WGS84 ellipsoid
% to cartesian WGS84 coordinates X,Y,Z [m]
%
% Usage:    X = geod2cartesian(lat,lon,H,EARTH)
%
% Input:    lat = latitude    [rad]
%           lon = lonitude   [rad]
%           H   = height      [m]
%
%          EARTH (see file simgui.m, function Calculate_pushbutton_Callback)
%
% Output:  Vector X = [Xcoord Ycoord Zcoord] [m] 
%

f      = EARTH.F_ELLIPSOID;        % WGS84
eEarth2= (2.0*f-f*f);
a      = EARTH.A_ELLIPSOID;       % (m) major semiaxis WGS84

X = [0.0, 0.0, 0.0];
if ( lat == 0.0 )
    l = a;    
else
    l = a/sqrt(1.0-eEarth2*sqr(sin(lat)));
end
rho = (l+H)*cos(lat);
z = (l*(1.0-eEarth2)+H)*sin(lat);
X = [ rho*cos(lon), rho*sin(lon), z ];
