function [lat,lon,H] = cartesian2geod(X,EARTH)
% cartesian2geod converts cartesian WGS84 coordinates X,Y,Z [m]
% to geodetic WGS84 coordinates latitude, longitude [rad] and
% height [m] above WGS84 ellipsoid
%
% Usage:  [lat,lon,H] = cartesian2geod(X,EARTH)
%
% Input:  Vector X    = [Xcoord Ycoord Zcoord]
%
%         EARTH         (see file simgui.m, function Calculate_pushbutton_Callback)
%
% Output: lat   =    latitude    [rad]
%         lon   =    longitude   [rad]
%         H     =    height      [m]
%

DEG = pi/180;
EPS = 1E-8;
MAXIT = 10;
f     = EARTH.F_ELLIPSOID;    % WGS84 
eEarth2= (2.0*f-f*f);
a      = EARTH.A_ELLIPSOID;   %(m) major semiaxis WGS84

rb     = a*(1-f);   
% rb is only used to calculate an initial value for 'H'

rho = norm([X(1) X(2)]);
if (rho ~= 0.0)
    tanPHI = X(3)/rho;
else
    tanPHI = 10E20;
end
r = norm(X);
H = r - rb;
b = asin(X(3)/r);  
k = a*(1-eEarth2);
i = 0;
b0 = 10E10;
H0 = 10E10; 
while(((abs(b-b0) >= EPS) | (abs(H-H0) >= EPS)) && (i <= MAXIT));
    b0 = b;
    H0 = H; 
    l = a/sqrt(1.0-eEarth2*sqr(sin(b0)));
    if(abs(b) < EPS)
        b = 0.0;
        H = rho - a;       
    else
        b = atan2(tanPHI*(l+H),(1-eEarth2)*l+H0);
        H = X(3)/sin(b) - k/sqrt(1.0- eEarth2*sqr(sin(b)));
    end
    i = i+1;
end    
if(i > MAXIT)
    lat = 0.0;
    lon = 0.0;
    H = 0.0;
else
    lat = b; 
    lon = atan2(X(2),X(1));
end
