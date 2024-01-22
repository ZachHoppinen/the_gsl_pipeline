function [X,Y,Z] = geod2cartesianForArrays(lat,lon,H,EARTH)
%  GEOD2CARTESIANFORARRAYS converts geodetic WGS84 coordinates
%	latitude, longitude [rad] and height [m] above WGS84 ellipsoid
%	to cartesian WGS84 coordinates X,Y,Z [m]. The function works
%	for arrays, too, as long as the vectors lat, lon and H have the
%	same length.
%
%	Usage:	[X,Y,Z] = geod2cartesianForArrays(lat,lon,H,EARTH)
%
%	Input:	lat = latitude	[rad]
%           	lon = lonitude	[rad]
%           	H   = height	[m]
%
%          	EARTH (see file simgui.m, function Calculate_pushbutton_Callback)
%
%	Output:	[X,Y,Z]		[m] 
% 
%  
%  	Created: 15 Feb. 2005 by Othmar Frey <ofrey@geo.unizh.ch>
% 

f      = EARTH.F_ELLIPSOID;        % WGS84
eEarth2= (2.0*f-f*f);
a      = EARTH.A_ELLIPSOID;       % (m) major semiaxis WGS84


if ( lat == 0.0 )
    l = a;    
else
    l = a./sqrt(1.0-eEarth2.*sqr(sin(lat)));
end
rho = (l+H).*cos(lat);

X = rho.*cos(lon);
Y = rho.*sin(lon);
Z = (l.*(1.0-eEarth2)+H).*sin(lat);
