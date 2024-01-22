function [dX,dY,dZ] = topocentric2geocentric(dnorthing,deasting,dup,lat,lon)
%  TOPOCENTRIC2GEOCENTRIC transforms a difference vector [dnorthing, deasting, dup]
%  	given in the topocentric-Cartesian coordinates [dX,dY,dZ] of the
%  	horizon system to the geocentric-Cartesian coordinates of an equator
%  	system (e.g. WGS84). The function works for arrays, too, as long as
%  	the lengths of dnorthing, deasting and dup remain the same.
%
%  
%	Usage:	[dX,dY,dZ] = topocentric2geocentric(dnorthing,deasting,dup,lat,lon)
%  
%	Input:	dnorthing,deasting,dup	[m]
%  		lat,lon			[rad]
%   
%  	Output:	dX,dY,dZ		[m]
%  
%  	Reference:
%  		Scriptum: Einfuehrung in die hoehere Geodaesie (H.G. Kahle) p.42
%	
%  
%  	Created: 15 Feb. 2005 by Othmar Frey <ofrey@geo.unizh.ch>
% 


    
    dX = (-sin(lat).*cos(lon)).*dnorthing	+ -sin(lon).*deasting	+ (cos(lat).*cos(lon)).*dup;
    dY = (-sin(lat).*sin(lon)).*dnorthing	+ cos(lon).*deasting	+ (cos(lat).*sin(lon)).*dup;
    dZ = cos(lat).*dnorthing						+ sin(lat).*dup;