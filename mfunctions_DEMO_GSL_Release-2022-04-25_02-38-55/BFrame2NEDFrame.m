function [north, east, down] = BFrame2NEDFrame(xbody,ybody,zbody,roll,pitch,heading)
%  BFRAME2NEDFRAME transforms the 3D vector [xbody, ybody, zbody] from the B-Frame
%  	to the NED-Frame [north, east, down] using the Euler angles roll,pitch
%	and heading. The function works for arrays, too, as long as
%  	the length of the array is the same for all input parameters or
%  	the length of the array is the same for roll, pitch and heading and only one
%  	element is given for each of xbody, ybody and zbody.
%  	xbody,ybody and zbody are given in [m];
%  	roll,pitch and heading in [rad].
%
%  	Usage: [north, east, down] = BFrame2NEDFrame(xbody,ybody,zbody,roll,pitch,heading)
%  
% 	Definition of coordinate systems:
%		B-frame (Body frame):
%			x = longitudinal (roll) axis of the aircraft
%  			y = lateral (pitch) axis of the aircraft (pointing to the right when
%  			    looking from the rear)
%  			z = normal (yaw, heading) axis (pointing downwards)
%  
%  		NED-frame (North-East-Down frame) (equiv. to topocentric frame, only difference:
%  					 	North-East-Down instead of East-North-Up):
%  			xe = axis pointing to North
%  			ye = axis pointing to East
%  			ze = axis pointing down in nadir direction
%  
%	See also:
%		http://www.globalsecurity.org/space/library/report/1997/basicnav.pdf
%  
%	Created: 15 Feb. 2005 by Othmar Frey <ofrey@geo.unizh.ch>
%
%   Copyright: 2017 Gamma Remote Sensing AG
%              Othmar Frey <frey@gamma-rs.ch>
%  

	a11 =  cos(heading);
	a12 = (-sin(heading));
	a21 =  sin(heading);
	a22 =  cos(heading);
	
	b11 =  cos(pitch);
	b13 =  sin(pitch);
	b31 = (-sin(pitch));
	b33 =  cos(pitch);
	
	c22 =  cos(roll);
	c23 = (-sin(roll));
	c32 =  sin(roll);
	c33 =  cos(roll);
	
	
	north = a11.*b11.*xbody + (a12.*c22 + a11.*b13.*c32).*ybody + (a12.*c23 + a11.*b13.*c33).*zbody;
	east  = a21.*b11.*xbody + (a22.*c22 + a21.*b13.*c32).*ybody + (a22.*c23 + a21.*b13.*c33).*zbody;
	down  = b31.*xbody + (b33.*c32).*ybody + (b33.*c33).*zbody;
	
	
	
	
%	The following implementation is more illustrative than the one above.
%  	However, it does not work for input ARRAYS:
%  
%  	Mheading = [[cos(heading) -sin(heading) 0.0]; ...
%  		[sin(heading) cos(heading) 0.0]; ...
%  		[0.0 0.0 1.0]];
%  
%  	Mpitch = [[cos(pitch) 0.0 sin(pitch)]; ...
%  		  [0.0 1.0 0.0]; ...
%  		  [-sin(pitch) 0.0 cos(pitch)]];
%  
%  	Mroll = [[1.0 0.0 0.0]; ...
%  		[0.0 cos(roll) -sin(roll)]; ...
%  		[0.0 sin(roll) cos(roll)]];
%  	
%  	
%  	bodyCoord = [xbody; ybody; zbody];
%  	
%  	NEDCoord = Mheading*Mpitch*Mroll*bodyCoord;
%  	
%  	north = NEDCoord(1);
%  	east  = NEDCoord(2);
%  	down  = NEDCoord(3);
%  	