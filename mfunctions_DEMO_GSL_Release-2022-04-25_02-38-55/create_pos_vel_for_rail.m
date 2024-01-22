function [railParam, pos_array, vel_array] = create_pos_vel_for_rail(inRawParFilename,outPosFile,outVelFile)
% create_pos_vel_for_rail calculates the position and velocity at each pulse
%  measured along the rail and writes the the position and velocity to 
%  ascii files to be used with the TDBP processor
%
%   USAGE:
%       [railParam, pos_array, vel_array] = create_pos_vel_for_rail(inRawParFilename,outPosFile,outVelFile)
%
%
%   SEE ALSO:
%       create_slc_par_for_gsl_rail.m, create_polar_rec_grid_for_gsl_rail.m 
%
%
%   Created:    2018-07-09 by Othmar Frey <frey@gamma-rs.ch>
%
%   Copyright:  2019 Gamma Remote Sensing AG
%               Othmar Frey <frey@gamma-rs.ch>
%

verbose = 0;

% Read GS-L raw parameter file
rawParStruct  = readGammaParFile(inRawParFilename);


start_pos = rawParStruct.RAIL_start_pos;
end_pos = rawParStruct.RAIL_end_pos; 
vel = rawParStruct.RAIL_velocity; 
acc = rawParStruct.RAIL_acceleration; 

s_tot= end_pos - start_pos;	% total length of the distance traveled
if(s_tot < 0.0)
    error(['ERROR RAIL motion params: end rail position is less than start rail position';
           'start position mm: %d   end position mm: %d',start_pos,end_pos])
end
if(vel <= 0.0)
    error('ERROR RAIL motion params: rail target velocity <= 0.0')
end
if(acc <= 0.0)
    error('ERROR RAIL motion params: rail acceleration <= 0.0')
end
  
tacc = vel/acc;             % time to accelerate to target velocity
sacc = 0.5*acc*tacc*tacc;   % distance traveled until target velocity is reached
scv = s_tot - 2.0*sacc;     % distance traveled at the target velocity
if (scv >= 0.0)
    tcv = scv/vel;          % time spent traveling at the target velocity
    t_tot = 2.0*tacc + tcv;
else
    printf('*** WARNING RAIL motion params : rail car does not reach target velocity ***')
    printf('start position mm: %d   end position mm: %d   velocity mm/s: %.3f  acceleration mm/s**2: %.3f', start_pos, end_pos, vel, acc);
    scv = 0.0;
    tcv = 0.0;
    sacc = s_tot/2.0;          % assume that acceleration and deceleration distances are equal
    tacc = sqrt(2.0*s0/acc);
    t_tot = 2.0*tacc;          % time to accelerate and decelerate
end
  
railParam.start_pos   = start_pos;
railParam.end_pos     = end_pos;
railParam.vel         = vel;
railParam.acc         = acc;
railParam.sacc        = sacc;
railParam.tacc        = tacc;
railParam.scv         = scv;
railParam.tcv         = tcv;
railParam.s_tot       = s_tot;
railParam.t_tot       = t_tot;

 
% Note : rawParStruct.RAW_IPP   is NOT the interpulse period,
% Do not use it!
% Use the PRI instead
PRI = rawParStruct.RAW_samp_IPP/rawParStruct.ADC_sample_rate;

% time and position array for the accelerating part
t1 = (0:1:floor(tacc/PRI))* PRI;
pos_acc = 0.5*acc*t1.^2;
vel_acc = acc*t1;
  
% time and position array for the constant velocity part
t2 = (1:1:floor(tcv/PRI))*PRI + (t1(end)-tacc);
pos_const_vel = vel.*t2 + sacc;
vel_const_vel = vel.*ones(size(pos_const_vel));

% time and position array for the decelerating part 
t3 = (1:1:floor(tacc/PRI)) *PRI + (t2(end)-tcv);
pos_dec =  -0.5*acc*(t3-tacc).^2 + 2*sacc + scv;
vel_dec =      -acc*(t3-tacc); 

time_vec = [t1 (t2+tacc) (t3+tacc+tcv)];
pos_azimuth = [pos_acc pos_const_vel pos_dec];
vel_azimuth = [vel_acc vel_const_vel vel_dec];

% center azimuth coordinate around average value
pos_azimuth = pos_azimuth - mean(pos_azimuth);

% create arrays of zeroes for the other components of Pos and Vel
pos_range = zeros(size(pos_azimuth));
pos_elev  = zeros(size(pos_azimuth));
vel_range = zeros(size(pos_azimuth));
vel_elev = zeros(size(pos_azimuth));

%whos
if(verbose)
    figure
    plot(t1,pos_acc)
    figure
    plot(t2+tacc,pos_const_vel)
    figure
    plot(t3+tacc+tcv,pos_dec)
    figure
    plot(time_vec,pos_azimuth)
    figure
    plot(diff(time_vec))
    figure
    plot(diff(pos_azimuth))
end

pos_array = ([pos_range.' pos_azimuth.' pos_elev.'])./1000.0;
vel_array = ([vel_range.' vel_azimuth.' vel_elev.'])./1000.0;

% write positioning file
if(exist(outPosFile, 'file')==2)
    error('File aleady exists! Abort writing pos file (sensor position).');
else
    dlmwrite(outPosFile, pos_array, 'delimiter', ',', 'precision', '%.12f');
end

% write velocity file
if(exist(outVelFile, 'file')==2)
    error('File aleady exists! Abort writing vel file (sensor velocity).');
else
    dlmwrite(outVelFile, vel_array, 'delimiter', ',', 'precision', '%.12f');
end
