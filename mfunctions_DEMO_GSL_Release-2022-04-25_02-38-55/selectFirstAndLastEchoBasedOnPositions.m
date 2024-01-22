function [arrayOfIndices] = selectFirstAndLastEchoBasedOnPositions(approxStartPos,approxEndPos,posFile,velFile, minVelThreshold, subaperture_length,outTxtFileWithIndices,verbose)
% selectFirstAndLastEchoBasedOnPositions is used in context of
%  carborne (or airborne) SAR focusing and repeat-pass InSAR.
%  It solves the problem of termining  a first and a last azimuth echo
%  index (and optionally the subaperture indices as calculated from equisdistant
%  curvilinear subsections of the sensor trajectory) where the problem arises
%  that (especially when driving on the road in traffic) the position of the
%  first and last echo, respectively are typically not at the same location
%  for repeat-pass SAR data acquisitions.
%
%  The function selectFirstAndLastEchoBasedOnPositions allows to determine 
%  - Finding the start and end azimuth index on a sensor trajectory defined
%    by two given positions as e.g. read from a map/Google Earth etc.
%  - subapertures within that track of a given length
%  - (approximately) common repeat-pass InSAR tracks (and spatially common sections of subapertures)
%  based on:
%   
%   1. the position and velocity vector files (pos and vel for each echo) of a sensor track
%       and given lat, lon , height coordinates of an approximate first and last position close to the track. 
%   2. a minimum velocity threshold to avoid false selections
%   3. a given subaperture length 
%
% The function calculates the shortest distance of the given start (and end) point position
%      to the normal planes (defined by the velocity vector)
%      at each position on the sensor trajectory
% The FIRST (and the LAST) valid azimuth echo indices are found by searching for the
%      shortest distance from the normal plane while ruling out any potential false
%      matches by evaluating only those array elements that have at the same time 
%      a velocity above a given threshold.
%
% Then, the indices for the subapertures are calculated
%    based on the cumulative distance along the track and the given
%    subaperture length.
%
% Eventually the array of indices are written to a text file:
%     - first and the last echo in the case of one full track,
%     - all indices corresponding to beginning and end of subapertures, otherwise.
%
%
% Usage:
%   [arrayOfIndices] = selectFirstAndLastEchoBasedOnPositions(approxStartPos,approxEndPos,posFile,velFile, minVelThreshold, subaperture_length,outTxtFileWithIndices,verbose)
%
%   where:
%       approxStartPos         : 3x1 array with lat(deg) lon(deg) ell_H as
%                                e.g. obtained from inspection the pos_kml track files using a
%                                manual selection of the approximate location where the SAR focusing
%                                should start.
%       approxEndPos           : likewise: 3x1 array with lat(deg) lon(deg) ell_H for the selected end position of the track
%       posFile                : position file (ascii) with WGS84 ECEF X,Y,Z coord. for each radar echo (as required by az_proc_tdbp_gpu)
%       velFile                : velocity file (ascii) with WGS84 ECEF v_X,v_Y,v_Z coord. for each radar echo (as required by az_proc_tdbp_gpu)
%       minVelThreshold        : minumum velocity threshold (e.g. 1 (m/s)) to rule out  
%       subaperture_length     : subaperture length in m (if set to 0.0 then only one full aperture is considered)
%       outTxtFileWithIndices  : text file that will contain the echo indices for the first and last echo (and inbetween subapertures if applicable)
%       verbose                : optional verbose parameter. If set to verbose = 1, abs. velocity and velocity components along azimuth are plotted to figures.
%                                verbose = 0 or omission of the verbose parameter will skip the plots.
%
% Created:  2020-01-16 frey@gamma-rs.ch
% Modified: 2020-02-24 frey@gamma-rs.ch
%                       Minor update of help text.
% Modified: 2021-03-18 frey@gamma-rs.ch
%                       Added check for lastEcho
% Modified: 2021-03-19 frey@gamma-rs.ch
%                       Added optional verbose parameter to plot (or not
%                       plot) velocity components.
%
% Modified: 2021-06-07 frey@gamma-rs.ch
%                       Added additional test to check whether velocity vectors
%                       and the vector between the given start and end position
%                       point to the same half space.
%


if(~exist('verbose','var'))
    verbose = 0;
else 
    if(verbose~=0 && verbose~=1)
        error('verbose must have values 0 or 1.')
    end
end 
    
EARTH.A_ELLIPSOID = 6378137.0;
EARTH.F_ELLIPSOID = 1.0/298.257223560;
EARTH.MU          = 3.986004418e14;

start_lat = approxStartPos(1)*(pi/180);
start_lon = approxStartPos(2)*(pi/180);
start_H   = approxStartPos(3);

end_lat = approxEndPos(1)*(pi/180);
end_lon = approxEndPos(2)*(pi/180);
end_H   = approxEndPos(3);


startPos = geod2cartesian(start_lat,start_lon,start_H,EARTH);
endPos = geod2cartesian(end_lat,end_lon,end_H,EARTH);

% Load position and velocity files
pos = load(posFile);
vel = load(velFile);

% Check consistency of pos and vel data
[numEntriesPos, numOfDimensionsPos] = size(pos);
[numEntriesVel, numOfDimensionsVel] = size(vel);
if(numEntriesPos~=numEntriesVel)    
    error('Number of entries in position and velocity files are not equal.');
end
if(numOfDimensionsPos~=3)    
    error('The dimension of the coordinate space of the position data must equal 3 (x,y,z).');
end
if(numOfDimensionsVel~=3)    
    error('The dimension of the coordinate space of the position data must equal 3 (x,y,z).');
end


% Plot the velocity over time for reference and layer
abs_vel = sqrt(vel(:,1).*vel(:,1) + vel(:,2).*vel(:,2) + vel(:,3).*vel(:,3));

if(verbose)
    figure;
    plot(abs_vel), title('Abs velocity along sensor trajectory');
    xlabel('Azimuth index');
    ylabel('Velocity [m/s]');
    grid on;
    hold off;
    
    figure;
    plot(vel), title('Velocity components along sensor trajectory');
    xlabel('Azimuth index');
    ylabel('Velocity [m/s]');
    grid on;
    hold off;
end

% 1.a) Calculate indices where the abs velocity is above a given treshold
ind_where_vel_above_tresh = abs_vel>minVelThreshold;

% 1.b) Calculate indices for which the velocity vectors points into the same
%      halfspace as the vector between given start and end position.
%      This is needed to avoid wrong selections if the data acquisions includes
%       u-turns or other highly non-linear features
vel_proj =  vel(:,1)./abs_vel.*(endPos(1)-startPos(1)) + ...
            vel(:,2)./abs_vel.*(endPos(2)-startPos(2)) + ...
            vel(:,3)./abs_vel.*(endPos(3)-startPos(3));

ind_same_directon = vel_proj>0;

% 2.a) calculate the shortest distance of the given start point position
%      to the normal planes (defined by the velocity vector)
%      at each position on the sensor trajectory
dist_from_plane =  vel(:,1)./abs_vel.*(startPos(1)-pos(:,1)) + ...
                    vel(:,2)./abs_vel.*(startPos(2)-pos(:,2)) + ...
                    vel(:,3)./abs_vel.*(startPos(3)-pos(:,3));

% 2.b) find the FIRST valid azimuth echo index by searching for the
%      shortest distance from the plan while ruling out any potential false
%      matches by evaluating only those array elements that have at the same time 
%      a velocity above a given threshold
%[value, firstEcho] = min(abs(dist_from_plane))
[firstEcho, ~] = find(abs(dist_from_plane)==min(abs(dist_from_plane(ind_where_vel_above_tresh & ind_same_directon))));

% 3.a) calculate the shortest distance of the given end point position to the normal planes (defined by the velocity vector)
%      at each position on the sensor trajectory
dist_from_plane =  vel(:,1)./abs_vel.*(endPos(1)-pos(:,1)) + ...
                   vel(:,2)./abs_vel.*(endPos(2)-pos(:,2)) + ...
                   vel(:,3)./abs_vel.*(endPos(3)-pos(:,3));
               
% 3.b) find the LAST valid azimuth echo index by searching for the
%      shortest distance from the plan while ruling out any potential false
%      matches by evaluating only those array elements that have at the same time 
%      a velocity above a given threshold
%[value, lastEcho] = min(abs(dist_from_plane))
[lastEcho, ~] = find(abs(dist_from_plane)==min(abs(dist_from_plane(ind_where_vel_above_tresh & ind_same_directon))));
if(lastEcho==numEntriesPos)
    lastEcho=lastEcho-1;
end

% 4. Calculate distance between sensor positions at each echo
%    and the cumulative distances 
%    between the selected indices firstEcho and last echo
diff_pos_vec = diff(pos);
dist_betw_echo_pos = sqrt(diff_pos_vec(firstEcho:lastEcho,1).^2 + diff_pos_vec(firstEcho:lastEcho,2).^2 + diff_pos_vec(firstEcho:lastEcho,3).^2);
cumulative_dist_along_track = cumsum(dist_betw_echo_pos);

% 5. Calculate number of subapertures based on total track distance 
%    and the given subaperture length in meter.
if(subaperture_length==0.0)
    numOfSubapertures = 1;
else 
    numOfSubapertures = floor(cumulative_dist_along_track(end)/subaperture_length)
    if(numOfSubapertures<2)
        error('The resulting number of subapertures is lower than 2. Decrease the input subaperture length!');
    end
end

% 6. Calculate the indices for the subapertures
%    based on the cumulative distance along the track and the given
%    subaperture length
if(subaperture_length==0.0)
     arrayOfIndices = [firstEcho,lastEcho];
else
    subaperture_number = floor(cumulative_dist_along_track./subaperture_length);
    [arrayOfIndices, ~] = find(diff(subaperture_number)==1);
    arrayOfIndices = [firstEcho; firstEcho+arrayOfIndices(:)];
end

% 7. write the array of indices to a text file:
%     - first and the last echo in the case of one full track,
%     - all indices corresponding to beginning and end of subapertures, otherwise.
dlmwrite(outTxtFileWithIndices, arrayOfIndices(:),'delimiter', ',', 'precision', '%d');   

