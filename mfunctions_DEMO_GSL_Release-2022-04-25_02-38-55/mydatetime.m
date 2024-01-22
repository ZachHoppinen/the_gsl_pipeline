function t = mydatetime(timestring)
% MYDATETIME mimicks the following usage of the Matlab function datetime:
%   t = datetime(nameAndValue{2},'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSXXX', ...
%                             'TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss.SSSSSS');
%   MYDATETIME is used in readGammaParFile for compatibility with Octave.
%
%   USAGE:
%       t = mydatetime(timestring)
%   
%   Created:   frey@gamma-rs.ch, 10. Sep 2019
%   Modified:   frey@gamma-rs.ch., 20. Sep 2021
%               Added test to distinguish between the 
%                "new" format YYYY-MM-ddTHH:mm:ss.SSSSSSXXX as used in the GAMMA L-band SAR
%               and the legacy time format YYYYY-MM-dd HH:mm:ss.SSSSSSXXX used
%               for the GPRI-II Ku-band SAR mode
%
%   Copyright: 2021 Gamma Remote Sensing AG
%              Othmar Frey <frey@gamma-rs.ch>
%

if(strfind(timestring,'T'))  % test to distinguish between the new and the legacy time-and-date format
    val = strsplit(timestring,'T');
    else
     val = strsplit(timestring,' ');
end

val2 = strsplit(val{2},'+');
if(length(val2)==2)
    date_and_time = [strtrim(val{1}),' ',strtrim(val2{1})];
    val3 = strsplit(val2{2},':');
    hoursOffSetFromUTC = str2num(val3{1}) + str2num(val3{2})/60.0;
    t = datestr(datenum(date_and_time)-hoursOffSetFromUTC/24,0);
else
    val2 = strsplit(val{2},'-');
    if(length(val2)==2) 
        date_and_time = [strtrim(val{1}),' ',strtrim(val2{1})];
        val3 = strsplit(val2{2},':');
        hoursOffSetFromUTC = str2num(val3{1}) + str2num(val3{2})/60.0;
        t = datestr(datenum(date_and_time)+hoursOffSetFromUTC/24,0);
    else
        date_and_time = [val{1} ' ' val2{1}];
        t = datestr(datenum(date_and_time));
    end
end