function mykmlwriteline(kmlFileName, lat,lon,H,geoidUndulation,lineWidth,colHEXBin)
% MYKMLWRITELINE writes a kml file (line, open polygon) from (1-D arrays of)
% lat (in deg) , lon (in deg) , H,
% geoid undulation (one value or 1-D array), line width (integer number)
% and HEXBIN color code for the line color.
%
% Examples for line color HEXBIN code are
%       red: ff0000ff
%
% The range of values for any one color is 0 to 255 (00 to ff).
% The order of expression is aabbggrr, where:
% aa=alpha (00 to ff);
% bb=blue (00 to ff);
% gg=green (00 to ff);
% rr=red (00 to ff).
% For alpha, 00 is fully transparent and
%            ff is fully opaque.
% Usage:
%   mykmlwriteline(kmlFileName, lat,lon,H,geoidUndulation,lineWidth,colHEXBin)
%
% SEE ALSO:
%   -
%
% Created:  Othmar Frey <frey@gamma-rs.ch>, 01. Jan 2020
%

    fid=fopen(kmlFileName,'wt');
    message = ferror(fid);
    if(~isempty(message))
        error(message);
        return;
    end

    fprintf(fid,'<?xml version="1.0" encoding="utf-8"?>\n');
    fprintf(fid,'<kml xmlns="http://www.opengis.net/kml/2.2">\n');
    fprintf(fid,'   <Document>\n');
    fprintf(fid,'      <name>%s</name>\n',kmlFileName);
    fprintf(fid,'      <Placemark>\n');
    fprintf(fid,'         <Snippet maxLines="0"> </Snippet>\n');
    fprintf(fid,'         <description> </description>\n');
    fprintf(fid,'         <name>Line 1</name>\n');
    fprintf(fid,'         <Style>\n');
    fprintf(fid,'            <LineStyle>\n');
    fprintf(fid,'               <width>%d</width>\n',lineWidth);
    fprintf(fid,'               <color>%s</color>\n',colHEXBin);
    fprintf(fid,'            </LineStyle>\n');
    fprintf(fid,'         </Style>\n');
    fprintf(fid,'         <LineString>\n');
    fprintf(fid,'            <altitudeMode>absolute</altitudeMode>\n');
    fprintf(fid,'            <coordinates>');
    if(length(geoidUndulation)==1)
        geoidUndulation = geoidUndulation*ones(length(lat),1);
    end
    for i=1:length(lat)
       fprintf(fid,' %.12f,%.12f,%.12f',lon(i),lat(i),H(i)-geoidUndulation(i));
    end
    fprintf(fid,'            </coordinates>\n');
    fprintf(fid,'         </LineString>\n');
    fprintf(fid,'      </Placemark>\n');
    fprintf(fid,'   </Document>\n');
    fprintf(fid,'</kml>\n');

    fclose(fid);
end
