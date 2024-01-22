function gammaParStructure = readGammaParFile(parFileName)
% READGAMMAPARFILE reads GAMMA style text parameter values from
%  the file parFileName, which are text parameter files such as:
%  
%  PROC_par
%  SAR_par
%  SLC_par
%  baseline_par
%  OFF_par
%  DIFF_par  
% 
%  The structure with all keywords and values is stored in the structure variable
%  gammaParStructure, which is the return value of the function.
%
% USAGE:
%  gammaParStructure = readGammaParFile(parFileName)
%
% SEE ALSO:
%  initPROCpar.m, writePROCpar.m  initSARpar.m, writeSARpar.m, initSLCpar.m, writeSLCpar.m
%  initBaselinepar.m, writeBaselinepar.m, initOFFpar.m, writeOFFpar.m, initDIFFpar.m, writeDIFFpar.m
%
%   Created:   frey@gamma-rs.ch, 27. Nov 2012
%   Modified:  frey@gamma-rs.ch, 14. Mar 2013 to also support:
%                                            - baseline_par
%                                            - OFF_par
%                                            - DIFF_par baseline parameter files
%
%   Modified:  frey@gamma-rs.ch, 03. Jun 2018 to also support the
%                                           Gamma-SAR L-band radar GS-L
%                                           raw_par format. in particular
%                                           the new ISO time format:
%                                           'yyyy-MM-dd''T''HH:mm:ss+00:00'
%
%   Modified:  frey@gamma-rs.ch, 03. Apr 2019 changed from 
%                   t = datetime(nameAndValue{2},'InputFormat','yyyy-MM-dd''T''HH:mm:ssXXX','TimeZone','UTC',...
%                    'Format','yyyy-MM-dd HH:mm:ss');
%                 to:
%                   t = datetime(nameAndValue{2},'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSXXX','TimeZone','UTC',...
%                    'Format','yyyy-MM-dd HH:mm:ss.SSSSSS');
%
%   Modified:  frey@gamma-rs.ch, 09. Sep 2019 changed from                  
%                   t = datetime(nameAndValue{2},'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSXXX','TimeZone','UTC',...
%                    'Format','yyyy-MM-dd HH:mm:ss.SSSSSS');
%                 to:
%                   t = mydatetime(nameAndValue{2});
%                 for compatibility with Octave.
%
%   Modified:  frey@gamma-rs.ch, 20. Sep 2021 included section 
%                 elseif(strcmp(fieldName,'time_start'))
%                  [nameAndValue] = regexp(textline,'time_start:','split')
%                  t = mydatetime(strtrim(nameAndValue{2}))
%                  gammaParStructure.(fieldName) = t;
%               to handle the legacy time and data format used in the GPRI-II raw_par file
%
%   Copyright: 2021 Gamma Remote Sensing AG
%              Othmar Frey <frey@gamma-rs.ch>
%

[fid, message]=fopen(parFileName,'rt');
  if(~isempty(message))
    error(sprintf('%s\n%s%s',parFileName,message));
    return;
  end

  % Read text file line by line until the end of file (EOF) is reached.
  % Each line is parsed and dynamically assigned to structure variable with
  % field name corresponding to the keywords in the structured text file.
  while(~feof(fid)) % Repeat the following command sequence until EOF is reached
    textline=fgetl(fid);
    if(~isempty(strfind(textline,':')))
      [nameAndValue] = regexp(textline,':','split');
      fieldName = strtrim(nameAndValue{1});
      if(isempty(fieldName) || ~isletter(fieldName(1))) % This check allows one to also read the parameters from SnowScat combined header/binary data
	    break
      end
      % Handle new date and time format 
      if(strcmp(fieldName,'RAW_start_time'))
          [nameAndValue] = regexp(textline,'RAW_start_time:','split');
          t = mydatetime(nameAndValue{2})
          gammaParStructure.(fieldName) = t;
      elseif(strcmp(fieldName,'time_start'))
          [nameAndValue] = regexp(textline,'time_start:','split')
          t = mydatetime(strtrim(nameAndValue{2}));
          gammaParStructure.(fieldName) = t;
      else
          if(length(fieldName)<=63 && length(fieldName)>=2)
            if(isletter(fieldName(1)) && isletter(fieldName(2)))
              if(~isempty(fieldName))
                [dummy] = regexp(fieldName,'\(','split'); % Handles the special case of baseline keyword initial_baseline(TCN) and  precision_baseline(TCN)' set to field: initial_baseline and precision_baseline
                fieldName = dummy{1};
                [temp] = regexp(fieldName,'/','split'); % Handles the special case of keyword channel/mode (instead of channel_mode)
                if(length(temp)>1)
                  fieldName = [strtrim(temp{1}) '_' strtrim(temp{2})];
                end
                for ii=1:length(fieldName)
                    if(isalpha_num(fieldName(ii)) ~= logical(1))
                     return; 
                  end    
                end    
                beforeSplit = strtrim(nameAndValue{2});
                [values] = regexp(beforeSplit,' ','split');
                if(isnan(str2double(values{1})))
                  % if the first token of the splitted right-hand side
                  % of the column is not a number then the field is treated as TEXT
                  gammaParStructure.(fieldName) = beforeSplit;               
                else
                  % otherwise it's a number or a set of numbers with potentially some text (entities like Hz, m/s) at the end
                  counter = 1;
                  foundValue = 0;
                  lengthOfvalues = length(values); 
                  for(i=1:length(values) )
                    if(~isnan(str2double(values{i})))
                      arrayOfNumbers(counter) =str2double(values{i});
                      counter=counter+1;
                    end
                  end
                  gammaParStructure.(fieldName) = arrayOfNumbers;               
                  clear arrayOfNumbers;
                  clear values;
                end
              end
            end
         end
      end
    end
  end
%Close the file
fclose(fid);
