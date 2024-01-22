#!/usr/bin/env bash
#{
/usr/bin/env octave-cli -q $0 $@
exit
#}

honeywellNaviTXTFile = '20230802/nav/ATLANS-20230802-161326_OutC_POSTPROCESSING-replay.xpf.txt'
filename = '20230802/nav/ATLANS-20230802-161326_OutC_POSTPROCESSING-replay.xpf.txt_plus18sec'

disp('no adjusted nav file. creating...')
nav_data = load(honeywellNaviTXTFile);
nav_data(:,1) = nav_data(:,1)+18;
dlmwrite(fullfile(filename),nav_data,'delimiter',',','precision','%.12f');

