%% detect saturations
close all
if isempty ( dir('*merged_B1_denoise_rethm.con'))
    [file] =  dir('*denoise.con'); % 
else
    [file] =  dir('*merged_B1_denoise_rethm.con');
end
    pathname = file.folder;
    filename = file.name;

clear sat;
[sat] = mq_detect_saturations_withRef(pathname,filename,0.01,'child','string_comp');

%if exist('TSPCA', 'dir')
%       rmdir('TSPCA')
%end
    
savepath = [pathname]

save([savepath,'\sat.mat'], 'sat');