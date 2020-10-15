clc
clear all
close all

% add fieldtrip in path
restoredefaultpath
%addpath D:\YananSun\MEG110\fieldtrip-20180825; % change path if necessary
addpath C:\Work\fieldtrip-20180825; % change path if necessary
ft_defaults
fprintf('\nAdd Fieldtrip into path.\n')   

% add scripts full path
%scriptpath='D:\YananSun\MMN\scripts_MMN_YSun';
scriptpath='C:\Work\scripts_MMN_YSun';
addpath(genpath(scriptpath));% change path if necessary
fprintf('\nAdd scripts into path.\n') 

% define the path of your data 
curPath = pwd %'E:\YananSun\MMN\MMN_analysed_Wei_YSun\MMN_ICA_done\'  
%addpath(pwd)
%% BATCH STARTS
%%if 'subject' is not defined, the script will automatically get the subjects from the
%%current folder.

%subject_rethm ={'3370',	'3374',	'3378',	'3380',	'3394',	'3402',	'3407',	'3409',	'3411',	'3412',	'3416',	'3418',	'3419',	'3421',	'3422',	'3423',	'3426',	'3427',	'3428',	'3429',	'3430',	'3433',	'3435',	'3438',	'3439',	'3440',	'3441',	'3443',	'3445',	'3448',	'3493',	'3505',	'3508',	'3573',	'3612',	'3664'}%36 subjects have ReTHM
subject_rethm = {'3372',	'3383',	'3391',	'3437',	'3494',	'3496',	'3531',	'3634',	'3652',	'3678',	'3679',	'3680'};% subjects with autism

subject = [subject_rethm];% Totally 36 subjects

%Subj_correctedTrig = {'3448'};
Subj_correctedTrig = {'3678''3372' };%  '3678'

bad_subject = {};%subjects have more than 1 block.
%bad_subject = {'3647'};%subjects with autism but without rethm.
if exist('subject')
 [logic,ind]=ismember (bad_subject,subject);
    act_subject= subject;

    if logic==1
        act_subject(ind)=[];
    end   
else 
y=dirdir(cd);
%act_subject = {y.name};
end
act_subject = {'3372'  '3678'} % '3428','3448'        '3380' '3416' '3440'
for i = 2:length(act_subject)
    cd([curPath,'\', act_subject{i}])
    %ME125_Phase2_step0 % Automatically detect saturations.
    %ME125_Phase2_Step1 % Automatically do filtering.
    %ME125_Phase2_Step2 % Manually reject artifacts and bad channels.
    %ME125_Phase2_Step3 % Automatically run ICA.
    %ME125_Phase2_Step4 % Manually remove ICs related to EOG.
    ME125_Phase2_Step5 % Automatically epoch the data.
end

%% scripts for stats at group level
%ME125_Master_SensorLevel_ClusterBasedStats_YSun
%ME125_Master_SensorLevel_ClusterBasedStats_Planar_YSun
