% Change orthogonal selection to happen in the main script so that
% dependent on the conditions being compared. And, so that if only one
% condition we won't apply a pre-selection.

% Updated process data to downsample to 100 hz after power transformation

% change the progress bar to in matlab text


% Move trial rejection into this script so that it'll be based on condition

% VERTEX_ALL will be the same across all freqs now, right, since doing the
% orthogonal pre-selection within this script? If yes, clean up in
% ProcessData.

% changed the y axis label to correct freq

% added in part to test interaction between factors and plot accordingly
% added in subplots for anterior, middle and posterior STG
% added in super title for the entire plot

clear; clc; close all
tic

% Setup paths - CHANGE
% dropboxpath = 'C:\Users\djbrang.UMROOT\Dropbox\'; % DB Windows
% % % % %dropboxpath = [filesep,'Users',filesep,'db',filesep,'Dropbox',filesep]; % DB MAC
dropboxpath = '/Users/zhewei/Dropbox (University of Michigan)'; % Cody MAC3

% Add ECoG Analysis path to scripts
addpath([dropboxpath,filesep,'ECoG_Analysis_Scripts']); % Add path for analysis scripts
addpath([dropboxpath,filesep,'ECoG_Analysis_Scripts',filesep,'dependencies']); % Add path for analysis scripts

% % % % % projectPathDir = [dropboxpath,filesep,'LME_AVSpeech',filesep]; % David
projectPathDir = [dropboxpath,filesep,'_Projects/Brang_LME_AVSpeech',filesep];% Cody

electRegPath = [dropboxpath,filesep,'Electrode_Registration',filesep]; %% Dropbox path to electrode labels
bipDataDir =[projectPathDir,'Preprocessed',filesep];  %% Directory that contains bipolar data.
behavDir = [projectPathDir,'Behavioral-Data',filesep];
outputDir =[projectPathDir,'Processed',filesep];  %% Directory that contains bipolar data.
codeDir =[projectPathDir,'Matlab_codes',filesep];  %% Directory that contains matlab code

% Analyze and plot group data on cvs_avg35_inMNI152
FS_Subject = 'cvs_avg35_inMNI152';
FS_Subject_Dir = [dropboxpath,filesep,'Electrode_Registration',filesep,FS_Subject,filesep]; % Created a copy of the fold so that FS not needed to be installed


% Setup directories - Should not need to change
cd(projectPathDir);
addpath(codeDir)
SubjPath = [electRegPath FS_Subject filesep];
[SUCCESS,MESSAGE,MESSAGEID] = mkdir(outputDir);
[SUCCESS,MESSAGE,MESSAGEID] = mkdir([outputDir,'/images']);



% Pre-Process Data?
PreProcess_Data = 0;
ReProcess_Vertices = 0; % not working

% Analysis Type:
% Analysis_Type = 'Voxel-Wise'; % Voxel-wise analysis
Analysis_Type = 'Time-Series'; % Time-series analysis


%%

%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   Statistics to be run    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%


% Karthik's LME paper used 1 and 4 mostly
LME_Comparison = 1; % Main effect of Condition
% LME_Comparison = 2; % Interaction of Condition x Frequency Range
% LME_Comparison = 3; % Interaction of Condition x Frequency Range x ROI (Time-series only)
% LME_Comparison = 4; % Interaction of Condition x Frequency Range x ROI x Time-points (Time-series only)
% LME_Comparison = 5; % Interaction of Condition x Time-points
% LME_Comparison = 6; % Interaction of Condition x Frequency Range x Time-points

% Comparison to be made
% Cond list - Follow the same structure as ECoG_Plotter
% SPECS.Conditions = {'Aud-Noise-McG-A','Aud-Clear-McG-B','Aud-Noise-Non-A','Aud-Clear-Non-B','Cong-Noise-McG-A','Cong-Clear-McG-B','Cong-Noise-Non-A','Cong-Clear-Non-B','Incong-Noise-McG-A','Incong-Clear-McG-B','Incong-Noise-Non-A','Incong-Clear-Non-B','VIS-Noise-McG-A','VIS-Clear-McG-B','VIS-Noise-Non-A','VIS-Clear-Non-B'};
% SPECS.CondNums = [1:16];

% LME comparing McG vs non-McG? Look for interaction. Differences based on response (Aud, Vis, mcg)
% Need to code comparisons that use behavior

% Conditions = {
% 'Aud-Noise-McG-A'
% 'Aud-Clear-McG-B'
% 'Aud-Noise-Non-A'
% 'Aud-Clear-Non-B'
% 'Cong-Noise-McG-A'
% 'Cong-Clear-McG-B'
% 'Cong-Noise-Non-A'
% 'Cong-Clear-Non-B'
% 'Incong-Noise-McG-A'
% 'Incong-Clear-McG-B'
% 'Incong-Noise-Non-A'
% 'Incong-Clear-Non-B'
% 'VIS-Noise-McG-A'
% 'VIS-Clear-McG-B'
% 'VIS-Noise-Non-A'
% 'VIS-Clear-Non-B'};
% SPECS.CondNums = [1:16];

% Stat_Comp = [5:8;1:4]; Stat_Comp_Label = 'V'; % Vis Alone
Stat_Comp = [5:8;1:4]; Stat_Labels = {'Cong','Aud'};Stat_Comp_Label = 'Cong_vs_Aud'; % Aud vs Cong
% Stat_Comp = [9:12;1:4]; Stat_Labels = {'InCong','Aud'};Stat_Comp_Label = 'Incong_vs_Aud'; % Cong vs Incong
% Stat_Comp = [5:8;9:12]; Stat_Labels = {'Cong','Incong'};Stat_Comp_Label = 'Cong_vs_Incong'; % Cong vs Incong
% Stat_Comp = [1:2:12;2:2:12]; Stat_Comp_Label = 'Noise_vs_Clear'; % Noise vs Clear
% Stat_Comp = [1:4 1:4;5:12]; Stat_Comp_Label = 'A_vs_AV';  % A vs AV (cong + incong)
% Stat_Comp = [9 10;11 12]; Stat_Comp_Label = 'Incong_McG_vs_NonMcG';  % A vs AV (cong + incong)
% Stat_Comp = [1:4;13:16]; Stat_Comp_Label = 'A_vs_V';  % A vs Vis alone
% Stat_Comp = [5 7;1 3]; Stat_Comp_Label = 'Cong_vs_Aud_Noise'; 
% Stat_Comp = [6 8;2 4]; Stat_Comp_Label = 'Cong_vs_Aud_Clear'; 


% Stat_Comp = [1,3;2,4;5,7;6,8]; Stat_Labels = {'A_Noise','A_Clear','C_Noise','C_Clear'}; Stat_Comp_Label = 'AV_Noise'; % 
% Stat_Comp = [1,3;5,7]; Stat_Labels = {'A_Noise','C_Noise'}; Stat_Comp_Label = 'AvsV_HighNoise'; % 
% Stat_Comp = [5:8;1:4]; Stat_Labels = {'Cong','Aud'}; Stat_Comp_Label = 'Cong_vs_Aud'; % 
% Stat_Comp = [5:8;1:4;9:12]; Stat_Labels = {'Cong','Aud','Incong'}; Stat_Comp_Label = 'Cong_vs_Aud_vs_Incong'; % 
% Stat_Comp = [13 15;14 16]; Stat_Labels = {'VisNoise','VisClear'}; Stat_Comp_Label = 'VisN_vs_VisC'; % 
%Stat_Comp = [13 15;14 16;1 3;2 4]; Stat_Labels = {'VisNoise','VisClear','AudNoise','AudClear'}; Stat_Comp_Label = 'VisNoiseClear_vs_AudNoiseClear'; % 


switch Stat_Comp_Label
    case'Cong_vs_Aud'
        % Use this to delete conds
        % tbl([find(contains(tbl.Trial_Cond,'vis'));find(contains(tbl.Trial_Cond,'incong'))],:) = [];

        % Relabel conds to numbers so that order for LME tests is
        % consistent
        % * Since there was a validity bias, there could be
        % pre-stim differences

        cond_0_name = 'Cong_vs_Null';
        cond_1_name = 'Aud_vs_Null';
        cond_all_name = 'Cong_vs_Aud';
        legend_array = {'Cong','Aud'};

    case'Cong_vs_Incong'
        cond_0_name = 'Cong_vs_Null';
        cond_1_name = 'Incong_vs_Null';
        cond_all_name = 'Cong_vs_Incong';
        legend_array = {'Cong','Incong'};

    case'Aud_BF_vs_Aud_DG'
        % Use this to delete conds
        % tbl([find(contains(tbl.Trial_Cond,'vis'));find(contains(tbl.Trial_Cond,'incong'))],:) = [];

        % Relabel conds to numbers so that order for LME tests is
        % consistent
        % * Since there was a validity bias, there could be
        % pre-stim differences

        cond_0_name = 'Aud_BF_vs_Null';
        cond_1_name = '_Aud_DG_vs_Null';
        cond_all_name = 'Aud_BF_vs_Aud_DG';
        
    case'Incong_McG_vs_NonMcG'
        % Use this to delete conds
        % tbl([find(contains(tbl.Trial_Cond,'vis'));find(contains(tbl.Trial_Cond,'incong'))],:) = [];

        % Relabel conds to numbers so that order for LME tests is
        % consistent
        % * Since there was a validity bias, there could be
        % pre-stim differences

        cond_0_name = 'Incong_McG_vs_Null';
        cond_1_name = 'NonMcG_vs_Null';
        cond_all_name = 'Incong_McG_vs_NonMcG';
        
     case'AV_McGurk_Interaction'
        % Use this to delete conds
        % tbl([find(contains(tbl.Trial_Cond,'vis'));find(contains(tbl.Trial_Cond,'incong'))],:) = [];

        % Relabel conds to numbers so that order for LME tests is
        % consistent
        % * Since there was a validity bias, there could be
        % pre-stim differences

        cond_0_name = 'Aud_vs_Null';
        cond_1_name = 'AV_vs_Null';
        cond_all_name = 'Aud_vs_AV';
        
      case'AV_Noise_Interaction'
        % Use this to delete conds
        % tbl([find(contains(tbl.Trial_Cond,'vis'));find(contains(tbl.Trial_Cond,'incong'))],:) = [];

        % Relabel conds to numbers so that order for LME tests is
        % consistent
        % * Since there was a validity bias, there could be
        % pre-stim differences

        cond_0_name = 'Aud_vs_Null';
        cond_1_name = 'AV_vs_Null';
        cond_all_name = 'McGurk_vs_AV';
end


% condRunningSeq_idx = 2;
% 
% if  condRunningSeq_idx == 1
%     Stat_Comp = [1,2;3,4;9,10; 11,12]; Stat_Labels = {'AudMcG','AudNonMcG','AVMcG','AVNonMcG'}; 
%     Stat_Comp_Label = 'AV_McGurk_Interaction';  % A vs AV (cong + incong)
% elseif condRunningSeq_idx == 2
%     Stat_Comp = [1,3;2,4;5,7;6,8]; Stat_Labels = {'AudNoise','AudClear','CongNoise','CongClear'}; 
%     Stat_Comp_Label = 'AV_Noise_Interaction'; % A vs AV (cong + incong)
% end

% Frequency Range
% 1 = HG, 2 = Theta, 3 = Beta
Freq_Options = {'HG','Theta','Beta'};

% Only include one value if testing main effect of condition (unless want to merge across conds)
Freq_Test = 2; % Include multiple numbers in brackets to either average or compare
% Freq_Test = [1:3]; % Include multiple numbers in brackets to either average or compare



% LME titles
LME_Options_Title = {'MainEffect','CondXFreq','CondXFreqXROI','CondXFreqXROIXTime','CondXTime','CondXFreqXTime'};

% Create title
Analysis_Title = LME_Options_Title{LME_Comparison};
for freqIdx=1:length(Freq_Test)
    Analysis_Title = [Analysis_Title '_' Freq_Options{Freq_Test(freqIdx)}];
end

% Check if more than 1 frequency included
if ismember(LME_Comparison,[2:3 6]) && length(Freq_Test)==1
    error('Must include more than 1 frequency to test interaction')
end
if ismember(LME_Comparison,[3 4]) 
    if Analysis_Type==1
        error('LME options 3 and 4 only work on Time-series analysis')
    end
end
%%
%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Orthogonal Analyses Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: Only runs for 2+ conds

orthogPthresh = .01; % set to -1 to skip
orthogTimeRng = [0 .5];
orthogPolarity = [1 0 0]; % Positive polarity for HG, both for Theta, Beta
%
% 
% Setup Orthogonal time range to look for significance
resampledRate = 100;
timeStart = -2;
timeEnd = 2;
orthogTime = timeStart:1/resampledRate:timeEnd; % Can recreate time-array since all data resampled
% nearest point search
% Idx of samples that are used for orthogonal analysis
orthogSampleIdx = dsearchn(orthogTime',orthogTimeRng(1)):dsearchn(orthogTime',orthogTimeRng(2));
%%%%%%%%%%%%%%%%%%%%% Orthogonal Analyses Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Select Surface to plot on %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: Script will generate options for all listed surfaces to skip 
% re-generating the data later if you change resolution

surfResList = {'','.2mm','.4mm','.10mm','.20mm'};
surfResListNames = {'x_1mm','x_2mm','x_4mm','x_10mm','x_20mm'};


surfResIdx=4;
surfRes = surfResList{surfResIdx};
surfResName = surfResListNames{surfResIdx};
[SUCCESS,MESSAGE,MESSAGEID] = mkdir([outputDir,'/images/',Stat_Comp_Label,'_',surfResName(3:end)]);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Behavior Files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sublist = {...
    '1162HF','ECOG_AV_NoisySpeech_v4_1162HF_7-27_15-15.mat';...
    '1164UM','ECOG_AV_NoisySpeech_v1_1164UM_7-23_15-37.mat';...
    '1165HF','ECOG_AV_NoisySpeech_v4_1165HF_8-1_11-52.mat';...
    '1167HF','ECOG_AV_NoisySpeech_v8_1167HF_8-10_14-51.mat';...
    '1168UM','ECOG_AV_NoisySpeech_v12_1168UM_9-7_14-47.mat';...
    '1171UC','ECOG_AV_NoisySpeech_v12_1171UC_9-26_15-49.mat';...
    '1181UM','ECOG_AV_NoisySpeech_v13_1181UM_12-14_14-14.mat';...
    '1182HF','ECOG_AV_NoisySpeech_v13_1182_HF_1-9_15-41.mat';...
    '1183HF','ECOG_AV_NoisySpeech_v13_1183HF_1-25_15-36.mat';...
    '1187HF','ECOG_AV_NoisySpeech_v13_1187HF_3-20_11-36.mat';...
    '1188UM','ECOG_AV_NoisySpeech_v13_1188UM_4-5_14-58.mat';...
    '1193UC','ECOG_AV_NoisySpeech_v15_1193UC_6-5_16-43.mat';...
    '1208UM','ECOG_AV_NoisySpeech_v15_1208UM_11-18_12-26.mat';...
    '1210UM','ECOG_AV_NoisySpeech_v15_1210UM_12-18_13-16.mat';...
    '1214HF','ECOG_AV_NoisySpeech_v15_1214HF_2-5_15-43.mat';...
    '1219UM','ECOG_AV_NoisySpeech_v15_1219UM_3-14_12-19.mat';...
    '1223UM','ECOG_AV_NoisySpeech_v15_1223UM_4-16_13-9.mat';...
    '1226UM','ECOG_AV_NoisySpeech_v15_1226UM_6-11_13-53.mat';...
    '1230UM','ECOG_AV_NoisySpeech_v15_1230UM_8-5_14-34.mat';...
    '1243UM','ECOG_AV_NoisySpeech_v15_1243UM_11-26_14-26.mat';...
    '1245UM','ECOG_AV_NoisySpeech_v15_1245UM_1-13_14-2.mat';...
    };

subjIDs = sublist(:,1)';
behavFiles = join([repmat({behavDir},length(sublist(:,2)'),1),sublist(:,2)],'')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Electrode Regions to Include %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

regions_to_include = [...
                        "superiortemporal",...
                        "middletemporal",...
                        "transversetemporal",...
                        "bankssts",...
                        "supramarginal",...
                        "lateraloccipital",...
                        "inferiorparietal",...
                        "inferiortemporal"...
                        ];           %% Regions to analyze
EL_radius = 10;
num_subs_req = 2;   %% Minimum number of subs at a vertex to model effect
reference = 2;  %% Bipolar

% Update freq options to set analyses
Freq_Options = Freq_Options(Freq_Test);
orthogPolarity = orthogPolarity(Freq_Test);
ROI_names = {'Anterior','Middle','Posterior'};
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%     Run plotter           %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: 
% Only need to run once to generate .mat files
% This gets us the bread and butter to run the vertex based or timeseries
% based analysis below

if PreProcess_Data==1
    for Freq_X=1:length(Freq_Options)
        Freq_Test = Freq_Options{Freq_X};
        LME_AVSpeech_ProcessData(subjIDs,behavFiles,SubjPath,electRegPath,bipDataDir,outputDir,Freq_Test,regions_to_include,EL_radius,surfResList,surfResListNames,FS_Subject_Dir)
    end
end


if ReProcess_Vertices==1
    for Freq_X=1:length(Freq_Options)
        LME_AVSpeech_ProcessVertices(subjIDs,SubjPath,electRegPath,outputDir,Freq_Test,EL_radius,surfResList,surfResListNames,reference)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%     Run plotter           %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Find vertex with the most participants and electrodes
% VERTEX-WISE ANALYSIS
if strcmp(Analysis_Type,'Voxel-Wise') 
    
    clear DATA_Subs DATA_Subs_Behavior
    for su=1:length(subjIDs)
        for Freq_X=1:length(Freq_Options)
            % Load sub info
            subjid = subjIDs{su};
            load([outputDir,'/',Freq_Options{Freq_X},'_',subjIDs{su},'_DATA'],'SUBX','DATA')
            DATA_Subs.(genvarname([Freq_Options{Freq_X},'_',subjIDs{su}])).DATA = DATA;
        end
        
        % Load behavioral data
        load([outputDir,'/',subjIDs{su},'_BEHAV'],'DATA_Behav')
        DATA_Subs_Behavior.(genvarname(['x',subjIDs{su}])) = DATA_Behav;
    end
    
    % Load vertices information
    % VERTEX.(genvarname(surfResListNames{surfRes}))
    %iterate through 3 freqs
    VERTEX_ALL = {};
    for su = 1:length(subjIDs)
        subjid = subjIDs{su};
        
        VERTEX_ALL_temp = {};
        for Freq_X=1:length(Freq_Options)
            
            % Load sub info
            load([outputDir,'/',Freq_Options{Freq_X},'_',subjIDs{su},'_VERTEX'],'SUBX','VERTEX')
            
            VERTEX_ALL_temp = [VERTEX_ALL_temp;VERTEX.(genvarname(surfResName)),num2cell(zeros(size(VERTEX.(genvarname(surfResName)),1),1)+Freq_X)];
        end
                
        % Populate all vertices
        VERTEX_ALL = [VERTEX_ALL;VERTEX_ALL_temp];
    end
   
    
    
    % Check orthogonality of data when 2+ conds
    % Remove electrodes that do not have a significant response
    % Will not run if only 1 condition because then circular analysis
    %     
    % Create array of subs and electrodes that should be excluded from
    % analyses. Apply this to the VERTEX_ALL structure

    
    %     n_conds = size(Stat_Comp,1);
%     if n_conds>=2 && orthogPthresh>0
%         combined_conds = Stat_Comp(:);
%         
%         % Create combined VERTEX_ALL variable to identify indices
%         temp=cellfun(@num2str,VERTEX_ALL(:,4),'un',0);
%         VERTEX_ALL_combined = join([VERTEX_ALL(:,1),VERTEX_ALL(:,2),temp],'_');
%         
%         %iterate through freqs being analyzed
%         VERTEX_ALL_rm_idx = [];
%         for Freq_X=1:length(Freq_Options)
%             curr_Freq = Freq_Options{Freq_X};
%             
%             % Iterate through subjects
%             for su = 1:length(subjIDs)
%                 subjid = subjIDs{su};
%                 
%                 % Pull data for combined conds.
%                 currSubData = DATA_Subs.([curr_Freq,'_',subjid]).DATA;
%                 Combined_ECoG = currSubData.ECOG(:,orthogSampleIdx,find(ismember(currSubData.Stim_Array, combined_conds)));
%                 
%                 % RM non-analyzed electrodes
%                 Combined_ECoG = Combined_ECoG(currSubData.Analyzed_Channels,:,:);
%                 temp_Sub_Data_chans = currSubData.Channels(currSubData.Analyzed_Channels);
%                 
%                 % check each electrode
%                 for elec_x = 1:size(Combined_ECoG,1)
%                     % pull data, rm nans
%                     Combined_ECoG_el = squeeze(Combined_ECoG(elec_x,:,:)); % time x trial
%                     Combined_ECoG_el(:,find(isnan(Combined_ECoG_el(1,:))))=[];
%                     
%                     % run stats
%                     [~,P,~,STATS] = ttest(Combined_ECoG_el,0,'dim',2);
%                     
%                     % if polarity requirement is set, alter p values
%                     if orthogPolarity(Freq_X)==1 % e.g., HGp
%                         temp = find(STATS.tstat<0);
%                         P(temp) = 1;
%                     end
%                     if orthogPolarity(Freq_X)==-1 % uncommon
%                         temp = find(STATS.tstat>0);
%                         P(temp) = 1;
%                     end
%                     
%                     % Find if at least 1 fdr-corrected time-point survives
%                     [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(P);
%                     adj_p = min(adj_p);
%                     
%                     if adj_p>=orthogPthresh
%                         % Find index to remove electrode
%                         temp_elec = temp_Sub_Data_chans{elec_x};
%                         temp_elec = temp_elec(1:end-2); % rm _b
%                         temp_name = [subjid,'_',temp_elec,'_',num2str(Freq_X)];
%                         VERTEX_ALL_rm_idx = [VERTEX_ALL_rm_idx,find(strcmp(temp_name,VERTEX_ALL_combined))'];
%                     end
%                 end
%                 
%             end
%         end
%         
%         % remove electrodes that are not sig on orthogonal comparison
%         VERTEX_ALL(VERTEX_ALL_rm_idx,:)=[];
%         
%     end
%     
    % This code counts the number of independent subjects who contribute data
    % to each vertex (after removing data if applied above)
    VERTEX_ALL_Numbers = cell2mat(VERTEX_ALL(:,3)); % convert cell array of vertices to numbers
    Unique_Vertices = unique(VERTEX_ALL_Numbers); % Get unique vertices
    Num_Subs_Per_Vertex = zeros(1,length(Unique_Vertices)); % Dummy array
    for i=1:length(Unique_Vertices) % Iterate through each vertex to count subs
        temp = find(VERTEX_ALL_Numbers==Unique_Vertices(i)); % Indentify the indices for each vertex
        Num_Subs_Per_Vertex(i) = size(unique(VERTEX_ALL(temp,1)),1); % Count unique subs at those indices
    end
    
    % % % % % % Find Vertex present in most subs
    % % % % % temp_vertices_index = find(temp_vertices_subs==max(temp_vertices_subs),1);
    % % % % % temp_vertices_hit = temp_vertices(temp_vertices_index);
    % % % % %
    % % % % % Analysis_electrode_Idx = find(temp_vertices_numbers==temp_vertices_hit); % 79608
    % % % % % Analysis_electrode_Data = VERTEX_ALL(Analysis_electrode_Idx,:);
    % % % % % Analysis_electrode_Data_Subs = unique(Analysis_electrode_Data(:,1));
    % % % % %
    % % % % % % Iterate through subjects, load in data into structure
    % % % % % clear DATA_Subs
    % % % % % for i=1:length(Analysis_electrode_Data_Subs)
    % % % % %     load([outputDir,'/',Freq_Test,'_',Analysis_electrode_Data_Subs{i},'_DATA'],'SUBX','DATA')
    % % % % %     DATA_Subs.(genvarname(['X_',Analysis_electrode_Data_Subs{i}])).DATA = DATA;
    % % % % % end
    
    
            

    
    
    %% Iterate through vertices with more than 2 subjects, w interaction
    
    % Identify vertices that have more than the required number of subs
    % In this case, we only include subs with 2 or more subs
    Included_Vertices = Unique_Vertices(find(Num_Subs_Per_Vertex>=num_subs_req));
    Included_Vertices_Subs = Num_Subs_Per_Vertex(find(Num_Subs_Per_Vertex>=num_subs_req));
    
    


    % Input times to be analyzed (e.g., from -.1 to 0 will avg that rng)
    time_rngs=[];
    time_rngs(:,1) = -1:.1:.4;
    time_rngs(:,2) = -.9:.1:.5;
    
    % Testing
    % time_rngs = [-.2 -.1;-.1 0;0 .1;.1 .2];

    % Setup Stat arrays
    STATS.P = nan(size(time_rngs,1),length(Included_Vertices));
    STATS.Est = nan(size(time_rngs,1),length(Included_Vertices));
    
    % %     %%%%% NEED TO UPDATE textprogressbar. BREAKS IF RUNNING AFTER ERROR
    %             BECAUSE THEY'RE USING PERSISTANT VARIABLES
    textprogressbar('Creating table for each Vertex: ');
    
    for roiIdx = 1:length(Included_Vertices)
        
        %waitbar(roiIdx/length(Included_Vertices),f)
        textprogressbar(round(100*roiIdx/length(Included_Vertices),0));
        
        Analysis_electrode_Idx = find(VERTEX_ALL_Numbers==Included_Vertices(roiIdx)); % 79608
        Analysis_electrode_Data = VERTEX_ALL(Analysis_electrode_Idx,:);
        
        % Use this to build up subject density plot
        Analysis_electrode_Data_Subs = unique(Analysis_electrode_Data(:,1));
        
        % Iterate through subs, create matrix of values, averaged across
        % electrodes
        grouped_sub_p = nan(size(Analysis_electrode_Data_Subs,1),size(time_rngs,1));
        for sub_x=1:size(Analysis_electrode_Data_Subs,1)
            currSub = Analysis_electrode_Data_Subs{sub_x,1};
            currSubFreq = Freq_Options{Analysis_electrode_Data{sub_x,4}};
            currSubData = DATA_Subs.([currSubFreq,'_',currSub]).DATA;
            
            % Find matching electrodes
            currElecName = Analysis_electrode_Data(strcmp(Analysis_electrode_Data(:,1),currSub),2);
            
            % Find electrode indices
            temp_Elec_Idx = nan(1,length(currElecName));
            for j=1:length(currElecName)
                temp_Elec_Idx(j) = find(strcmp(currSubData.Channels,[currElecName{j},'_b'])); % Add back in bipolar
            end
            
            % Pull data and average across selected electrodes
            temp_Sub_ECOG = squeeze(nanmean(currSubData.ECOG(temp_Elec_Idx,:,:),1))'; % Trials x time                
 
            % Limit data to analyzed trials
            % Moved this step down from earlier in the code since indexing
            % a 3D matrix takes much longer
            currSubData.Analyzed_Trials = find(~isnan(temp_Sub_ECOG(:,1)));
            temp_Sub_ECOG = temp_Sub_ECOG(currSubData.Analyzed_Trials,:);
            currSubData.Stim_Array = currSubData.Stim_Array(currSubData.Analyzed_Trials);
            
            % Average across time windows
  
            temp_Sub_ECOG_Avg = nan(size(temp_Sub_ECOG,1),size(time_rngs,1)); % trials x time
            for time_x=1:size(time_rngs,1)
                time_rngs_val = dsearchn(orthogTime',time_rngs(time_x,1)):dsearchn(orthogTime',time_rngs(time_x,2));
                temp_Sub_ECOG_Avg(:,time_x) = nanmean(temp_Sub_ECOG(:,time_rngs_val),2); % trial by time: ECOG
            end
            
            % Calculate statistics for this sub/vertex
            % One-sample for one condition
            if n_conds==1
                temp_cond_idx = Stat_Comp(1,:);
                temp_cond_idx = find(ismember(currSubData.Stim_Array,temp_cond_idx));
                [~,grouped_sub_p(sub_x,:)] = ttest(temp_Sub_ECOG_Avg(temp_cond_idx,:));
                grouped_sub_p(sub_x,:) = grouped_sub_p(sub_x,:).*sign(nanmean(temp_Sub_ECOG_Avg(temp_cond_idx,:),1));
            end
            
            % If 2 conditions, run an unpaired t-test
            if n_conds==2
                temp_cond_idx1 = Stat_Comp(1,:);
                temp_cond_idx1 = find(ismember(currSubData.Stim_Array,temp_cond_idx1));
                temp_cond_idx2 = Stat_Comp(2,:);
                temp_cond_idx2 = find(ismember(currSubData.Stim_Array,temp_cond_idx2));
                [~,grouped_sub_p(sub_x,:)] = ttest2(temp_Sub_ECOG_Avg(temp_cond_idx1,:),temp_Sub_ECOG_Avg(temp_cond_idx2,:));
                grouped_sub_p(sub_x,:) = grouped_sub_p(sub_x,:).*sign(nanmean(temp_Sub_ECOG_Avg(temp_cond_idx1,:),1)-nanmean(temp_Sub_ECOG_Avg(temp_cond_idx2,:),1));
            end
            
            % If more than 2 conditions, run univariate anova
            if n_conds>2
                % trial x time
                tempECoG = [];
                tempgroup = [];
                for Condx=1:n_conds % 1:4
                    % go in the stimulus array and find the corresponding
                    % trials for condition (e.g.) 1 and 2
                    tempgoup_idx = find(ismember(currSubData.Stim_Array,Stat_Comp(Condx,:)));
                    tempgroup = [tempgroup zeros(1,length(tempgoup_idx))+Condx];
                    tempECoG = [tempECoG;temp_Sub_ECOG_Avg(tempgoup_idx,:)];
                end
                
                % Run anova on each time-point
                for time_x=1:size(time_rngs,1)
                    tempdata = tempECoG(:,time_x)';
                    tempP=anova1(tempdata,tempgroup,'off');
                    grouped_sub_p(sub_x,time_x)= tempP;
                end
                
%                 % REPEATED from line 500?
%                 for Condx=1:n_conds % 1:4
%                     % go in the stimulus array and find the corresponding
%                     % trials for condition (e.g.) 1
%                     tempgoup_idx = find(ismember(currSubData.Stim_Array,Stat_Comp(Condx,:)));
%                     tempgroup = [tempgroup zeros(1,length(tempgoup_idx))+Condx]; % creating the factors
%                     tempECoG = [tempECoG;temp_Sub_ECOG_Avg(tempgoup_idx,:)]; % trial x time 
%                 end 
                subjectData.Conditions = tempgroup;
                subjectData.ECoG = tempECoG;
                
                % tell the table what is non-McG(0) and McG (1)
                
                if strcmp(Stat_Comp_Label,'AV_McGurk_Interaction')
                    subjectData.McGurk= nan(size(tempgroup));
                    subjectData.Audiovisual= nan(size(tempgroup));
                    

                    for trial_idx=1:length(tempgroup)
                        switch tempgroup(trial_idx)
                            case 1
                                subjectData.McGurk(trial_idx) = 1;
                                subjectData.Audiovisual(trial_idx) =0;
                            case 2
                                subjectData.McGurk(trial_idx) = 0;
                                subjectData.Audiovisual(trial_idx) =0;
                            case 3
                                subjectData.McGurk(trial_idx) = 1;
                                subjectData.Audiovisual(trial_idx) =1;
                            case 4
                                subjectData.McGurk(trial_idx) = 0;
                                subjectData.Audiovisual(trial_idx) =1;
                        end
                    end
                    
                            % here we create a table 
                            % this is what it should look like at the end of the table building
                                % SubjID McG AV t1 t2 t3 t4.....t15
                                %        1    1
                                %        1    2
                                %        2    1
                                %        2    2

% % % %                 %THIS SHOULD BE IN THE IF BLOCK
% % % %                 % NOT USING TABLE?
                       subInteractionTable = table(tempgroup',subjectData.McGurk',subjectData.Audiovisual',subjectData.ECoG,...
                           'VariableNames',{'Group','McGurk','Audiovisual','ECoG'});

                       % run a 2x2 ANOVA for each time point

                       % % % % %           % USE THIS INSTEAD -- gives same result but cleaner
                       % % % % %           anovan(subInteractionTable.ECoG(:,10),{subInteractionTable.Audiovisual subInteractionTable.Noise},'model','interaction','varnames',{'AV_Type','Noise_Type'})
                       
                       % stouffer might only work well for tests since
                       % anovas and interactions don't have negative stats
                       % Maybe treating as one-sample t-test works?

                       % in the end, we should get a p value from the 2x2 anova for
                       % each timepoint for each subject at each vertex
                       interaction_pval_arr=nan(1,size(tempECoG,2));

                       % find cond trial index (again)
                        % McG 0, AV 0
                       nonMcG_nonAV_trial_idx = find(subjectData.Conditions==2);
                        % McG 0, AV 1
                       nonMcG_AV_trial_idx = find(subjectData.Conditions==4);
                        % McG 1, AV 0
                       McG_nonAV_trial_idx = find(subjectData.Conditions==1);
                        % McG 1, AV 1
                       McG_AV_trial_idx = find(subjectData.Conditions==3);

                       for time_idx = 1:size(tempECoG,2)   
                           % fill in the table
                           time_ECoG_Value=tempECoG(:,time_idx);

                           nonMcG_ECoG=[time_ECoG_Value(nonMcG_nonAV_trial_idx),time_ECoG_Value(nonMcG_AV_trial_idx)];
                           McG_ECoG=[time_ECoG_Value(McG_nonAV_trial_idx),time_ECoG_Value(McG_AV_trial_idx)];

                           time_ECoG_Table = [nonMcG_ECoG;McG_ECoG];

                           [p,interaction_tbl] = anova2(time_ECoG_Table,length(subjectData.Conditions)/4,'off'); % turning off the anova table produced
                           % go into table to find where the interaction value is
                           time_interaction_val = cell2mat(interaction_tbl(4,6));

                           interaction_pval_arr(time_idx)= time_interaction_val;
                       end
               
                elseif strcmp(Stat_Comp_Label,'AV_Noise_Interaction')
                      subjectData.Noise= nan(size(tempgroup));
                      subjectData.Audiovisual= nan(size(tempgroup));
                    

                    for trial_idx=1:length(tempgroup)
                        switch tempgroup(trial_idx)
                            case 1
                                subjectData.Noise(trial_idx) = 1;
                                subjectData.Audiovisual(trial_idx) =0;
                            case 2
                                subjectData.Noise(trial_idx) = 0;
                                subjectData.Audiovisual(trial_idx) =0;
                            case 3
                                subjectData.Noise(trial_idx) = 1;
                                subjectData.Audiovisual(trial_idx) =1;
                            case 4
                                subjectData.Noise(trial_idx) = 0;
                                subjectData.Audiovisual(trial_idx) =1;
                        end
                    end

                    % % % %                 %THIS SHOULD BE IN THE IF BLOCK
                    % % % %                 % NOT USING TABLE?
                    subInteractionTable = table(tempgroup',subjectData.Noise',subjectData.Audiovisual',subjectData.ECoG,...
                        'VariableNames',{'Group','Noise','Audiovisual','ECoG'});
                    
                    % % % % %           % USE THIS INSTEAD -- gives same result but cleaner
                    % % % % %           anovan(subInteractionTable.ECoG(:,10),{subInteractionTable.Audiovisual subInteractionTable.Noise},'model','interaction','varnames',{'AV_Type','Noise_Type'})

                       % run a 2x2 ANOVA for each time point

                       % in the end, we should get a p value from the 2x2 anova for
                       % each timepoint for each subject at each vertex
                       interaction_pval_arr=nan(1,size(tempECoG,2));

                       % find cond trial index (again)
                        % Noise 0, AV 0
                       clear_nonAV_trial_idx = find(subjectData.Conditions==2);
                        % Noise 0, AV 1
                       clear_AV_trial_idx = find(subjectData.Conditions==4);
                        %  Noise 1, AV 0
                       noise_nonAV_trial_idx = find(subjectData.Conditions==1);
                        %  Noise 1, AV 1
                       noise_AV_trial_idx = find(subjectData.Conditions==3);

                       for time_idx = 1:size(tempECoG,2)   
                           % fill in the table
                           time_ECoG_Value=tempECoG(:,time_idx);

                           clear_ECoG=[time_ECoG_Value(clear_nonAV_trial_idx),time_ECoG_Value(clear_AV_trial_idx)];
                           noise_ECoG=[time_ECoG_Value(noise_nonAV_trial_idx),time_ECoG_Value(noise_AV_trial_idx)];

                           time_ECoG_Table = [clear_ECoG;noise_ECoG];

                           [p,interaction_tbl] = anova2(time_ECoG_Table,length(subjectData.Conditions)/4,'off'); % turning off the anova table produced
                           % go into table to find where the interaction value is
                           time_interaction_val = cell2mat(interaction_tbl(4,6));

                           interaction_pval_arr(time_idx)= time_interaction_val;
                       end
                end
                grouped_sub_p(sub_x,:)=interaction_pval_arr;
            end
        end
  
        
        % Update P
        for time_x=1:size(time_rngs,1)
            STATS.P(time_x,roiIdx) = stouffer(grouped_sub_p(:,time_x), 2);
                
            % add sign info
            % STATS.Est(time_x,roiIdx) = sign(mean(grouped_sub_p(:,time_x)));            
            STATS.Est(time_x,roiIdx) = sign(STATS.P(time_x,roiIdx));            
        end

    end
    %close(f)
        
    % Fix any pcomb values of 0
    % Matlab will output pcomb < 1e-323 as 0
    temp = find(STATS.P==0);
    if ~isempty(temp)
        STATS.P(temp) = 1.0e-323*sign(STATS.Est);
    end
    
    % Plot vertices
    % put time_x in 3d plot; fdr for time & space
    plot_elects = 0;
    LME_AVSpeech_3DBrainPlot(plot_elects,STATS,Included_Vertices,Included_Vertices_Subs,SubjPath,surfRes,surfResName,regions_to_include,time_rngs,[outputDir,'/images/',Stat_Comp_Label,'_',surfResName(3:end),'/',Analysis_Title,'_',],FS_Subject_Dir)
    
    % Save stats for easy replotting
    save([outputDir,'/images/',Stat_Comp_Label,'_',Analysis_Title,'_',surfResName(3:end),'.mat'],'STATS','time_rngs','Included_Vertices','SubjPath','surfRes','outputDir','Stat_Comp_Label','Freq_Test','surfResName','regions_to_include','LME_Comparison','Analysis_Title','surfResIdx')
    
end











%%
% Process 3 STG ROI Time-Series

if strcmp(Analysis_Type,'Time-Series') 
    
    % Already have the vertices next to each electrode
    % Match these vertices in the 1mm data against the 3 ROI labels' vertices
    
    % Use 1mm
    surfResIdx=1;
    surfRes = surfResList{surfResIdx};
    surfResName = surfResListNames{surfResIdx};
    
    % Load stgdivisions.csv
    % Add 1 since label file starts at 0, whereas matlab idx starts at 1
    stgdivisions = readmatrix('stgdivisions.csv');
    stgdivisions(:,1) = stgdivisions(:,1)+1;
    ROI_Num = length(unique(stgdivisions(:,2)));
    

    
    %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% load the vertices of the subject IDs %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Note:  vertex_all is 4 columns
    % subjID, electrode Name, the vertice number, and Freq thats being
    % analyzed


    VERTEX_ALL = {};
    for subjIdx = 1:length(subjIDs)
        subjid = subjIDs{subjIdx};
        for Freq_X=1:length(Freq_Options)       
            % Load sub info
            subjVertex = load([outputDir,'/',Freq_Options{Freq_X},'_',subjIDs{subjIdx},'_VERTEX'],'SUBX','VERTEX');
            subjVertexData = subjVertex.VERTEX.(genvarname(surfResName));
            freqCol = num2cell(zeros(size(subjVertex.VERTEX.(genvarname(surfResName)),1),1)+Freq_X);
            % Populate all vertices
            % Since only 3 ROIs, don't need to exlude 'vertices' w/o enough data
            VERTEX_ALL = [VERTEX_ALL;...
                subjVertexData,...
                freqCol];
            clear subjVertex subjVertexData
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%    Load all subjects' data %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Iterate through subjects, load in data into structure
    clear DATA_Subs
    %textprogressbar('Loading subject data: ');
    for su=1:length(subjIDs)
        %textprogressbar(round(100*su/length(subjIDs),0));
        for Freq_X=1:length(Freq_Options)
            % Load sub info
            subjid = subjIDs{su};
            load([outputDir,'/',Freq_Options{Freq_X},'_',subjIDs{su},'_DATA'],'SUBX','DATA')
            DATA_Subs.(genvarname([Freq_Options{Freq_X},'_',subjIDs{su}])).DATA = DATA;
        end
    end
    %textprogressbar('Done');


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Check orthogonality of data when 2+ conds %%%%%%%%%%%%%%
    %%%%% Remove electrodes that do not have a significant response %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Note: Will not run if only 1 condition because then circular analysis
    

    % Currently setup to require that electrodes be signifcant in each
    % frequency for inclusion. LME would take care of this missing data.
    % Do separately for different freqs
 

    % Create array of subs and electrodes that should be excluded from
    % analyses. Apply this to the VERTEX_ALL structure

    % for example, STAT_COMP here is 
    % 5,6,7,8
    % 1,2,3,4
    % there are two conditions(Cong and Aud), and
    % each condition has 4 corresponding exp index

    n_conds = size(Stat_Comp,1);
%     if n_conds>=2 && orthogPthresh>0
%         combined_conds = Stat_Comp(:);
%         
%         % Create combined VERTEX_ALL variable to identify indices
%         temp=cellfun(@num2str,VERTEX_ALL(:,4),'un',0);
%         VERTEX_ALL_combined = join([VERTEX_ALL(:,1),VERTEX_ALL(:,2),temp],'_');
%         
%         %iterate through freqs being analyzed
%         VERTEX_ALL_rm_idx = [];
%         for Freq_X=1:length(Freq_Options)
%             curr_Freq = Freq_Options{Freq_X};
%             
%             % Iterate through subjects
%             for su = 1:length(subjIDs)
%                 subjid = subjIDs{su};
%                 
%                 % Pull data for combined conds.
%                 
%                 % channel x time x trials
%                 currSubData = DATA_Subs.([curr_Freq,'_',subjid]).DATA;
%                 
%                 %  only keep trials of what we are interested in
%                 Combined_ECoG = currSubData.ECOG(:,orthogSampleIdx,find(ismember(currSubData.Stim_Array, combined_conds)));
%                 
%                 % Remove non-analyzed electrodes
%                 Combined_ECoG = Combined_ECoG(currSubData.Analyzed_Channels,:,:);
%                 % get names of the analyzed channels
%                 temp_Sub_Data_chans = currSubData.Channels(currSubData.Analyzed_Channels);
%                 
%                 % check each electrode
%                 for elec_x = 1:size(Combined_ECoG,1)
%                     % pull data from each elec, squeeze to time x trial
%                     Combined_ECoG_el = squeeze(Combined_ECoG(elec_x,:,:)); 
%                     % remove nans
%                     Combined_ECoG_el(:,find(isnan(Combined_ECoG_el(1,:))))=[];
%                     
%                     % run stats
%                     % 1 tail or 2 tail should depend on polarity req
%                     % CHANGE
%                     [~,P,~,STATS] = ttest(Combined_ECoG_el,0,'dim',2);
%                     
%                     % if polarity requirement is set, alter p values
%                     if orthogPolarity(Freq_X)==1 % e.g., HGp
%                         temp = find(STATS.tstat<0);
%                         P(temp) = 1;
%                     end
%                     
%                     if orthogPolarity(Freq_X)==-1 % uncommon
%                         temp = find(STATS.tstat>0);
%                         P(temp) = 1;
%                     end
%                     
%                     % Find if at least 1 fdr-corrected time-point survives
%                     [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(P);
%                     adj_p = min(adj_p);
%                     
%                     if adj_p>=orthogPthresh
%                         % Find index to remove electrode
%                         temp_elec = temp_Sub_Data_chans{elec_x};
%                         temp_elec = temp_elec(1:end-2); % rm _b
%                         temp_name = [subjid,'_',temp_elec,'_',num2str(Freq_X)];
%                         VERTEX_ALL_rm_idx = [VERTEX_ALL_rm_idx,find(strcmp(temp_name,VERTEX_ALL_combined))'];
%                     end
%                 end
%                 
%             end
%         end
%         
%         % remove electrodes that are not sig on orthogonal comparison
%         VERTEX_ALL(VERTEX_ALL_rm_idx,:)=[];
%     end
     % find all unique electrodes 
     % Assuming VERTEX_ALL is a cell array with strings in the first two columns

    % Concatenate the first two string columns into one
    combinedStrings = strcat(VERTEX_ALL(:,1), '_', VERTEX_ALL(:,2));
    
    % Find unique rows based on the combined string, preserving original order
    [uniqueCombinedStrings, ia] = unique(combinedStrings, 'stable');
    
    % Extract the unique rows from the original cell array
    uniqueElectrodes = VERTEX_ALL(ia,1:2);
    includedElecDir = '/Users/zhewei/Dropbox (University of Michigan)/_Projects/Brang_LME_AVSpeech/Processed/ThetaTest/';
    % now we have all the electrodes before the SigList update
    for subjIdx = 1:length(subjIDs)
        subjid = [subjIDs{subjIdx}(5:6) subjIDs{subjIdx}(1:4)];
        % significant electrode Name is stored at
        subjIncludedElec = load([includedElecDir subjid '_Updated.mat']);
        subjIncludedElec = subjIncludedElec.SigListUpdatedElecNames;
    end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Check for replicate electrodes %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%

    % now that we narrowed down to just the responsive electrodes
    % Go through each sub/elec, find if any intersection
    allSubjElecsName = join([VERTEX_ALL(:,1) VERTEX_ALL(:,2)],2);
    % get all the unique electrodes from everybody 
    allSubjElecsNameUniq = unique(allSubjElecsName);
    
    %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%
    %%%%%%%%% sort the elec into diff regions of the STG  %%%%%%%%% 
    for i=1:ROI_Num
        % create a space for elecs in each ROI
        ROI_Elecs.(['x',num2str(i)]) = {};
    end

    
    for elecIdx = 1:length(allSubjElecsNameUniq)
        % Find vertices for elec
        verticeIndices = find(strcmp(allSubjElecsNameUniq{elecIdx},allSubjElecsName));
        % are these vals the voxel index?
        verticeValues = cell2mat(VERTEX_ALL(verticeIndices,3));
        % the freq no.s of these electrodes
        electFrqs = unique(cell2mat(VERTEX_ALL(verticeIndices,4)));
        % get the sub and the electrode
        temp_idx_subelec = [VERTEX_ALL(verticeIndices(1),1) VERTEX_ALL(verticeIndices(1),2)];
        
        % iterate through each of the ROIs
        % find the proportion of vertices for each ROI w length(find(stgdivisions(:,2)==1))
        % Instead, just count vertices
        temp_hits = [0 0 0];
        for ROI_x = 1:ROI_Num
            temp_ROI_vertices = stgdivisions(stgdivisions(:,2)==ROI_x,1);
            [C,ia,ib] = intersect(temp_ROI_vertices,verticeValues);
            temp_hits(ROI_x) = length(C); 
        end
        
        % Find max vertices
        if sum(temp_hits)>0
            [~,Idx] = max(temp_hits);
            ROI_Elecs.(['x',num2str(Idx)]) = [ROI_Elecs.(['x',num2str(Idx)]);...
                temp_idx_subelec electFrqs];
        end
    end
    
   % at the end of this step, we have assigned each electrode to the three
   % subregions of the STG. This informatin is stored in the ROI_Elecs
   % varaible, x1, x2, x3 corresponds to the three parts 
    %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%
  

    % Set time
    time = timeStart:1/resampledRate:timeEnd;
    time_rngs = [-1 .5];
    time_rngs = dsearchn(time',time_rngs(1)):1:dsearchn(time',time_rngs(2));

%         time_rngs = dsearchn(time',time_rngs(1)):20:dsearchn(time',time_rngs(2));
% % % % %         time_rngs = dsearchn(time',time_rngs(1)):1:dsearchn(time',time_rngs(2));
    
% % % % %     if ismember(LME_Comparison,1:3)
% % % % %         time_rngs = dsearchn(time',time_rngs(1)):1:dsearchn(time',time_rngs(2));
% % % % %     elseif ismember(LME_Comparison,4)
% % % % %         % Need to reduce the size of the data; computer had memory issues
% % % % %         time_rngs = dsearchn(time',time_rngs(1)):5:dsearchn(time',time_rngs(2));
% % % % %     end
    


    for i=1:ROI_Num
        % create a space for elecs in each ROI
        ROI_Elecs_New.(['x',num2str(i)]) = {};
    end

    for roiUdx = 1:ROI_Num
        elecsCurrROI = ROI_Elecs.(['x',num2str(roiUdx)]);
        for subjIdx = 1:length(subjIDs)
            subjid = subjIDs{subjIdx};
            subjid2 = [subjIDs{subjIdx}(5:6) subjIDs{subjIdx}(1:4)];
            % significant electrode Name is stored at

            subjIncludedElec = load([includedElecDir  subjid2 '_Updated.mat']);
            subjIncludedElec = subjIncludedElec.SigListUpdatedElecNames;
            subjIncludedElec = cellfun(@(x) strrep(x, '_b', ''), subjIncludedElec, 'UniformOutput', false);
            for elecIdx = 1: length(subjIncludedElec)
                elecName = subjIncludedElec{elecIdx};
                currElecIdx = find(strcmp(elecsCurrROI(:,1),subjid) & strcmp(elecsCurrROI(:,2),elecName));
                if ~isempty(currElecIdx)
                    ROI_Elecs_New.(['x',num2str(roiUdx)])(end+1,:) = elecsCurrROI(currElecIdx,:);
                    clear idx
                end
            end
        end
          
    % Remove 'RTP4_b' 'RTP5_b' from 1164 due to noise
    % this is no longer needed because neither was included
%         idx1 = find(strcmp(elecsCurrROI(:,1),'1164UM') & strcmp(elecsCurrROI(:,2),'RTP4'));
%         idx2 = find(strcmp(elecsCurrROI(:,1),'1164UM') & strcmp(elecsCurrROI(:,2),'RTP5'));
%         ROI_Elecs.(['x',num2str(roiUdx)])([idx1 idx2],:)=[];
    end

    ROI_Elecs = ROI_Elecs_New;

    %% Iterate through times and ROIs, create matrix of values
    
    DATA_ALL.Text = {}; % SUB Cond Elec Trial
    % DATA_ALL.Conds = fieldnames(DATA);
    
    DATA_ALL.elecData = {}; % Fill in data values
    counter = 1;
    %textprogressbar('Creating Table: ');
    for roiIdx = 1:ROI_Num
        %textprogressbar(round(100*roiIdx/ROI_Num,0));
        Analysis_electrode_Data = ROI_Elecs.(['x',num2str(roiIdx)]); % Pull relevance Sub/Elec Names
               
        for roiElecIdx=1:size(Analysis_electrode_Data,1)
            currSub = Analysis_electrode_Data{roiElecIdx,1};
            currElecName = [Analysis_electrode_Data{roiElecIdx,2},'_b']; % Add back in bipolar
            currSubFreqs = Freq_Options(Analysis_electrode_Data{roiElecIdx,3});
            
            for Freq_X=1:length(currSubFreqs)
                % identify frequency
                currSubFreq = currSubFreqs{Freq_X};
                currSubData = DATA_Subs.([currSubFreq,'_',currSub]).DATA;
                
                for j=1:n_conds
                    
                    % currCond = DATA_ALL.Conds{j};
                    currCond = ['c',num2str(j)];
                    
                    % Find trials
                    trial_idx = find(ismember(currSubData.Stim_Array,Stat_Comp(j,:)));
                    el_idx = find(strcmp(currSubData.Channels,currElecName));
                    
                    % Data
                    % DATA_ALL.elecData(counter) = {[currSubData.(currCond).(currElecName).DATA]}; % Trials x time
                    DATA_ALL.elecData(counter) = {squeeze(currSubData.ECOG(el_idx,:,trial_idx))'}; % Trials x time
                    
                    % Text
                    % temp_trials = [DATA_Subs.([currSubFreq,'_',currSub]).DATA.(currCond).(currElecName).TRIAL_NUMBER]';
                    temp_trials = [trial_idx]';
                    DATA_ALL.Text = [DATA_ALL.Text;...
                        repmat({currSub},[length(temp_trials) 1]),...
                        repmat({currCond},[length(temp_trials) 1]),...
                        repmat({currElecName},[length(temp_trials) 1]),...
                        repmat({currSubFreq},[length(temp_trials) 1]),...
                        repmat({['x',num2str(roiIdx)]},[length(temp_trials) 1]),...
                        num2cell(temp_trials)];
                    
                    counter = counter+1;
                end
            end
        end
    end
    clear DATA_Subs DATA

    % Put into rows
    DATA_ALL.Values = [cat(1,DATA_ALL.elecData{:})]; 
        
    % Update DATA_ALL.Values to only the relevant time-points to reduce
    % memory pressure before next analyses
    DATA_ALL.Values = DATA_ALL.Values(:,time_rngs);    
    
    % Setup the final array dimensions
    if ismember(LME_Comparison,1:2) % Main effect of Cond or Cond x Freq
        DIM_time = length(time_rngs);
        DIM_ROI = ROI_Num;
    elseif ismember(LME_Comparison,3) % Cond x Freq x ROI
        DIM_time = length(time_rngs);
        DIM_ROI = 1;
    elseif ismember(LME_Comparison,4) % Cond x Freq x ROI x Time
        DIM_time = 1;
        DIM_ROI = 1;
    end
    
    %% Setup the final arrays

    tbl_array_main = cell(DIM_time,DIM_ROI); % Create array to fill. time x ROI
    clear tbl_array_conds

    for condx=1:n_conds
        tbl_array_conds.(['x',num2str(condx)]) = cell(DIM_time,DIM_ROI);
    end
    
    % Setup Stat arrays
    if LME_Comparison==1,stats_dim = 1+n_conds;else, stats_dim = 1;end
    STATS.P = nan(DIM_time,DIM_ROI,stats_dim);
    STATS.Est = nan(DIM_time,DIM_ROI,stats_dim); % Code for polarity
    STATS.SE = nan(DIM_time,DIM_ROI,stats_dim); % Code for polarity
    STATS.T = nan(DIM_time,DIM_ROI,stats_dim); % Code for polarity
    STATS.DF = nan(DIM_time,DIM_ROI,stats_dim); % Code for polarity
    STATS.ANOVA = cell(DIM_time,DIM_ROI,stats_dim);
    
    if strcmp(Stat_Comp_Label, 'AV_Noise_Interaction')
         STATS.P_AV = nan(DIM_time,DIM_ROI,stats_dim);
         STATS.P_Noise = nan(DIM_time,DIM_ROI,stats_dim);
         
    elseif strcmp(Stat_Comp_Label, 'AV_McGurk_Interaction')
         STATS.P_AV = nan(DIM_time,DIM_ROI,stats_dim);
         STATS.P_McGurk = nan(DIM_time,DIM_ROI,stats_dim);
    end
        
    % Create tables for LME_Comparisons 1-3 that consider time independent
    % Iterate through time and ROI as necessary
    if ismember(LME_Comparison,1:3)
%         f = waitbar(0,'Creating Table');
        for time_x=1:DIM_time
            for ROI_x=1:DIM_ROI
%                 waitbar((((time_x-1)*DIM_ROI)+ROI_x)/(DIM_time*DIM_ROI),f) % Update waitbar
                
                DATA_ALL.Combined_Temp = [DATA_ALL.Text repmat({time_x},[size(DATA_ALL.Text,1) 1]) num2cell(mean(DATA_ALL.Values(:,time_x),2))];
                
                if ismember(LME_Comparison,1:2) % Select only current ROI
                    temp = find(~contains(DATA_ALL.Combined_Temp(:,5),['x' num2str(ROI_x)]));
                    DATA_ALL.Combined_Temp(temp,:) = [];
                end
                
                
%                 % Zscore values within sub (regardless of cond)
%                 grouped_sub_array = unique(DATA_ALL.Combined_Temp(:,1));
%                 for subx=1:length(grouped_sub_array)
%                     temp_dat = find(strcmp(DATA_ALL.Combined_Temp(:,1),grouped_sub_array(subx)));
%                     DATA_ALL.Combined_Temp(temp_dat,8) = num2cell(zscore(cell2mat(DATA_ALL.Combined_Temp(temp_dat,8))));
%                 end
               
                tbl = cell2table(DATA_ALL.Combined_Temp,'VariableNames',{'Subject_ID','Trial_Cond','Electrode_name','FrequencyBand','ROI','Trial_number','Time','ECoG_value'});
                %                 tbl_cond1 = tbl;
                %                 tbl_cond2 = tbl;
                
                if strcmp(Stat_Comp_Label, 'AV_Noise_Interaction')
                    NoiseFactor = nan(size(tbl,1),1);
                    AVFactor = nan(size(tbl,1),1);
                    for tbl_item_idx = 1: size(tbl,1)
                        switch tbl.Trial_Cond{tbl_item_idx}
                            case 'c1'
                               NoiseFactor(tbl_item_idx) = 1;
                               AVFactor(tbl_item_idx) = 0;
                            case 'c2'
                               NoiseFactor(tbl_item_idx) = 0;
                               AVFactor(tbl_item_idx) = 0;
                            case 'c3'
                               NoiseFactor(tbl_item_idx) = 1;
                               AVFactor(tbl_item_idx) = 1;
                            case 'c4'
                               NoiseFactor(tbl_item_idx) = 0;
                               AVFactor(tbl_item_idx) = 1;
                        end
                    end
                    tbl = addvars(tbl,NoiseFactor,'Before','FrequencyBand');
                    tbl = addvars(tbl,AVFactor,'Before','FrequencyBand');
                    tbl.NoiseFactor = categorical(tbl.NoiseFactor);
                    tbl.AVFactor = categorical(tbl.AVFactor);
                    clear NoiseFactor AVFactor
                elseif strcmp(Stat_Comp_Label, 'AV_McGurk_Interaction')
                    McGurkFactor = nan(size(tbl,1),1);
                    AVFactor = nan(size(tbl,1),1);
                    for tbl_item_idx = 1: size(tbl,1)
                        switch tbl.Trial_Cond{tbl_item_idx}
                            case 'c1'
                               McGurkFactor(tbl_item_idx) = 1;
                               AVFactor(tbl_item_idx) = 0;
                            case 'c2'
                               McGurkFactor(tbl_item_idx) = 0;
                               AVFactor(tbl_item_idx) = 0;
                            case 'c3'
                               McGurkFactor(tbl_item_idx) = 1;
                               AVFactor(tbl_item_idx) = 1;
                            case 'c4'
                               McGurkFactor(tbl_item_idx) = 0;
                               AVFactor(tbl_item_idx) = 1;
                        end
                    end
                    tbl = addvars(tbl,McGurkFactor,'Before','FrequencyBand');
                    tbl = addvars(tbl,AVFactor,'Before','FrequencyBand');
                    tbl.McGurkFactor = categorical(tbl.McGurkFactor);
                    tbl.AVFactor = categorical(tbl.AVFactor);
                    clear McGurkFactor AVFactor
                end
                
                % Create separate tables for each effect alone 
                for condx=1:n_conds
                    tbl_conds.(['x',num2str(condx)]) = tbl;
                    rm_cond = unique(tbl.Trial_Cond)';
                    rm_cond(condx)=[];
                    tbl_conds.(['x',num2str(condx)])(find(contains(tbl_conds.(['x',num2str(condx)]).Trial_Cond,rm_cond)),:) = [];
                    
                    % Change into categorical variables
                    tbl_conds.(['x',num2str(condx)]).Subject_ID = categorical(tbl_conds.(['x',num2str(condx)]).Subject_ID);
                    tbl_conds.(['x',num2str(condx)]).Trial_Cond = categorical(tbl_conds.(['x',num2str(condx)]).Trial_Cond);
                    tbl_conds.(['x',num2str(condx)]).Electrode_name = categorical(tbl_conds.(['x',num2str(condx)]).Electrode_name);
                    tbl_conds.(['x',num2str(condx)]).FrequencyBand = categorical(tbl_conds.(['x',num2str(condx)]).FrequencyBand);
                    tbl_conds.(['x',num2str(condx)]).Trial_number = categorical(tbl_conds.(['x',num2str(condx)]).Trial_number);
                    tbl_conds.(['x',num2str(condx)]).ROI = categorical(tbl_conds.(['x',num2str(condx)]).ROI);

                    % Average multiple electrodes within a single subject
                      
                    
                    if strcmp(Stat_Comp_Label, 'AV_Noise_Interaction')
                       tbl_conds.(['x',num2str(condx)])  = varfun(@mean,tbl_conds.(['x',num2str(condx)]),'InputVariables','ECoG_value','GroupingVariables',{'Subject_ID','Trial_number','Trial_Cond','AVFactor','NoiseFactor','FrequencyBand','ROI'});
                  
                    elseif strcmp(Stat_Comp_Label, 'AV_McGurk_Interaction')
                        tbl_conds.(['x',num2str(condx)])  = varfun(@mean,tbl_conds.(['x',num2str(condx)]),'InputVariables','ECoG_value','GroupingVariables',{'Subject_ID','Trial_number','Trial_Cond','AVFactor','McGurkFactor','FrequencyBand','ROI'});   
                    else
                        tbl_conds.(['x',num2str(condx)])  = varfun(@mean,tbl_conds.(['x',num2str(condx)]),'InputVariables','ECoG_value','GroupingVariables',{'Subject_ID','Trial_number','Trial_Cond','FrequencyBand','ROI'});
                    end
                    
                    % Create cell arrays of tbl for each vertex
                    tbl_array_conds.(['x',num2str(condx)])(time_x,ROI_x) = {tbl_conds.(['x',num2str(condx)])};
                end

                % Change into categorical variables
                tbl.Subject_ID = categorical(tbl.Subject_ID);
                tbl.Trial_Cond = categorical(tbl.Trial_Cond);
                tbl.Electrode_name = categorical(tbl.Electrode_name);
                tbl.FrequencyBand = categorical(tbl.FrequencyBand);
                tbl.Trial_number = categorical(tbl.Trial_number);
                tbl.ROI = categorical(tbl.ROI);
                
                % Average multiple electrodes within a single subject
                % Match trials
                if strcmp(Stat_Comp_Label, 'AV_Noise_Interaction')
                   tbl  = varfun(@mean,tbl,'InputVariables','ECoG_value','GroupingVariables',{'Subject_ID','Trial_number','Trial_Cond','AVFactor','NoiseFactor','FrequencyBand','ROI'});
                elseif strcmp(Stat_Comp_Label, 'AV_McGurk_Interaction')
                   tbl  = varfun(@mean,tbl,'InputVariables','ECoG_value','GroupingVariables',{'Subject_ID','Trial_number','Trial_Cond','AVFactor','McGurkFactor','FrequencyBand','ROI'});
                else 
                   tbl  = varfun(@mean,tbl,'InputVariables','ECoG_value','GroupingVariables',{'Subject_ID','Trial_number','Trial_Cond','FrequencyBand','ROI'});
                   % change 1/30/24
%                    tbl  = varfun(@mean,tbl,'InputVariables','ECoG_value','GroupingVariables',{'Subject_ID','Electrode_name', 'Trial_Cond','FrequencyBand','ROI'});                
                end
                
%                 if strcmp(Stat_Comp_Label, 'AV_Noise_Interaction')
%                     NoiseFactor = nan(size(tbl,1),1);
%                     AVFactor = nan(size(tbl,1),1);
%                     for tbl_item_idx = 1: size(tbl,1)
%                         switch tbl.Trial_Cond(tbl_item_idx)
%                             case 'c1'
%                                NoiseFactor(tbl_item_idx) = 1;
%                                AVFactor(tbl_item_idx) = 0;
%                             case 'c2'
%                                NoiseFactor(tbl_item_idx) = 0;
%                                AVFactor(tbl_item_idx) = 0;
%                             case 'c3'
%                                NoiseFactor(tbl_item_idx) = 1;
%                                AVFactor(tbl_item_idx) = 1;
%                             case 'c4'
%                                NoiseFactor(tbl_item_idx) = 0;
%                                AVFactor(tbl_item_idx) = 1;
%                         end
%                     end
%                     tbl = addvars(tbl,NoiseFactor,'Before','FrequencyBand');
%                     tbl = addvars(tbl,AVFactor,'Before','FrequencyBand');
%                     tbl.NoiseFactor = categorical(tbl.NoiseFactor);
%                     tbl.AVFactor = categorical(tbl.AVFactor);
%                     clear NoiseFactor AVFactor
%                 elseif strcmp(Stat_Comp_Label, 'AV_McGurk_Interaction')
%                     McGurkFactor = nan(size(tbl,1),1);
%                     AVFactor = nan(size(tbl,1),1);
%                     for tbl_item_idx = 1: size(tbl,1)
%                         switch tbl.Trial_Cond(tbl_item_idx)
%                             case 'c1'
%                                McGurkFactor(tbl_item_idx) = 1;
%                                AVFactor(tbl_item_idx) = 0;
%                             case 'c2'
%                                McGurkFactor(tbl_item_idx) = 0;
%                                AVFactor(tbl_item_idx) = 0;
%                             case 'c3'
%                                McGurkFactor(tbl_item_idx) = 1;
%                                AVFactor(tbl_item_idx) = 1;
%                             case 'c4'
%                                McGurkFactor(tbl_item_idx) = 0;
%                                AVFactor(tbl_item_idx) = 1;
%                         end
%                     end
%                     tbl = addvars(tbl,McGurkFactor,'Before','FrequencyBand');
%                     tbl = addvars(tbl,AVFactor,'Before','FrequencyBand');
%                     tbl.McGurkFactor = categorical(tbl.McGurkFactor);
%                     tbl.AVFactor = categorical(tbl.AVFactor);
%                     clear McGurkFactor AVFactor
%                 end
%                 
% % % % % %                 error('a')
% % % % %                 % Zscore values within sub (regardless of cond)
% % % % %                 grouped_sub_array = unique(tbl.Subject_ID);
% % % % %                 for subx=1:length(grouped_sub_array)
% % % % %                     temp_dat = find(tbl.Subject_ID==grouped_sub_array(subx));
% % % % %                     tbl.mean_ECoG_value(temp_dat) = zscore(tbl.mean_ECoG_value(temp_dat));
% % % % %                 end
                
                
                
                % Create cell arrays of tbl for each vertex
                % each time and vertex has 
                tbl_array_main(time_x,ROI_x) = {tbl};
            end
        end
%         close(f)
    end
%     save(['Power_' Stat_Comp_Label '_byElec.mat'],'tbl_array_main')
%     return
%{    
    %%
    
    % Create tables for LME_Comparison 4
    % One table for everything
    if ismember(LME_Comparison,4)
        DATA_ALL.Combined_Temp_Array = cell(1,length(time_rngs));
        for time_x=1:length(time_rngs)
            % Populate the array with each timepoint
            DATA_ALL.Combined_Temp_Array(time_x) = {[DATA_ALL.Text repmat({time_x},[size(DATA_ALL.Text,1) 1]) num2cell(mean(DATA_ALL.Values(:,time_x),2))]};
        end
        % Merge data from array
        DATA_ALL.Combined_Temp = [cat(1,DATA_ALL.Combined_Temp_Array{:})]; % Transpose because 2x faster to fill matrix by columns
        
        % Create Table
        tbl = cell2table(DATA_ALL.Combined_Temp,'VariableNames',{'Subject_ID','Trial_Cond','Electrode_name','FrequencyBand','ROI','Trial_number','Time','ECoG_value'});
        if strcmp(Stat_Comp,'C_vs_A')
            % Pull Aud and Cong
            tbl([find(contains(tbl.Trial_Cond,'vis'));find(contains(tbl.Trial_Cond,'incong'))],:) = [];
            tbl(find(contains(tbl.Trial_Cond,'aud')),2) ={'0'};
            tbl(find(contains(tbl.Trial_Cond,'cong')),2) = {'1'};  % Pos Vals Cong > Aud
        end
        
        % Change into categorical variables
        tbl.Subject_ID = categorical(tbl.Subject_ID);
        tbl.Trial_Cond = categorical(tbl.Trial_Cond);
        tbl.Electrode_name = categorical(tbl.Electrode_name);
        tbl.FrequencyBand = categorical(tbl.FrequencyBand);
        tbl.Trial_number = categorical(tbl.Trial_number);
        tbl.ROI = categorical(tbl.ROI);
        
        % Average multiple electrodes within a single subject
        % Match trials
        tbl = varfun(@mean,tbl,'InputVariables','ECoG_value','GroupingVariables',{'Subject_ID','Trial_number','Trial_Cond','FrequencyBand','ROI','Time'});
        
%         % Zscore values within sub (regardless of cond)
%         grouped_sub_array = unique(grouped_data.Subject_ID);
%         for subx=1:length(grouped_sub_array)
%             temp_dat = find(grouped_data.Subject_ID==grouped_sub_array(subx));
%             grouped_data.mean_ECoG_value(temp_dat) = zscore(grouped_data.mean_ECoG_value(temp_dat));
%         end

        
        % Create cell arrays of tbl for each vertex
        tbl_array_main = {tbl};
    end
%}   
    %%

    % Clear large variables to free up memory
    clear DATA_ALL
   
    DATA_OUTPUT = cell(DIM_time,DIM_ROI);
    DATA_OUTPUT_INTERACTION = cell(DIM_time,1);
    clear DATA_OUTPUT_CONDS
    
    for condx=1:n_conds
        DATA_OUTPUT_CONDS.(['x',num2str(condx)]) = cell(DIM_time,DIM_ROI);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% RUN THE LME MODEL %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nRow,nCol] = size(tbl_array_main);
    for rowIdx = 1: nRow
        for colIdx = 1: nCol
            tbl_array_main{rowIdx,colIdx}.ROI = categorical(tbl_array_main{rowIdx,colIdx}.ROI);
            tbl_array_main{rowIdx,colIdx}.Trial_Cond = categorical(tbl_array_main{rowIdx,colIdx}.Trial_Cond);
        end
    end
    if ismember(LME_Comparison,1:3)
        % Start up par before parfor so that the progressbar works
        p = gcp;
        % Iterate through time and ROI as necessary
        for ROI_x=1:DIM_ROI
%             f = parfor_progressbar(DIM_time,['ROI ' num2str(ROI_x) ' of ' num2str(DIM_ROI) ', Running LME']); %create the progress bar
            parfor time_x=1:DIM_time
                if LME_Comparison==1                    
                    if strcmp(Stat_Comp_Label, 'AV_Noise_Interaction')
                        DATA_OUTPUT(time_x,ROI_x) = {fitlme(tbl_array_main{time_x,ROI_x},'mean_ECoG_value~AVFactor*NoiseFactor + (1|Subject_ID)','FitMethod','REML')};
                    elseif strcmp(Stat_Comp_Label, 'AV_McGurk_Interaction')
                        DATA_OUTPUT(time_x,ROI_x) = {fitlme(tbl_array_main{time_x,ROI_x},'mean_ECoG_value~AVFactor*McGurkFactor + (1|Subject_ID)','FitMethod','REML')};
                    else
                    %                 DATA_OUTPUT(roiIdx,1) = {fitlme(grouped_data,'mean_ECoG_value~Trial_Cond + (Trial_Cond|Subject_ID) ','FitMethod','REML')};
                    %                 DATA_OUTPUT(time_x,ROI_x) = {fitlme(tbl_array_main{time_x,ROI_x},'mean_ECoG_value~Trial_Cond + (1|Subject_ID) + (1|Subject_ID:Trial_number) + (1|Subject_ID:Electrode_name)','FitMethod','REML')};
                    %                 DATA_OUTPUT_cond1(time_x,ROI_x) = {fitlme(tbl_array_cond1{time_x,ROI_x},'ECoG_value~ 1 + (1|Subject_ID) + (1|Subject_ID:Trial_number) + (1|Subject_ID:Electrode_name)','FitMethod','REML')};
                    %                 DATA_OUTPUT_cond2(time_x,ROI_x) = {fitlme(tbl_array_cond2{time_x,ROI_x},'ECoG_value~ 1 + (1|Subject_ID) + (1|Subject_ID:Trial_number) + (1|Subject_ID:Electrode_name)','FitMethod','REML')};
                    
                    % here we should have everything needed to compute ci
                        DATA_OUTPUT(time_x,ROI_x) = {fitlme(tbl_array_main{time_x,ROI_x},...
                            'mean_ECoG_value~Trial_Cond + (Trial_Cond|Subject_ID)','FitMethod','REML')};
                    end
                    
                elseif LME_Comparison==2
                    DATA_OUTPUT(time_x,ROI_x) = {anova(fitlme(tbl_array_main{time_x,ROI_x},'mean_ECoG_value~Trial_Cond*FrequencyBand + (Trial_Cond|Subject_ID)','FitMethod','REML'))};
                elseif LME_Comparison==3
                    %                 DATA_OUTPUT(time_x,ROI_x) = {anova(fitlme(tbl_array_main{time_x,ROI_x},'ECoG_value~Trial_Cond*FrequencyBand*ROI + (1|Subject_ID) + (1|Subject_ID:Trial_number) + (1|Subject_ID:Electrode_name)','FitMethod','REML'))};
                    DATA_OUTPUT(time_x,ROI_x) = {anova(fitlme(tbl_array_main{time_x,ROI_x},'mean_ECoG_value~Trial_Cond*FrequencyBand*ROI + (Trial_Cond|Subject_ID)','FitMethod','REML'))};
                end
            end
        end

        % run interaction
        % first combine across ROIs
        numTimepoints = size(tbl_array_main, 1); % Assuming 151 timepoints
        tbl_array_combined = cell(numTimepoints, 1); % Preallocate a new cell array
        
        for time_x = 1:numTimepoints
            % Extract all tables for this timepoint
            tables_for_timepoint = tbl_array_main(time_x, :);
            
            % Concatenate tables across the ROI dimension
            combined_table = vertcat(tables_for_timepoint{:});
            
            % Store the combined table in the new cell array
            tbl_array_combined{time_x} = combined_table;
        end

        DATA_OUTPUT_INTERACTION = cell(numTimepoints, 1); % Preallocate cell array for results

        parfor time_x = 1:numTimepoints
            % Fit the LME model to the combined table for this timepoint
            DATA_OUTPUT_INTERACTION{time_x} = anova(fitlme(tbl_array_combined{time_x}, ...
                                'mean_ECoG_value ~ Trial_Cond*ROI + (Trial_Cond|Subject_ID)', ...
                                'FitMethod', 'REML'));
        end

        % apply FDR correction to the interaction significance
        % Step 1: Extract all interaction p-values
        interaction_p_values = zeros(numTimepoints, 1);
        for i = 1:numTimepoints
            % Assuming the interaction p-value is always in the same position (5th column of the 4th row)
            interaction_p_values(i) = DATA_OUTPUT_INTERACTION{i}.pValue(4);
        end
        
        % Step 2: Apply FDR correction
        [~, ~, ~, adjusted_p_values] = fdr_bh(interaction_p_values);
        
        % Step 3: Assign corrected p-values back to the cell array
        for i = 1:numTimepoints
            DATA_OUTPUT_INTERACTION{i}.pValue(4) = adjusted_p_values(i);
        end

        % Assuming 'adjusted_p_values' is a vector of FDR-corrected p-values for the interaction term
        timepoints = linspace(-1, 0.5, length(adjusted_p_values)); % Create a timepoints vector if needed
        significant_threshold = 0.05; % Define a threshold for significance
        
        % Create a logical vector where true corresponds to significant p-values
        is_significant = adjusted_p_values < significant_threshold;
        
        % Visualization
        figure; 
        hold on;
        plot(timepoints, adjusted_p_values, 'o-'); % Plot all corrected p-values
        plot(timepoints(is_significant), adjusted_p_values(is_significant), 'ro'); % Emphasize significant p-values
        line([min(timepoints), max(timepoints)], [significant_threshold significant_threshold], 'Color', 'red', 'LineStyle', '--'); % Significance threshold line
        xlabel('Time')
        ylabel('FDR Corrected p-value');
        title(['ROI x Condition Interaction '  strrep(Stat_Comp_Label,'_',' ')]);
        legend({'Corrected p-values', 'Significant p-values'}, 'Location', 'northoutside');
        hold off;


        % dont need this
        % calculate CI from the between cond
        % change 

        % do 2x2 repeated measures anova for ITPC and Power
        % power, averaged across trials and electrodes, one value for each
        % cond each subject (120-125ms)

        % itpc, averaged across electrodes, one value for each cond each
        % subject

        % to demonstrate itpc and power effects come from two different
        % populations 



        % now run the LME model for every single condition 
        if LME_Comparison==1
            for condx=1:n_conds
                DATA_OUTPUT_CURR_COND = DATA_OUTPUT_CONDS.(['x',num2str(condx)]);
                for ROI_x=1:DIM_ROI
                    parfor time_x=1:DIM_time

                         if strcmp(Stat_Comp_Label, 'AV_Noise_Interaction')
                            DATA_OUTPUT_CURR_COND(time_x,ROI_x) = {fitlme(tbl_array_conds.(['x',num2str(condx)]){time_x,ROI_x},'mean_ECoG_value~ 1 + (AVFactor*NoiseFactor|Subject_ID)','FitMethod','REML')};
                         
                         elseif strcmp(Stat_Comp_Label, 'AV_McGurk_Interaction')
                            DATA_OUTPUT_CURR_COND(time_x,ROI_x) = {fitlme(tbl_array_conds.(['x',num2str(condx)]){time_x,ROI_x},'mean_ECoG_value~ 1 + (AVFactor*McGurkFactor|Subject_ID)','FitMethod','REML')};
                         
                         else
                            DATA_OUTPUT_CURR_COND(time_x,ROI_x) = {fitlme(tbl_array_conds.(['x',num2str(condx)]){time_x,ROI_x},'mean_ECoG_value~ 1 + (Trial_Cond|Subject_ID)','FitMethod','REML')};
                        end
                    end   
                end
                DATA_OUTPUT_CONDS.(['x',num2str(condx)])=DATA_OUTPUT_CURR_COND;
            end
        end



    end 
    %% Pull stats and save 
    for time_x=1:DIM_time
        for ROI_x=1:DIM_ROI    
            if LME_Comparison==1
                ROI_timepoint_Coef = dataset2cell(DATA_OUTPUT{time_x,ROI_x}.Coefficients);
                STATS.ANOVA(time_x,ROI_x,1)={ROI_timepoint_Coef};
                if strcmp(Stat_Comp_Label, 'AV_McGurk_Interaction')
                    STATS.P_AV(time_x,ROI_x,1) = ROI_timepoint_Coef{4,6};
                    STATS.P_McGurk(time_x,ROI_x,1) = ROI_timepoint_Coef{4,6};
                    STATS.P_interact(time_x,ROI_x,1) = ROI_timepoint_Coef{5,6};
                elseif strcmp(Stat_Comp_Label, 'AV_Noise_Interaction')
                    STATS.P_AV(time_x,ROI_x,1) = ROI_timepoint_Coef{4,6};
                    STATS.P_Noise(time_x,ROI_x,1) = ROI_timepoint_Coef{4,6};
                    STATS.P_interact(time_x,ROI_x,1) = ROI_timepoint_Coef{5,6};
                else
                    % so in the first column...
                    % 'Trial_Cond_c2''s p value 

                     STATS.P(time_x,ROI_x,1) = ROI_timepoint_Coef{3,6};
                end
                % check hard coded indices here
                % the first column is for the overall one
                STATS.Est(time_x,ROI_x,1) = ROI_timepoint_Coef{3,2};
                STATS.SE(time_x,ROI_x,1) = ROI_timepoint_Coef{3,3};
                STATS.T(time_x,ROI_x,1) = ROI_timepoint_Coef{3,4};
                STATS.DF(time_x,ROI_x,1) = ROI_timepoint_Coef{3,5};
                

%                 the second and third columns are for other conds
                for condx=1:n_conds
                    ROI_timepoint_Coef_byCond = dataset2cell(DATA_OUTPUT_CONDS.(['x',num2str(condx)]){time_x,ROI_x}.Coefficients);
                    STATS.ANOVA(time_x,ROI_x,condx+1)={ROI_timepoint_Coef_byCond};
                    STATS.P(time_x,ROI_x,condx+1) = ROI_timepoint_Coef_byCond{2,6};       
                    STATS.Est(time_x,ROI_x,condx+1) = ROI_timepoint_Coef_byCond{2,2};
                    STATS.SE(time_x,ROI_x,condx+1) = ROI_timepoint_Coef_byCond{2,3};
                    STATS.T(time_x,ROI_x,condx+1) = ROI_timepoint_Coef_byCond{2,4};
                    STATS.DF(time_x,ROI_x,condx+1) = ROI_timepoint_Coef_byCond{2,5};
                end

            elseif ismember(LME_Comparison,2:4)
                ROI_timepoint_Coef = anova(DATA_OUTPUT{time_x,ROI_x});
                % ROI_timepoint_Coef = DATA_OUTPUT{time_x,ROI_x};
                STATS.ANOVA(time_x,ROI_x)={ROI_timepoint_Coef};
                STATS.P(time_x,ROI_x) = ROI_timepoint_Coef.pValue(end);
                STATS.Est(time_x,ROI_x) = NaN;
                STATS.SE(time_x,ROI_x) = NaN;
                STATS.T(time_x,ROI_x) = ROI_timepoint_Coef.FStat(end);
                STATS.DF(time_x,ROI_x) = ROI_timepoint_Coef.DF2(end);
            end
        end
    end
    
    % Save stats
%     save([outputDir,Stat_Comp_Label,'_',Analysis_Title,'_TimeSeries.mat'],'ROI_Num','orthogTime','STATS','time_rngs','Stat_Comp','Freq_Test','LME_Comparison','Analysis_Title','-v7.3')
    
    VarInteract = STATS.Est(:,:,2:3);
    [numTimepoints, numROIs, numConditions] = size(VarInteract);
    
    % Create vectors for each dimension
    timepoints = repmat((1:numTimepoints)', [numROIs*numConditions, 1]);
    rois = repmat(repelem(1:numROIs, numTimepoints), [1, numConditions]);
    conditions = repmat(repelem(1:numConditions, numTimepoints*numROIs), [1, 1]);
    
    % Reshape VarInteract into a column vector
    values = reshape(VarInteract, [], 1);
    
    % Create a table
    data4Interaction = table(timepoints, rois', conditions', values, 'VariableNames', {'Timepoint', 'ROI', 'Condition', 'Value'});
    
    % Two-way ANOVA with interaction
    [p, tbl, stats] = anovan(data4Interaction.Value, ...
                             {data4Interaction.ROI, data4Interaction.Condition}, ...
                             'model', 'interaction', ...
                             'varnames', {'ROI', 'Condition'});
    
    % Post-hoc analysis if interaction is significant
    if p(3) < 0.05 % Assuming the third p-value corresponds to the interaction
        results = multcompare(stats, 'Dimension', [1 2]);
    end
    
    % Interaction plot
    figure;
    interactionplot(data4Interaction.Value, ...
                    {data4Interaction.ROI, data4Interaction.Condition}, ...
                    'varnames', {'ROI', 'Condition'});
    sgtitle(strrep(Stat_Comp_Label,'_',' '))
    %%
    % Plot if LME_Comparison == 1
    if LME_Comparison==1
        % Plot vertices
        legend_handles = [];
        % 1=red, 2=blue, 3=green, 4=cyan, 5=yellow, 6=magenta, 7=orange, 8=blue-green, light blue,fuscia,lime-green,purple, grey,black
        colorarray=[190/255 32/255 38/255;39/255 126/255 183/255;59/255 181/255 74/255; 251/255 191/255 21/255; .188 .188 .188;1 0 1; 1 .5 0; 0 1 .5;0 .5 1;1 0 .5;.5 1 0;.5 0 1;.5 .5 .5; 0 0 0];
        color_conds = [2 1 3:length(colorarray)];
        LineThickness = 2.5;
        % legend_array = {'Congruent','Auditory-Alone'};
        legend_array=[Stat_Labels];
        
        
%         legend_handles = [];
        close all
        f1=figure;
        %%%%%%%%%%%%%%%%%% MAKING THE THREE SUBPLOTS START %%%%%%%%%%%%%%%%%%
        for roiIdx = 1:ROI_Num
            subplot(1,ROI_Num,4-roiIdx) % 4-roiIdx to put anterior on left
            title(ROI_names{roiIdx}) % add the name of ROI to subtitle to make it clear
            
            currROI_est = squeeze(STATS.Est(:,roiIdx,:));
            currROI_p = squeeze(STATS.P(:,roiIdx,:));
            
            if strcmp(Stat_Comp_Label, 'AV_Noise_Interaction')
                temp_p_av = squeeze(STATS.P_AV(:,roiIdx,:));
                temp_p_noise = squeeze(STATS.P_Noise(:,roiIdx,:));
                temp_p_interact = squeeze(STATS.P_interact(:,roiIdx,:));
                
            elseif strcmp(Stat_Comp_Label, 'AV_McGurk_Interaction')
                temp_p_av = squeeze(STATS.P_AV(:,roiIdx,:));
                temp_p_mcgurk = squeeze(STATS.P_McGurk(:,roiIdx,:));
                temp_p_interact = squeeze(STATS.P_interact(:,roiIdx,:));
            end
            
            currROI_SE = squeeze(STATS.SE(:,roiIdx,:))*1.96/2; % convert to CI
            
            
            %% Now we have pulled all the data for plotting

            figHolder = []; 

            % Initialize an array to store line handles
            mainLineHandles = gobjects(n_conds, 1);

            
             % Plotting main lines and adding their handles to legend_handles
             % it starts from column 2 because column 1 is the stat for overall regardless of conditions 
            for condx = 2:n_conds+1
                hx = shadedErrorBar(time(time_rngs),...
                    currROI_est(:,condx),...% change 3/4 to the first column, which is the standard error 
                    currROI_SE(:,1),...% change 3/4 to the first column
                    'k',...
                    .7);
                set(hx.mainLine, ...
                    'LineWidth',LineThickness,...
                    'Color', colorarray(color_conds(condx-1),:));
                set(hx.patch, ...
                    'facecolor', colorarray(color_conds(condx-1),:),...
                    'facealpha', .25);
                delete(hx.edge);

                 % Store the handle of the main line
                mainLineHandles(condx-1) = hx.mainLine;

                hold on;
            end
           
                        
            figHolder = [figHolder hx.mainLine];

            if Freq_Test==1
                ylim([-.1 .35])
                yAVBar = -0.05;
                yMcGurkBar = -0.07;
                yNoiseBar = -0.07;
                yInteractBar = -0.09;
                
            elseif Freq_Test==2
                ylim([-.05 .25])
                yAVBar = -0.1;
                yMcGurkBar = -0.12;
                yNoiseBar = -0.12;
                yInteractBar = -0.14;
            elseif Freq_Test==3
                ylim([-.2 .15])
                yAVBar = -0.15;
                yMcGurkBar = -0.17;
                yNoiseBar = -0.17;
                yInteractBar = -0.19;
            end
            

            %%%%%%%%%%%%%%%% adding in zero lines, horizontal and vertical %%%%%%%%%%%%%%%% 

            % Zero bar, horizontal
            temptime = time(time_rngs);
            tempLim = get(gca,'ylim');
            LineColor = [.5 .5 .5]+.1;
            tempx = [temptime(1) temptime(end)];
            tempy = [0 0];
            hy = line(tempx,tempy,'Color',LineColor);
            set(hy, 'LineStyle', '--')
            set(hy, 'LineWidth', 2)
            
            % Zero bar, vertical
            tempmin = tempLim(1);
            tempmax = tempLim(2);
            tempy = [tempmin tempmax];
            tempx = [0 0];
            hy = line(tempx,tempy,'Color',LineColor);
            set(hy, 'LineStyle', '--');
            set(hy, 'LineWidth', 2)
                        
            % Add text
            %     h=text(-.10,30,'Sound Onset');
            %     set(h,'Rotation',90);
            
            xlabel('Time (Seconds)');
            ylabel(['a.u.' ]);
            
            %legend(figHolder,legend_array,'Location','Northwest');
           
            % ylim([tempy]);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           










            %%%%%%%%%%%%%%%%%% MAKING THE STATS BOXES Start %%%%%%%%%%%%%%%%%%
             % STATS
            % run on just -.5 to -.5
            analysis_rng = dsearchn(time(time_rngs)',-.5):1:dsearchn(time(time_rngs)',.5);
            
            if strcmp(Stat_Comp_Label, 'AV_Noise_Interaction')
                adj_p_av = nan(1,length(temp_p_av(:,1)));
                [h_av, crit_p_av, adj_ci_cvrg_av, adj_p_av(analysis_rng)] = fdr_bh(temp_p_av(analysis_rng,1));
                tempx_av = zeros(1,length(adj_p_av));tempx_av(adj_p_av<.05)=1;
                k_av=find(tempx_av==1);
                
                adj_p_noise = nan(1,length(temp_p_noise(:,1)));
                [h_noise, crit_p_noise, adj_ci_cvrg_noise, adj_p_noise(analysis_rng)] = fdr_bh(temp_p_noise(analysis_rng,1));
                tempx_noise = zeros(1,length(adj_p_noise));tempx_noise(adj_p_noise<.05)=1;
                k_noise=find(tempx_noise==1);
                
                adj_p_interact = nan(1,length(temp_p_interact(:,1)));
                [h_interact, crit_p_interact, adj_ci_cvrg_interact, adj_p_interact(analysis_rng)] = fdr_bh(temp_p_interact(analysis_rng,1));
                tempx_interact = zeros(1,length(adj_p_interact));tempx_interact(adj_p_interact<.05)=1;
                k_interact=find(tempx_interact==1);
                % = [k_av,k_noise,k_interact];
                
            elseif strcmp(Stat_Comp_Label, 'AV_McGurk_Interaction')
                adj_p_av = nan(1,length(temp_p_av(:,1)));
                [h_av, crit_p_av, adj_ci_cvrg_av, adj_p_av(analysis_rng)] = fdr_bh(temp_p_av(analysis_rng,1));
                tempx_av = zeros(1,length(adj_p_av));tempx_av(adj_p_av<.05)=1;
                k_av=find(tempx_av==1);
               
                
                adj_p_mcgurk = nan(1,length(temp_p_mcgurk(:,1)));
                [h_mcgurk, crit_p_mcgurk, adj_ci_cvrg_mcgurk, adj_p_mcgurk(analysis_rng)] = fdr_bh(temp_p_mcgurk(analysis_rng,1));
                tempx_mcgurk = zeros(1,length(adj_p_mcgurk));tempx_noise(adj_p_mcgurk<.05)=1;
                k_mcgurk=find(tempx_mcgurk==1);
                
                adj_p_interact = nan(1,length(temp_p_interact(:,1)));
                [h_interact, crit_p_interact, adj_ci_cvrg_interact, adj_p_interact(analysis_rng)] = fdr_bh(temp_p_interact(analysis_rng,1));
                tempx_interact = zeros(1,length(adj_p_interact));tempx_interact(adj_p_interact<.05)=1;
                k_interact=find(tempx_interact==1);
                % = [k_av,k_mcgurk,k_interact];
            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% for all the other statistical tests %%%%
                %%%% i.e. the non-interaction ones %%%%
                % we want to run fdr on the existing statistics
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                adj_p = nan(1,length(currROI_p(:,1)));
                % run fdr to
                [h, crit_p, adj_ci_cvrg, adj_p(analysis_rng)] = fdr_bh(currROI_p(analysis_rng,1));
                tempx = zeros(1,length(adj_p));tempx(adj_p<.05)=1;
                fdrSigSampleIdx=find(tempx==1);
            end

            %%%%%%%%%%%%%%%%%% AV Noise Start %%%%%%%%%%%%%%%%%%
             if strcmp(Stat_Comp_Label, 'AV_Noise_Interaction')
                 
                 hy_av =[];hy_av_temp =[];
                if ~isempty(k_av)         
                    box_hits_av = [];
                    counter_av = 1; % Row
                    box_hits_av(1,1) = k_av(1);
                    nexthit_av = 0;
                    
                    for i=1:length(k_av)
                        if i<length(k_av)              
                            % open group if last iteration closed it
                            if nexthit_av == 1
                                nexthit_av = 0;
                                box_hits_av(counter_av,1) = k_av(i);
                            end

                            % close group
                            if k_av(i)~=k_av(i+1)-1
                                box_hits_av(counter_av,2) = k_av(i);
                                counter_av = counter_av+1;
                                nexthit_av = 1;
                            end
                        else
                            if nexthit_av==0
                                % close group
                                box_hits_av(counter_av,2) = k_av(i);
                            else
                                box_hits_av(counter_av,1) = k_av(i);
                                box_hits_av(counter_av,2) = k_av(i);
                            end
                        end
                    end
                    
                    
                    temptime = time(time_rngs);
                    box_hits_av = temptime(box_hits_av);
                    
                    

                    for i=1:size(box_hits_av,1)
                        % xstart ystart xlength ylength
                        % hy = rectangle('Position',[box_hits(i,1),tempmin,box_hits(i,2)-box_hits(i,1),abs(tempmin*.5)],'FaceColor',[0.1 .1 .1]);

                        % try with patch
                        %tempx_av = [box_hits_av(i,1) box_hits_av(i,2) box_hits_av(i,2) box_hits_av(i,1)];
                        tempx_av = [box_hits_av(i,1) box_hits_av(i,2)];
                        tempy_av = [tempmin tempmin tempmin+abs(tempmin*.5) tempmin+abs(tempmin*.5)];
                      % hy_av = patch(tempx_av,tempy_av,'red');
                        hy_av_temp = line(tempx_av, [yAVBar yAVBar],'Color','#A2142F','Linewidth',5);
                        % c = [0 1 1 0];
                        % hy = patch(tempx,tempy,c);
                        % caxis([0 .5])
                        uistack(hy_av_temp,'bottom');
                        hy_av = [hy_av;hy_av_temp];
                    end
%                     ah1=axes('position',get(gca,'position'),'visible','off');
%                     legend(ah1,hy_av,'AV','Location','southeast','Orientation','horizontal');
                end    
                
               
                    
                 hy_noise =[]; hy_noise_temp =[];
                 
                    if ~isempty(k_noise)         
                        box_hits_noise = [];
                        counter_noise = 1; % Row
                        box_hits_noise(1,1) = k_noise(1);
                        nexthit_noise = 0;

                        for i=1:length(k_noise)
                            if i<length(k_noise)              
                                % open group if last iteration closed it
                                if nexthit_noise == 1
                                    nexthit_noise = 0;
                                    box_hits_noise(counter_noise,1) = k_noise(i);
                                end

                                % close group
                                if k_noise(i)~=k_noise(i+1)-1
                                    box_hits_noise(counter_noise,2) = k_noise(i);
                                    counter_noise = counter_noise+1;
                                    nexthit_noise = 1;
                                end
                            else
                                if nexthit_noise==0
                                    % close group
                                    box_hits_noise(counter_noise,2) = k_noise(i);
                                else
                                    box_hits_noise(counter_noise,1) = k_noise(i);
                                    box_hits_noise(counter_noise,2) = k_noise(i);
                                end
                            end
                        end
                        temptime = time(time_rngs);
                        box_hits_noise = temptime(box_hits_noise);

                        for i=1:size(box_hits_noise,1)
                            % xstart ystart xlength ylength
                            % hy = rectangle('Position',[box_hits(i,1),tempmin,box_hits(i,2)-box_hits(i,1),abs(tempmin*.5)],'FaceColor',[0.1 .1 .1]);

                            % try with patch
                            %tempx_noise = [box_hits_noise(i,1) box_hits_noise(i,2) box_hits_noise(i,2) box_hits_noise(i,1)];
                            tempx_noise = [box_hits_noise(i,1) box_hits_noise(i,2)];
                            tempy_noise = [tempmin tempmin tempmin+abs(tempmin*.5) tempmin+abs(tempmin*.5)];
                            %hy_mcgurk = patch(tempx_mcgurk,tempy_mcgurk,'green');
                            hy_noise_temp = line(tempx_av, [yNoiseBar yNoiseBar],'Color','#77AC30','Linewidth',5);      
                            % c = [0 1 1 0];
                            % hy = patch(tempx,tempy,c);
                            % caxis([0 .5])
                            uistack(hy_noise_temp,'bottom');
                            hy_noise = [hy_noise;hy_noise_temp];
                        end
                        
%                         if isempty(k_av)
%                             ah1=axes('position',get(gca,'position'),'visible','off');
%                             legend(ah1,hy_noise,'Noise','Location','southeast','Orientation','horizontal')
%                         else 
%                             legappend('Noise')
%                         end
                    end     

                    hy_interact =[];hy_interact_temp =[];
                    
                    if ~isempty(k_interact)         
                        box_hits_interact = [];
                        counter_interact = 1; % Row
                        box_hits_interact(1,1) = k_interact(1);
                        nexthit_interact = 0;

                        for i=1:length(k_interact)
                            if i<length(k_interact)              
                                % open group if last iteration closed it
                                if nexthit_interact == 1
                                    nexthit_interact = 0;
                                    box_hits_interact(counter_interact,1) = k_interact(i);
                                end

                                % close group
                                if k_interact(i)~=k_interact(i+1)-1
                                    box_hits_interact(counter_interact,2) = k_interact(i);
                                    counter_interact = counter_interact+1;
                                    nexthit_interact = 1;
                                end
                            else
                                if nexthit_interact==0
                                    % close group
                                    box_hits_interact(counter_interact,2) = k_interact(i);
                                else
                                    box_hits_interact(counter_interact,1) = k_interact(i);
                                    box_hits_interact(counter_interact,2) = k_interact(i);
                                end
                            end
                        end
                        temptime = time(time_rngs);
                        box_hits_interact = temptime(box_hits_interact);

                        for i=1:size(box_hits_interact,1)
                            % xstart ystart xlength ylength
                            % hy = rectangle('Position',[box_hits(i,1),tempmin,box_hits(i,2)-box_hits(i,1),abs(tempmin*.5)],'FaceColor',[0.1 .1 .1]);

                            % try with patch
                            %tempx_interact = [box_hits_interact(i,1) box_hits_interact(i,2) box_hits_interact(i,2) box_hits_interact(i,1)];
                            tempx_interact = [box_hits_interact(i,1) box_hits_interact(i,2)];
                            tempy_interact = [tempmin tempmin tempmin+abs(tempmin*.5) tempmin+abs(tempmin*.5)];
                            % hy_interact = patch(tempx_interact,tempy_interact*0.25,'blue');
                            % hy_interact = line(tempx_interact,-0.08,'Color','#0072BD','DisplayName','Interaction');
                            hy_interact_temp = line(tempx_av, [yInteractBar yInteractBar],'Color','#0072BD','Linewidth',5);
                        
                            % c = [0 1 1 0];
                            % hy = patch(tempx,tempy,c);
                            % caxis([0 .5])
                            uistack(hy_interact_temp,'bottom');
                            hy_interact = [hy_interact;hy_interact_temp];
                        end
              
%                         if isempty(k_av) && isempty(k_noise)
%                             ah1=axes('position',get(gca,'position'),'visible','off');
%                             legend(ah1,hy_interact,'Interaction','Location','southeast','Orientation','horizontal')
%                         else 
%                             legappend('Interaction')
%                         end
                    end
%                     legend(figHolder,legend_array,'Location','Northwest','Orientation','vertical');
                     %%%%%%%%%%%%%%%%%% AV Noise End %%%%%%%%%%%%%%%%%%
                      
                     
                     %%%%%%%%%%%%%%%%%% AV McGurk Start %%%%%%%%%%%%%%%%%%
                 
             elseif strcmp(Stat_Comp_Label, 'AV_McGurk_Interaction')
                if ~isempty(k_av)         
                    box_hits_av = [];
                    counter_av = 1; % Row
                    box_hits_av(1,1) = k_av(1);
                    nexthit_av = 0;
                    
                    for i=1:length(k_av)
                        if i<length(k_av)              
                            % open group if last iteration closed it
                            if nexthit_av == 1
                                nexthit_av = 0;
                                box_hits_av(counter_av,1) = k_av(i);
                            end

                            % close group
                            if k_av(i)~=k_av(i+1)-1
                                box_hits_av(counter_av,2) = k_av(i);
                                counter_av = counter_av+1;
                                nexthit_av = 1;
                            end
                        else
                            if nexthit_av==0
                                % close group
                                box_hits_av(counter_av,2) = k_av(i);
                            else
                                box_hits_av(counter_av,1) = k_av(i);
                                box_hits_av(counter_av,2) = k_av(i);
                            end
                        end
                    end
                    
                    
                    temptime = time(time_rngs);
                    box_hits_av = temptime(box_hits_av);
                    hy_av =[];hy_av_temp =[];
                    

                    for i=1:size(box_hits_av,1)
                        % xstart ystart xlength ylength
                        % hy = rectangle('Position',[box_hits(i,1),tempmin,box_hits(i,2)-box_hits(i,1),abs(tempmin*.5)],'FaceColor',[0.1 .1 .1]);

                        % try with patch
                        %tempx_av = [box_hits_av(i,1) box_hits_av(i,2) box_hits_av(i,2) box_hits_av(i,1)];
                        tempx_av = [box_hits_av(i,1) box_hits_av(i,2)];
                        % tempy_av = [tempmin tempmin tempmin+abs(tempmin*.5) tempmin+abs(tempmin*.5)];
                        tempy_av = -0.04;
                        % hy_av = patch(tempx_av,tempy_av,'red');
                        hy_av_temp = line(tempx_av, [yAVBar yAVBar],'Color','#A2142F','Linewidth',5,'DisplayName','AV');
                        legend_handles(end + 1) = hy_av_temp;
                        % c = [0 1 1 0];
                        % hy = patch(tempx,tempy,c);
                        % caxis([0 .5])
                        uistack(hy_av_temp,'bottom');
                        hy_av = [hy_av;hy_av_temp];
                    end
%                     ah1=axes('position',get(gca,'position'),'visible','off');
%                     legend(ah1,hy_av,'Location','southeast','Orientation','horizontal');  
%                     
                    
                end       
                
                    hy_mcgurk =[];hy_mcgurk_temp =[];
                    
                    if ~isempty(k_mcgurk)         
                        box_hits_mcgurk = [];
                        counter_mcgurk = 1; % Row
                        box_hits_mcgurk(1,1) = k_mcgurk(1);
                        nexthit_mcgurk = 0;

                        for i=1:length(k_mcgurk)
                            if i<length(k_mcgurk)              
                                % open group if last iteration closed it
                                if nexthit_mcgurk == 1
                                    nexthit_mcgurk = 0;
                                    box_hits_mcgurk(counter_mcgurk,1) = k_mcgurk(i);
                                end

                                % close group
                                if k_mcgurk(i)~=k_mcgurk(i+1)-1
                                    box_hits_mcgurk(counter_mcgurk,2) = k_mcgurk(i);
                                    counter_mcgurk = counter_mcgurk+1;
                                    nexthit_mcgurk = 1;
                                end
                            else
                                if nexthit_mcgurk==0
                                    % close group
                                    box_hits_mcgurk(counter_mcgurk,2) = k_mcgurk(i);
                                else
                                    box_hits_mcgurk(counter_mcgurk,1) = k_mcgurk(i);
                                    box_hits_mcgurk(counter_noise,2) = k_mcgurk(i);
                                end
                            end
                        end
                        temptime = time(time_rngs);
                        box_hits_mcgurk = temptime(box_hits_mcgurk);

                        for i=1:size(box_hits_mcgurk,1)
                            % xstart ystart xlength ylength
                            % hy = rectangle('Position',[box_hits(i,1),tempmin,box_hits(i,2)-box_hits(i,1),abs(tempmin*.5)],'FaceColor',[0.1 .1 .1]);

                            % try with patch
                            % tempx_mcgurk = [box_hits_mcgurk(i,1) box_hits_mcgurk(i,2) box_hits_mcgurk(i,2) box_hits_mcgurk(i,1)];
                            tempx_mcgurk = [box_hits_mcgurk(i,1) box_hits_mcgurk(i,2)];     
                            tempy_mcgurk = [tempmin tempmin tempmin+abs(tempmin*.5) tempmin+abs(tempmin*.5)];
                            %hy_mcgurk = patch(tempx_mcgurk,tempy_mcgurk,'green');
                            hy_mcgurk_temp = line(tempx_mcgurk, [yMcGurkBar yMcGurkBar],'Color','#77AC30','Linewidth',5,'DisplayName','McGurk');      
                            % c = [0 1 1 0];
                            % hy = patch(tempx,tempy,c);
                            % caxis([0 .5])
                            uistack(hy_mcgurk_temp,'bottom');
                            hy_mcgurk = [hy_mcgurk;hy_mcgurk_temp];
                        end
    
%                        
%                         ah1=axes('position',get(gca,'position'),'visible','off');
%                         legend(ah1,hy_mcgurk,'McGurk','Location','southeast','Orientation','horizontal')
%                       
                    end
                   
                    
                    hy_interact =[];hy_interact_temp =[];
                    
                    if ~isempty(k_interact)         
                        box_hits_interact = [];
                        counter_interact = 1; % Row
                        box_hits_interact(1,1) = k_interact(1);
                        nexthit_interact = 0;

                        for i=1:length(k_interact)
                            if i<length(k_interact)              
                                % open group if last iteration closed it
                                if nexthit_interact == 1
                                    nexthit_interact = 0;
                                    box_hits_interact(counter_interact,1) = k_interact(i);
                                end

                                % close group
                                if k_interact(i)~=k_interact(i+1)-1
                                    box_hits_interact(counter_interact,2) = k_interact(i);
                                    counter_interact = counter_interact+1;
                                    nexthit_interact = 1;
                                end
                            else
                                if nexthit_interact==0
                                    % close group
                                    box_hits_interact(counter_interact,2) = k_interact(i);
                                else
                                    box_hits_interact(counter_interact,1) = k_interact(i);
                                    box_hits_interact(counter_interact,2) = k_interact(i);
                                end
                            end
                        end
                        temptime = time(time_rngs);
                        box_hits_interact = temptime(box_hits_interact);

                        for i=1:size(box_hits_interact,1)
                            % xstart ystart xlength ylength
                            % hy = rectangle('Position',[box_hits(i,1),tempmin,box_hits(i,2)-box_hits(i,1),abs(tempmin*.5)],'FaceColor',[0.1 .1 .1]);

                            % try with patch
                            %tempx_interact = [box_hits_interact(i,1) box_hits_interact(i,2) box_hits_interact(i,2) box_hits_interact(i,1)];
                            tempx_interact = [box_hits_interact(i,1) box_hits_interact(i,2)];
                  
                            tempy_interact = [tempmin tempmin tempmin+abs(tempmin*.5) tempmin+abs(tempmin*.5)];
                            % hy_interact = patch(tempx_interact,tempy_interact*0.25,'blue');
                            % hy_interact = line(tempx_interact,-0.08,'Color','#0072BD','DisplayName','Interaction');
                            hy_interact_temp = line(tempx_av, [yInteractBar yInteractBar],'Color','#0072BD','Linewidth',5,'DisplayName','Interact');
                        
                            % c = [0 1 1 0];
                            % hy = patch(tempx,tempy,c);
                            % caxis([0 .5])
                            uistack(hy_interact_temp,'bottom');
                            hy_interact = [hy_interact;hy_interact_temp];
                        end
                
                        %ah1=axes('position',get(gca,'position'),'visible','off');        
                        %legend(ah1,hy_interact,'Interaction','Location','southeast','Orientation','horizontal')
                    end
                    
%                     if ~isempty(hy_av) && isempty(hy_mcgurk) && isempty(hy_interact)
%                         legend(hy_av,'AV','Location','southeast','Orientation','horizontal')
%                     elseif isempty(hy_av) && ~isempty(hy_mcgurk) && isempty(hy_interact)
%                         legend('McGurk','Location','southeast','Orientation','horizontal')
%                     elseif isempty(hy_av) && isempty(hy_mcgurk) && ~isempty(hy_interact)
%                         legend('AV McGurk Interact','Location','southeast','Orientation','horizontal')
%                         
%                     elseif ~isempty(hy_av) && ~isempty(hy_mcgurk) && isempty(hy_interact)
%                         legend('AV','McGurk','Location','southeast','Orientation','horizontal')
%                     elseif ~isempty(hy_av) && isempty(hy_mcgurk) && ~isempty(hy_interact)
%                         legend('AV','AV McGurk Interact','Location','southeast','Orientation','horizontal')
%                     elseif isempty(hy_av) && ~isempty(hy_mcgurk) && ~isempty(hy_interact)
%                         legend('McGurk','AV McGurk Interact','Location','southeast','Orientation','horizontal')
%                     elseif ~isempty(hy_av) && ~isempty(hy_mcgurk) && ~isempty(hy_interact)
%                         legend('AV','McGurk','AV McGurk Interact','Location','southeast','Orientation','horizontal')
%                     end
%                     
                    %%%%%%%%%%%%%%%%%% AV McGurk End %%%%%%%%%%%%%%%%%%
                    
             else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%% Plot the statistics %%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%% i.e. btm blackbox of fig %%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % remove middle hits
                if ~isempty(fdrSigSampleIdx)
                    box_hits = [];
                    counterA = 1; % Rows
                    box_hits(1,1) = fdrSigSampleIdx(1);
                    nexthit = 0;
                    for i=1:length(fdrSigSampleIdx)
                        if i<length(fdrSigSampleIdx)              
                            % open group if last iteration closed it
                            if nexthit == 1
                                nexthit = 0;
                                box_hits(counterA,1) = fdrSigSampleIdx(i);
                            end

                            % close group
                            if fdrSigSampleIdx(i)~=fdrSigSampleIdx(i+1)-1
                                box_hits(counterA,2) = fdrSigSampleIdx(i);
                                counterA = counterA+1;
                                nexthit = 1;
                            end
                        else
                            if nexthit==0
                                % close group
                                box_hits(counterA,2) = fdrSigSampleIdx(i);
                            else
                                box_hits(counterA,1) = fdrSigSampleIdx(i);
                                box_hits(counterA,2) = fdrSigSampleIdx(i);
                            end
                        end
                    end
                    temptime = time(time_rngs);
                    box_hits = temptime(box_hits);

                    for i=1:size(box_hits,1)
                        % xstart ystart xlength ylength
                        % hy = rectangle('Position',[box_hits(i,1),tempmin,box_hits(i,2)-box_hits(i,1),abs(tempmin*.5)],'FaceColor',[0.1 .1 .1]);

                        % try with patch
                        tempx = [box_hits(i,1) box_hits(i,2) box_hits(i,2) box_hits(i,1)];
                        tempy = [tempmin tempmin tempmin+abs(tempmin*.5) tempmin+abs(tempmin*.5)];
                        hy = patch(tempx,tempy,[0.1 .1 .1]);
                        % c = [0 1 1 0];
                        % hy = patch(tempx,tempy,c);
                        % caxis([0 .5])
                        uistack(hy,'bottom');
                    end
                end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%% Plot the statistics %%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%% i.e. btm blackbox of fig %%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             end          
             % Add legend outside the loop
            legend(mainLineHandles, legend_array, 'Location', 'Northwest');
            save(['/Users/zhewei/Dropbox (University of Michigan)/_Projects/Brang_LME_AVSpeech/_ScriptsPaper/',Stat_Comp_Label,'_ROI_',num2str(roiIdx),'_fdr_pValue.mat'],'adj_p')
        end
         %%%%%%%%%%%%%%%%%% MAKING THE THREE SUBPLOTS END %%%%%%%%%%%%%%%%%%
         
         ax = gcf;
         exportgraphics(ax,[Stat_Comp_Label '_ThetaPower.eps'],'ContentType','vector')
         
        % add main title OVER the main plots
%         legend(figHolder,legend_array,'Location','Northwest');
%         sgtitle(strrep([Stat_Comp_Label,' ',Analysis_Title],'_',' '));
        % Save Image
        plot_title = strrep([Stat_Comp_Label,' ',Analysis_Title],'_',' ');
        % suptitle(plot_title)
        print(f1,'-depsc','-painters',[outputDir,'/images/',Stat_Comp_Label,'_',Analysis_Title,'_TimeSeries.eps'])

    end
end

% end
toc




