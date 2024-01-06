function aa_EpisodeCodingAtTaskSwitches_220622()
% Automatic Analysis User master script for DMN Context-decoding During Task Switch
% fMRI experiment, by Ashley Zhou and Daniel Mitchell, Nov 2022
%
% Basd on example aa_user_fmri (aa version 5.*.*)

%% INITIALISE
dbstop if error
addpath /imaging/local/software/AA/release-5.4.0_202008 % (last supported version at CBU)
%addpath /imaging/local/software/spm_cbu_svn/releases/spm12_latest
aa_ver5

here=fileparts(which(mfilename));
addpath(here)

%% DEFINE SPECIFIC PARAMETERS
%  Create recipe, load defaults
aap=aarecipe('aap_parameters_defaults_CBSU.xml','aap_tasklist_fmri_ECATS_220622.xml');
%{
Example tasklist contained:
aamod_autoidentifyseries_timtrio
aamod_get_dicom_structural
aamod_get_dicom_epi
aamod_get_dicom_fieldmap
aamod_convert_structural
aamod_convert_epis
aamod_convert_fieldmaps
aamod_fieldmap2VDM
aamod_tsdiffana
aamod_realignunwarp
aamod_tsdiffana
aamod_slicetiming
aamod_tsdiffana
aamod_coreg_extended_1
aamod_segment8
aamod_coreg_extended_2epi
aamod_coreg_extended_2meanepi
aamod_norm_write
aamod_norm_write_meanepi
aamod_smooth

aamod_firstlevel_model
aamod_firstlevel_contrasts
aamod_firstlevel_threshold
aamod_secondlevel_model
aamod_secondlevel_contrasts
aamod_secondlevel_threshold
%}

% ensure SPM is on the path
%aap.directory_conventions.toolbox(1).dir= '/imaging/local/software/spm_cbu_svn/releases/spm12_latest';
aap.directory_conventions.toolbox(1).dir= '/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7771';
spmdir = aap.directory_conventions.toolbox(strcmp({aap.directory_conventions.toolbox.name},'spm')).dir;
spmhit = which('spm_spm');
%remove the last '/'
if any(spmhit)
    %spmdir='/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7771';
    assert(strcmp(fileparts(spmhit), spmdir), 'spm on path differs from that in aap.directory_conventions');
else
    fprintf('adding spmdir to path: %s\n', spmdir);
    SPMtool = spmClass(spmdir);
    SPMtool.load;
end

% Modify standard recipe module selection here
aap.acq_details.root = here; % project folder
aap.acq_details.numdummies=8; % about 10s worth

aap.options.wheretoprocess = 'localsingle'; % queuing system			% OPTIONS: 'localsingle'|'qsub' for aa engine, typical value 'qsub'
aap.options.autoidentifyfieldmaps=true;  							% typical value 1
aap.options.autoidentifystructural_chooselast = 1; aap.options.autoidentifystructural_average=0;
%aap.options.autoidentifystructural_chooselast = 0; aap.options.autoidentifystructural_average=1;
aap.options.NIFTI4D = 1;										% typical value 1
aap.options.email='ashley.zhou@mrc-cbu.cam.ac.uk';
%aap.options.garbagecollection=0; % causing some odd error?
aap.options.aaworkerroot='/group/duncan-lab/users/az01'; % otherwise defaults to home space

% may need to change these to autoidentify series:
aap.directory_conventions.protocol_structural='T1w_MPR'; % Moataz's HCP T1; note, prescan normalize should have been enabled, and both images may be saved, so use 2nd or average
%aap.directory_conventions.protocol_fieldmap='FieldMap_AP'; % Moataz's HCP fieldmap; there's also a PA one; I can't manage to process either with aa/spm
aap.directory_conventions.protocol_fieldmap='fieldmap_gre'; % regular fieldmap

aap.directory_conventions.analysisid = 'aa5_analysis_220622'; % analysis folder
aap.directory_conventions.subject_directory_format=3; % 1=based on CBUID; 2=subjectnumber; 3=manual
aap.directory_conventions.subjectoutputformat='%s';
aap.directory_conventions.seriesoutputformat='Series%03d*'; % what it will look for in rawdata

%% Add data

% sessions
aap = aas_addsession(aap,'Run1');
aap = aas_addsession(aap,'Run2');
aap = aas_addsession(aap,'Run3');
aap = aas_addsession(aap,'Run4');
aap.acq_details.selected_sessions=1:4;

% subjects

%Pilots
%aap = aas_addsubject(aap,'CBU220400',mri_findvol(aap,'CBU220400_MR22003'),'functional',{12,14,16,18});
% pilot; scenes incorrect
%aap = aas_addsubject(aap,'CBU220408',mri_findvol(aap,'CBU220408_MR22003'),'functional',{18,20}); % pilot
%aap = aas_addsubject(aap,'CBU220425',mri_findvol(aap,'CBU220425_MR22003'),'functional',{12,14,18,20});
% look for the indexes for the functionals
%
%aap = aas_addsubject(aap,'CBU220428',mri_findvol(aap,'CBU220428_MR22003'),'functional',{12,14,18,20});
%aap = aas_addsubject(aap,'CBU220428',mri_findvol(aap,'CBU220428_MR22003'),'functional',{12,14,18,20});
%aap = aas_addsubject(aap,'CBU220469',mri_findvol(aap,'CBU220469_MR22003'),'functional',{10,12,16,18});

%Participants
aap = aas_addsubject(aap,'CBU220465',mri_findvol(aap,'CBU220465_MR22003'),'functional',{12,18,22,24});  %behavfile = 5
aap = aas_addsubject(aap,'CBU220469',mri_findvol(aap,'CBU220469_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220471',mri_findvol(aap,'CBU220471_MR22003'),'functional',{8,10,14,16}); %Charlotte's, has T1w and Tw2 from other folder copied over
aap = aas_addsubject(aap,'CBU220473',mri_findvol(aap,'CBU220473_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220474',mri_findvol(aap,'CBU220474_MR22003'),'functional',{8,10,14,16});
aap = aas_addsubject(aap,'CBU220475',mri_findvol(aap,'CBU220475_MR22003'),'functional',{10,12,16,18});
%aap =aas_addsubject(aap,'CBU220479',mri_findvol(aap,'CBU220479_MR22003'),'functional',{10,12,16,18});
%subject did not complete 2nd run
aap = aas_addsubject(aap,'CBU220480',mri_findvol(aap,'CBU220480_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220481',mri_findvol(aap,'CBU220481_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220482',mri_findvol(aap,'CBU220482_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220485',mri_findvol(aap,'CBU220485_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220491',mri_findvol(aap,'CBU220491_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220492',mri_findvol(aap,'CBU220492_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220493',mri_findvol(aap,'CBU220493_MR22003'),'functional',{10,12,20,22});
aap = aas_addsubject(aap,'CBU220494',mri_findvol(aap,'CBU220494_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220495',mri_findvol(aap,'CBU220495_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220496',mri_findvol(aap,'CBU220496_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220499',mri_findvol(aap,'CBU220499_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220503',mri_findvol(aap,'CBU220503_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220506',mri_findvol(aap,'CBU220506_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220509',mri_findvol(aap,'CBU220509_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220510',mri_findvol(aap,'CBU220510_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220512',mri_findvol(aap,'CBU220512_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220513',mri_findvol(aap,'CBU220513_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220514',mri_findvol(aap,'CBU220514_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220519',mri_findvol(aap,'CBU220519_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220523',mri_findvol(aap,'CBU220523_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220524',mri_findvol(aap,'CBU220524_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220525',mri_findvol(aap,'CBU220525_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220526',mri_findvol(aap,'CBU220526_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220533',mri_findvol(aap,'CBU220533_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220535',mri_findvol(aap,'CBU220535_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220536',mri_findvol(aap,'CBU220536_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220542',mri_findvol(aap,'CBU220542_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220538',mri_findvol(aap,'CBU220538_MR22003'),'functional',{10,12,16,18});
aap = aas_addsubject(aap,'CBU220539',mri_findvol(aap,'CBU220539_MR22003'),'functional',{10,12,16,18});
%6
%% task settings
aap.tasksettings.aamod_slicetiming.autodetectSO = 1; % automatically detect slice order
aap.tasksettings.aamod_slicetiming.refslice = 1; % with autodetectSO, this should be (index of?) time (ms) of reference slice?

aap.tasksettings.aamod_norm_write.vox = [2 2 2];
aap.tasksettings.aamod_norm_write_meanepi.vox = [2 2 2];

aap.tasksettings.aamod_smooth.FWHM = 10;

for model=1:12
aap.tasksettings.aamod_firstlevel_model(model).xBF.name = 'hrf';       %'hrf (with time and dispersion derivatives)';
aap.tasksettings.aamod_firstlevel_model(model).xBF.UNITS = 'secs';    	% OPTIONS: 'scans'|'secs' for onsets and durations, typical value 'secs'
aap.tasksettings.aamod_firstlevel_model(model).xBF.T0=[]; % empty to default to refslice, and # microtime bins per TR = # slices (SPM says must do this if have done slicetiming?)
aap.tasksettings.aamod_firstlevel_model(model).includemovementpars = 1;% Include/exclude Moco params
aap.tasksettings.aamod_firstlevel_model(model).TR=[];
% aap.tasksettings.aamod_firstlevel_threshold(model).threshold.correction='none';
% aap.tasksettings.aamod_firstlevel_threshold(model).threshold.p=0.025;
end
aap.tasksettings.aamod_firstlevel_threshold(5).threshold.correction='FDR';
aap.tasksettings.aamod_firstlevel_threshold(5).threshold.p=0.05;

aap.tasksettings.aamod_firstlevel_threshold(6).threshold.correction='FDR';
aap.tasksettings.aamod_firstlevel_threshold(6).threshold.p=0.05;

%setting decoding model's input directories and regressor names
decoding_dirs = {'DecodeContext','DecodeNovelty','DecodeSide','DecodeAnswer','DecodeSideSwitch_dur','DecodeSideSwitch_delta'};
model_dirs = [5:8,10,11];
regressor_names={{'.*_beach','.*_cafe'},{'.*_novel','.*_repeat'},{'.*_left','.*_right'},{'.*_yes','.*_no'},{'.*_same','.*_different'},{'.*_same','.*_different'}};

for index=1:6
    decodingdir=sprintf(decoding_dirs{index});
    reg = regressor_names{index};
    aap.tasksettings.aamod_decoding4_az(1).outdir=decodingdir;
    aap.tasksettings.aamod_decoding4_az(1).modeldir=['aamod_firstlevel_model_0000' num2str(model_dirs(index))];
    aap.tasksettings.aamod_decoding4_az(1).outfileprefix=decodingdir;
    aap.tasksettings.aamod_decoding4_az(1).doneflagsuffix=decodingdir;
    aap.tasksettings.aamod_decoding4_az(1).regressorfilt=reg;
    aap.tasksettings.aamod_decoding4_az(1).labelfilt=reg(1);
end

[aap.tasksettings.aamod_decoding4_az([1:4,11:12]).testsets]=deal({'task_stay','within_domain','between_domain',' rest_','restart','extended_rest'}); % to train across groups, but test separately on these groups (should be non overlapping)
[aap.tasksettings.aamod_decoding4_az(:).testtype]=deal('LxROCV'); % 'LxROCV', (don't use 'L1ROCV' - this isn't set up properly for splitting test sets)
[aap.tasksettings.aamod_decoding4_az(:).decodingsoftware]=deal('liblinear');
[aap.tasksettings.aamod_decoding4_az(:).scalemethod]=deal('none'); %deal('none');%
[aap.tasksettings.aamod_decoding4_az(:).scaleestimation]=deal('none'); % none, across, all, all_percondition, perstep
[aap.tasksettings.aamod_decoding4_az(:).masks]=deal(fullfile(aap.acq_details.root,'DMN_visual_motor_MD_ROIs')); % this will do an ROI analysis on ALL image files in the directory; if empty, will use a searchlight.
[aap.tasksettings.aamod_decoding4_az(:).subsampletrainitems]=deal(Inf);
[aap.tasksettings.aamod_decoding4_az(:).subsampleclassifiers]=deal(Inf);
%%%%%%%

%%  Define model

% Obtain TR from the first session
h = spm_dicom_headers(mri_finddcm(aap, 'CBU220400_MR22003',12)); % should be 1.2080
TR = h{1}.RepetitionTime/1000; % in seconds
n=2;

%retrieves result data from each run of each subject 
for sub = {aap.acq_details.subjects.subjname}
    for sess=aap.acq_details.selected_sessions
        subname = sub{1};
        subname = subname(4:end);
        
        %This participant had data file recorded in atypical format
        if strcmp(subname,'220465')
            subname = '5';
        end
        load(fullfile(fileparts(aap.acq_details.root),'scripts',['exp_pilot_' subname '_run_' num2str(sess) '.mat']))
        
        %exp_starttime = result(1).frame_onset - 1 - n*1.208;
        %%%%
        exp_starttime = result(1).exp_start_time + aap.acq_details.numdummies*TR;
        %%%%
        
        %Model is defined in the order of context, novelty, side, answer,
        %and taskside switch
        scene = {'beach', 'cafe'};
        novelty = {'novel', 'repeat'};
        side = {'left', 'right'};
        answer = {'yes','no'};
        switch_side={'same','different'};
        cs_dt={'context_switch','dummy_trial'};
        
        %% Defining models
        for index = 1:2
            for cond = {'task_stay','within_domain','between_domain','rest','restart','extended_rest'}
                %rts.([cond '_' scene])=[];
                mod1_onset.([cond{1} '_' scene{index}])=[];
                mod1_dur.([cond{1} '_' scene{index}])=[];
                
                mod2_onset.([cond{1} '_' novelty{index}])=[];
                mod2_dur.([cond{1} '_' novelty{index}])=[];
                
                mod3_onset.([cond{1} '_' side{index}])=[];
                mod3_dur.([cond{1} '_' side{index}])=[];
                
                mod4_onset.([cond{1} '_' answer{index}])=[];
                mod4_dur.([cond{1} '_' answer{index}])=[];
                
                mod_switchside_onset.([cond{1} '_' switch_side{index}])=[];
                mod_switchside_dur.([cond{1} '_' switch_side{index}])=[];
            end
            
            for cs= cs_dt
                mod_switchside_onset.([cs{1} '_' switch_side{index}])=[];
                mod_switchside_dur.([cs{1} '_' switch_side{index}])=[];
            end
        end
        mod_onset.dummy_trial = [];
        mod_onset.context_switch = [];
        
        mod_dur.dummy_trial= [];
        mod_dur.context_switch  = [];
        
        task_types = {'a1','a2','b1','b2'};
        for type = task_types
            mod_task_onset.(type{1})=[];
            mod_task_dur.(type{1})=[];
        end
        
        for first_trial = {'a1','a2','b1','b2','r'}
            for second_trial = {'a1','a2','b1','b2','r'}
                reg_name = [first_trial{1} '_' second_trial{1}];
                model_switch_by_domain_onset.(reg_name) = [];
                model_switch_by_domain_dur.(reg_name) = [];
                
            end
        end
            %initialise the struct with each switch (task pair)
            %struct will have fields a1_a1, a2_a1, b1_a1....
            
            
        %calculating outlier threshold for reaction times based on third
        %derivative past mean
        ok=~cellfun(@isempty,{result.rt}); rts=[result(ok).rt]-[result(ok).frame_onset]; slowoutlier=mean(rts(~isnan(rts)))+3*std(rts(~isnan(rts)));
        ok=~cellfun(@isempty,{result.rt}); rts=[result(ok).rt]-[result(ok).context_onset]; context_slowoutlier=mean(rts(~isnan(rts)))+3*std(rts(~isnan(rts)));
        ok=~cellfun(@isempty,{result.rt}); rts=[result(ok).rt]-[result(ok).stim_onset]; stim_slowoutlier=mean(rts(~isnan(rts)))+3*std(rts(~isnan(rts)));
        
        for trial = 1:length(result)
            cond = result(trial).switch_type;
            scene = result(trial).scene;

            if ~isempty(result(trial).rt) && ~isnan(result(trial).rt)
                rt = result(trial).rt - result(trial).frame_onset;  %rt is calculated from frame onset to response made
                context_rt = result(trial).rt - result(trial).context_onset;  %context_rt is calculated from context/background onset to response made
                stim_rt = result(trial).rt - result(trial).stim_onset;  %rt is calculated from stimulus onset to response made
            else
                %if rt is empty in result data, ie. if there wasn't a response, take diff between next trial fixation
                %If it's the last trial, use average response time in rest
                if trial == length(result)
                    t = find(strcmp({result.switch_type},'rest'));
                    rt = result(t(end-1)+1).frame_onset - 1.5 - result(t(end-1)).frame_onset;
                    context_rt = result(t(end-1)+1).frame_onset - 1.5 - result(t(end-1)).context_onset;
                    stim_rt = result(t(end-1)+1).frame_onset - 1.5 - result(t(end-1)).stim_onset;
                else
                    rt = result(trial+1).frame_onset - 1.5 - result(trial).frame_onset;
                    context_rt = result(trial+1).frame_onset - 1.5 - result(trial).context_onset;
                    stim_rt = result(trial+1).frame_onset - 1.5 - result(trial).stim_onset;
                end
            end
            assert(~isnan(rt));
            assert(~isnan(context_rt));
            assert(~isnan(stim_rt));
            
            %Take the smaller of the third der outlier or rt
            rt =min(rt,slowoutlier);
            context_rt=min(context_rt,context_slowoutlier);
            stim_rt=min(stim_rt,stim_slowoutlier);
            
            %onset is calculated from frame showing up
            onset = result(trial).frame_onset-exp_starttime;
            
            
            
            if strcmp(cond,'dummy_trial') ||strcmp(cond,'context_switch')
                
                %record onsets and durations for dummy trials and context
                %switch (jungle detection) trials separately
                mod_onset.(cond) = [mod_onset.(cond) onset];
                mod_dur.(cond) = [mod_dur.(cond) rt];
                
                %store the duration and onset of sideswitch trials,
                %including dummy trials and context switch trials
                if trial==1
                    continue;
                elseif result(trial).side==result(trial-1).side
                    mod_switchside_dur.([cond '_same']) = [mod_switchside_dur.([cond '_same']) rt];
                    mod_switchside_onset.([cond '_same']) = [mod_switchside_onset.([cond '_same']) onset];
                else
                    
                    mod_switchside_dur.([cond '_different']) = [mod_switchside_dur.([cond '_different']) rt];
                    mod_switchside_onset.([cond '_different']) = [mod_switchside_onset.([cond '_different']) onset];
                end
            else
                pre_type = result(trial-1).type;
                this_type = result(trial).type;
                
                if contains(pre_type,'r')
                    pre_type = 'r';
                end
                if contains(this_type,'r')
                    this_type='r';
                end
                this_trial_switch = [pre_type '_' this_type];
               
                model_switch_by_domain_onset.(this_trial_switch) = [model_switch_by_domain_onset.(this_trial_switch) onset];
                model_switch_by_domain_dur.(this_trial_switch) = [model_switch_by_domain_dur.(this_trial_switch) stim_rt];

                
                %store duration and onsets by category
                
                %beach vs cafe
                mod1_dur.([cond '_' scene]) = [mod1_dur.([cond '_' scene]) context_rt];
                mod1_onset.([cond '_' scene]) = [mod1_onset.([cond '_' scene]) onset];
                
                %task type
                if ~contains(result(trial).type,'r')
                    mod_task_onset.(result(trial).type) = [mod_task_onset.(result(trial).type) onset];
                    mod_task_dur.(result(trial).type) = [mod_task_dur.(result(trial).type) rt];
                end
                
                %novel vs repeat
                if ~isempty(result(trial).scene_novel)
                    mod2_dur.([cond '_novel']) = [mod2_dur.([cond '_novel']) context_rt];
                    mod2_onset.([cond '_novel']) = [mod2_onset.([cond '_novel']) onset];
                else
                    mod2_dur.([cond '_repeat']) = [mod2_dur.([cond '_repeat']) context_rt];
                    mod2_onset.([cond '_repeat']) = [mod2_onset.([cond '_repeat']) onset];
                end
                
                %left vs right
                if result(trial).side
                    mod3_onset.([cond '_right'])=[mod3_onset.([cond '_right']) onset];
                    mod3_dur.([cond '_right'])=[mod3_dur.([cond '_right']) rt];
                else
                    mod3_onset.([cond '_left'])=[mod3_onset.([cond '_left']) onset];
                    mod3_dur.([cond '_left'])=[mod3_dur.([cond '_left']) rt];
                end
                
                %answer yes vs no (excluding rest trials)
                if ~isempty(result(trial).response) && ~(strcmp(result(trial).response,'context'))
                    mod4_onset.([cond '_' result(trial).response]) = [mod4_onset.([cond '_' result(trial).response]) onset];
                    mod4_dur.([cond '_' result(trial).response]) = [mod4_dur.([cond '_' result(trial).response]) stim_rt];
                end
                
                %same vs diff sides from previous trial
                if trial==1
                    continue;
                elseif result(trial).side==result(trial-1).side
                    mod_switchside_dur.([cond '_same']) = [mod_switchside_dur.([cond '_same']) rt];
                    mod_switchside_onset.([cond '_same']) = [mod_switchside_onset.([cond '_same']) onset];
                else
                    
                    mod_switchside_dur.([cond '_different']) = [mod_switchside_dur.([cond '_different']) rt];
                    mod_switchside_onset.([cond '_different']) = [mod_switchside_onset.([cond '_different']) onset];
                end
                
            end
        end
        
        
        %there is 13 models already (~?)
        
        %for model 14, add onsets and durs
        for name =fieldnames(model_switch_by_domain_dur)'
            aap=aas_addevent(aap,'aamod_firstlevel_model_00012',sub{1},aap.acq_details.sessions(sess).name,...
                name{1},model_switch_by_domain_onset.(name{1}),model_switch_by_domain_dur.(name{1}));
           
        end
        
        
        %add models to aa pipeline with onsets and durations of trials
        model_durs = {mod1_dur,mod2_dur,mod3_dur,mod4_dur,mod_task_dur,mod_switchside_dur};
        model_onsets ={mod1_onset,mod2_onset,mod3_onset,mod4_onset,mod_task_onset,mod_switchside_onset};
        
        %two models for each (smoothed and unsmoothed)
        model_dirs = {[1,5],[2,6],[3,7],[4,8],9,[10,11]};
        
        
        %make model where prev task and current task are regressors
        
        
        for models=1:6
            model_dir = model_dirs{models};
            model_dur = model_durs{models};
            model_onset = model_onsets{models};
            
            %retrieve name for each model
            for model=model_dir
                model_name = 'aamod_firstlevel_model';
                if model >1
                    if model<10
                        model_name =[model_name '_0000' num2str(model)];
                    else
                        model_name =[model_name '_000' num2str(model)];
                    end
                end
                
                %For task type model, store as real task names instead
                if model ==9
                    task_names_full={'living','shoebox','letter_A','letter_I'};
                    ind =1;
                    for name =fieldnames(mod_task_dur)'
                        aap=aas_addevent(aap,'aamod_firstlevel_model_00009',sub{1},aap.acq_details.sessions(sess).name,...
                            task_names_full{ind},model_onset.(name{1}),model_dur.(name{1}));
                        ind=ind+1;
                    end
                    
                %for second sideswitch model use delta function instead
                elseif model == 11
                    for name =fieldnames(model_dur)'
                        aap=aas_addevent(aap,model_name,sub{1},aap.acq_details.sessions(sess).name,...
                            name{1},model_onset.(name{1}),0);
                    end
                else
                    for name =fieldnames(model_dur)'
                        aap=aas_addevent(aap,model_name,sub{1},aap.acq_details.sessions(sess).name,...
                            name{1},model_onset.(name{1}),model_dur.(name{1}));
                    end
                end
                
                %for all models also add dummy trials and context switch
                %trials as regressors
                    
                aap=aas_addevent(aap,model_name,sub{1},aap.acq_details.sessions(sess).name,...
                    'dummy_trial',mod_onset.dummy_trial,mod_dur.dummy_trial);
                aap=aas_addevent(aap,model_name,sub{1},aap.acq_details.sessions(sess).name,...
                    'context_switch',mod_onset.context_switch,mod_dur.context_switch);
            end
        end
    end
    
end
% end

%% Define contrasts

%contrast model for all the switch task pairs 
doms = {'a','b'};
for f = 1:2
    
    ts_con=[];
    wd_con=[];
    bd_con=[];
    restart_con=[];
    rest_con=[];
    for s =1:2
        
        switchname = [doms{f} num2str(s) '_' doms{f} num2str(s)];
        ts_con = [ts_con sprintf('+2x%s|',upper(switchname))];
        
        wdswitchname = [doms{f} num2str(s) '_' doms{f} num2str(3-s)];
        wd_con = [wd_con sprintf('+2x%s|',upper(wdswitchname))];
        
        bdswitchname1 = [doms{f} num2str(s) '_' doms{3-f} num2str(3-s)];
        bdswitchname2 = [doms{f} num2str(s) '_' doms{3-f} num2str(s)];
        bd_con=[bd_con sprintf('+1x%s|',upper(bdswitchname1)) sprintf('+1x%s|',upper(bdswitchname2))];

        restart = ['r_' doms{f} num2str(s)];
        restart_con =[restart_con sprintf('+2x%s|',upper(restart))];
        
        rest = [doms{f} num2str(s) '_r' ];
        rest_con =[restart_con sprintf('+2x%s|',upper(rest))];
    end
    
    if f ==1
            %a1-a1,a2-a2, sem ts   
            aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00009','*','sameforallsessions',ts_con(1:length(ts_con)-1),'sem_ts','T');
            
            %a1-a2,a2,a1,sem wd
            aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00011','*','sameforallsessions',wd_con(1:length(wd_con)-1),'sem_wd','T');
            
            %a1-b1,a1-b2, 
            %a2-b1,a2-b2
            aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00013','*','sameforallsessions',bd_con(1:length(bd_con)-1),'sem_lex_bd','T');
            
            %r-a1,r-a2
            aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00015','*','sameforallsessions',restart_con(1:length(restart_con)-1),'sem_restart','T');

            %a1-r,a2-r
            aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00017','*','sameforallsessions',rest_con(1:length(rest_con)-1),'sem_rest','T');

            
        else
            %b1-b1,b2-b2, lex ts
            aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00010','*','sameforallsessions',ts_con(1:length(ts_con)-1),'lex_ts','T');
            
            %b1-b2,b2-b1,lex wd
            aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00012','*','sameforallsessions',wd_con(1:length(wd_con)-1),'lex_wd','T');
            aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00014','*','sameforallsessions',bd_con(1:length(bd_con)-1),'lex_sem_bd','T');
    
            %r-b1,r-b2
            aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00016','*','sameforallsessions',restart_con(1:length(restart_con)-1),'lex_restart','T');
            %b1-r,b2-r
            aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00018','*','sameforallsessions',rest_con(1:length(rest_con)-1),'lex_rest','T');


    end
        
end

rest_con=sprintf('+2x%s|',upper(restart));
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00019','*','sameforallsessions',rest_con(1:length(rest_con)-1),'rest_rest','T');



%For each model the 12 regressors are stored by 6 conditions x 2 groups(ex. scene x conditions, dummy, context switch)
%Ex. contrast for beach (+1) vs cafe (-1)

regressor_names={{'beach','cafe'},{'novel','repeat'},{'left','right'},{'yes','no'},{'same','different'},{'same','different'}};
con_names = {fieldnames(mod1_dur)',fieldnames(mod2_dur)',fieldnames(mod3_dur)',fieldnames(mod4_dur)',fieldnames(mod_switchside_dur)',fieldnames(mod_switchside_dur)'};

%Assign +1 for all conditions in group A (ex. beach) and -1 for all
%conditions in group B (ex. cafe)
%Flip for opposite contrast

for contrast = 1:6
    con=[];
    opp_con=[];
    names=con_names{contrast};
    for ind = 1:12
        if ind<7
            if contrast ==4 %yes/no contrasts we exclude rest conditions (4,6,10,12)
                if ind ==4 || ind ==6
                    continue;
                end
            end
            con = [con sprintf('+1x%s|',upper(names{ind}))];
            opp_con = [opp_con sprintf('-1x%s|',upper(names{ind}))];
        else
            if contrast ==4 %yes/no contrasts we exclude rest conditions (4,6,10,12)
                if ind ==10 || ind ==12
                    continue;
                end
            end
            con = [con sprintf('-1x%s|',upper(names{ind}))];
            opp_con = [opp_con sprintf('+1x%s|',upper(names{ind}))];
        end
    end
    con = con(1:length(con)-1);
    opp_con = opp_con(1:length(opp_con)-1);
    
    contrast_name = 'aamod_firstlevel_contrasts';
    contrast_num = contrast;
    %switchside contrasts are number 7,8
    if contrast_num>4
        contrast_num=contrast_num+2;
    end
    if contrast >1
        contrast_name = [contrast_name '_0000' num2str(contrast_num)];
    end
    condition_names = regressor_names{contrast};
    aap = aas_addcontrast(aap,contrast_name,'*','sameforallsessions',con,[condition_names{1} '_' condition_names{2} '_conditions'],'T');
    aap = aas_addcontrast(aap,contrast_name,'*','sameforallsessions',opp_con,[condition_names{2} '_' condition_names{1} '_conditions'],'T');
end

%for task type contrasts, each combination-pair is contrasted
con_names = {'living','shoebox','letter_A'};
opp_con_names = {{'shoebox','letters_A','letter_I'},{'letters_A','letter_I'},{'letter_I'}};
decode_num=5;

for num = 1:3
   
    for oppcon_num = opp_con_names{num}
        task_con = [sprintf('+1x%s|',upper(con_names{num})) sprintf('-1x%s',upper(oppcon_num{1}))];
        
        %Each combination-pair is also defined as regressors in decoding
        %model for a total of 6 decoding models
        decodingdir=sprintf([con_names{num} '_vs_' oppcon_num{1}]); %
        aap.tasksettings.aamod_decoding4_az(decode_num).modeldir='aamod_firstlevel_model_00009';
        aap.tasksettings.aamod_decoding4_az(decode_num).outdir=decodingdir;
        aap.tasksettings.aamod_decoding4_az(decode_num).outfileprefix=decodingdir;
        aap.tasksettings.aamod_decoding4_az(decode_num).doneflagsuffix=decodingdir;
        aap.tasksettings.aamod_decoding4_az(decode_num).regressorfilt={con_names{num},oppcon_num{1}};
        aap.tasksettings.aamod_decoding4_az(decode_num).labelfilt=con_names(num);

        aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00006','*','sameforallsessions',task_con,[con_names{num} '_vs_' oppcon_num{1}],'T');
        decode_num=decode_num+1;
        
    end
end

%add contrasts for conditions, only did this for side and switchside models
rest_con=[];
for cond = {'within_domain','between_domain','rest','restart','extended_rest','task_stay'}
    
    %added another contrast to compare rest vs task
    if strcmp(cond,'rest') || strcmp(cond,'extended_rest')
        rest_con = [rest_con sprintf('+2x%s|',upper([cond{1} '_same']))];
        rest_con = [rest_con sprintf('+2x%s|',upper([cond{1} '_different']))];
    else
        rest_con = [rest_con sprintf('-1x%s|',upper([cond{1} '_same']))];
        rest_con = [rest_con sprintf('-1x%s|',upper([cond{1} '_different']))];
    end
    
    for cond2 = {'task_stay','within_domain','between_domain','rest','restart','extended_rest'}
        %for each combination pair of taskswitch conditions
        con=[];  %contrast for just comparing each condition
        delta_con=[];  %contrast for comparing each condition in delta model
        same_switch_con=[];  %contrast for comparing all conditions in the same side
        diff_switch_con=[];  %contrast for comparing all conditions switching to a different side
        
        if strcmp(cond{1},cond2{1})
            continue;
        end
        
        pos_index = find(contains(side_names,[cond{1} '_'])-contains(side_names,['_' cond{1}]));
        neg_index = find(contains(side_names,[cond2{1} '_'])-contains(side_names,['_' cond2{1}]));
        
        dpos_index = find(contains(switchside_names,[cond{1} '_'])-contains(switchside_names,['_' cond{1}]));
        dneg_index = find(contains(switchside_names,[cond2{1} '_'])-contains(switchside_names,['_' cond2{1}]));
        
        same_switch_con = [sprintf('+1x%s|',upper([cond{1} '_same'])) sprintf('-1x%s|',upper([cond2{1} '_same']))];
        diff_switch_con = [sprintf('+1x%s|',upper([cond{1} '_different'])) sprintf('-1x%s|',upper([cond2{1} '_different']))];
        
        for pos =1:length(pos_index)
            con = [con sprintf('+1x%s|',upper(side_names{pos_index(pos)}))];
            delta_con = [delta_con sprintf('+1x%s|',upper(switchside_names{dpos_index(pos)}))];
        end
        for neg =1:length(neg_index)
            con = [con sprintf('-1x%s|',upper(side_names{neg_index(neg)}))];
            delta_con = [delta_con sprintf('-1x%s|',upper(switchside_names{dneg_index(neg)}))];
        end
        con = con(1:length(con)-1);
        delta_con = delta_con(1:length(delta_con)-1);
        aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00003','*','sameforallsessions',con,[cond{1} '_minus_' cond2{1}],'T');
        %unsmoothed conditions comparison
        aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00005','*','sameforallsessions',con,[cond{1} '_minus_' cond2{1}],'T');
        aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00007','*','sameforallsessions',same_switch_con,['same_' cond{1} '_minus_' cond2{1}],'T');
        aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00007','*','sameforallsessions',diff_switch_con,['diff_' cond{1} '_minus_' cond2{1}],'T');
        aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00007','*','sameforallsessions',delta_con,'delta_conditions','T');
        aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00008','*','sameforallsessions',same_switch_con,['same_' cond{1} '_minus_' cond2{1}],'T');
        aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00008','*','sameforallsessions',diff_switch_con,['diff_' cond{1} '_minus_' cond2{1}],'T');
    end
end

%contrast all dummy trials with task-stay trials
dummy_con = '+1xDUMMY_TRIAL_SAME|+1xDUMMY_TRIAL_DIFFERENT|-1xTASK_STAY_SAME|-1xTASK_STAY_DIFFERENT';
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00007','*','sameforallsessions',dummy_con,'dummy_vs_task_stay','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00008','*','sameforallsessions',dummy_con,'dummy_vs_task_stay','T');

%contrast rest with task, looks like this:
%rest_task_con = '+1xrest_same|+1xextended_rest_same|+1xrest_different|+1xextended_rest_different|-1xtask_stay_same|-1xtask_stay_different|-1xwithin_domain_same|-1xwithin_domain_different|-1xbetween_domain_same|-1xbetween_domain_different|-1xrestart_same|-1xrestart_different';
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00007','*','sameforallsessions',rest_con,'rest_task_con','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00008','*','sameforallsessions',rest_con,'rest_task_con','T');


%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));

return;