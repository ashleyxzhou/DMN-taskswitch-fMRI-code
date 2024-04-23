function data_analysis()
%Script for generating results and figures for manuscript 
%This script should be in the same folder as all data files from the Github
%in order to be run:
%exp_data folder (contains behav results), pilotsub_mem_task folder
%(contains memory task results), and two mat files for univariate and
%multivariate data
%Note that certain functions included at end of script are from toolboxes
%that have copyrights (disclaimers are included)
%v0904

dbstop if error

%% Behavioural results
behav_data.RT_regular=[];
behav_data.RT_cswitch=[];
behav_data.accuracy_regular=[];
behav_data.accuracy_cswitch=[];

subs= {220513,220525,5,220469,220471,220473,220474,220475,220480,220481,220482,220485,220491,220492,220493,220494,220495,220496,220499,220503,220506,220509,220510,220512,220514,220519,220523,220524,220526,220533,220535,220536,220538,220539,220542};

nsub = length(subs);

first_half_sub_avg_rt = nan(1,nsub);
second_half_sub_avg_rt = nan(1,nsub);
allsub_accuracy=[];
allsub_context_accuracy=[];
allsub_context_rt=[];
allsub_rt=[];
dprime_allsubs = nan(length(subs),5);
se_allsubs = nan(length(subs),6);
subnum=0;
sub_mean.regular_RT = nan(length(subs),1);
sub_mean.regular_accuracy = nan(length(subs),1);
sub_mean.cswitch_RT = nan(length(subs),1);
sub_mean.cswitch_accuracy = nan(length(subs),1);
alloutlier_num=[];
all_stds=[];
all_means=[];
all_outliers=[];
all_rest=[];
sub_dir=fullfile(pwd,'exp_data'); % extract the participant's data from zip file if downloaded from Github
mem_dir=fullfile(pwd,'pilotsub_mem_task');
for sub = subs
    
    subnum=subnum+1;
    conditions.within_domain=[];
    conditions.between_domain=[];
    conditions.task_stay=[];
    conditions.rest=[];
    conditions.restart=[];
    conditions.extended_rest=[];
    dat.RT=struct('regular',conditions,'cswitch',conditions);
    dat.accuracy=struct('regular',conditions,'cswitch',conditions);
    
    first_half=[];
    second_half=[];
    
    %calculate mean +3std of all trials from subject as threshold for slow
    %outliers
    sub_outliers=[];
    sub_stds=[];
    submeans=[];
    sub_outlier=0;
    sub_rest=[];
    for run = 1:4
        load(fullfile(sub_dir,['exp_pilot_' num2str(sub{1}) '_run_' num2str(run) '.mat']),'result');
        ok=~cellfun(@isempty,{result.rt}); rts=[result(ok).rt]-[result(ok).stim_onset]; slowoutlier=mean(rts(~isnan(rts)))+3*std(rts(~isnan(rts)));
        sub_stds = [sub_stds std(rts(~isnan(rts)))];
        submeans = [submeans mean(rts(~isnan(rts)))];
        sub_outliers=[sub_outliers slowoutlier];
        %rest
        restind=find(strcmp({result.switch_type},'rest'));
        sub_rest = [sub_rest ([result(restind+1).stim_onset]-[result(restind).stim_onset]-1.5)];
    end
    slowoutlier=mean(sub_outliers);
    all_outliers=[all_outliers mean(sub_outliers)];
    sub_outliers=[];
    
    all_stds = [all_stds mean(sub_stds)];
    all_means = [all_means mean(submeans)];
    all_rest = [all_rest mean(sub_rest)];
    
    for run=1:4
        load(fullfile(sub_dir,['exp_pilot_' num2str(sub{1}) '_run_' num2str(run) '.mat']),'result');
        
        for trial=1:length(result)
            
            switch_type=regexprep(result(trial).switch_type,'-','_'); % ensure suitable for field name
            if isempty(result(trial).rt) % response time not stored for rest trials without response, but not needed here
                responsetime=nan;
            else
                responsetime = result(trial).rt;
            end
            
            if strcmp(switch_type,'dummy_trial')
                continue; % because switch type is not defined
            elseif strcmp(switch_type,'context_switch')
                % work out the focal task switch type
                task_type = result(trial).type;
                prev_task_type = result(trial-1).type;
                if strcmp(task_type(1),'r')
                    if strcmp(prev_task_type(1),'r')
                        cond='extended_rest';
                    else
                        cond='rest';
                    end
                    
                elseif strcmp(prev_task_type,'r')
                    cond='restart';
                elseif strcmp(task_type,prev_task_type)
                    cond='task_stay';
                elseif task_type(1)==prev_task_type(1)
                    cond='within_domain';
                else
                    cond='between_domain';
                end
                %task_type = result(trial).switch_type;
                
                if isnan(responsetime) % if they DID NOT respond
                    dat.RT.cswitch.(cond) = [dat.RT.cswitch.(cond), nan]; % RT won't be counted anyway
                    dat.accuracy.cswitch.(cond) = [dat.accuracy.cswitch.(cond) 0];
                    allsub_context_accuracy = [allsub_context_accuracy 0];
                else
                    if strcmp(result(trial).response,'context')
                        dat.RT.cswitch.(cond) = [dat.RT.cswitch.(cond), responsetime-result(trial).context_onset];
                        dat.accuracy.cswitch.(cond) = [dat.accuracy.cswitch.(cond) 1];
                        allsub_context_accuracy = [allsub_context_accuracy 1];
                        allsub_context_rt = [allsub_context_rt responsetime-result(trial).context_onset];
                    else % they responded to the regular task
                        dat.RT.cswitch.(cond) = [dat.RT.cswitch.(cond), nan];
                        dat.accuracy.cswitch.(cond) = [dat.accuracy.cswitch.(cond) 0];
                        allsub_context_accuracy = [allsub_context_accuracy 0];
                    end
                end
            else % for regular responses
                if ~isnan(responsetime)
                    %omit the rt it's larger than
                    %three std devs for this subject
                    if (responsetime-result(trial).stim_onset)>slowoutlier
                        sub_outlier = sub_outlier+1;
                        continue;
                    end
                    dat.RT.regular.(switch_type) = [dat.RT.regular.(switch_type), responsetime-result(trial).stim_onset];
                    allsub_rt = [allsub_rt responsetime-result(trial).stim_onset];
                    if run<3
                        first_half = [first_half responsetime-result(trial).stim_onset];
                    else
                        second_half = [second_half responsetime-result(trial).stim_onset];
                    end
                    if ~isempty(result(trial).accuracy)
                        allsub_accuracy = [allsub_accuracy result(trial).accuracy];
                        task_type = result(trial).type;
                        if ~strcmp(task_type,'r')
                            dat.accuracy.regular.(switch_type) = [dat.accuracy.regular.(switch_type) result(trial).accuracy];
                        end
                    end
                end
                
            end
            
        end % next trial
        
        %% calculate the average number of trials between the same type of switches
        run_switches={result.switch_type};
        if run ==1
            for switches = unique(run_switches)
                run_trials.(switches{1}) = [];
                run_time.(switches{1}) = [];
            end
        end
        for switches = unique(run_switches)
            
            indexes=find(contains(run_switches,switches{1}));
            run_trials.(switches{1}) = [run_trials.(switches{1}) mean(diff(indexes))];
            run_time.(switches{1}) = [run_time.(switches{1}) mean(diff([result(indexes).stim_onset]))];
        end
    end % next run
    
    for switches = unique(run_switches)
        avg_trials_between(subnum).(switches{1})=mean(run_trials.(switches{1}));
        avg_time_between(subnum).(switches{1})=mean(run_time.(switches{1}));
    end
    alloutlier_num =[alloutlier_num sub_outlier];
    first_half_sub_avg_rt(subnum)=mean(first_half);
    second_half_sub_avg_rt(subnum)=mean(second_half);
    
    for responsetype={'regular','cswitch'}
        for measure={'RT','accuracy'}
            
            switch measure{1}
                case 'RT'
                    wd = dat.RT.(responsetype{1}).within_domain(dat.accuracy.(responsetype{1}).within_domain==1);
                    bd = dat.RT.(responsetype{1}).between_domain(dat.accuracy.(responsetype{1}).between_domain==1);
                    ts = dat.RT.(responsetype{1}).task_stay(dat.accuracy.(responsetype{1}).task_stay==1);
                    restart = dat.RT.(responsetype{1}).restart(dat.accuracy.(responsetype{1}).restart==1);
                    rest = dat.RT.(responsetype{1}).rest(dat.accuracy.(responsetype{1}).rest==1);
                    extended_rest = dat.RT.(responsetype{1}).extended_rest(dat.accuracy.(responsetype{1}).extended_rest==1);
                    
                case 'accuracy'
                    wd = dat.accuracy.(responsetype{1}).within_domain;
                    bd = dat.accuracy.(responsetype{1}).between_domain;
                    ts = dat.accuracy.(responsetype{1}).task_stay;
                    restart = dat.accuracy.(responsetype{1}).restart;
                    rest = dat.accuracy.(responsetype{1}).rest;
                    extended_rest = dat.accuracy.(responsetype{1}).extended_rest;
            end
            
            %mean for each condition per sub
            sub_means = [nanmean(ts) nanmean(wd) nanmean(bd) nanmean(restart) nanmean(rest) nanmean(extended_rest)];
            
            sub_mean.([responsetype{1} '_' measure{1}])(subnum) = mean([ts wd bd restart rest extended_rest ]);
            behav_data.([measure{1} '_' responsetype{1}]) = [behav_data.([measure{1} '_' responsetype{1}]);sub_means];
            eachsubmean=mean(sub_means);
        end % next measure
    end % next response type
    
    %% memory tasks...
    load(fullfile(mem_dir,['pilotsub_' num2str(sub{1}) 'mem_task.mat']),'result');
    
    if(isempty(result))
        continue;
    end
    
    dprime=nan(1,length(result));
    sorterror=nan(1,length(result));
    
    for trial =1:length(result)
        
        targets=false(1,16);
        targets(result(trial).correct_repeat_position)=true;
        targets(result(trial).correct_novel_position)=true;
        selection=false(1,16);
        selection(result(trial).select_pics)=true;
        hitrate=sum(targets & selection)/sum(targets);
        farate=sum(~targets & selection)/sum(~targets);
        hitrate=min( max(hitrate,0.5/sum(targets)), (sum(targets)-0.5)/sum(targets) );
        farate=min( max(farate,0.5/sum(~targets)), (sum(~targets)-0.5)/sum(~targets) );
        dprime(trial)=norminv(hitrate)-norminv(farate);
        sorterror(trial)=sum(abs(result(trial).correct_order-result(trial).sort_order));
        
    end % next trial
    
    DP(1)=mean(dprime(strcmp('task-stay',[result.condition])));
    DP(2)=mean(dprime(strcmp('within-domain',[result.condition])));
    DP(3)=mean(dprime(strcmp('between-domain',[result.condition])));
    DP(4)=mean(dprime(strcmp('restart',[result.condition])));
    DP(5)=mean([dprime(strcmp('rest',[result.condition]))  dprime(strcmp('extended-rest',[result.condition]))]);
    dprime_allsubs(subnum,:)=DP;
    
    SE(1)=mean(sorterror(strcmp('task-stay',[result.condition])));
    SE(2)=mean(sorterror(strcmp('within-domain',[result.condition])));
    SE(3)=mean(sorterror(strcmp('between-domain',[result.condition])));
    SE(4)=mean(sorterror(strcmp('restart',[result.condition])));
    SE(5)=mean(sorterror(strcmp('rest',[result.condition])));
    SE(6)=mean(sorterror(strcmp('extended-rest',[result.condition])));
    
    se_allsubs(subnum,:)=SE;
    
end

means.regular_accuracy = mean(allsub_accuracy);
means.cswitch_accuracy = mean(allsub_context_accuracy);
means.regular_RT = mean(allsub_rt);
means.cswitch_RT = mean(allsub_context_rt);
alltasks_mean=mean(eachsubmean);
for responsetype={'regular','cswitch'}
    measure={'RT'};
    
    sub_data = behav_data.([measure{1} '_' responsetype{1}]);
    
    if strcmp(responsetype{1},'cswitch')
        sub_data = sub_data([1,4:end],:); %for context switch trials 3 participants are omitted
        sub_data(:,5)=mean(sub_data(:,5:6),2);
        [h,p,ci,stats]=ttest(sub_data(:,1),mean(sub_data(:,2:4),2)); %t-test of task repeats against tas switches averaged
        fprintf('\nContext-switch trial rt:\ntask switch combined versus task repeats\nt(%d)=%.3f, p=%.3f, BF=%.3f',stats.df,stats.tstat,p,t1smpbf(stats.tstat,stats.df+1));
        [tbl,rm]=simple_mixed_anova(sub_data(:,2:4),[],{'Conditions'}); %anova between task switch conditions
        fprintf('\nOne-way ANOVA for regular trial rt:\n');
        disp(tbl);
        F=tbl.F(3);x=tbl.DF(3);y=tbl.DF(3+1); % roi
        BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nCondition BF=%.3e\n',BF10);
                
        [h,p,ci,stats]=ttest(mean(sub_data(:,5),2),mean(sub_data(:,1:4),2)); %ttest between tasks and rest (djm: removed condition 6 from first input, as already combined above)
        brest=t1smpbf(stats.tstat,stats.df+1);
        fprintf('\nContext-switch trial rt:\ntask switch combined versus rest\nt(%d)=%.3f, p=%.3f, BF=%.3e',stats.df,stats.tstat,p,brest);
        
        [~,pvalue,CI,stats_rt_cs] = ttest(sub_data);
        error = abs(CI-repmat(nanmean(sub_data,1),2,1));
        sub_means=nanmean(sub_data);
        
        csvwrite('rt_cs.csv',sub_means(1:5));
        csvwrite('allsubs_rt_cs.csv',sub_data(:,1:5));
        csvwrite('rt_cs_error.csv',error(1,1:5));
        
    else
        [~,p,ci,stats]=ttest(sub_data(:,1),mean(sub_data(:,2:4),2)); %task switch combined versus task repeats
        fprintf('Regular trial rt:\ntask switch combined versus task repeats\nt(%d)=%.3f, p=%.3f, BF=%.3e',stats.df,stats.tstat,p,t1smpbf(stats.tstat,stats.df+1));
        csvwrite('alltaskrt_reg.csv',mean(sub_data(:,1:5),1));
        csvwrite('allsubs_rt.csv',sub_data(:,1:5));
        rtregmean=mean(mean(sub_data(:,1:4),1));
        
        [tbl,rm]=simple_mixed_anova(sub_data(:,2:4),[],{'Conditions'});
        fprintf('\nOne-way ANOVA for regular trial rt:\n');
        disp(tbl);
        F=tbl.F(3);x=tbl.DF(3);y=tbl.DF(3+1); % 
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('BF=%.3e',BF10);
    
        [~,~,ci,stats]=ttest(sub_data(:,1:5));
        error=abs(ci-repmat(mean(sub_data(:,1:5)),2,1));
        csvwrite('alltaskrt_reg_error.csv',error(1,2:5));
        
        [~,p,ci,stats]=ttest(sub_data(:,2),sub_data(:,3)); %task switch combined versus task repeats
        fprintf('\nwithin-domain versus between-domain\nt(%d)=%.3f, p=%.3f, BF =%.3f',stats.df,stats.tstat,p,t1smpbf(stats.tstat,stats.df+1));
        
        [~,p,ci,stats]=ttest(sub_data(:,2),sub_data(:,4)); %task switch combined versus task repeats
        fprintf('\nwithin-domain versus restart\nt(%d)=%.3f, p=%.3f, BF =%.3f',stats.df,stats.tstat,p,t1smpbf(stats.tstat,stats.df+1));
        
        [~,p,ci,stats]=ttest(sub_data(:,3),sub_data(:,4)); %task switch combined versus task repeats
        fprintf('\nbetween-domain versus restart\nt(%d)=%.3f, p=%.3f, BF =%.3f\n',stats.df,stats.tstat,p,t1smpbf(stats.tstat,stats.df+1));
        
        temp=sub_data(:,1:5)-repmat(sub_data(:,1),1,5);
        [H,p,CI,stats_rt_reg]=ttest(temp);
        error=abs(CI-repmat(mean(temp),2,1));
        csvwrite('reg_rt_error.csv',error(1,2:5));
        
    end
    % next measure
end % next response type



%% load univariate results
%/group/duncan-lab/users for linux
load(fullfile(pwd,'others_minus_task_stay211022.mat'));
%data is 5 conditions (wd,bd,rest,restart,ex-rest against ts) x 6 rois (MD,
%visual,somatosensory, core, dmpfc, mtl) x 35 subjects

%ANOVA of core and mtl rois
%first, average the rest and extended rest
data(6,:,:) = mean(data([3,5],:,:),1);
unv_data = data([1,2,4],[4,6],:); %wd, bd, restart x core,mtl x 35 subs
anova_matrix = permute(unv_data,[3,1,2]);
tbl=simple_mixed_anova(anova_matrix,[],{'Conditions','ROIs'});
disp('Univariate data two-way ANOVA with condition and ROI');
disp(tbl);
F=tbl.F(5);
x=tbl.DF(5); % DF of effect
y=tbl.DF(5+1); % residual DF of error
BF10=rmANOVAbf_FB23(F,x,y);
disp(['BF for ROIs: ' num2str(BF10)]);
F=tbl.F(7);x=tbl.DF(7);y=tbl.DF(7+1); % interaction
BF10=rmANOVAbf_FB23(F,x,y);
disp(['BF for interaction of condition/ROIs: ' num2str(BF10)]);

%4=core, 6=mtl, 1=md, 5=dmpfc
roi_names={'core_dmn','mtl','md','dmpfc'};
r=1;
for roi=[4,6,1,5]
    mean_data=mean(data,3);
    mean_data = mean_data([1,2,4,3,5],roi);
    rest_ext_rest = squeeze(data(6,roi,:));
    mean_data(4) = mean(rest_ext_rest);
    
    [~,~,CI,~] = ttest(rest_ext_rest);
    rest_error= abs(CI-mean_data(4));
    mean_data = mean_data(1:4);

        csvwrite(sprintf('%s_unv.csv',roi_names{r}),mean_data);
        csvwrite(sprintf('%s_unv_allsubs.csv',roi_names{r}),squeeze(data([1,2,4,6],roi,:)));
 
    
    ROI_name = split(dataspecs.rois{roi},'.');
    
    %plot error bar
    for cond = 1:5
        [~,p_ts,CI,stats] = ttest(data(cond,roi,:));
        mean_data2=mean(data,3);
        abs_value = abs(CI-mean_data2(cond,roi));
        error(1,cond)=abs_value(:,:,1);
        error(2,cond)=abs_value(:,:,2);
        
        if cond ==1
            p_values_ts.(ROI_name{1})=p_ts;
            t_values_ts.(ROI_name{1})=stats.tstat;
        else
            p_values_ts.(ROI_name{1})=[p_values_ts.(ROI_name{1}) p_ts];
            t_values_ts.(ROI_name{1})=[t_values_ts.(ROI_name{1}) stats.tstat];
        end
        %do pair-wise t-tests between each condition and make a tbl
        for cond2=1:5
            if cond>=cond2
                continue;
            end
            [~,p,~,stats]=ttest(data(cond,roi,:),data(cond2,roi,:));
            bf_ttest=t1smpbf(stats.tstat,stats.df+1);
            first_cond =  split(dataspecs.contrasts{cond},'_');
            first_cond =first_cond{1};
            second_cond =  split(dataspecs.contrasts{cond2},'_');
            second_cond=second_cond{1};
            pair_ttest.([first_cond '_' second_cond])=[p stats.tstat stats.df bf_ttest];
            
        end
        p_values.(ROI_name{1})=pair_ttest;
    end
    hold on
    e=errorbar(1:4,mean_data,[error(1,[1,2,4]) rest_error(1)],[error(2,[1,2,4]) rest_error(1)]);
    e.LineStyle='none';
    
    csvwrite(sprintf('%s_unv_error_bar.csv',roi_names{r}),[error(1,[1,2,4]) rest_error(1)]);
    r=r+1;
    
    fprintf('\nROI:%s\n',ROI_name{1});
    all_switch_avg = mean(squeeze(data([1,2,4],roi,:)));
    [~,p,~,vs_ts_tstats] = ttest(all_switch_avg');
    disp('Avg task condition against task stay: ');
    fprintf('t(%d)=%.3f, p=%.3f, BF=%.3e\n',vs_ts_tstats.df,vs_ts_tstats.tstat,p,t1smpbf(vs_ts_tstats.tstat,vs_ts_tstats.df+1)); 
    
    rest_avg = squeeze(data(6,roi,:));
    [~,p,~,vs_rest_tstats] = ttest(all_switch_avg'-rest_avg);
    disp('Avg task condition against rest: ');
    fprintf('t(%d)=%.3f, p=%.3f\n, BF = %.3e\n',vs_rest_tstats.df,vs_rest_tstats.tstat,p, t1smpbf(vs_rest_tstats.tstat,vs_rest_tstats.df+1));
    
    %One way anova of conditions
    [tbl,rm]=simple_mixed_anova(squeeze(data([1,2,4],roi,:))',[],{'Conditions'},{});
    disp('One way ANOVA of types of task transitions');
    disp(tbl);
    F=tbl.F(3);x=tbl.DF(3);y=tbl.DF(3+1); % 
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('BF=%.3e',BF10);
    
    disp('T-tests, in order of p values, tstats and df');
    disp(p_values.(ROI_name{1}));
end

%% memory task
[~,~,CI] = ttest(dprime_allsubs(2:end,:));
error = abs(CI-repmat(nanmean(dprime_allsubs(2:end,:)),2,1));

[~,p,CI,stats] = ttest(mean(dprime_allsubs(2:end,:),2));
fprintf('\nMemory task recognition accuracy:\nAll conditions combined versus 0:\nt(%d)=%.3f, p=%.3f, BF=%.3e',stats.df,stats.tstat,p,t1smpbf(stats.tstat,stats.df+1));
[~,p,CI,stats] = ttest(dprime_allsubs(2:end,1),mean(dprime_allsubs(2:end,2:4),2));
fprintf('\nTransitions in focal task combined versus task repeats:\nt(%d)=%.3f, p=%.3f, BF=%.3f',stats.df,stats.tstat,p,t1smpbf(stats.tstat,stats.df+1));
[~,p,CI,stats] = ttest(dprime_allsubs(2:end,5),mean(dprime_allsubs(2:end,1:4),2)); %all tasks vs rest
fprintf('\nTransitions in focal task combined versus rest:\nt(%d)=%.3f, p=%.3f, BF=%.3f',stats.df,stats.tstat,p,t1smpbf(stats.tstat,stats.df+1));

[tbl,rm]=simple_mixed_anova(dprime_allsubs(2:end,2:4),[],{'Conditions'});
F=tbl.F(3);x=tbl.DF(3);y=tbl.DF(3+1); % 
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('BF=%.3e',BF10);
fprintf('\nOne-way ANOVA for conditions:\n');
disp(tbl);

csvwrite('recognition_accuracy.csv',nanmean(dprime_allsubs(2:end,:)));
csvwrite('recognition_error.csv',error(1,:));
csvwrite('allsubs_rec_accuracy.csv',dprime_allsubs(2:end,:));

%% load decoding results

load('decoding_data.mat');

% context_dir = 'aa5_analysis_220622unsmoothed_GLMcontext_conditions\aamod_decoding4_az_00001';
% novelty_dir = 'aa5_analysis_220622unsmoothed_GLMnovelty_conditions\aamod_decoding4_az_00002';
% 
% conditions_dir = {context_dir, novelty_dir};
% conditions = {'Context','Novelty'};
% switch_types =  {'task_stay','within_domain','between_domain','rest','restart','extended_rest'};

decode_results = struct('Context',[],'Novelty',[]);
%for each decoding model
for num = 1:2
%     
%     decode_data = fullfile(datadir,conditions_dir{num});
%     names = dir(decode_data);
%     subjects = 0;
%     anova_matrix = [];
    
    %for each subject
%     for k = 1:length(names)
%         subname = names(k).name;
%         if strcmp(subname,'.')||strcmp(subname,'..')
%             continue;
%         end
%         subjects = subjects +1;
%         toload = fullfile(decode_data,subname);
%         toload = fullfile(toload,['Decode' conditions{num}]);
%         
%         %for each condition
%         for cond = 1:6
%             file_toload = ['Decode' conditions{num} '_dprime_set000' num2str(cond) '.mat'];
%             load(fullfile(toload,file_toload),'results');
%             %subject's data for this condition in this model
%             anova_matrix(subjects,:,cond) = results.dprime.output';
%         end
%     end
    
    
    anova_matrix=decoding_data.(['Yeo_' conditions{num}]);
    all_cond_networks_avg=mean(anova_matrix,3); %combine all conditions
    [sorted_networks,index]=sort(mean(all_cond_networks_avg,1),2,'descend');
    
    [H,p,CI,~] = ttest(all_cond_networks_avg(:,index));
    
    decode_results.(conditions{num}).sig_fdr=fdr_bh(p);
    decode_results.(conditions{num}).index=index;
    decode_results.(conditions{num}).sig_p=H;
    decode_results.(conditions{num}).sorted=sorted_networks;
    decode_results.(conditions{num}).error=abs(CI-repmat(mean(all_cond_networks_avg(:,index),1),2,1));
    
    %first find significant parcellations
    network_uncorrected_sig = index(H==1);
    network_to_plot = index(fdr_bh(p)==1);
    
    %average all the significant parcellations
    comb_sig = squeeze(mean(anova_matrix(:,network_to_plot,:),2));
    [H,p,CI,~] = ttest(comb_sig);
    [h,p,ci,stats]=ttest(mean(comb_sig(:,[2,3,5]),2),comb_sig(:,1)); %ttest between sig combined net's task transitions and ts
    
    fprintf('\nDecoding %s\nFor significant FDR-corrected Yeo parcellations:\nAll transitions combined versus task-stay:\nt(%d)=%.3f, p=%.3f, BF=%.3f',conditions{num},stats.df,stats.tstat,p,t1smpbf(stats.tstat,stats.df+1));
    
    [h,p,ci,stats]=ttest(mean(comb_sig(:,[1,2,3,5]),2),mean(comb_sig(:,[4,6]),2)); %ttest between all tasks and rest
    fprintf('\nAll transitions combined versus rest:\nt(%d)=%.3f, p=%.3f, BF=%.3f',stats.df,stats.tstat,p,t1smpbf(stats.tstat,stats.df+1));
    
    [tbl,rm]=simple_mixed_anova(comb_sig(:,[2,3,5]),[],{'Conditions'});
    fprintf('\nOne-way ANOVA for conditions:\n');
    disp(tbl);
    F=tbl.F(3);x=tbl.DF(3);y=tbl.DF(3+1); % 
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('BF=%.3e',BF10);
    decode_results.(conditions{num}).combined = mean(comb_sig);
    decode_results.(conditions{num}).all =comb_sig;
    decode_results.(conditions{num}).combined_error = abs(CI-repmat(mean(comb_sig),2,1));
    
    %% Core and MTL decoding
    anova_matrix=decoding_data.(['mtl_core_' conditions{num}]);
    
    %MTL and core
    [h,p,CI,~]=ttest(squeeze(anova_matrix(:,1:2,:)));
    CI=CI-repmat(mean(anova_matrix(:,1:2,:),1),2,1);
    
    csvwrite([lower(conditions{num}) '_mtl_core_error.csv'],CI);
    csvwrite(['mtl_core_' lower(conditions{num}) '_decoding.csv'],squeeze(mean(anova_matrix(:,1:2,:)))); %2 roi x 6 condition
    csvwrite(['allsubs_mtl_core_' lower(conditions{num}) '_decoding.csv'],anova_matrix(:,1:2,:)); %2 roi x 6 condition

    comb_transitions = mean(anova_matrix(:,1:2,[2,3,5]),3);
    nanova=zeros(35,2,2);
    nanova(:,:,1)=squeeze(comb_transitions);
    nanova(:,:,2)=squeeze(anova_matrix(:,1:2,1));
    [tbl,rm]=simple_mixed_anova(nanova,[],{'ROI','Conditions'},{});
    fprintf('\nDecoding of %s for Core and MTL\nTwo-way ANOVA for conditions(task repeat and averaged transitions) and ROI:\n',conditions{num});
    disp(tbl);
    F=tbl.F(3);x=tbl.DF(3);y=tbl.DF(3+1); % roi
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nROI BF=%.3e\n',BF10);
        F=tbl.F(5);x=tbl.DF(5);y=tbl.DF(5+1); % cond
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nCondition BF=%.3e\n',BF10);
        F=tbl.F(7);x=tbl.DF(7);y=tbl.DF(7+1); % interaction
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nInteraction BF=%.3e\n',BF10);
    nanova(:,:,1)=squeeze(mean(anova_matrix(:,1:2,[1,2,3,5]),3)); % djm added (to compare to all tasks, not just switches)
    nanova(:,:,2)=squeeze(mean(anova_matrix(:,1:2,[4,6]),3));
    [tbl,rm]=simple_mixed_anova(nanova,[],{'ROI','Conditions'},{});
    fprintf('\nDecoding of %s for Core and MTL\nTwo-way ANOVA for conditions(rest and averaged transitions) and ROI:\n',conditions{num});
    disp(tbl);
    F=tbl.F(3);x=tbl.DF(3);y=tbl.DF(3+1); % roi
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nROI BF=%.3e\n',BF10);
        F=tbl.F(5);x=tbl.DF(5);y=tbl.DF(5+1); % cond
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nCondition BF=%.3e\n',BF10);
        F=tbl.F(7);x=tbl.DF(7);y=tbl.DF(7+1); % interaction
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nInteraction BF=%.3e\n',BF10);
    [tbl,rm]=simple_mixed_anova(anova_matrix(:,1:2,[2,3,5]),[],{'ROI','Conditions'},{});
    fprintf('\nDecoding of %s for Core and MTL\nTwo-way ANOVA for types of transitions and ROI:\n',conditions{num});
    disp(tbl);
    F=tbl.F(3);x=tbl.DF(3);y=tbl.DF(3+1); % roi
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nROI BF=%.3e\n',BF10);
        F=tbl.F(5);x=tbl.DF(5);y=tbl.DF(5+1); % cond
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nCondition BF=%.3e\n',BF10);
        F=tbl.F(7);x=tbl.DF(7);y=tbl.DF(7+1); % interaction
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nInteraction BF=%.3e\n',BF10);
    
end

filename = 'yeoDecodingResults.json';
jsontxt = jsonencode(decode_results);
fid= fopen(fullfile(pwd,filename),'w');
fprintf(fid,jsontxt);
fclose(fid);


return;

function bf10 = t1smpbf(t,n,r)
%
% bf10 = t1smpbf(t,n,[r=0.707])
%
% Calculates JZS Bayes Factor for a one-sample t-test given t and sample size n.
% The optional input r is the scale factor which defaults to 0.707.
% This quantifies the evidence in favour of the alternative hypothesis. 
% See Rouder et al, 2009, Psychon Bull Rev for details.
%

% Default scale factor
if nargin < 3
    r = 0.707;
end

% Function to be integrated
F = @(g,t,n,r) (1+n.*g.*r.^2).^(-1./2) .* (1 + t.^2./((1+n.*g.*r.^2).*(n-1))).^(-n./2) .* (2.*pi).^(-1./2) .* g.^(-3./2) .* exp(-1./(2.*g));

% Bayes factor calculation
bf01=nan(size(t));
for i=1:numel(t)
    bf01(i) = (1 + t(i).^2./(n-1)).^(-n/2) ./ integral(@(g) F(g,t(i),n,r),0,Inf);
end

% Invert Bayes Factor
bf10 = 1 ./ bf01;

function [tbl,rm] = simple_mixed_anova(datamat, varargin)

% tbl = simple_mixed_anova(datamat, between_factors, within_factor_names,
% between_factor_names)
%
% Repeated-measures or mixed ANOVA with any number of factors. 
%
% Function built on top of existing MATLAB functions in an attempt to
% simplify user inputs and manipulations while still being able to perform
% all common ANOVA designs.
%
% DATAMAT is a numerical matrix containing the values associated with each 
% level of each within-subject factor in each subject (responses). The 
% subjects should always be the first dimension of the matrix. Each 
% within-subject factor is another dimension.
%
% BETWEEN_FACTORS is a numerical matrix where each row represents a subject 
% and each column represents a between-subjects factor. A given value
% represents the level of the column's factor associated with the row's 
% subject. Optional.
%
% WITHIN_FACTOR_NAMES is a cell array of strings indicating the name of
% each within-subject factor (one for each dimension of the datamat
% variable except the first, in order). Optional.
%
% BETWEEN_FACTOR_NAMES is a cell array of strings indicating the name of
% each between-subjects factor (one for each column of the between_factors
% matrix, in order). These factors are assumed to be categorical (groups). 
% Optional.
%
% TABLE is a table indicating the F statistics, p-values and other
% statistics associated with each term of the model. The "intercept" terms
% can be ignored (e.g. "(Intercept):WS01" indicates the main effect of 
% WS01).
%
% RM is a structure with the repeated measures model parameters and
% statistics.
%
% Does not support covariates or partial models (without all interactions)
% for now.
%
%
% EXAMPLE
%
% A design with 24 subjects and 2 within-subject factors, the first one 
% having 3 levels (time: pre-test, 1st post-test, 2nd post-test) and the 
% second one 4 levels (experimental condition: A, B, C, D). The subjects 
% are grouped in 4 groups: 2 variables with 2 levels each (gender: male or
% female; age group: young or old).
%
% The input datamat should be a 24 x 3 x 4 matrix with each row
% corresponding to a subject. So the element datamat(1,1,1) will correspond
% to the response of the subject #1 in the pre-test in experimental
% condition A, the element datamat(2,3,2) will correspond to the response
% of the subject #2 in the 2nd post-test in experimental condition B, and
% so on.
%
% The input between_factors will be a 24 x 2 matrix with each row
% corresponding to a subject and each column to a between-subjects factor
% (gender and age). Each column will be filled with 1s and 2s, or
% 0s and 1s, or other numbers, indicating the gender/age group of the 
% respective subject.
%
% tbl = simple_mixed_anova(datamat, between_factors, {'Time', 'Exp_cond'},
% {'Gender', 'Age_group'})
%
% Copyright 2017, Laurent Caplette
% https://www.researchgate.net/profile/Laurent_Caplette

% Check if correct number of inputs
narginchk(1,4)

% Assign inputs to variables; if none, will be empty array
between_factors = [];
within_factor_names = [];
between_factor_names = [];
if nargin>1    
    between_factors = varargin{1};
    if nargin>2
        within_factor_names = varargin{2};
        if nargin>3
            between_factor_names = varargin{3};
        end
    end
end

% Determine numbers of variables and measures
nWithin = ndims(datamat)-1;
nBetween = size(between_factors,2);
nVars = size(datamat);
nVars = nVars(2:end); % don't use the nb of subjects on the first dim
nMeas = prod(nVars);

% Check if dimensions of matrices are ok
if size(datamat,1)<2
    error('There must be more than one subject.')
end
if ~isempty(between_factors)
    if size(between_factors,1)~=size(datamat,1)
        error('Both input matrices must have the same nb of subjects.')
    end
end

% Check if there is more than one unique value
if length(unique(datamat))<2
    error('The data matrix must contain more than one unique value.')
end
for ii = 1:size(between_factors,2)
    if length(unique(between_factors(:,ii)))<2
        error('Each between-subjects factor must contain more than one unique value.')
    end
end

% Error if more variable names than variables as input
if length(between_factor_names)>nBetween
    error('Too many between-subject factor names or not enough between-subject variables as input.')
end
if length(within_factor_names)>nWithin
    error('Too many within-subject factor names or not enough within-subject variables as input.')
end

% Check validity of variable names
for ii = 1:length(between_factor_names)
    if ~isvarname(between_factor_names{ii})
        error('Variable names must be continuous strings starting with a letter and without symbols.')
    end
end
for ii = 1:length(within_factor_names)
    if ~isvarname(within_factor_names{ii})
        error('Variable names must be continuous strings starting with a letter and without symbols.')
    end
end

% Assign variable names if not enough or empty
if length(between_factor_names)<nBetween
    nMissing = nBetween - length(between_factor_names);
    BS = repmat('BS', [nMissing 1]); % list of 'BS'
    missing_factor_names = cellstr([BS num2str([1:nMissing]', '%02.0f')]);
    between_factor_names = [between_factor_names missing_factor_names];
end
if length(within_factor_names)<nWithin
    nMissing = nWithin - length(within_factor_names);
    WS = repmat('WS', [nMissing 1]); % list of 'WS'
    missing_factor_names = cellstr([WS num2str([1:nMissing]', '%02.0f')]);
    within_factor_names = [within_factor_names missing_factor_names];
end

% Create table detailing within-subject design
withinVarLevels = fullfact(nVars); % all level combinations
within_table = array2table(withinVarLevels, 'VariableNames', within_factor_names);
for ii = 1:nWithin % ensure that each within-subject factor is categorical (levels==discrete)
    evalc(sprintf('within_table.%s = categorical(within_table.%s)', within_factor_names{ii}, within_factor_names{ii}));
end

% Vectorize all dimensions after first one of the data matrix
y = datamat(:,:);

% Create data table
yList = repmat('Y', [nMeas 1]); % list of 'Y'
numList = num2str([1:nMeas]', '%03.0f'); % support up to 999 measures
measureNames = cellstr([yList numList]); % create names for every measure
for ii = 1:nBetween % add between-subject factors
    measureNames{nMeas+ii} = between_factor_names{ii};
end
total_table = array2table([y between_factors],'VariableNames', measureNames);
for ii = 1:nBetween % ensure that each between-subject factor is categorical (levels/groups==discrete)
    evalc(sprintf('total_table.%s = categorical(total_table.%s)', between_factor_names{ii}, between_factor_names{ii}));
end

% Create between-subjects model using Wilkinson notation
betweenModel = '';
for ii = 1:nBetween
    betweenModel = [betweenModel,measureNames{nMeas+ii},'*'];
end
betweenModel = betweenModel(1:end-1); % remove last star
if isempty(betweenModel)
    betweenModel = '1'; % if no between-subjects factor, put constant term (usually implicit)
end

% Create within-subject model using Wilkinson notation
withinModel = '';
for ii = 1:nWithin
    withinModel = [withinModel,within_factor_names{ii},'*']; % stars for full model (all interactions)
end
withinModel = withinModel(1:end-1); % remove last star

% Fit repeated measures model
rm = fitrm(total_table, sprintf('%s-%s~%s', measureNames{1}, measureNames{nMeas}, betweenModel),...
    'WithinDesign', within_table);

% Run ANOVA
tbl = ranova(rm, 'WithinModel', withinModel);


% Copyright (c) 2015, David Groppe
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report)

% djm: to ignore nans
sz=size(pvals); 
pvals=pvals(:);
notnans=find(~isnan(pvals));
pvals=pvals(notnans);

if nargin<1,
    error('You need to provide a vector or matrix of p-values.');
else
    if ~isempty(find(pvals<0,1)),
        error('Some p-values are less than 0.');
    elseif ~isempty(find(pvals>1,1)),
        error('Some p-values are greater than 1.');
    end
end

if nargin<2,
    q=.05;
end

if nargin<3,
    method='pdep';
end

if nargin<4,
    report='no';
end

s=size(pvals);
if (length(s)>2) || s(1)>1,
    [p_sorted, sort_ids]=sort(reshape(pvals,1,prod(s)));
else
    %p-values are already a row vector
    [p_sorted, sort_ids]=sort(pvals);
end
[dummy, unsort_ids]=sort(sort_ids); %indexes to return p_sorted to pvals order
m=length(p_sorted); %number of tests

if strcmpi(method,'pdep'),
    %BH procedure for independence or positive dependence
    thresh=(1:m)*q/m;
    wtd_p=m*p_sorted./(1:m);
    
elseif strcmpi(method,'dep')
    %BH procedure for any dependency structure
    denom=m*sum(1./(1:m));
    thresh=(1:m)*q/denom;
    wtd_p=denom*p_sorted./(1:m);
    %Note, it can produce adjusted p-values greater than 1!
    %compute adjusted p-values
else
    error('Argument ''method'' needs to be ''pdep'' or ''dep''.');
end

if nargout>3,
    %compute adjusted p-values; This can be a bit computationally intensive
    adj_p=zeros(1,m)*NaN;
    [wtd_p_sorted, wtd_p_sindex] = sort( wtd_p );
    nextfill = 1;
    for k = 1 : m
        if wtd_p_sindex(k)>=nextfill
            adj_p(nextfill:wtd_p_sindex(k)) = wtd_p_sorted(k);
            nextfill = wtd_p_sindex(k)+1;
            if nextfill>m
                break;
            end;
        end;
    end;
    adj_p=reshape(adj_p(unsort_ids),s);
end

rej=p_sorted<=thresh;
max_id=find(rej,1,'last'); %find greatest significant pvalue
if isempty(max_id),
    crit_p=0;
    h=pvals*0;
    adj_ci_cvrg=NaN;
else
    crit_p=p_sorted(max_id);
    h=pvals<=crit_p;
    adj_ci_cvrg=1-thresh(max_id);
end

%%%% djm: to ignore nans
tested=h;
h=nan(prod(sz),1);
h(notnans)=tested;
h=reshape(h(:),sz);

if nargout>3
    tested=adj_p;
    adj_p=nan(prod(sz),1);
    adj_p(notnans)=tested;
    adj_p=reshape(adj_p(:),sz);
end
%%%%%

if strcmpi(report,'yes'),
    n_sig=sum(p_sorted<=crit_p);
    if n_sig==1,
        fprintf('Out of %d tests, %d is significant using a false discovery rate of %f.\n',m,n_sig,q);
    else
        fprintf('Out of %d tests, %d are significant using a false discovery rate of %f.\n',m,n_sig,q);
    end
    if strcmpi(method,'pdep'),
        fprintf('FDR/FCR procedure used is guaranteed valid for independent or positively dependent tests.\n');
    else
        fprintf('FDR/FCR procedure used is guaranteed valid for independent or dependent tests.\n');
    end
end
