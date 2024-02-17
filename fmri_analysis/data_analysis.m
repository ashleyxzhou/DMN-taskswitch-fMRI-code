function data_analysis()
%Script for generating results and figures for manuscript

dbstop if error
addpath Z:\Duncan-lab\users\dm01\MoreTools;


%% load univariate results
%/group/duncan-lab/users for linux
load('Z:\Duncan-lab\users\az01\task_switch\DataAnalysis\aa5_analysis_220622unsmoothed_GLMside_conditions\others_minus_task_stay211022.mat');
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

%4=core, 6=mtl, 1=md
for roi=[4,6,1]
    mean_data=mean(data,3);
    mean_data = mean_data([1,2,4,3,5],roi);
    rest_ext_rest = squeeze(data(6,roi,:));
    mean_data(4) = mean(rest_ext_rest);
    
    [~,~,CI,~] = ttest(rest_ext_rest);
    rest_error= abs(CI-mean_data(4));
    mean_data = mean_data(1:4);
    
    if roi==4
        csvwrite('core_dmn_unv.csv',mean_data);
    elseif roi==6
        csvwrite('mtl_unv.csv',mean_data);
    end
    
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
            first_cond =  split(dataspecs.contrasts{cond},'_');
            first_cond =first_cond{1};
            second_cond =  split(dataspecs.contrasts{cond2},'_');
            second_cond=second_cond{1};
            pair_ttest.([first_cond '_' second_cond])=[p stats.tstat stats.df];
            
        end
        p_values.(ROI_name{1})=pair_ttest;
    end
    hold on
    e=errorbar(1:4,mean_data,[error(1,[1,2,4]) rest_error(1)],[error(2,[1,2,4]) rest_error(1)]);
    e.LineStyle='none';
    
    if roi ==4
        csvwrite('core_dmn_unv_error_bar.csv',[error(1,[1,2,4]) rest_error(1)]);
    elseif roi==6
        csvwrite('mtl_unv_error_bar.csv',[error(1,[1,2,4]) rest_error(1)]);
    end
    
    fprintf('\nROI:%s\n',ROI_name{1});
    all_switch_avg = mean(squeeze(data([1,2,4],roi,:)));
    [~,p,~,vs_ts_tstats] = ttest(all_switch_avg');
    disp('Avg task condition against task stay: ');
    fprintf('t(%d)=%.3f, p=%.3f, BF=%.3e\n',vs_ts_tstats.df,vs_ts_tstats.tstat,p,bf.bfFromT(vs_ts_tstats.tstat,vs_ts_tstats.df)); 
    
    rest_avg = squeeze(data(6,roi,:));
    [~,p,~,vs_rest_tstats] = ttest(all_switch_avg'-rest_avg);
    disp('Avg task condition against rest: ');
    fprintf('t(%d)=%.3f, p=%.3f\n, BF = %.3f',vs_rest_tstats.df,vs_rest_tstats.tstat,p, bfFromT(vs_rest_tstats.tstat,vs_rest_tstats.df));
    
    %One way anova of conditions
    [tbl,rm]=simple_mixed_anova(squeeze(data([1,2,4],roi,:))',[],{'Conditions'},{});
    disp('One way ANOVA of types of task transitions');
    disp(tbl);
    F=tbl.F(3);x=tbl.DF(3);y=tbl.DF(3+1); % 
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('BF=%.3e',BF10);
    
    disp('T-tests, in order of p values, tstats and df');
    disp(p_values.(ROI_name{1}));
end

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
all_outliers=[];

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
    addpath Z:\Duncan-lab\users\az01\task_switch\scripts;
    
    %calculate mean +3std of all trials from subject as threshold for slow
    %outliers
    for run = 1:4
        load(['exp_pilot_' num2str(sub{1}) '_run_' num2str(run) '.mat'],'result');
        ok=~cellfun(@isempty,{result.rt}); rts=[result(ok).rt]-[result(ok).stim_onset]; slowoutlier=mean(rts(~isnan(rts)))+3*std(rts(~isnan(rts)));
        
        all_outliers=[all_outliers slowoutlier];
    end
    slowoutlier=mean(all_outliers);
    all_outliers=[];
    
    for run=1:4
        load(['exp_pilot_' num2str(sub{1}) '_run_' num2str(run) '.mat'],'result');
        
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
    end % next run
    
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
            
        end % next measure
    end % next response type
    
    %% memory tasks...
    load(['pilotsub_' num2str(sub{1}) 'mem_task.mat'],'result');
    
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

for responsetype={'regular','cswitch'}
    measure={'RT'};
    
    sub_data = behav_data.([measure{1} '_' responsetype{1}]);
    
    if strcmp(responsetype{1},'cswitch')
        sub_data = sub_data([1,4:end],:); %for context switch trials 3 participants are omitted
        sub_data(:,5)=mean(sub_data(:,5:6),2);
        [h,p,ci,stats]=ttest(sub_data(:,1),mean(sub_data(:,2:4),2)); %t-test of task repeats against tas switches averaged
        fprintf('\nContext-switch trial rt:\ntask switch combined versus task repeats\nt(%d)=%.3f, p=%.3f, BF=%.3f',stats.df,stats.tstat,p,bfFromT(stats.tstat,stats.df));
        [tbl,rm]=simple_mixed_anova(sub_data(:,2:4),[],{'Conditions'}); %anova between task switch conditions
        fprintf('\nOne-way ANOVA for regular trial rt:\n');
        disp(tbl);
        F=tbl.F(3);x=tbl.DF(3);y=tbl.DF(3+1); % roi
        BF10=rmANOVAbf_FB23(F,x,y);fprintf('\nCondition BF=%.3e\n',BF10);
                
        [h,p,ci,stats]=ttest(mean(sub_data(:,5),2),mean(sub_data(:,1:4),2)); %ttest between tasks and rest (djm: removed condition 6 from first input, as already combined above)
        brest=bfFromT(stats.tstat,stats.df);
        fprintf('\nContext-switch trial rt:\ntask switch combined versus rest\nt(%d)=%.3f, p=%.3f, BF=%.3e',stats.df,stats.tstat,p,brest);
        
        [~,pvalue,CI,stats_rt_cs] = ttest(sub_data);
        error = abs(CI-repmat(nanmean(sub_data,1),2,1));
        sub_means=nanmean(sub_data);
        
        csvwrite('rt_cs.csv',sub_means(1:5));
        csvwrite('rt_cs_error.csv',error(1,1:5));
        
    else
        [~,p,ci,stats]=ttest(sub_data(:,1),mean(sub_data(:,2:4),2)); %task switch combined versus task repeats
        fprintf('Regular trial rt:\ntask switch combined versus task repeats\nt(%d)=%.3f, p=%.3f, BF=%.3e',stats.df,stats.tstat,p,bfFromT(stats.tstat,stats.df));
        csvwrite('alltaskrt_reg.csv',mean(sub_data(:,1:5),1));
        
        [tbl,rm]=simple_mixed_anova(sub_data(:,2:4),[],{'Conditions'});
        fprintf('\nOne-way ANOVA for regular trial rt:\n');
        disp(tbl);
        F=tbl.F(3);x=tbl.DF(3);y=tbl.DF(3+1); % 
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('BF=%.3e',BF10);
    
        [~,~,ci,stats]=ttest(sub_data(:,1:5));
        error=abs(ci-repmat(mean(sub_data(:,1:5)),2,1));
        csvwrite('alltaskrt_reg_error.csv',error(1,2:5));
        
        [~,p,ci,stats]=ttest(sub_data(:,2),sub_data(:,3)); %task switch combined versus task repeats
        fprintf('\nwithin-domain versus between-domain\nt(%d)=%.3f, p=%.3f, BF =%.3f',stats.df,stats.tstat,p,bfFromT(stats.tstat,stats.df));
        
        [~,p,ci,stats]=ttest(sub_data(:,2),sub_data(:,4)); %task switch combined versus task repeats
        fprintf('\nwithin-domain versus restart\nt(%d)=%.3f, p=%.3f, BF =%.3f',stats.df,stats.tstat,p,bfFromT(stats.tstat,stats.df));
        
        [~,p,ci,stats]=ttest(sub_data(:,3),sub_data(:,4)); %task switch combined versus task repeats
        fprintf('\nbetween-domain versus restart\nt(%d)=%.3f, p=%.3f, BF =%.3f\n',stats.df,stats.tstat,p,bfFromT(stats.tstat,stats.df));
        
        temp=sub_data(:,1:5)-repmat(sub_data(:,1),1,5);
        [H,p,CI,stats_rt_reg]=ttest(temp);
        error=abs(CI-repmat(mean(temp),2,1));
        csvwrite('reg_rt_error.csv',error(1,2:5));
        
    end
    % next measure
end % next response type

[~,~,CI] = ttest(dprime_allsubs(2:end,:));
error = abs(CI-repmat(nanmean(dprime_allsubs(2:end,:)),2,1));

[~,p,CI,stats] = ttest(mean(dprime_allsubs(2:end,:),2));
fprintf('\nMemory task recognition accuracy:\nAll conditions combined versus 0:\nt(%d)=%.3f, p=%.3f, BF=%.3e',stats.df,stats.tstat,p,bfFromT(stats.tstat,stats.df));
[~,p,CI,stats] = ttest(dprime_allsubs(2:end,1),mean(dprime_allsubs(2:end,2:4),2));
fprintf('\nTransitions in focal task combined versus task repeats:\nt(%d)=%.3f, p=%.3f, BF=%.3f',stats.df,stats.tstat,p,bfFromT(stats.tstat,stats.df));
[~,p,CI,stats] = ttest(dprime_allsubs(2:end,5),mean(dprime_allsubs(2:end,1:4),2)); %all tasks vs rest
fprintf('\nTransitions in focal task combined versus rest:\nt(%d)=%.3f, p=%.3f, BF=%.3f',stats.df,stats.tstat,p,bfFromT(stats.tstat,stats.df));

[tbl,rm]=simple_mixed_anova(dprime_allsubs(2:end,2:4),[],{'Conditions'});
F=tbl.F(3);x=tbl.DF(3);y=tbl.DF(3+1); % 
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('BF=%.3e',BF10);
fprintf('\nOne-way ANOVA for conditions:\n');
disp(tbl);

csvwrite('recognition_accuracy.csv',nanmean(dprime_allsubs(2:end,:)));
csvwrite('recognition_error.csv',error(1,:));


%% load decoding results

datadir='Z:\Duncan-lab\users\az01\task_switch\DataAnalysis';

context_dir = 'aa5_analysis_220622unsmoothed_GLMcontext_conditions\aamod_decoding4_az_00001';
novelty_dir = 'aa5_analysis_220622unsmoothed_GLMnovelty_conditions\aamod_decoding4_az_00002';

conditions_dir = {context_dir, novelty_dir};
conditions = {'Context','Novelty'};
switch_types =  {'task_stay','within_domain','between_domain','rest','restart','extended_rest'};
decode_results = struct('Context',[],'Novelty',[]);
%for each decoding model
for num = 1:2
    
    decode_data = fullfile(datadir,conditions_dir{num});
    names = dir(decode_data);
    subjects = 0;
    anova_matrix = [];
    
    %for each subject
    for k = 1:length(names)
        subname = names(k).name;
        if strcmp(subname,'.')||strcmp(subname,'..')
            continue;
        end
        subjects = subjects +1;
        toload = fullfile(decode_data,subname);
        toload = fullfile(toload,['Decode' conditions{num}]);
        
        %for each condition
        for cond = 1:6
            file_toload = ['Decode' conditions{num} '_dprime_set000' num2str(cond) '.mat'];
            load(fullfile(toload,file_toload),'results');
            %subject's data for this condition in this model
            anova_matrix(subjects,:,cond) = results.dprime.output';
        end
    end
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
    
    fprintf('\nDecoding %s\nFor significant FDR-corrected Yeo parcellations:\nAll transitions combined versus task-stay:\nt(%d)=%.3f, p=%.3f, BF=%.3f',conditions{num},stats.df,stats.tstat,p,bfFromT(stats.tstat,stats.df));
    
    [h,p,ci,stats]=ttest(mean(comb_sig(:,[1,2,3,5]),2),mean(comb_sig(:,[4,6]),2)); %ttest between all tasks and rest
    fprintf('\nAll transitions combined versus rest:\nt(%d)=%.3f, p=%.3f, BF=%.3f',stats.df,stats.tstat,p,bfFromT(stats.tstat,stats.df));
    
    [tbl,rm]=simple_mixed_anova(comb_sig(:,[2,3,5]),[],{'Conditions'});
    fprintf('\nOne-way ANOVA for conditions:\n');
    disp(tbl);
    F=tbl.F(3);x=tbl.DF(3);y=tbl.DF(3+1); % 
    BF10=rmANOVAbf_FB23(F,x,y);fprintf('BF=%.3e',BF10);
    decode_results.(conditions{num}).combined = mean(comb_sig);
    decode_results.(conditions{num}).combined_error = abs(CI-repmat(mean(comb_sig),2,1));
    
end

filename = 'yeoDecodingResults.json';
jsontxt = jsonencode(decode_results);
fid= fopen(fullfile('Z:\Duncan-lab\users\az01\task_switch\scripts',filename),'w');
fprintf(fid,jsontxt);
fclose(fid);

%for mtl and core rois
context_dir = 'original_aa5_analysis_220622unsmoothed_GLMscontext_conditions\aamod_decoding4_az_00001';
novelty_dir = 'original_aa5_analysis_220622unsmoothed_GLMsnovelty_conditions\aamod_decoding4_az_00002';
conditions_dir={context_dir,novelty_dir};
for num = 1:2
    measures = {'mean_decision_value_djm','dprime','AUC_minus_chance','accuracy_minus_chance'};
    
    data = fullfile(datadir,conditions_dir{num});
    names = dir(data);
    subjects = 0;
    anova_matrix = [];
    
    %for each subject
    for k = 1:length(names)
        subname = names(k).name;
        if strcmp(subname,'.')||strcmp(subname,'..')
            continue;
        end
        subjects = subjects +1;
        toload = fullfile(data,subname);
        toload = fullfile(toload,['Decode' conditions{num}]);
        
        %for each measure
        i=2;
        %for i = 1:length(measures)
        
        for cond = 1:6
            file_toload = ['Decode' conditions{num} '_' measures{i} '_set000' num2str(cond) '.mat'];
            load(fullfile(toload,file_toload),'results');
            
            %subject's data for this condition in this model
            sub_data = results.(measures{i}).output';
            anova_matrix(subjects,:,cond)=sub_data;
            
        end
    end
    
    %MTL and core
    [h,p,CI,~]=ttest(squeeze(anova_matrix(:,1:2,:)));
    CI=CI-repmat(mean(anova_matrix(:,1:2,:),1),2,1);
    csvwrite([lower(conditions{num}) '_mtl_core_error.csv'],CI);
    csvwrite(['mtl_core_' lower(conditions{num}) '_decoding.csv'],squeeze(mean(anova_matrix(:,1:2,:)))); %2 roi x 6 condition
    
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


return;