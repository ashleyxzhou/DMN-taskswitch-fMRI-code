function poster_figures()
dbstop if error
addpath Z:\Duncan-lab\users\dm01\MoreTools;

%initialise figures
index=1;
column=2;
figure(11); clf(11)
pos=fitplots2([2 1],'',2);%regular trial RT, core and mtl univariate
figure(20); clf(20) %RT for context switch, context decoding, novelty decoding, memory test d'
pos2=fitplots2([2 2],'',2);
try
    set(10,'WindowStyle','normal', 'WindowState','maximized','color','w');
catch
end
ax=nan(6,1);

%% load univariate results
%/group/duncan-lab/users for linux
load('Z:\Duncan-lab\users\az01\task_switch\DataAnalysis\aa5_analysis_220622unsmoothed_GLMside_conditions\others_minus_task_stay211022.mat');
%load('delta_others_minus_task_stay241122.mat');
%data is 5 conditions (wd,bd,rest,restart,ex-rest against ts) x 6 rois (MD,
%visual,somatosensory, core, dmpfc, mtl) x 35 subjects

%ANOVA of core and mtl rois
%first, average the rest and extended rest
data(6,:,:) = mean(data([3,5],:,:),1);
unv_data = data([1,2,4],[4,6],:); %wd, bd, restart x core,mtl x 35 subs
anova_matrix = permute(unv_data,[3,1,2]);
tbl=simple_mixed_anova(anova_matrix,[],{'Conditions','ROIs'});
F=tbl.F(5);
x=tbl.DF(5); % DF of effect
y=tbl.DF(5+1); % residual DF of error
BF10=rmANOVAbf_FB23(F,x,y);
disp(['BF for ROIs: ' num2str(BF10)]);
F=tbl.F(7);x=tbl.DF(7);y=tbl.DF(7+1); % interaction
BF10=rmANOVAbf_FB23(F,x,y);
disp(['BF for interaction of condition/ROIs: ' num2str(BF10)]);

for roi=[4,6]
    mean_data=mean(data,3);
    mean_data = mean_data([1,2,4,3,5],roi);
    rest_ext_rest = squeeze(data(6,4,:));
    mean_data(4) = mean(rest_ext_rest);
    
    [~,~,CI,~] = ttest(rest_ext_rest);
    rest_error= abs(CI-mean_data(4));
    mean_data = mean_data(1:4);
    
    if roi==4
        csvwrite('mtl_unv.csv',mean_data);
    else
        csvwrite('core_dmn_unv.csv',mean_data);
    end
    color_bar = {[1 157/255 27/255],[251/255 187/255 16/255],[233/255 77/255 54/255],[103/255 192/255 77/255],[190/255 43/255 187/255],[0 167/255 136/255]};
    colornum=1;
    
    ROI_name = split(dataspecs.rois{roi},'.');
    
    
    ax(column,index)=subplot('position',pos{column,index});
    bar(mean_data,'FaceColor',color_bar{colornum},'EdgeColor','none');
    box off;
    colornum=colornum+1;
    conditions_names = regexprep(dataspecs.contrasts([1,2,4,3]),'_',' ');
    
    set(gca, 'xticklabel',{'within-domain','between-domain','restart','rest'},'FontSize',13);
    ylabel('BOLD Response','FontSize',18);
    %ylim([0 1.2]);
    %     xtickangle(30);
    title('Core DMN response for regular trials, compared to task-stay','interpreter','none','FontSize',18);
    %legend({'Core DMN');
    %mean_data=mean(data,3);
    
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
        csvwrite('mtl_unv_error_bar.csv',[error(1,[1,2,4]) rest_error(1)]);
    else
        csvwrite('core_dmn_unv_error_bar.csv',[error(1,[1,2,4]) rest_error(1)]);
        
    end
    
    all_switch_avg = mean(squeeze(data([1,2,4],4,:)));
    [~,p,~,vs_ts_tstats] = ttest(all_switch_avg);
    %Bayesfactor_switch_vs_ts = t1smpbf(vs_ts_tstats.tstat,length(all_switch_avg));
    rest_avg = mean([squeeze(data(3,4,:)) squeeze(data(5,4,:))],2)';
    [~,p,~,vs_rest_tstats] = ttest(all_switch_avg-rest_avg);
    %Bayesfactor_switch_vs_rest = t1smpbf(vs_rest_tstats.tstat,length(all_switch_avg));
    simple_mixed_anova(squeeze(data([1,2,4],1,:))',[],{'Conditions'},{});
    
end

behav_data.RT_regular=[];
behav_data.RT_cswitch=[];
behav_data.accuracy_regular=[];
behav_data.accuracy_cswitch=[];

subs= {220469,220471,220473,220474,220475,220480,220481,220482,220485,220491,220492,220493,220494,220495,220496,220499,220503,220506,220509,220510,220512,220514,220519,220524,220526,220533,220535,220536,220538,220539,220542};
%subs = {220513,220539,220525};
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

%calculate each sub mean for all trials, and grand mean nfor all subs
%for each condition mean, minus sub mean and add grand mean
%calculate std and mean for all condition sub means
%error = std/square(length(sub))


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
    for run = 1:4
        
        load(['exp_pilot_' num2str(sub{1}) '_run_' num2str(run) '.mat'],'result');
        
        %if sub{1}==5
        %slowoutlier=5; %pilot subject has results file in a different format
        %else
        ok=~cellfun(@isempty,{result.rt}); rts=[result(ok).rt]-[result(ok).stim_onset]; slowoutlier=mean(rts(~isnan(rts)))+3*std(rts(~isnan(rts)));
        %end
        
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
                    %instead of omitting the rts, test if it's larger than
                    %three std devs than this subjects
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
    
    %%% memory tasks...
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

column=1;
for responsetype={'regular','cswitch'}
    measure={'RT'};
    
    sub_data = behav_data.([measure{1} '_' responsetype{1}]);
    
    if strcmp(responsetype{1},'cswitch')
        sub_data = sub_data(3:end,:); %for context switch trials 3 participants are omitted
        sub_data(:,5)=mean(sub_data(3:end,5:6),2);
        [h,p,ci,tstat]=ttest(sub_data(:,1),mean(sub_data(:,2:4),2)) %t-test of task repeats against tas switches averaged
        simple_mixed_anova(sub_data(:,2:4),[],{'Conditions'}); %anova between task switch conditions
        BF_anova_between_ts = rmANOVAbf_FB23(0.0073691,2,62);
        
        [h,p,ci,tstat]=ttest(mean(sub_data(:,5:6),2),mean(sub_data(:,1:4),2)); %ttest between tasks and rest
        bfFromT(-5.5456,31);
        
        adj=sub_data-repmat(nanmean(sub_data,2),1,6)+nanmean(sub_data(:));
        [~,pvalue,CI,stats_rt_cs] = ttest(adj)
        error = abs(CI-repmat(nanmean(adj,1),2,1));
        ax(column,index)=subplot('position',pos2{column,index});
        sub_means=nanmean(adj);
        bar(sub_means(1:5),'FaceColor',color_bar{colornum},'EdgeColor','none');
        
        hold on
        e=errorbar(1:5,sub_means(1:5),error(1,1:5),error(2,1:5));
        %set(gca,'XTick',[]);
        title('RT for context switch (jungle detection)','FontSize',18);
        ylabel('Reaction time (s)','FontSize',15);
        ylim([0.5 1]);
        set(gca, 'xticklabel',{'task-stay','within-domain','between-domain','restart','rest'},'FontSize',13);
        xtickangle(30);
        
        csvwrite('rt_cs.csv',sub_means(1:5));
        csvwrite('rt_cs_error.csv',error(1,1:5));
        
    else
        sub_data(:,5)=mean(sub_data(:,5:6),2);
        ax(column,index)=subplot('position',pos{column,index});
        temp=sub_data(:,1:5)-repmat(sub_data(:,1),1,5);
        [H,p,CI,stats_rt_reg]=ttest(temp);
        [h,pp,cc,statswdbd]=ttest(temp(:,3),temp(:,2)); %wd bd diff
        
        csvwrite('rt_reg.csv',temp(2:5));
        csvwrite('rt_reg_error.csv',error(1,2:5));
    end
    % next measure
end % next response type

adj=dprime_allsubs-repmat(nanmean(dprime_allsubs,2),1,5)+nanmean(dprime_allsubs(:));
[~,~,CI] = ttest(adj);
error = abs(CI-repmat(nanmean(dprime_allsubs),2,1));

csvwrite('recognition_accuracy.csv',nanmean(dprime_allsubs));
csvwrite('recognition_error.csv',error(1,:));


datadir='Z:\Duncan-lab\users\az01\task_switch\DataAnalysis';

context_dir = 'aa5_analysis_220622unsmoothed_GLMcontext_conditions\aamod_decoding4_az_00001';
novelty_dir = 'aa5_analysis_220622unsmoothed_GLMnovelty_conditions\aamod_decoding4_az_00002';

conditions_dir = {context_dir, novelty_dir};
conditions = {'Context','Novelty'};
switch_types =  {'task_stay','within_domain','between_domain','rest','restart','extended_rest'};

%for each decoding model
index = 2;
column=1;
anova_pvalues = {};
for num = 1:2
    
    sub_avg_dprime = struct('task_stay',zeros(1,17),'within_domain',zeros(1,17),'between_domain',zeros(1,17),'rest',zeros(1,17),'restart',zeros(1,17),'extended_rest',zeros(1,17));
    
    sub_avg_array ={sub_avg_dprime};
    networks_avg = [];
    measures = {'dprime'};
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
        
        %for each measure
        for i = 1:length(measures)
            
            sub_networks=[];
            
            %for each condition
            for cond = 1:6
                file_toload = ['Decode' conditions{num} '_' measures{i} '_set000' num2str(cond) '.mat'];
                load(fullfile(toload,file_toload),'results');
                
                %subject's data for this condition in this model
                networks = results.(measures{i}).output;
                
                sub_networks =[sub_networks networks];
                networks = networks';
                if cond==4 %rest
                    rest_network = networks;
                    
                elseif cond ==6 %extended_rest
                    networks = mean([rest_network; networks]);
                    sub_avg_array{i}.(switch_types{4}) = sub_avg_array{i}.(switch_types{4})+networks;
                else
                    sub_avg_array{i}.(switch_types{cond}) = sub_avg_array{i}.(switch_types{cond})+networks;
                end
                
                anova_matrix(subjects,:,cond) = networks;
                
            end
            
            sub_network_avg = mean(sub_networks,2);
            
            networks_avg = [networks_avg sub_network_avg];
            
        end
        
    end
    
    all_cond_networks_avg=mean(networks_avg,2);
    [sorted,net_index]=sort(all_cond_networks_avg,1,'descend');
    
    network_to_plot = net_index(1:3);
    
    ax(column,index)=subplot('position',pos2{column,index});
    column=num+1;
    
    decode_results = [];
    for taskswitch = 1:6
        allsubs_combo_three=[];
        for net = network_to_plot'
            allsubs_combo_three = [allsubs_combo_three;anova_matrix(:,net,taskswitch)];
        end
        decode_results = [decode_results nanmean(allsubs_combo_three)];
        [~,~,CI]=ttest(allsubs_combo_three);
        error(num,taskswitch)=abs(CI(1)-mean(allsubs_combo_three));
    end
    
    decode_results =decode_results([1 2 3 5 4]);
    error(num,1:5)=error(num,[1 2 3 5 4]);
    
    toexport=[lower(conditions{num}) '_decoding2.csv'];
    csvwrite(toexport,decode_results);
    
    toexport=[lower(conditions{num}) '_error2.csv'];
    csvwrite(toexport,error(num,1:5));
end

return;