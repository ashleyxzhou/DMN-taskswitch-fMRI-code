function taskswitch_data_analysis_masterscript()
dbstop if error
addpath Z:\Duncan-lab\users\dm01\MoreTools;

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
%tbl_core = simple_mixed_anova(anova_matrix(:,:,1),[],{'Conditions'});
F=tbl.F(5);
x=tbl.DF(5); % DF of effect
y=tbl.DF(5+1); % residual DF of error
BF10=rmANOVAbf_FB23(F,x,y);
disp(['BF for ROIs: ' num2str(BF10)]);
F=tbl.F(7);x=tbl.DF(7);y=tbl.DF(7+1); % interaction
BF10=rmANOVAbf_FB23(F,x,y);
disp(['BF for interaction of condition/ROIs: ' num2str(BF10)]);

%%%load behaviour results
subs= {220513,220525,5,220469,220471,220473,220474,220475,220480,220481,220482,220485,220491,220492,220493,220494,220495,220496,220499,220503,220506,220509,220510,220512,220514,220519,220523,220524,220526,220533,220535,220536,220538,220539,220542};

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
    network_to_plot = index(fdr_bh(p)==1);
    
    %average all the significant parcellations
    comb_sig = squeeze(mean(anova_matrix(:,network_to_plot,:),2));
    [H,p,CI,~] = ttest(comb_sig);
    [h,p,ci,tstat]=ttest(mean(comb_sig(:,[2,3,5]),2),comb_sig(:,1)); %ttest between sig combined net's task transitions and ts
    [h,p,ci,tstat]=ttest(mean(comb_sig(:,[2,3,5]),2),mean(comb_sig(:,[4,6]),2));
    tbl=simple_mixed_anova(comb_sig(:,[2,3,5]),[],{'Conditions'});
    decode_results.(conditions{num}).combined = mean(comb_sig);
    decode_results.(conditions{num}).combined_error = abs(CI-repmat(mean(comb_sig),2,1));
    
end

filename = 'yeoDecodingResults.json';
jsontxt = jsonencode(decode_results);
fid= fopen(fullfile('Z:\Duncan-lab\users\az01\task_switch\scripts',filename),'w');
fprintf(fid,jsontxt);
fclose(fid);

return;