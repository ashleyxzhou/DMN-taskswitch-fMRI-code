function mri_pilot

dbstop if error
commandwindow;

TR=1.208;

%scansync('reset', TR)

subinfo = inputdlg({'Subject ID:','run'});
%set up .mat var to record data
resultname = ['exp_pilot_' subinfo{1} '_run_' subinfo{2} '.mat'];

while exist(resultname,'file')
    subinfo = inputdlg({'Subject ID:','run'});
    resultname = ['exp_pilot_' subinfo{1} '_run_' subinfo{2} '.mat'];
end

result={};

here=fileparts(mfilename('fullpath'));
stimuli_folder_dir=fullfile(fileparts(here),'stimuli');
save(fullfile(here,resultname),'result');


% Setup PTB with some default values
PsychDefaultSetup(2);

% Seed the random number generator.
%id x 10 + run
rng(str2double(subinfo{1})*10+str2double(subinfo{2}));
[series, switch_seq] = random_seq_gen();

% Skip sync tests for demo purposes only
Screen('Preference', 'SkipSyncTests', 1);

HideCursor;

%----------------------------------------------------------------------
%                       Screen setup
%----------------------------------------------------------------------

% Set the screen number to the external secondary monitor if there is one
% connected
screenNumber = max(Screen('Screens'));

% Define black, white and grey
white = WhiteIndex(screenNumber);
grey = white / 2;
black = BlackIndex(screenNumber);

% Open the screen
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [], 32, 2);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Flip to clear
Screen('Flip', window);

% Set the text size
Screen('TextSize', window, 30);

% Query the maximum priority level
topPriorityLevel = MaxPriority(window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Set the blend funciton for the screen
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

DrawFormattedText(window, 'Please wait while experiment is loading', 'center', 'center', black);
Screen('Flip', window);



%----------------------------------------------------------------------
%                       Keyboard information
%----------------------------------------------------------------------

% Keybpard setup
spaceKey = KbName('space');
escapeKey = KbName('ESCAPE');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');
RestrictKeysForKbCheck([spaceKey leftKey rightKey escapeKey]);



%----------------------------------------------------------------------
%                      Experimental Image List
%----------------------------------------------------------------------

% % Get the image files for the experiment
% liv_small_imageFolder = 'Z:\users\az01\task_switch\stimuli\Living_Small\';
% imgList = dir(fullfile(liv_small_imageFolder, '*.bmp'));
% imgList = {imgList(:).name};
% numImages = length(imgList);
%
% liv_large_imageFolder = 'Z:\users\az01\task_switch\stimuli\Living_Large\';
% non_liv_large_imageFolder = 'Z:\users\az01\task_switch\stimuli\Man_Large\';
% non_liv_small_imageFolder = 'Z:\users\az01\task_switch\stimuli\Man_Small\';

liv_small_imageFolder = [stimuli_folder_dir '\Living_Small\'];
imgList = dir(fullfile(liv_small_imageFolder, '*.bmp'));
imgList = {imgList(:).name};
numImages = length(imgList);

liv_large_imageFolder = [stimuli_folder_dir '\Living_Large\'];
non_liv_large_imageFolder = [stimuli_folder_dir '\Man_large\'];
non_liv_small_imageFolder = [stimuli_folder_dir '\Man_Small\'];

objects=[];
objects.nonliving_small = [];
objects.nonliving_large = [];
objects.living_small = [];
objects.living_large = [];

%Load all the files into textures
%40 objects in each category:
%nonliving_small,nonliving_large,living_small,living_large
for i =1:numImages
    imageName = imgList{i};
    
    %load images, make textures (vectors) in before loop
    theImage = imread([non_liv_large_imageFolder imageName]);
    non_l = Screen('MakeTexture', window, theImage);
    objects.nonliving_large = [objects.nonliving_large non_l];
    
    theImage = imread([non_liv_small_imageFolder imageName]);
    non_s = Screen('MakeTexture', window, theImage);
    objects.nonliving_small = [objects.nonliving_small non_s];
    
    theImage = imread([liv_large_imageFolder imageName]);
    liv_l = Screen('MakeTexture', window, theImage);
    objects.living_large = [objects.living_large liv_l];
    
    theImage = imread([liv_small_imageFolder imageName]);
    liv_s = Screen('MakeTexture', window, theImage);
    objects.living_small = [objects.living_small liv_s];
    
    DrawFormattedText(window, sprintf('Please wait while experiment is loading%s',repmat('.',1,i)), 'center', 'center', black);
    Screen('Flip', window);
end


imgList = dir(fullfile([stimuli_folder_dir '\dummy_cafe'], '*.jpg'));
imgList = {imgList(:).name};
numImages = length(imgList);
dummy_cafe = nan(1,numImages);
for i =1:numImages
    imageName = imgList{i};
    theImage = imread([stimuli_folder_dir '\dummy_cafe\' imageName]);
    dummy_cafe(i) = Screen('MakeTexture', window, theImage);
end

imgList = dir(fullfile([stimuli_folder_dir '\dummy_beach\'], '*.jpg'));
imgList = {imgList(:).name};
numImages = length(imgList);
dummy_beach = nan(1,numImages);
for i =1:numImages
    imageName = imgList{i};
    theImage = imread([stimuli_folder_dir '\dummy_beach\' imageName]);
    dummy_beach(i) = Screen('MakeTexture', window, theImage);
end


text_type_names = {'LettersA' 'LettersI' 'LettersBoth' 'LettersNeither'};
obj_type_names = {'living_small' 'nonliving_small' 'living_large' 'nonliving_large'};

beach_imageFolder = [stimuli_folder_dir '\run' subinfo{2} '\scenes_beach\'];
cafe_imageFolder = [stimuli_folder_dir '\run' subinfo{2} '\scenes_cafe\'];
jungle_imageFolder = [stimuli_folder_dir '\scenes_jungle\'];

context_imgList = dir(fullfile(cafe_imageFolder, '*.jpg'));
context_imgList = {context_imgList(:).name};

jcontext_imgList = dir(fullfile(jungle_imageFolder, '*.jpg'));
jcontext_imgList = {jcontext_imgList(:).name};

num_context = length(context_imgList);
all_beach =nan(1,num_context);
all_cafe =nan(1,num_context);
all_jungle =nan(1,18);

for index = 1:18
    contextName = jcontext_imgList{index};
    jungle_image = imread([jungle_imageFolder contextName]);
    all_jungle(index) = Screen('MakeTexture', window, jungle_image);
    clear jungle_image;
end
all_jungle = Shuffle(all_jungle);

for index = 1:num_context
    contextName = context_imgList{index};
    beach_image = imread([beach_imageFolder contextName]);
    all_beach(index) = Screen('MakeTexture', window, beach_image);
    clear beach_image;
    
    cafe_image = imread([cafe_imageFolder contextName]);
    all_cafe(index) = Screen('MakeTexture', window, cafe_image);
    clear cafe_image;
    
    DrawFormattedText(window, sprintf('Please wait while experiment is loading%s',repmat('.',1,i+index)), 'center', 'center', black);
    Screen('Flip', window);
end

beach_index = all_beach;
cafe_index = all_cafe;
all_beach = Shuffle(all_beach);
all_cafe = Shuffle(all_cafe);

fileid{1} = fopen(fullfile(stimuli_folder_dir, '\Letters\Both4.txt'), 'r');
fileid{2} = fopen(fullfile(stimuli_folder_dir, '\Letters\A_only4.txt'), 'r');
fileid{3} = fopen(fullfile(stimuli_folder_dir, '\Letters\I_only4.txt'), 'r');
fileid{4} = fopen(fullfile(stimuli_folder_dir, '\Letters\Neither4.txt'), 'r');

for n = 1:4
    data = textscan(fileid{n},'%s');
    
    fclose(fileid{n});
    
    switch n
        case 1
            LettersBoth = data{1};
        case 2
            LettersA = data{1};
        case 3
            LettersI = data{1};
        case 4
            LettersNeither = data{1};
    end
end



%----------------------------------------------------------------------
%                        Condition Matrix
%----------------------------------------------------------------------

%[mask,switch_seq,series] = context_switch_gen(series,switch_seq);
[mask,switch_seq,series] = context_switch_gen_v5d_djm(series,switch_seq);

% Make the matrix for each condition types' context pictures
% For each condition, there is an array of 12 pics, 1 repeated + 6 novel

% 1* 1 2 1 3 1 4 1 5 1 6 1

cond_names = {'task_stay','between_domain','within_domain','rest','extended_rest','restart'};

% context_matrix.beach={};
% context_matrix.cafe={};
array_repeated_pics_beach =nan(1,6);
array_repeated_pics_cafe =nan(1,6);

array_novel_pics_cafe = nan(1,6);
array_novel_pics_beach = nan(1,6);
for cond = 1: 6
    %choose random repeated picture number, remove from pool
    
    repeated_pic_beach = all_beach(1);
    repeat_beach_num = find(beach_index==repeated_pic_beach);
    all_beach(1)=[];
    array_repeated_pics_beach(cond) = repeat_beach_num;
    
    repeated_pic_cafe = all_cafe(1);
    repeat_cafe_num = find(cafe_index==repeated_pic_cafe);
    array_repeated_pics_cafe(cond)= repeat_cafe_num;
    all_cafe(1)=[];
    
    %     %choose 7 random novel pics
    %     context_matrix.beach(cond)=[];
    %     context_matrix.cafe(cond)=[];
    context_matrix(cond).beach = [repeated_pic_beach repeated_pic_beach];
    context_matrix(cond).cafe = [repeated_pic_cafe repeated_pic_cafe];
    for ind = 2:6
        
        if ind == 2
            array_novel_pics_cafe(cond)=find(cafe_index==all_cafe(1));
            array_novel_pics_beach(cond)=find(beach_index==all_beach(1));
        end
        
        
        context_matrix(cond).beach = [context_matrix(cond).beach all_beach(1) repeated_pic_beach ];
        all_beach(1)=[];
        
        context_matrix(cond).cafe = [context_matrix(cond).cafe all_cafe(1)  repeated_pic_cafe];
        all_cafe(1)=[];
        
    end
    context_matrix(cond).beach_counter = 1;
    context_matrix(cond).cafe_counter = 1;
    
end

array_repeated_pics_beach=[array_repeated_pics_beach array_novel_pics_beach];
array_repeated_pics_cafe=[array_repeated_pics_cafe array_novel_pics_cafe];
repeated_beach = ['sub' subinfo{1} '_run_' subinfo{2} 'repeatBeach.mat'];
save(repeated_beach,'array_repeated_pics_beach');

repeated_cafe = ['sub' subinfo{1} '_run_' subinfo{2} 'repeatCafe.mat'];
save(repeated_cafe,'array_repeated_pics_cafe');

rest_avg = 2;
rt_avg = rest_avg;

sides = zeros(1, length(switch_seq));
answers = zeros(1,length(switch_seq));
%counterbalancing the sides
for i = 1:6
    all_index = find(ismember(switch_seq,cond_names(i)));
    mask_index = mask(all_index);
    beach_indexes = Shuffle(find(mask_index == 1));
    right_side = all_index(beach_indexes(1:6));
    sides(right_side)=1;
    left_side = all_index(beach_indexes(7:12));
    answers(randsample(right_side,3))=1;
    answers(randsample(left_side,3))=1;
    
    cafe_indexes = Shuffle(find(mask_index == 0));
    right_side = all_index(cafe_indexes(1:6));
    sides(right_side)=1;
    left_side = all_index(cafe_indexes(7:12));
    answers(randsample(right_side,3))=1;
    answers(randsample(left_side,3))=1;
end
%----------------------------------------------------------------------
%                        Fixation Cross
%----------------------------------------------------------------------

% Screen Y fraction for fixation cross
crossFrac = 0.0167;

% Here we set the size of the arms of our fixation cross
fixCrossDimPix = windowRect(4) * crossFrac;

% Now we set the coordinates (these are all relative to zero we will let
% the drawing routine center the cross in the center of our monitor for us)
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];

% Set the line width for our fixation cross
lineWidthPix = 4;


%----------------------------------------------------------------------
%                      Experimental Loop
%----------------------------------------------------------------------
try
    % Start screen
    DrawFormattedText(window, 'Loaded successfully. Please wait for experiment to begin.', 'center', 'center', black);
    Screen('Flip', window);
    KbWait;
    
    % wait for scanner & dummy scans
    DrawFormattedText(window,'Waiting for scanner...','center','center',[0 0 0]);
    Screen('Flip',window);
    expstarttime=[];
    ndummies=10;
    pulsetime=nan(1,ndummies);
    for dummy=1:ndummies
        [pulsetime(dummy), pulsenum, daqstate]=scansync(1,Inf);
        if isempty(expstarttime), expstarttime=GetSecs; result(1).exp_start_time = expstarttime; end
        fprintf('\nPulse %d: %2.3f.',dummy,pulsetime(dummy));
        DrawFormattedText(window,strcat('Waiting for scanner...',num2str(11-dummy)),'center','center',[0 0 0]);
        Screen('Flip',window);
    end
    fprintf('\nMeasured TR = %2.3fs\n',mean(diff(pulsetime)));
    WaitSecs(1);
    
    %stores the sequence of trial
    
    %a1, b1, a2, b2,
    %is it living, 'A', is it larger than shoebox, 'I'
    Screen('TextSize', window, 60);
    for trial = 1:length(series)
        %
        %         if isempty(switches)==0 && trial == switches(1)
        %            %if this was a context_switch trial
        %            switches(1) =[];
        %            which = all_jungle(1);
        %            all_jungle(1) = 0;
        %         end
        
        scene_is_beach =mask(trial);
        %saves trial info
        type = series{trial};
        
        if mask(trial) ==2
            result(trial).scene = 'jungle';
        end
        
        %a = semantic; b = lexical
        result(trial).domain =type(1);
        result(trial).switch_type = switch_seq{trial};
        result(trial).trial = trial;
        result(trial).type = series{trial};
        
        %Screen('DrawTexture', window,which, [], CenterRectOnPointd([0 0 screenXpixels screenYpixels], xCenter, yCenter), 0);
        % Draw a fixation cross for the start of the trial
        Screen('FillRect', window, grey);
        Screen('DrawLines', window, allCoords,...
            lineWidthPix, white, [xCenter yCenter], 2);
        [~, trialstart]=Screen('Flip', window);
        %WaitSecs(1.5);
        % Draw the fixation cross in white, set it to the center of our screen
        
        %randomly pick which side of screen to put picture
        r = randi([0,1], 1);
        answer  = randi([0,1], 1);
        
        respMade = 0;
        Priority(topPriorityLevel);
        rest = strcmp(switch_seq{trial},'rest') || strcmp(switch_seq{trial},'extended-rest');
        
        
        %change dummy background to jungle? + dummy trial with scene that
        %fits with next sequence of backgrounds
        
        
        if strcmp(switch_seq{trial}, 'dummy_trial')
            if scene_is_beach == 1
                result(trial).scene = 'beach';
                dummy_background = randsample(dummy_beach,1);
                result(trial).scene_pic_num = find(dummy_beach==dummy_background);
            else
                result(trial).scene = 'cafe';
                dummy_background = randsample(dummy_cafe,1);
                result(trial).scene_pic_num = find(dummy_cafe==dummy_background);
            end
            which = dummy_background;
            type = series{trial};
        elseif strcmp(switch_seq{trial}, 'context_switch')
            which = all_jungle(1);
            all_jungle(1)=[];
            type = series{trial};
        else
            r = sides(trial);
            answer = answers(trial);
            cond_index =  find(ismember(cond_names,switch_seq{trial}));
            %checks which picture to show for each condition
            if scene_is_beach == 1
                next_scene = context_matrix(cond_index).beach_counter;
                context_matrix(cond_index).beach_counter = context_matrix(cond_index).beach_counter +1;
                which = context_matrix(cond_index).beach(next_scene);
                
                %                 if rem(length(context_matrix(cond_index).beach),2) == 0
                %                     result(trial).scene_repeat = 1;
                %                 else
                %                     result(trial).scene_novel = 1;
                %                 end
                %context_matrix(cond_index).beach(1)=[];
                result(trial).scene = 'beach';
                result(trial).scene_pic_num = find(beach_index==which);
                
            else
                next_scene = context_matrix(cond_index).cafe_counter;
                context_matrix(cond_index).cafe_counter =  context_matrix(cond_index).cafe_counter +1;
                which = context_matrix(cond_index).cafe(next_scene);
                %context_matrix(cond_index).cafe(1)=[];
                result(trial).scene = 'cafe';
                result(trial).scene_pic_num = find(cafe_index==which);
                
            end
            if next_scene ==1
                result(trial).scene_novel = 1;
            elseif rem(next_scene,2) == 0
                result(trial).scene_repeat = 1;
            else
                result(trial).scene_novel = floor(next_scene/2)+1;
            end
            
        end
        
        result(trial).correct_answer = answer;
        result(trial).side = r;
        
        Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.75, yCenter));
        Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.3, yCenter));
        if type(1) == 'a'
            if r == 0
                
                f = type(2);
                if strcmp(f,'1')
                    Screen('FrameRect', window, [1 0 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                elseif strcmp(f,'2')
                    Screen('FrameRect', window, [0 1 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                    
                end
                
            else
                
                f = type(2);
                if strcmp(f,'1')
                    Screen('FrameRect', window, [1 0 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                elseif strcmp(f,'2')
                    Screen('FrameRect', window, [0 1 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                end
            end
        elseif type(1) == 'b'
            %if text (b) is the domain
            %and if the next trial text should be on the right
            if r == 1
                f = type(2);
                if strcmp(f,'1')
                    Screen('FrameRect', window, [0 0 1], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                elseif strcmp(f,'2')
                    Screen('FrameRect', window, [225 225 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                    
                end
            else
                
                f = type(2);
                if strcmp(f,'1')
                    Screen('FrameRect', window, [0 0 1], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                elseif strcmp(f,'2')
                    Screen('FrameRect', window, [225 225 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                    
                end
            end
        else
            
            if r == 1
                Screen('FrameRect', window, black, CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                
            else
                Screen('FrameRect', window, black, CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                
            end
            
        end
        
        [~, frame_onset] = Screen('Flip', window, trialstart+1.5);
        result(trial).frame_onset = frame_onset;
        
        Screen('DrawTexture', window,which, [], CenterRectOnPointd([0 0 screenXpixels screenYpixels], xCenter, yCenter), 0);
        Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.75, yCenter));
        Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.3, yCenter));
        if type(1) == 'a'
            if r == 0
                
                f = type(2);
                if strcmp(f,'1')
                    Screen('FrameRect', window, [1 0 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                elseif strcmp(f,'2')
                    Screen('FrameRect', window, [0 1 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                    
                end
                
            else
                
                f = type(2);
                if strcmp(f,'1')
                    Screen('FrameRect', window, [1 0 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                elseif strcmp(f,'2')
                    Screen('FrameRect', window, [0 1 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                end
            end
        elseif type(1) == 'b'
            %if text (b) is the domain
            %and if the next trial text should be on the right
            if r == 1
                f = type(2);
                if strcmp(f,'1')
                    Screen('FrameRect', window, [0 0 1], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                elseif strcmp(f,'2')
                    Screen('FrameRect', window, [225 225 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                    
                end
            else
                
                f = type(2);
                if strcmp(f,'1')
                    Screen('FrameRect', window, [0 0 1], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                elseif strcmp(f,'2')
                    Screen('FrameRect', window, [225 225 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                    
                end
            end
        else
            
            if r == 1
                Screen('FrameRect', window, black, CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                
            else
                Screen('FrameRect', window, black, CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                
            end
            
        end
        
        [~,context_onset] = Screen('Flip', window,trialstart + 1.5 + 0.5);
        result(trial).context_onset = context_onset;
        
        
        if strcmp(switch_seq{trial}, 'context_switch')
            while GetSecs < (context_onset + 0.5)
                [keyIsDown, secs, keyCode] = KbCheck(-1);
                [secs, ~, daqstate]=scansync(2,0.001);
                if keyCode(KbName('space')) == 1 || daqstate.previousflags(2)
                    result(trial).response = 'context';
                    result(trial).rt = secs;
                    respMade=1;
                    %break;
                    %                 elseif keyCode(KbName('ESCAPE')) == 1
                    %                     sca;
                    %                     disp('*** Experiment terminated ***');
                    %                     return
                end
            end
            
            %         else
            %             WaitSecs(0.5);
        end
        
        if type(1) == 'a'
            
            if type(2) == '1'
                %if required answer is yes to question 'is it
                %living?'
                if answer
                    %random choose from living pool
                    if randi([1,2],1) == 1
                        pool = objects.living_large;
                        result(trial).stimuli_type = 'living_large';
                    else
                        pool = objects.living_small;
                        result(trial).stimuli_type = 'living_small';
                    end
                    
                else
                    
                    if randi([1,2],1) == 1
                        pool = objects.nonliving_large;
                        result(trial).stimuli_type = 'nonliving_large';
                    else
                        pool = objects.nonliving_small;
                        result(trial).stimuli_type = 'nonliving_small';
                    end
                end
            else
                %if required answer is yes to question 'does it fit into a
                %shoebox?'
                %TO-DO: change everything
                
                if answer==0
                    
                    if randi([1,2],1) == 1
                        pool = objects.nonliving_large;
                        result(trial).stimuli_type = 'nonliving_large';
                    else
                        pool = objects.living_large;
                        result(trial).stimuli_type = 'living_large';
                    end
                    
                else
                    if randi([1,2],1) == 1
                        pool = objects.living_small;
                        result(trial).stimuli_type = 'living_small';
                    else
                        pool = objects.nonliving_small;
                        result(trial).stimuli_type = 'nonliving_small';
                    end
                end
            end
            
            %object_to_show is randomly picked from pool
            object_to_show = randsample(pool,1);
            object_index = find(pool == object_to_show);
            result(trial).stimuli = object_index;
            
            %make texture here?
            
            %text is random
            random_text_type = randi([1,4],1);
            result(trial).irrelevant_stim_type = text_type_names{random_text_type};
            
            textpool = {LettersA LettersI LettersBoth LettersNeither};
            text_to_show = randsample(textpool{random_text_type},1);
        elseif type(1) == 'b'
            if type(2) == '1'
                %if required answer is yes to question 'does A fit in
                %to make a word?'
                
                if answer
                    if randi([1,2],1) ==1
                        pool = LettersA;
                        result(trial).stimuli_type = 'LettersA';
                    else
                        pool = LettersBoth;
                        result(trial).stimuli_type = 'LettersBoth';
                    end
                    %random choose from this pool
                else
                    if randi([1,2],1) ==1
                        pool = LettersI;
                        result(trial).stimuli_type = 'LettersI';
                    else
                        pool = LettersNeither;
                        result(trial).stimuli_type = 'LettersNeither';
                    end
                    
                end
            else
                %if required answer is yes to question 'Does I fit in
                %to make a word?'
                if answer
                    if randi([1,2],1) ==1
                        pool = LettersI;
                        result(trial).stimuli_type = 'LettersI';
                    else
                        pool = LettersBoth;
                        result(trial).stimuli_type = 'LettersBoth';
                    end
                    
                else
                    if randi([1,2],1) ==1
                        pool = LettersA;
                        result(trial).stimuli_type = 'LettersA';
                    else
                        pool = LettersNeither;
                        result(trial).stimuli_type = 'LettersNeither';
                    end
                    
                end
            end
            
            %object is random
            %TO-DO: record which irrelvant sample was selected
            %pick text from pool
            text_to_show = randsample(pool,1);
            result(trial).stimuli = text_to_show;
            
            random_obj_type = randi([1,4],1);
            %put before loop
            result(trial).irrelevant_stim_type = obj_type_names{random_obj_type};
            
            objectpool = {objects.living_small,objects.nonliving_small,objects.living_large,objects.nonliving_large};
            object_to_show = randsample(objectpool{random_obj_type},1);
            
            %if trial is rest, randomly pick text and object
        else
            %text_to_show = randsample(pool,1);
            %object_to_show = randsample(pool,1);
            random_obj_type = randi([1,4],1);
            
            %result(trial).rest_irrelevant_obj_type = obj_type_names{random_obj_type};
            
            objectpool = {objects.living_small,objects.nonliving_small,objects.living_large,objects.nonliving_large};
            object_to_show = randsample(objectpool{random_obj_type},1);
            
            random_text_type = randi([1,4],1);
            %result(trial).rest_irrelevant_text_type = text_type_names{random_text_type};
            
            textpool = {LettersA LettersI LettersBoth LettersNeither};
            text_to_show = randsample(textpool{random_text_type},1);
            
        end
        rest_rand = randi([0,1],1);
        
        
        Screen('DrawTexture', window,which, [], CenterRectOnPointd([0 0 screenXpixels screenYpixels], xCenter, yCenter), 0);
        Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.75, yCenter));
        Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.3, yCenter));
        %checks which stimuli to put, and which side
        
        if type(1) == 'a'
            %if object (a) is the domain
            %and if the next trial object should be on the left
            if r == 0
                Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.75, yCenter));
                Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.3, yCenter));
                Screen('DrawTexture', window, object_to_show, [], CenterRectOnPointd([0 0 300 300], screenXpixels * 0.3, yCenter), 0);
                %Screen('FillRect', window, white, CenterRectOnPointd([0 0 200 100], screenXpixels * 0.75 + 100, yCenter));
                DrawFormattedText(window, text_to_show{1}, screenXpixels * 0.7, 'center', black);
                
                f = type(2);
                if strcmp(f,'1')
                    Screen('FrameRect', window, [1 0 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                elseif strcmp(f,'2')
                    Screen('FrameRect', window, [0 1 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                    
                end
                
                
                
                % if next trial should be on the right
            else
                Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.75, yCenter));
                Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.3, yCenter));
                Screen('DrawTexture', window, object_to_show, [], CenterRectOnPointd([0 0 300 300], screenXpixels * 0.75, yCenter), 0);
                
                %Screen('FillRect', window, white, CenterRectOnPointd([0 0 200 100], screenXpixels * 0.25  + 100, yCenter));
                DrawFormattedText(window, text_to_show{1}, screenXpixels * 0.25, 'center', black);
                
                f = type(2);
                if strcmp(f,'1')
                    Screen('FrameRect', window, [1 0 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                elseif strcmp(f,'2')
                    Screen('FrameRect', window, [0 1 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                    
                end
                
                
            end
        elseif type(1) == 'b'
            
            
            %if text (b) is the domain
            %and if the next trial text should be on the right
            if r == 1
                Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.75, yCenter));
                Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.3, yCenter));
                Screen('DrawTexture', window, object_to_show, [], CenterRectOnPointd([0 0 300 300], screenXpixels * 0.3, yCenter), 0);
                %Screen('FillRect', window, white, CenterRectOnPointd([0 0 200 100], screenXpixels * 0.75 + 100, yCenter));
                DrawFormattedText(window, text_to_show{1}, screenXpixels * 0.7, 'center', black);
                
                f = type(2);
                if strcmp(f,'1')
                    Screen('FrameRect', window, [0 0 1], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                elseif strcmp(f,'2')
                    Screen('FrameRect', window, [225 225 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                    
                end
                
                
                
                % if next trial should be on the left
            else
                Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.75, yCenter));
                Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.3, yCenter));
                Screen('DrawTexture', window, object_to_show, [], CenterRectOnPointd([0 0 300 300], screenXpixels * 0.75, yCenter), 0);
                
                %Screen('FillRect', window, white, CenterRectOnPointd([0 0 200 100], screenXpixels * 0.25  + 100, yCenter));
                DrawFormattedText(window,text_to_show{1}, screenXpixels * 0.25, 'center', black);
                
                f = type(2);
                if strcmp(f,'1')
                    Screen('FrameRect', window, [0 0 1], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                elseif strcmp(f,'2')
                    Screen('FrameRect', window, [225 225 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                    
                end
                %vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                
            end
        else
            if rest_rand
                Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.75, yCenter));
                Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.3, yCenter));
                Screen('DrawTexture', window, object_to_show, [], CenterRectOnPointd([0 0 300 300], screenXpixels * 0.75, yCenter), 0);
                
                %Screen('FillRect', window, white, CenterRectOnPointd([0 0 200 100], screenXpixels * 0.25  + 100, yCenter));
                DrawFormattedText(window,text_to_show{1}, screenXpixels * 0.25, 'center', black);
                if r == 1
                    Screen('FrameRect', window, black, CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                    result(trial).irrelevant_stim_type = text_type_names{random_text_type};
                    result(trial).stimuli_type = obj_type_names{random_obj_type};;
                else
                    Screen('FrameRect', window, black, CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                    result(trial).irrelevant_stim_type = obj_type_names{random_obj_type};;
                    result(trial).stimuli_type = text_type_names{random_text_type};
                end
            else
                Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.75, yCenter));
                Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.3, yCenter));
                Screen('DrawTexture', window, object_to_show, [], CenterRectOnPointd([0 0 300 300], screenXpixels * 0.3, yCenter), 0);
                
                %Screen('FillRect', window, white, CenterRectOnPointd([0 0 200 100], screenXpixels * 0.25  + 100, yCenter));
                DrawFormattedText(window,text_to_show{1}, screenXpixels * 0.7, 'center', black);
                if r == 1
                    Screen('FrameRect', window, black, CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                    result(trial).irrelevant_stim_type = obj_type_names{random_obj_type};;
                    result(trial).stimuli_type = text_type_names{random_text_type};
                else
                    Screen('FrameRect', window, black, CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                    result(trial).irrelevant_stim_type = text_type_names{random_text_type};
                    result(trial).stimuli_type = obj_type_names{random_obj_type};;
                end
                
            end
            
        end
        
        
        %TO-DO: check
        %[vbl,restart] = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        [~, stim_onset]= Screen('Flip', window,trialstart + 1.5 + 0.5 + 0.5);
        result(trial).stim_onset = stim_onset;
        
        while respMade == 0
            respMade = 0;
            
            %TO-DO: record if participant has pressed any buttons during
            %rest
            if type(1)=='r'
                
                if length(rt_avg) >1
                    rest_avg = mean(rt_avg);
                end
                while GetSecs < (stim_onset + rest_avg)
                    [keyIsDown, secs, keyCode] = KbCheck(-1);
                    [secs, ~, daqstate]=scansync(5,0.001);
                    if keyIsDown || isempty(daqstate.previousflags)==0
                        if keyCode(KbName('space')) == 1 || daqstate.previousflags(5)
                            respMade = 1;
                            result(trial).response = 'context';
                            result(trial).accuracy = 1;
                            result(trial).context_response = 1;
                            result(trial).rt = secs;
                            
                            if strcmp(switch_seq{trial},'context_switch')
                                
                                break;
                            end
                            result(trial).context_response = 1; %true if participants still responded during rest
                        end
                        
                    end
                end
                break;
                %after showing for 2 seconds, while loop should break, should
                %enter next for trial loop
                %if rest get out of bigger while loop
                
                
            elseif strcmp(switch_seq{trial},'context_switch')
                
                %while GetSecs < (onset_time + 2)
                [keyIsDown, secs, keyCode] = KbCheck(-1);
                [secs, ~, daqstate]=scansync(2:5,inf);
                
                if keyIsDown || isempty(daqstate.previousflags)==0
                    respMade = 1;
                    save(fullfile(here,resultname),'result');
                    if keyCode(KbName('space')) == 1|| daqstate.previousflags(5)
                        result(trial).rt = secs(4);
                        result(trial).response = 'context';
                        result(trial).accuracy = 1;
                        result(trial).context_response = 1; %true if participants still responded during rest
                    elseif keyCode(leftKey) == 1 || daqstate.previousflags(3)
                        result(trial).rt = secs(3);
                        result(trial).response = 'yes';
                        if answer ==1
                            result(trial).accuracy =1;
                        else
                            result(trial).accuracy =0;
                        end
                    elseif keyCode(rightKey) == 1 || daqstate.previousflags(2)
                        result(trial).rt = secs(2);
                        result(trial).response = 'no';
                        if answer ==0
                            result(trial).accuracy =1;
                        else
                            result(trial).accuracy =0;
                        end
                    elseif keyCode(KbName('ESCAPE')) == 1
                        resultname = ['exp_pilot_' subinfo{1} '_run_' subinfo{2} '.mat'];
                        save(fullfile(here,resultname),'result');
                        sca;
                        disp('*** Experiment terminated ***');
                        return
                    end
                end
            else
                %end
                %break;
                
                
                % Poll the keyboard for keys
                [keyIsDown, secs, keyCode] = KbCheck(-1);
                [secs, ~, daqstate]=scansync(2:5,inf);
                
                result(trial).rt = GetSecs;
                if keyCode(KbName('space')) == 1 || daqstate.previousflags(5)
                    %respMade = 1;
                    result(trial).rt = secs(4);
                    result(trial).context_response = 1; %true if participants still responded during rest
                end
                if keyCode(leftKey) == 1 || daqstate.previousflags(3)
                    result(trial).rt = secs(2);
                    to_add = secs(2)-stim_onset;
                    rt_avg = [rt_avg to_add];
                    respMade = 1;
                    result(trial).response = 'yes';
                    if answer ==1
                        result(trial).accuracy =1;
                    else
                        result(trial).accuracy =0;
                    end
                elseif keyCode(rightKey) == 1 || daqstate.previousflags(2)
                    respMade = 1;
                    result(trial).rt = secs(1);
                    to_add = secs(1)-stim_onset;
                    rt_avg = [rt_avg to_add];
                    result(trial).response = 'no';
                    if answer ==0
                        result(trial).accuracy =1;
                    else
                        result(trial).accuracy =0;
                    end
                elseif keyCode(KbName('ESCAPE')) == 1 %|| daqstate.previousflags(4)
                    resultname = ['exp_pilot_' subinfo{1} '_run_' subinfo{2} '.mat'];
                    save(resultname,'result');
                    sca;
                    disp('*** Experiment terminated ***');
                    return
                end
            end
        end
        %TO-DO: judgement of results, accuracy, trial info etc.
        
    end
    save(fullfile(here,resultname),'result');
    
    % Set the text size
    Screen('TextSize', window, 30);
    
    Screen('FillRect', window, grey);
    DrawFormattedText(window, 'Well done! That was the end of the run, please wait.', 'center', 'center', black);
    Screen('Flip', window);
    WaitSecs(10);
    
    % Close the onscreen window
    sca
    
catch e
    % This part of the code will be executed even if there is an error
    % above in the code.
    % Close the onscreen window
    sca
    % This is to display the caught error in the matlab command window, so
    % you'll see what was wrong.
    rethrow(e);
end

