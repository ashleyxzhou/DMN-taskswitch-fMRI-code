function pilot_practice

dbstop if error
commandwindow;

TR=1.208;
scansync('reset', TR)

here=fileparts(mfilename('fullpath'));
%one level up
stimuli_folder_dir=fullfile(fileparts(here),'stimuli');
results_folder_dir=fullfile(fileparts(here),'Results');

subinfo = inputdlg({'Subject ID:'});
%set up .mat var to record data
resultname = ['practice_sub' subinfo{1} '.mat'];

while exist(resultname,'file')
    subinfo = inputdlg({'Subject ID:'});
    resultname = ['practice_sub' subinfo{1} '.mat'];
end

result={};

save(fullfile(results_folder_dir,resultname),'result');

% Setup PTB with some default values
PsychDefaultSetup(2);

% Seed the random number generator.
%id x 10 + run
rng(str2double(subinfo{1})*10);
[series, switch_seq] = random_seq_gen();

% Skip sync tests for demo purposes only
Screen('Preference', 'SkipSyncTests', 1);


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
Screen('TextSize', window,30);

% Query the maximum priority level
topPriorityLevel = MaxPriority(window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Set the blend funciton for the screen
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

DrawFormattedText(window, 'Please wait while experiment is loading, usually takes 10 seconds...', 'center', 'center', black);
Screen('Flip', window);
%WaitSecs(10);



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

rootdir = stimuli_folder_dir;

% Get the image files for the experiment
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
    
end

text_type_names = {'LettersA' 'LettersI' 'LettersBoth' 'LettersNeither'};
obj_type_names = {'living_small' 'nonliving_small' 'living_large' 'nonliving_large'};

beach_imageFolder = [stimuli_folder_dir '\dummy_beach\'];
cafe_imageFolder =[stimuli_folder_dir '\dummy_cafe\'];
jungle_imageFolder =[stimuli_folder_dir '\scenes_jungle\'];


jcontext_imgList = dir(fullfile(jungle_imageFolder, '*.jpg'));
jcontext_imgList = {jcontext_imgList(:).name};

all_jungle =[];
for index = 1:18
    contextName = jcontext_imgList{index};
    jungle_image = imread([jungle_imageFolder contextName]);
    jungle_tex = Screen('MakeTexture', window, jungle_image);
    all_jungle = [all_jungle jungle_tex];
    clear jungle_image;
end
all_jungle = Shuffle(all_jungle);

imgList = dir(fullfile(cafe_imageFolder, '*.jpg'));
imgList = {imgList(:).name};
numImages = length(imgList);
all_cafe = nan(1,numImages);
for i =1:numImages
    imageName = imgList{i};
    theImage = imread([stimuli_folder_dir '\dummy_cafe\' imageName]);
    all_cafe(i) = Screen('MakeTexture', window, theImage);
end

imgList = dir(fullfile(beach_imageFolder, '*.jpg'));
imgList = {imgList(:).name};
numImages = length(imgList);
all_beach = nan(1,numImages);
for i =1:numImages
    imageName = imgList{i};
    theImage = imread([stimuli_folder_dir '\dummy_beach\' imageName]);
    all_beach(i) = Screen('MakeTexture', window, theImage);
end

all_beach = Shuffle(all_beach);
all_cafe = Shuffle(all_cafe);

fileid{1} = fopen(fullfile(rootdir, '\Letters\Both4.txt'), 'r');
fileid{2} = fopen(fullfile(rootdir, '\Letters\A_only4.txt'), 'r');
fileid{3} = fopen(fullfile(rootdir, '\Letters\I_only4.txt'), 'r');
fileid{4} = fopen(fullfile(rootdir, '\Letters\Neither4.txt'), 'r');

for n = 1:4
    data = textscan(fileid{n},'%s');
    
    %close file
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

[mask,series] = context_switch_gen_practice(series,switch_seq);
mask=[' ' mask];
series = ['endblock' series];
types = {'a1' 'a2' 'b1' 'b2'};
for j = 1:4
    for i = 1:10
        series = [types{j} series];
        if i<6
            mask =[0 mask];
        else
            mask = [1 mask];
        end
    end
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
    DrawFormattedText(window, 'Use right hand index finger for yes, middle finger for no.\nPress Space To Begin. ', 'center', 'center', black);
    Screen('Flip', window);
    KbWait;
    percentage = 0;
    
    Screen('TextSize', window,60);
    while percentage<0.8
        %stores the sequence of trials
        
        %a1, b1, a2, b2,
        %is it living, 'A', is it larger than shoebox, 'I'
        total_accuracy = [];
        
        for trial = 1:length(series)
            
            %saves trial info
            type = series{trial};
            result(trial).type = series{trial};
            
            if strcmp(type,'endblock')
                Screen('TextSize', window,30);
                DrawFormattedText(window, 'Well done! There is an extra rule: Press the leftmost key when you see a jungle in the background.\nYou do not need to respond to the tasks on the jungle trial.\nPress any key to continue.', 'center', 'center', black);
                Screen('Flip', window);
                KbWait;
                Screen('TextSize', window,60);
            else
                
                %Screen('DrawTexture', window,which_scene, [], CenterRectOnPointd([0 0 screenXpixels screenYpixels], xCenter, yCenter), 0);
                % Draw a fixation cross for the start of the trial
                Screen('FillRect', window, grey);
                Screen('DrawLines', window, allCoords,...
                    lineWidthPix, white, [xCenter yCenter], 2);
                Screen('Flip', window);
                WaitSecs(1.5);
                % Draw the fixation cross in white, set it to the center of our screen
                
                %randomly pick which_scene side of screen to put picture
                r = randi([0,1], 1);
                answer  = randi([0,1], 1);
                
                respMade = 0;
                Priority(topPriorityLevel);

                %change dummy background to jungle? + dummy trial with scene that
                %fits with next sequence of backgrounds
                
                if mask(trial)==2
                    which_scene = randsample(all_jungle,1);
                elseif mask(trial)==1
                    
                    %checks which_scene picture to show for each condition
                    
                    which_scene = randsample(all_beach,1);
                    result(trial).scene = 'beach';
                    
                else
                    which_scene = randsample(all_cafe,1);
                    result(trial).scene = 'cafe';
                    
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
                
                Screen('Flip', window);
                result(trial).frame_onset = GetSecs;
                WaitSecs(0.5);
                
                
                Screen('DrawTexture', window,which_scene, [], CenterRectOnPointd([0 0 screenXpixels screenYpixels], xCenter, yCenter), 0);
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
                
                Screen('Flip', window);
                result(trial).context_onset = GetSecs;
                onset_time = GetSecs;

                if mask(trial)==2
                    while GetSecs < (onset_time + 0.5)
                        [keyIsDown, secs, keyCode] = KbCheck(-1);
                        [secs, ~, daqstate]=scansync(2,0.001);
                        if keyCode(KbName('space')) == 1 || daqstate.previousflags(2)
                            result(trial).response = 'space';
                            result(trial).rt = secs;
                            respMade=1;
                           
                        end
                    end
                    
                else
                    WaitSecs(0.5);
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
                        %if required answer is yes to question 'does it fit
                        %into a shoebox?'
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
                    %TO-DO: record which_scene irrelvant sample was selected
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
                
                result(trial).stim_onset = GetSecs;
                
                while respMade == 0
                    respMade = 0;
                    
                    Screen('DrawTexture', window,which_scene, [], CenterRectOnPointd([0 0 screenXpixels screenYpixels], xCenter, yCenter), 0);
                    %checks which_scene stimuli to put, and which_scene side
                    Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.75, yCenter));
                    Screen('FillRect', window, white, CenterRectOnPointd([0 0 400 400], screenXpixels * 0.3, yCenter));
                    if type(1) == 'a'
                        %if object (a) is the domain
                        %and if the next trial object should be on the left
                        if r == 0
                            
                            Screen('DrawTexture', window, object_to_show, [], CenterRectOnPointd([0 0 300 300], screenXpixels * 0.3, yCenter), 0);
                            DrawFormattedText(window, text_to_show{1}, screenXpixels * 0.7, 'center', black);
                            
                            f = type(2);
                            if strcmp(f,'1')
                                Screen('FrameRect', window, [1 0 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                            elseif strcmp(f,'2')
                                Screen('FrameRect', window, [0 1 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                                
                            end
                            
                            
                            
                            % if next trial should be on the right
                        else
                            
                            Screen('DrawTexture', window, object_to_show, [], CenterRectOnPointd([0 0 300 300], screenXpixels * 0.75, yCenter), 0);
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
                            
                            Screen('DrawTexture', window, object_to_show, [], CenterRectOnPointd([0 0 300 300], screenXpixels * 0.3, yCenter), 0);
                           
                            DrawFormattedText(window, text_to_show{1}, screenXpixels * 0.7, 'center', black);
                            
                            f = type(2);
                            if strcmp(f,'1')
                                Screen('FrameRect', window, [0 0 1], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                            elseif strcmp(f,'2')
                                Screen('FrameRect', window, [225 225 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.75, yCenter), 8);
                                
                            end
                            
                            
                            
                            % if next trial should be on the left
                        else
                            
                            Screen('DrawTexture', window, object_to_show, [], CenterRectOnPointd([0 0 300 300], screenXpixels * 0.75, yCenter), 0);
                           
                            DrawFormattedText(window,text_to_show{1}, screenXpixels * 0.25, 'center', black);
                            
                            f = type(2);
                            if strcmp(f,'1')
                                Screen('FrameRect', window, [0 0 1], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                            elseif strcmp(f,'2')
                                Screen('FrameRect', window, [225 225 0], CenterRectOnPointd([0 0 350 350], screenXpixels * 0.3, yCenter), 8);
                                
                            end

                        end
                    else
                        if rest_rand
                            
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
                    
                    
                   
                    Screen('Flip', window);
                    onset_time = GetSecs;
                    

                    if type(1)=='r'
                        
                        while GetSecs < (onset_time + 2)
                            [keyIsDown, secs, keyCode] = KbCheck(-1);
                            %TODO
                            [secs, ~, daqstate]=scansync(5,0.001);
                            if keyIsDown  || isempty(daqstate.previousflags)==0
                                if keyCode(KbName('space')) == 1 || daqstate.previousflags(5)
                                    respMade = 1;
                                    
                                    if mask(trial)==2
                                        
                                        break;
                                    end
                                    result(trial).response = 'space'; %true if participants still responded during rest
                                end
                                result(trial).subres = 1; %true if participants still responded during rest
                                result(trial).rt = GetSecs;
                            end
                        end
                        break;
                        %after showing for 2 seconds, while loop should break, should
                        %enter next for trial loop
                        %if rest get out of bigger while loop
                        
                        
                    elseif mask(trial)==2
                        
                        %while GetSecs < (onset_time + 2)
                        [keyIsDown, secs, keyCode] = KbCheck(-1);
                        %TODO
                        [secs, ~, daqstate]=scansync(2:5,inf);
                        result(trial).rt = secs;
                        if keyIsDown || isempty(daqstate.previousflags)==0
                            respMade = 1;
                            if keyCode(KbName('space')) == 1 || daqstate.previousflags(5)
                                result(trial).context_response = 1; %true if participants still responded during rest
                                DrawFormattedText(window, 'Correct!', 'center', 'center', black);
                                Screen('Flip', window);
                                WaitSecs(1);
                                %total_accuracy = [total_accuracy 1];
                            elseif keyCode(leftKey) == 1 || daqstate.previousflags(3)
                                result(trial).response = 'yes';
                                DrawFormattedText(window, 'That was a jungle!', 'center', 'center', black);
                                Screen('Flip', window);
                                WaitSecs(1);
                            elseif keyCode(rightKey) == 1 || daqstate.previousflags(2)
                                result(trial).response = 'no';
                                DrawFormattedText(window, 'That was a jungle!', 'center', 'center', black);
                                Screen('Flip', window);
                                WaitSecs(1);
                              
                            %TODO: take off daqstate
                            elseif keyCode(KbName('ESCAPE')) == 1 || daqstate.previousflags(4)
                                resultname = ['practice_' subinfo{1} '_run_' subinfo{2} '.mat'];
                                save(resultname,'result');
                                sca;
                                disp('*** Experiment terminated ***');
                                return
                            end
                        end
                    end
                    %end
                    %break;
                    
                    
                    % Poll the keyboard for keys
                    [keyIsDown, secs, keyCode] = KbCheck(-1);
                    
                    %TODO
                    [secs, ~, daqstate]=scansync(2:5,inf);
                    result(trial).rt = secs;
                    if keyCode(KbName('space')) == 1 || daqstate.previousflags(5)
                        %respMade = 1;
                        result(trial).context_response = 1; % true if participants still responded during rest
                    end
                    if keyCode(leftKey) == 1 || daqstate.previousflags(3)
                        respMade = 1;
                        result(trial).response = 'yes';
                        if answer ==1
                            result(trial).accuracy =1;
                            DrawFormattedText(window, 'Correct!', 'center', 'center', black);
                            Screen('Flip', window);
                            WaitSecs(1);
                            
                            total_accuracy = [total_accuracy 1];
                        else
                            result(trial).accuracy =0;
                            DrawFormattedText(window, 'Incorrect', 'center', 'center', black);
                            Screen('Flip', window);
                            WaitSecs(1);
                            
                            total_accuracy = [total_accuracy 0];
                        end
                    elseif keyCode(rightKey) == 1  || daqstate.previousflags(2)
                        respMade = 1;
                        result(trial).response = 'no';
                        if answer ==0
                            result(trial).accuracy =1;
                            DrawFormattedText(window, 'Correct!', 'center', 'center', black);
                            Screen('Flip', window);
                            WaitSecs(1);
                            
                            total_accuracy = [total_accuracy 1];
                        else
                            result(trial).accuracy =0;
                            DrawFormattedText(window, 'Incorrect', 'center', 'center', black);
                            Screen('Flip', window);
                            WaitSecs(1);
                            
                            total_accuracy = [total_accuracy 0];
                        end
                    elseif keyCode(KbName('ESCAPE')) == 1
                        resultname = ['practice_sub' subinfo{1} '.mat'];
                        save(resultname,'result');
                        sca;
                        disp('*** Experiment terminated ***');
                        return
                    end
                end
                %TO-DO: judgement of results, accuracy, trial info etc.
            end
            save(fullfile(results_folder_dir,resultname),'result');
        end
        
        percentage = mean(total_accuracy);
        end_screen = ['Percentage correct: ' num2str(round(percentage*100),2) '%'];
        
        DrawFormattedText(window, end_screen, 'center', 'center', black);
        Screen('Flip', window);
        WaitSecs(3);
        
        Screen('TextSize', window,30);
        DrawFormattedText(window, 'In the actual experiment, you will not be provided feedback.', 'center', 'center', black);
        Screen('Flip', window);
        WaitSecs(3);
    end
    %save(resultname,'result');
    
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

