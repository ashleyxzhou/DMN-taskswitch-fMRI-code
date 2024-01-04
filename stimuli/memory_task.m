function memory_task

%TODO: save reaction times and selected/sort order, accuracy, pic number
dbstop if error
commandwindow;

subinfo = inputdlg({'Subject ID:'});
%set up .mat var to record data
resultname = ['pilotsub_' subinfo{1} 'mem_task.mat'];
%
% while exist(resultname,'file')
%     subinfo = inputdlg({'Subject ID:'});
%     resultname = ['pilotsub_' subinfo{1} 'mem_task.mat'];
% end

result={};

% save(resultname,'result');

PsychDefaultSetup(2);
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
%[window, windowRect] = PsychImaging('OpenWindow', 1, grey, [], 32, 2);
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [], 32, 2);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Flip to clear
Screen('Flip', window);

% Set the text size
Screen('TextSize', window, 30);


% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Set the blend funciton for the screen
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

DrawFormattedText(window, 'Please wait while experiment is loading...', 'center', 'center', black);
Screen('Flip', window);

here=fileparts(mfilename('fullpath'));
%one level up
stimuli_folder_dir=fullfile(fileparts(here),'stimuli');

run_context_cafe=zeros(4,6);
run_context_beach=zeros(4,6);

run_context_cafe_tex=zeros(4,6);
run_context_beach_tex=zeros(4,6);

novel_cafe_matrix = zeros(4,6);
novel_beach_matrix = zeros(4,6);
new_cafe_matrix = zeros(8,6);
new_beach_matrix = zeros(8,6);

novel_beach_imageFolder = [stimuli_folder_dir '\run5\scenes_beach\'];
novel_cafe_imageFolder =[stimuli_folder_dir '\run5\scenes_cafe\'];

all_novel_beach = nan(1,48);
all_novel_cafe = nan(1,48);
for num = 1:48
    if num<10
        contextName = ['0' num2str(num) '.jpg'];
    else
        contextName = [num2str(num) '.jpg'];
    end
    all_novel_cafe(num)= Screen('MakeTexture', window,imread([novel_cafe_imageFolder contextName]));
    clear cafe_image;
    all_novel_beach(num)= Screen('MakeTexture', window,imread([novel_beach_imageFolder contextName]));
    clear beach_image;
end

rand_beach = Shuffle(all_novel_beach);
rand_cafe = Shuffle(all_novel_cafe);

for c = 1:6
    new_beach_matrix(1:8,c) = rand_beach (1:8);
    rand_beach(1:8) = [];
    
    new_cafe_matrix(1:8,c) = rand_cafe(1:8);
    rand_cafe(1:8)=[];
end

%create a matrix where each row represents a run, each column represents a
%condition
cond_names = {'task-stay','between-domain','within-domain','rest','extended-rest','restart'};

for run = 1:4
    load(['sub' subinfo{1} '_run_' num2str(run) 'repeatCafe.mat'], 'array_repeated_pics_cafe');
    rep_cafe=array_repeated_pics_cafe;
    run_context_cafe(run,:)=rep_cafe(1:6);
    %result(1).correct_repeated_order = rep_cafe(1:6);
    run_novel_cafe(run,:)=rep_cafe(7:12);
    
    load(['sub' subinfo{1} '_run_' num2str(run) 'repeatBeach.mat'], 'array_repeated_pics_beach');
    rep_beach=array_repeated_pics_beach;
    run_context_beach(run,:)=rep_beach(1:6);
    run_novel_beach(run,:)=rep_beach(7:12);
    
    beach_imageFolder = [stimuli_folder_dir '\run' num2str(run) '\scenes_beach\'];
    cafe_imageFolder =[stimuli_folder_dir '\run' num2str(run) '\scenes_cafe\'];
    
    for pic = 1:6
        pic_index = rep_cafe(pic);
        
        if pic_index<10
            contextName = ['0' num2str(pic_index) '.jpg'];
        else
            contextName = [num2str(pic_index) '.jpg'];
        end
        
        run_context_cafe_tex(run,pic) = Screen('MakeTexture', window, imread([cafe_imageFolder contextName]));
        clear cafe_image;
        
        pic_index = rep_cafe(pic+6);
        
        if pic_index<10
            contextName = ['0' num2str(pic_index) '.jpg'];
        else
            contextName = [num2str(pic_index) '.jpg'];
        end
        novel_cafe_matrix(run,pic) =  Screen('MakeTexture', window, imread([cafe_imageFolder contextName]));
        
        pic_index = rep_beach(pic);
        
        if pic_index<10
            contextName = ['0' num2str(pic_index) '.jpg'];
        else
            contextName = [num2str(pic_index) '.jpg'];
        end
        run_context_beach_tex(run,pic) = Screen('MakeTexture', window, imread([beach_imageFolder contextName]));
        
        pic_index = rep_beach(pic+6);
        
        if pic_index<10
            contextName = ['0' num2str(pic_index) '.jpg'];
        else
            contextName = [num2str(pic_index) '.jpg'];
        end
        novel_beach_matrix(run,pic) = Screen('MakeTexture', window, imread([beach_imageFolder contextName]));
        
        clear beach_image;
        
        
    end
end

try
    oldtype=ShowCursor(0);
    
    cond_array = Shuffle(1:6);
    cond2 = flip(cond_array)+6;
    cond_array=[cond_array cond2];
    
    %For the cafe context, each of the six conditions
    for index = 1: 12
        cond = cond_array(index);
        if cond<7
            result(index).condition = cond_names(cond);
        else
            result(index).condition = cond_names(cond-6);
        end
        Screen('FillRect', window, white);
        DrawFormattedText(window, 'Click on 8 pictures that you think you have seen in the experiment.\n\nPress Space To Begin.', 'center', 'center', black);
        Screen('Flip', window);
        KbWait;
        
        if cond<7
            pics_to_show = [novel_cafe_matrix(1:4,cond)' new_cafe_matrix(1:8,cond)' run_context_cafe_tex(1:4,cond)'];
            result(index).novel_pics = run_novel_cafe(:,cond)';
            result(index).repeat_pics = run_context_cafe(:,cond)';
        else
            pics_to_show = [novel_beach_matrix(1:4,cond-6)' new_beach_matrix(1:8,cond-6)' run_context_beach_tex(1:4,cond-6)'];
            result(index).novel_pics = run_novel_beach(:,cond-6)';
            result(index).repeat_pics = run_context_beach(:,cond-6)';
        end
        sort_index = Shuffle(1:16);
        result(index).correct_novel_position = sort_index(1:4);
        result(index).correct_repeat_position = sort_index(13:16);
        to_show(sort_index)=pics_to_show;
        
        Screen('FillRect', window, white);

        allRects = nan(16,4);
        for ycor = 1:4
            for xcor = 1:4
                allRects(xcor+(ycor-1)*4,:) = CenterRectOnPointd([0 0 300 200], screenXpixels*0.2*xcor, screenYpixels*0.2*ycor);
                Screen('DrawTexture', window, to_show(xcor+(ycor-1)*4), [], allRects(xcor+(ycor-1)*4,:), 0);

            end
        end

        % Draw the rect to the screen
        [~, result(index).select_onset] = Screen('Flip', window, [], 1);
        
        result(index).select_rts =[];
        result(index).select_pics=[];
        mouseclicks=0;
        selected= zeros(1,16); %checks whether the picture is already selected
         [keyIsDown, secs, keyCode] = KbCheck(-1);
        if keyCode(KbName('ESCAPE')) == 1
            
%             save(resultname,'result');
            sca;
            disp('*** Experiment terminated ***');
            return
        end
        
        while mouseclicks<8
            [x,y,buttons]=GetMouse;
            
            while ~any(buttons)
                %waits for a key-press
                [x,y,buttons]=GetMouse;
                
                if any(buttons)
                    click_in_rect = 0;
                    for rect_num=1:16
                        if ~selected(rect_num)
                            inside = IsInRect(x,y,allRects(rect_num,:));
                            if inside
                                click_in_rect=1;
                                Screen('FillRect',window,0.8,allRects(rect_num,:));
                                [~,click_rt] = Screen('Flip', window, [], 1);
                                result(index).select_rts = [ result(index).select_rts click_rt];
                                result(index).select_pics=[ result(index).select_pics rect_num];
                                selected(rect_num)=1;
                                break;
                            end
                        end
                        
                    end
                    if ~click_in_rect
                        buttons=[];
                    end
                end
            end
            while any(buttons) % wait for release
                [x,y,buttons] = GetMouse;

            end
            
            mouseclicks=mouseclicks+1;

        end
        
        result(index).num_repeat_selected = length(intersect(result(index).correct_repeat_position,result(index).select_pics));
        result(index).num_novel_selected = length(intersect(result(index).correct_novel_position,result(index).select_pics));
        
        Screen('FillRect', window, white);
        DrawFormattedText(window, 'Click on the pictures in the order that you have seen them.\n\nPress Space To Begin.', 'center', 'center', black);
        Screen('Flip', window);
        KbWait;
        
        Screen('FillRect', window, white);
        order=Shuffle(1:4);
        result(index).correct_order = order;
        
        xpts=[xCenter xCenter xCenter/2 xCenter*1.5];
        ypts=[yCenter/2 yCenter*1.5 yCenter yCenter];
        
        four_rects = nan(4,4);
        for i = 1:4
            
            four_rects(i,:) = CenterRectOnPointd([0 0 300 300], xpts(i), ypts(i));
            if cond<7
                Screen('DrawTexture', window, run_context_cafe_tex(order(i),cond), [], four_rects(i,:), 0);
            else
                Screen('DrawTexture', window, run_context_beach_tex(order(i),cond-6), [], CenterRectOnPointd([0 0 300 300], xpts(i), ypts(i)), 0);
            end
            
        end
        
        
        [~, result(index).sort_onset] = Screen('Flip', window,[],1);
      
        sort_num=1;
        result(index).sort_rts =[];
        result(index).sort_order = [];
        selected=zeros(1,4);
        while sort_num<5
             [x,y,buttons]=GetMouse;
            
            while ~any(buttons)
                %waits for a key-press
                [x,y,buttons]=GetMouse;
                
                if any(buttons)
                    click_in_rect =0;
                    for rects = 1:4
                        if ~selected(rects)
                            if IsInRect(x,y,four_rects(rects,:))
                                Screen('FillRect',window,0.8,four_rects(rects,:));
                                [rectx,recty]=RectCenter(four_rects(rects,:));
                                DrawFormattedText(window, num2str(sort_num), rectx, recty, black);
                                [~,click_rt] = Screen('Flip', window, [], 1);
                                result(index).sort_rts = [result(index).sort_rts click_rt];
                                result(index).sort_order = [result(index).sort_order rects];
                                selected(rects)=1;
                                click_in_rect=1;
                                break;
                            end
                        end
                    end
                    if ~click_in_rect
                        buttons = [];
                    end
                end
            end
            while any(buttons) % wait for release
                [x,y,buttons] = GetMouse;
            end
            sort_num = sort_num+1;
            
        end
        WaitSecs(1);

        [keyIsDown, secs, keyCode] = KbCheck(-1);
        if keyCode(KbName('ESCAPE')) == 1
            
%             save(resultname,'result');
            sca;
            disp('*** Experiment terminated ***');
            return
        end
%         save(resultname,'result');
    end
 
    
    Screen('FillRect', window, white);
    DrawFormattedText(window, 'Thank you! Press any key to quit.', 'center', 'center', black);
    Screen('Flip', window);
    KbWait;
    
    sca
catch e
    sca
    rethrow(e);
end


