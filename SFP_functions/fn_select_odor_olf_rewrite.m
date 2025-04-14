function d = fn_select_odor_olf_rewrite(res)
% Rewritten from Cogent to Psychtoolbox.
% Preserves overall logic and trial flow of the original code, but uses
% PTB calls (Screen, DrawFormattedText, Flip, GetSecs, WaitSecs, GetMouse)
% for stimulus presentation, timing, and input instead of Cogent 2000.

% By default, PTB coordinates start at the top-left of the screen, but
% for text and shapes, we often rely on center coordinates. We'll query
% the screen size, then reference (xCenter, yCenter) for drawing.

    %% ============== SETUP & VARIABLES ==============
    Screen('Preference', 'SkipSyncTests', 1); % For quick debugging; remove in production
    screens = Screen('Screens');
    screenNumber = max(screens);
    backgroundColor = [0 0 0]; % black
    
    [window, windowRect] = Screen('OpenWindow', screenNumber, backgroundColor);
    [xCenter, yCenter]   = RectCenter(windowRect);
    % set text and defaults
    Screen('TextFont',  window, 'Arial');
    Screen('TextSize',  window, 30);
    Screen('TextColor', window, [255 255 255]);
    
    % For timing
    t0 = GetSecs;  % “time zero”
    
    % Variables from your original code
    odornames    = res.odornames;
    sniffDur     = 2000;                   % ms odor presentation
    odDelay      = res.delay;              % ms delay between odor line open & sniff
    pre_cue_dur  = 1500;                   % ms from odor open to sniff-cue onset
    Q1_dur       = 2600;                   % ms from question-1 onset to question-2 onset
    Q1_wait      = 800;                    % ms between detect question & descriptor scale
    cueDur       = sniffDur + odDelay + Q1_dur;
    nRep         = res.nRep;
    nStim        = res.nStim;
    ntrials      = nRep * nStim;
    
    %% Create descriptor lists based on session
    if res.sess == 1
        percept_list = {' Intensity';' Pleasantness';' Edible'; ' Familiar'; ...
                        ' Fishy';' Burnt';' Sour'; ' Decayed'; ' Musky'};
        pd_id = (1:1:nRep)';
    elseif res.sess == 2
        percept_list = {' Intensity';' Pleasantness';' Edible'; ' Familiar'; ...
                        ' Fishy';' Burnt';' Sour'; ' Decayed'; ' Musky'};
        pd_id = (1:9)';
        pd_id = pd_id(1:nRep);
    else
        percept_list = {' Fruity'; ' Sweaty'; ' Cool';' Floral';' Sweet',...
                        ' Warm';' Bakery-like'; ' Spicy';' Ammonia'};
        pd_temp = (10:18)';
        pd_id = pd_temp(1:nRep);
    end
    
    % Data array
    % [1: odor order (1-10)
    %  2: CID of odor
    %  3: onset of odor trigger (sec)
    %  4: onset of sniff cue (sec)
    %  5: descriptor ID
    %  6: detect rating
    %  7: button pressed (mouse left/right)
    %  8: detect RT (s)
    %  9: time at which detection was submitted
    % 10: descriptor rating
    % 11: descriptor RT (s)
    % 12: time at which descriptor rating was submitted
    % 13: (optional) initial starting point of scale
    %  -- You can expand these columns as needed.
    d = zeros(ntrials,13);
    
    % Pseudorandom stimulus order
    [~, order] = sort(rand(nStim, nRep)); 
    order = order(:)';
    d(:,1) = order;
    % Assign odor IDs
    CID_list   = cell2mat(odornames(1:end-1,1));
    d(:,2)     = CID_list(order);
    
    % Perceptual descriptor for each trial
    % We replicate the “trial_shuffler(des_id, order, 'mat')” logic
    des_id       = (1:1:nRep)';
    QID          = trial_shuffler_ptb(des_id, order);
    percept_mat  = percept_list(QID);
    d(:,5)       = pd_id(QID);
    
    % Generate ITIs from a distribution
    itis = trial_shuffler_ptb(res.pd, order);
    
    % Olfactometer init (if you have that hardware code in place)
    if res.run_odors
        daq = OlfConfigDaq;  %#ok<NASGU>  (User-supplied)
        OlfFlowControl(daq, odornames{end,4}(1), odornames{end,4}(2));
        OlfOpenLines(odornames{end,3}, daq, 16);
    end
    
    %% ============== WAIT FOR SCANNER/TRIGGER, ETC. ==============
    DrawFormattedText(window, 'Please wait.', 'center', 'center', [255 255 255]);
    Screen('Flip', window);
    
    % Example: wait for '5' key (scanner pulse) or user press 'p'
    validTriggerKeys = KbName({'5%', 'p'});
    FlushEvents; 
    while true
       [keyIsDown, ~, keyCode] = KbCheck;
       if keyIsDown && any(keyCode(validTriggerKeys))
           break;
       end
    end
    
    % small 2s delay
    WaitSecs(2.0); 
    t0 = GetSecs; % reset time zero after “trigger”
    
    if res.run_odors
       OlfOpenLines(odornames{end,3}, daq, 14); 
       WaitSecs(0.05);
       OlfOpenLines(odornames{end,3}, daq, 16);
    end

    % Precompute onsets for odor triggers
    % Ctrigger = cumsum([t0+0.5, sniffDur + itis]) in seconds, but let's do it in a loop:
    intervals = sniffDur/1000 + itis/1000; % convert ms -> s
    Ctrigger  = zeros(1, ntrials);
    Ctrigger(1) = t0 + 0.5; % 500 ms after t0
    for ii = 2:ntrials
        Ctrigger(ii) = Ctrigger(ii-1) + intervals(ii-1);
    end
    
    %% ============== MAIN TRIAL LOOP ==============
    for t = 1:ntrials
        
        % 1) Wait until next odor trigger time
        waitUntil_ptb(Ctrigger(t));
        
        if res.run_odors
            % Switch MFC flow
            OlfFlowControl(daq, odornames{d(t,1),4}(1), odornames{d(t,1),4}(2));
            OlfOpenLines(odornames{d(t,1),3}, daq, d(t,1));
        end
        tOdorTrigOns = GetSecs;
        d(t,3)       = tOdorTrigOns - t0;  % odor trigger onset (sec from t0)
        
        % Wait up to (odDelay - pre_cue_dur) ms
        WaitSecs( (odDelay - pre_cue_dur)/1000 );
        
        % 2) Show fixation cross in white
        DrawFormattedText(window, '+', 'center', 'center', [255 255 255]);
        Screen('Flip', window);
        WaitSecs(pre_cue_dur/1000);
        
        % 3) Sniff cue: change color to e.g. blue
        DrawFormattedText(window, '+', 'center', 'center', [0 0 255]);
        tsniffCueOns   = Screen('Flip', window);
        d(t,4)         = tsniffCueOns - t0;
        
        % Switch flow to clean
        if res.run_odors
            OlfFlowControl(daq, odornames{end,4}(1), odornames{end,4}(2));
            OlfOpenLines(odornames{end,3}, daq, 16);
        end
        
        WaitSecs(sniffDur/1000);
        
        % 4) Ask detectability
        [detect, rtDetect, sw] = detectability_ptb(window, t, Q1_dur-Q1_wait, res.mode);
        d(t,6)  = detect;
        d(t,7)  = sw;
        d(t,8)  = rtDetect;        % RT
        d(t,9)  = d(t,3) + cueDur/1000 + rtDetect; % time at which detection rating was done

        % 5) If subject said "yes" or if we are in certain sessions, show scale
        if (detect ~= 0 && res.sess ~= 4)
            WaitSecs(Q1_wait/1000);
            [ratingVal, ratingRT, scaleStart] = fn_scale_horizontal_ptb(window, ...
                                              percept_mat{t}, 4.0);
            d(t,10) = ratingVal;
            d(t,11) = ratingRT;
            d(t,13) = scaleStart;  % if you want to store that
            d(t,12) = d(t,3) + cueDur/1000 + ratingRT;
        else
            d(t,10) = nan;
            d(t,11) = nan;
            d(t,12) = d(t,3) + cueDur/1000; % or also NaN
        end
        
        % Optionally clear screen or show a blank
        Screen('Flip', window);
    end
    
    % small delay after final trial
    WaitSecs(2.0);
    
    %% ============== SAVE DATA ==============
    clock_int = fix(clock);
    time_str  = sprintf('%02.0f%02.0f', clock_int(4), clock_int(5));
    subject   = sprintf('Complete_%s_%s', res.subj, time_str);
    if ~exist('res','dir'), mkdir('res'); end
    save(fullfile('res', [subject '.mat']), 'd');
    
    %% ============== CLOSE ==============
    sca;  % closes all PTB windows
end


%% =========== SUBFUNCTIONS BELOW ===========

function Qout = trial_shuffler_ptb(des_id, order)
% A simple mock of your 'trial_shuffler' that returns a matrix of
% descriptor indices aligned with the trial order. Adjust as needed.
    nStim   = length(unique(order));
    nRep    = length(des_id);
    desGrid = repmat(des_id, 1, nStim);   % each column is a rep vector
    % Flatten in the same order as "order"
    Qtemp   = desGrid(:)';
    Qout    = Qtemp; % or more sophisticated random logic as needed
end

function waitUntil_ptb(tTarget)
% Busy-wait until GetSecs >= tTarget. 
% This is the PTB version of Cogent's waituntil(Ctrigger(t)).
    while GetSecs < tTarget
        % do nothing
    end
end


function [detect, rt, sw] = detectability_ptb(window, trialNo, t_total, modeFlag)
% Replaces detectability() from your Cogent code with PTB.
% "Smell anything?" with 2 clickable boxes: "Yes" / "No".
%
% Inputs:
%   window   - the PTB window pointer
%   trialNo  - for debugging text
%   t_total  - maximum time allowed in sec
%   modeFlag - if ==99, we display trialNo in the question for debugging
% Outputs:
%   detect   - 1 or 0
%   rt       - reaction time in sec
%   sw       - which mouse button (0 or 1), or -1 if none

    startT = GetSecs;
    rt     = NaN;
    detect = NaN;
    sw     = -1;
    
    yesBox   = CenterRectOnPoint([0 0 100 50],  xCenter(window)+60, yCenter(window));
    noBox    = CenterRectOnPoint([0 0 100 50],  xCenter(window)-60, yCenter(window));
    
    if modeFlag == 99
        questionStr = sprintf('Smell anything?_%d', trialNo);
    else
        questionStr = 'Smell anything?';
    end
    
    while (GetSecs - startT) < t_total
        % Draw question
        Screen('FillRect', window, [255 255 255], yesBox);
        Screen('FillRect', window, [255 255 255], noBox);
        
        DrawFormattedText(window, 'Yes',  xCenter(window)+60-15, yCenter(window)+10, [0 0 0]);
        DrawFormattedText(window, 'No',   xCenter(window)-60-15, yCenter(window)+10, [0 0 0]);
        
        DrawFormattedText(window, questionStr, 'center', yCenter(window)-100, [255 255 255]);
        
        Screen('Flip', window);
        
        % Check mouse
        [mx, my, buttons] = GetMouse(window);
        if any(buttons)
            clickT = GetSecs;
            rt     = clickT - startT;
            if IsInRect(mx, my, yesBox)
                detect = 1;
                sw     = 1;
            elseif IsInRect(mx, my, noBox)
                detect = 0;
                sw     = 0;
            else
                detect = NaN;
                sw     = -1;
            end
            break;
        end
    end
    
    % If timed out, detect remains NaN
    % If user clicked, we have detect & rt
    % Gray-out the boxes to indicate selection
    if sw > -1
        Screen('FillRect', window, [150 150 150], yesBox);
        Screen('FillRect', window, [150 150 150], noBox);
        DrawFormattedText(window, 'Yes', xCenter(window)+60-15, yCenter(window)+10, [0 0 0]);
        DrawFormattedText(window, 'No',  xCenter(window)-60-15, yCenter(window)+10, [0 0 0]);
        Screen('Flip', window);
        % You may optionally let it remain up for the leftover time
        leftover = t_total - (GetSecs - startT);
        if leftover > 0 && modeFlag == 1
            WaitSecs(leftover);
        end
    end
end


function [rating, rt, startPt] = fn_scale_horizontal_ptb(window, descriptor, timeLimit)
% Replaces fn_scale_horizontal from Cogent with a simpler PTB version.
% A vertical rating scale is drawn; user moves the mouse up/down and
% clicks to confirm. You can adapt for horizontal as well.
%
% Inputs:
%   window     - PTB window
%   descriptor - text label of the dimension being rated
%   timeLimit  - max time in seconds to collect response
%
% Outputs:
%   rating - numeric rating
%   rt     - how long until user clicked
%   startPt- initial position (if you want to store it)

    startT   = GetSecs;
    rating   = NaN;
    rt       = NaN;
    startPt  = 0;  % optional
    
    % scale range (in “screen” Y coordinates)
    max_y  = 120; 
    yRange = [-max_y, max_y];
    curY   = 0;  % start in the middle
    step   = 2;  % how quickly the marker moves per mouse increment
    
    % Optionally randomize scale direction
    coin_toss    = randi(2); 
    sign_flipper = (-1)^(coin_toss+1);  % +1 or -1
    
    % Map: for "Intensity" => from Very weak ... Very strong, etc.
    labels = {'Not at all','Somewhat','Moderate','Strong','Extremely'};
    if strcmpi(strtrim(descriptor),'Intensity')
        labels = {'Very weak','','','', 'Very strong'};
    elseif strcmpi(strtrim(descriptor),'Pleasantness')
        labels = {'Dislike','','Neutral','', 'Like very much'};
    end
    % ... you can do more conditions if you like
    
    while (GetSecs - startT) < timeLimit
        % Get mouse
        [mx, my, buttons] = GetMouse(window);
        
        % Move rating based on vertical mouse movement
        % The simplest approach: repeatedly read mouse Y, transform to range
        % If you want it incremental, track differences from a previous frame, etc.
        
        % (Example) we do a quick approach: read “my,” center it around yCenter => map
        % to scale. Alternatively track deltas. For a direct approach:
        curY = my - yCenter(window);
        % clamp
        if curY < yRange(1), curY = yRange(1); end
        if curY > yRange(2), curY = yRange(2); end
        
        % Draw the scale line
        Screen('DrawLine', window, [255 255 255], xCenter(window), yCenter(window)+yRange(1), ...
               xCenter(window), yCenter(window)+yRange(2), 3);
        
        % Draw scale labels
        labelPositions = linspace(yRange(1), yRange(2), length(labels));
        for iL = 1:length(labels)
            DrawFormattedText(window, labels{iL}, ...
                xCenter(window)+50, yCenter(window)+labelPositions(iL)-10, [255 255 255]);
        end
        
        % Draw descriptor at top
        DrawFormattedText(window, descriptor, 'center', yCenter(window)-200, [255 255 255]);
        
        % Draw marker
        Screen('FillRect', window, [255 0 0], ...
          CenterRectOnPoint([0 0 5 25], xCenter(window), yCenter(window)+curY));
        
        Screen('Flip', window);
        
        % check if user clicked
        if any(buttons)
            rt = GetSecs - startT;
            % compute rating in [-1..+1], for example
            rating = sign_flipper * (curY / max_y); 
            
            % optional: wait the remainder or break immediately
            % WaitSecs(timeLimit - rt); 
            break;
        end
    end
    
    % If no click, rating=NaN; rt=NaN.
    startPt = 0;  % or store the initial position if you want
    Screen('Flip', window);
end


%% Helper: return xCenter, yCenter from a window
function [xC, yC] = xCenter(window)
% “Trick” to get xCenter,yCenter easily inside subfunctions without passing them down repeatedly.
    global __PTB_windowRect
    if isempty(__PTB_windowRect)
        __PTB_windowRect = Screen('Rect', window);
    end
    [xC, yC] = RectCenter(__PTB_windowRect);
end
