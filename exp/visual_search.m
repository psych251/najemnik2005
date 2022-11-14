% Visual search task by Najemnik & Geisler (2005, Nature)
%   PSYCH251 replication project
%   Required: mgl toolbox (https://github.com/justingardner/mgl)

% TODO
% # sanity check the gamma table

% Change log
% # added: red cursor to indicate mouse locations
% # added: time stamps to help eye-tracker data analysis
% # fixed: bugs on indicating target's locations 

% Credit: Jiwon Yeon and Josh Wilson

%
function myscreen = visual_search(varargin)

mglClose; clear; close all; clc;

% Experiment parameters
eye_track                = 0;  % use eye-tracker
stimulus.noise.contrast  = .2; % background noise contrast
stimulus.gabor.contrasts = [.15, .2, .25, .35, .5]; % visibility map - based

%
myscreen.hideCursor         = 1;
myscreen.displayName        = 'vpixx';
myscreen.screenNumber       = 2;
myscreen.calibType          = 'Specify particular calibration';
myscreen.calibFilename      = '0001_dn0a221834_221005.mat';
myscreen.calibFullFilename  = '/Users/gru/proj/mgl/task/displays/0001_dn0a221834_221005';
myscreen.datadir            = '/Users/gru/proj/hwgu/data/geisler_search';
myscreen.saveData           = 1;

%
myscreen = initScreen(myscreen);
mglSetParam('abortedStimfilesDir', [myscreen.datadir, '/aborted'], 1);

% 
stimulus.nBlocks = 6;    
stimulus.cBlock  = 0;
stimulus.TrialsPerBlock = 32;
stimulus.responsekeys = 50;   
stimulus.noise.size   = 15;   
stimulus.gabor.size   = 1;      
stimulus.gabor.tilt   = 315;
stimulus.gabor.cycle  = 6;
stimulus = create_pink_filter(stimulus, myscreen);
stimulus.noise.noise_frame_pixel = va_to_pix(stimulus.noise.size+3, ...
    [myscreen.screenWidth, myscreen.screenHeight]);
stimulus = define_locations(stimulus); 

% 
task{1}{1}.waitForBacktick = 0;
task{1}{1}.segmin = [inf, .1, inf, inf, 2];  
task{1}{1}.segmax = [inf, .5, inf, inf, 2];
task{1}{1}.getResponse = [0 0 1 0 0];
task{1}{1}.numTrials = stimulus.nBlocks * stimulus.TrialsPerBlock;
task{1}{1}.randVars.block.gabor_location = 1:stimulus.gabor.nLoc;
task{1}{1}.randVars.block.gabor_contrast = stimulus.gabor.contrasts; 
task{1}{1}.randVars.calculated.detection_rt = nan;
task{1}{1}.randVars.calculated.decision_rt = nan;
task{1}{1}.randVars.calculated.framecount = nan;
task{1}{1}.randVars.calculated.mousePos = [nan nan];
task{1}{1}.randVars.calculated.response_offset = [nan nan];
task{1}{1}.randVars.calculated.t_mouse_onset = nan;
task{1}{1}.randVars.calculated.t_stim_onset  = nan;

%
myscreen = initStimulus('stimulus', myscreen);
mglClearScreen(.5); mglFlush;
mglStencilCreateBegin(1);

% 
x = linspace(-7.5, 7.5, 100000);
y = sqrt(7.5^2 - x.^2);
x = [x, fliplr(x)];
y = [y, -fliplr(y)];
mglPolygon(x,y,[1 1 1]);
mglStencilCreateEnd;
mglClearScreen(.5); mglFlush;
mglDisplayCursor(0);

% 
disp('Initializing task ...');
[task{1}{1}, myscreen] = initTask(task{1}{1}, myscreen,...
    @startSegmentCallback, @updateScreenCallback, @getResponseCallback, ...
    @startTrialCallback);

% 
if eye_track
    disp(' Calibrating eye-tracker ...')
    myscreen = eyeCalibDisp(myscreen);
end

% 
mglClearScreen(.5);
mglTextSet([],32,1);
mglTextDraw('Starting the experiment',[0,0]); mglFlush;
while 1
    k = mglGetKeys;
    if k(myscreen.keyboard.backtick)==1, break; end
end

% Start the task
while (task{1}{1}.trialnum <= task{1}{1}.numTrials) && ~myscreen.userHitEsc
    [task{1}, myscreen] = updateTask(task{1},myscreen,1);
    myscreen = tickScreen(myscreen, task{1});
end

% End the task
mglClearScreen(0.5);
mglTextSet([],32,1);
mglTextDraw('Experiment ends',[0, .7]);
mglTextDraw('Please wait..', [0, -.7]); mglFlush;
mglWaitSecs(3);
myscreen = endTask(myscreen,task);

end


%% callbacks
function [task, myscreen] = startTrialCallback(task, myscreen)
    global stimulus
    
    if mod(task.trialnum, stimulus.TrialsPerBlock) == 1
        stimulus.cBlock = stimulus.cBlock + 1;
    end
    
    % generate noise images
    stimulus = create_pink_noise(stimulus, myscreen);
    
    % generate gabor
    stimulus = create_gabor(stimulus,task);
    
    % convert gabor locations to pixels
    if task.trialnum == 1
        % convert visual angle of the locations to pixels
        displaySize = size(stimulus.noise.im);
        stimulus.gabor_locations_pix = va_to_pix(stimulus.gabor_locations_va, displaySize);
    end

    stimulus = combine_stimuli(stimulus,task);
    
    % create texture in advance
    stimulus.tex_target = mglCreateTexture(stimulus.final_im{1});
    stimulus.tex_nontarget = mglCreateTexture(stimulus.final_im{2});
    
    % reset framecount
    task.thistrial.framecount = 0;
end

%
function [task, myscreen] = startSegmentCallback(task, myscreen)
    global stimulus
    
    if task.thistrial.thisseg == 1
        % show how many blocks are left 
        if mod(task.trialnum, stimulus.TrialsPerBlock) == 1
            mglClearScreen(.5);
            mglTextSet([],32,1);
            mglTextDraw(sprintf('Starting block %d out of %d blocks', ...
                stimulus.cBlock, stimulus.nBlocks),[0,0]);
            mglFlush(); mglWaitSecs(2);
        end
        
        % show a fixation cross and wait for the button press
        mglClearScreen(.5);
        mglFillOval(0,0,[.2 .2],0);
        mglFlush();
        
        % waiting for the subject to start the trial
        while 1
            keycode = mglGetKeys;
            if any(keycode(stimulus.responsekeys)==1) % spacebar
                break
            end
        end        
        task = jumpSegment(task);
        
    elseif task.thistrial.thisseg == 2 
        % present an empty screen
        mglClearScreen(.5);    
        
    elseif task.thistrial.thisseg == 3
        disp(['current target:' num2str(task.thistrial.gabor_location)])

        % stimulus-locked time
        task.thistrial.t_stim_onset = mglGetSecs;

    elseif task.thistrial.thisseg == 4
        % set the mouse cursor at the center
        mglSetMousePosition(myscreen.screenWidth/2, myscreen.screenHeight/2, ...
            myscreen.screenNumber);
    
        % start response time recording
        stimulus.t0 = mglGetSecs;
        task.thistrial.t_mouse_onset = stimulus.t0;
        
    elseif task.thistrial.thisseg == 5
        mglClearScreen(stimulus.bg_color{2}/255*1);
        mglStencilSelect(1);
        mglBltTexture(stimulus.tex_nontarget,[0 0]);
        mglStencilSelect(0);
        mglFillOval(0,0,[.2 .2],0);     % fixation
    end
end

%
function [task, myscreen] = updateScreenCallback(task, myscreen)

    global stimulus
    if task.thistrial.thisseg == 3
        % stimulus presentation
        mglClearScreen(stimulus.bg_color{1}/255*1);
        mglStencilSelect(1);
        mglBltTexture(stimulus.tex_target,[0 0]);
        mglStencilSelect(0);

    elseif task.thistrial.thisseg == 4
        % update framecount
        task.thistrial.framecount = task.thistrial.framecount+1;
    
        % decision prompt
        mglClearScreen(stimulus.bg_color{2}/255*1);
        mglStencilSelect(1);
        mglBltTexture(stimulus.tex_nontarget,[0 0]);
        mglStencilSelect(0);
        mglFillOval(0,0,[.2 .2], 0);
        mglTextDraw('Click on the screen where the target appeared', [0,10]);
        
        % record mouse positions
        mInfo = mglGetMouse(myscreen.screenNumber);        
        x = (mInfo.x - myscreen.screenWidth/2) * myscreen.imageWidth/myscreen.screenWidth;
        y = (mInfo.y - myscreen.screenHeight/2) * myscreen.imageHeight/myscreen.screenHeight;
        mousePos(task.thistrial.framecount,:) = [x,y];
    
        % display a cursor
        mglMetalDots([x;y;0], [1;0;0;1], [0.2;0.2], 1, 1);
    
        if mInfo.buttons == 1
            % get decision RT
            task.thistrial.decision_rt = mglGetSecs(stimulus.t0);
    
            % save response info
            task.thistrial.mousePos = mousePos;     % in degrees
            task.thistrial.response_offset = task.thistrial.gabor_location - [mInfo.x, mInfo.y];
    
            task = jumpSegment(task);
        end
    end
end 

%
function [task, myscreen] = getResponseCallback(task, myscreen)
    global stimulus
    if task.thistrial.thisseg == 3
        % get detection response
        while 1
            keycode = mglGetKeys;
            if any(keycode(stimulus.responsekeys)==1)
                task.thistrial.detection_rt = task.thistrial.reactionTime;
                break
            end
        end
        task = jumpSegment(task);
    end
end


%% utils
function stimulus = create_pink_filter(stimulus, myscreen)
    w = myscreen.screenWidth;
    h = myscreen.screenHeight;
    sz = max(w,h);
    
    % make the odd size of the image
    if mod(sz,2)==0, sz = sz-1; end
    
    % make pink filter
    last_freq = ceil(sz/2);
    pink_filter = zeros(sz,sz);
    [x y] = meshgrid(-ceil(sz/2)+1:ceil(sz/2)-1, -ceil(sz/2)+1:ceil(sz/2)-1);
    index = sqrt(x.^2 + y.^2);
    for f = 1:last_freq
        pink_filter(index > f-1 & index < f+1) = 1/f;
    end
    stimulus.pink_filter = pink_filter;
end

%
function stimulus = create_pink_noise(stimulus, myscreen)

    noise_with_buffer = stimulus.noise.size + 3;    % visual angle
    noise_frame_pixel = va_to_pix(noise_with_buffer, ...
        [myscreen.screenWidth, myscreen.screenHeight]);
    stimulus.noise.noise_frame_pixel = noise_frame_pixel;
    % make the size of the image an odd number
    if mod(noise_frame_pixel,2)==0, noise_frame_pixel = noise_frame_pixel+1; end    
    
    filter_sz = size(stimulus.pink_filter);
    pink_filter = stimulus.pink_filter(...
        floor(filter_sz(1)/2)+1-(noise_frame_pixel-1)/2:floor(filter_sz(1)/2)+1+(noise_frame_pixel-1)/2, ...
        floor(filter_sz(2)/2)+1-(noise_frame_pixel-1)/2:floor(filter_sz(2)/2)+1+(noise_frame_pixel-1)/2);
    
    % fft on white noise
    white = randn(noise_frame_pixel, noise_frame_pixel);
    fwhite = fftshift(fft2(white));
    phase = angle(fwhite);
    
    % create new magnitude
    new_mag = fwhite .* pink_filter;
    new_Fourier = new_mag .* (cos(phase) + sqrt(-1)*sin(phase));
    im = ifft2(ifftshift(new_Fourier));
    
    % change contrast
    contrast = stimulus.noise.contrast;
    N = length(im(:));
    m_im = mean(im(:));
    coeff = sqrt((N*contrast^2) / sum((im(:)-m_im).^2));
    stimulus.noise.im = coeff .* im;
end

%
function stimulus = create_gabor(stimulus,task)

    grating = mglMakeGrating(stimulus.gabor.size, stimulus.gabor.size, ...
        stimulus.gabor.cycle, stimulus.gabor.tilt, 0);
    contrast = task.thistrial.gabor_contrast;
    grating = grating * contrast;
    gaussian = mglMakeGaussian(stimulus.gabor.size, stimulus.gabor.size, 1, 1);
    stimulus.gabor.im = grating.*gaussian;

end 

%
function stimulus = define_locations(stimulus)

    h_dist = 1.5; 
    v_dist = sqrt(h_dist^2-(h_dist/2)^2);
    y_lim = stimulus.noise.size - 1;
    
    trigrid = [];
    y_current = 0;
    xx = 0;
    displacement = 0;
    while y_current < y_lim
        if displacement == 0
            xx = [0:h_dist:y_lim]';
            yy = ones(length(xx), 1)*y_current;
            displacement = 1;
        else
            xx = [h_dist/2:h_dist:y_lim]';
            yy = ones(length(xx), 1)*y_current;
            displacement = 0;
        end
        trigrid = [trigrid;[xx, yy]];
        y_current = y_current+v_dist;
    end
    
    trigrid = trigrid - repmat(max(trigrid)./2,size(trigrid,1),1);
    inside = sqrt(trigrid(:,1).^2+trigrid(:,2).^2) <= y_lim/2;
    locations_va = trigrid(inside,:);
    locations = va_to_pix(locations_va, stimulus.noise.noise_frame_pixel);
    
    stimulus.gabor_locations_va = locations_va;
    stimulus.gabor_locations = locations;   % in pixels
    stimulus.gabor.nLoc = size(locations,1);

end
    
%
function stimulus = combine_stimuli(stimulus,task)
    noise = stimulus.noise.im;
    gabor = stimulus.gabor.im;
    location_num = task.thistrial.gabor_location;   
    location = stimulus.gabor_locations_pix(location_num,:);
    
    % make circular gabor patch
    radius_pixel = va_to_pix(stimulus.gabor.size/2); 
    gabor_sz     = size(gabor);
    [stencil_x, stencil_y] = meshgrid(-(gabor_sz(1)-1)/2:(gabor_sz(1)-1)/2, ...
        -(gabor_sz(2)-1)/2:(gabor_sz(2)-1)/2);
    stencil = (sqrt(stencil_x.^2 + stencil_y.^2) <= radius_pixel);
    gabor_circle = stencil' .* gabor;
    
    % determine the location to display
    gabor_position = zeros(size(noise,1), size(noise,2));
    gabor_position([location(1)-(gabor_sz(1)-1)/2:location(1)+(gabor_sz(1)-1)/2], ...
        [location(2)-(gabor_sz(2)-1)/2:location(2)+(gabor_sz(2)-1)/2]) = gabor_circle;
    
    % add gabor to the noise
    final_im = noise + gabor_position;
    
    % clip to max and min
    final_im(final_im > 1) = 1;
    final_im(final_im < -1) = -1;
    
    % scale it to [0 255], for both stimulus images
    final_im =  255 .* ((final_im + 1) ./ 2);
    stimulus.final_im{1} = final_im;
    stimulus.final_im{2} = 255 .* ((stimulus.noise.im + 1) ./ 2);    % noise-only image
    
    % decide background color
    for image = 1:2
        bg_color = mean(stimulus.final_im{image}(:));
        stimulus.bg_color{image} = bg_color;
    end
end

% 
function pix = va_to_pix(va, displaySize)

    if nargin < 2
        displaySize = [mglGetParam('screenWidth'), mglGetParam('screenHeight')];
    end
    if isscalar(displaySize)
        displaySize = [displaySize, displaySize];
    end
    
    if isvector(va)
        if size(va,1) == 2 
            va = va';
            needTranspose = 1;
        elseif all(size(va) > 2) || length(size(va)) > 2
            error('Visual angle must be a scalar or 2D coordinates')
        else
            needTranspose = 0;
        end
    else
        needTranspose = 0;
    end
    
    xDeg2Pix = mglGetParam('xDeviceToPixels');
    yDeg2Pix = mglGetParam('yDeviceToPixels');
    
    % change visual angle values to pixels
    if isscalar(va)
        pix = round(va .* max([xDeg2Pix,yDeg2Pix]));        
    else    
        pix = [va(:,1) .* yDeg2Pix, va(:,2) .* xDeg2Pix];    
        pix = round([pix(:,1) + floor(displaySize(1)/2), ...
            pix(:,2) + floor(displaySize(2)/2)]);     
    end
    
    if needTranspose
        pix = pix';
    end
end