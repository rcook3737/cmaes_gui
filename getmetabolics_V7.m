%% getmetabolics_V5 TESTING
% INPUTS
% BEN LOOK HERE
% test_mode
%     1 = saved video
%     2 = real time camera
% test_method
%     0 = not resizing the boxes
%     1 = click for box, no tracking
%     2 = click for box, tracking

%%% Camera Calibration
% test_mode = 2;
% test_method = 1;
% mins = .25;
% widthBOX = 300;
% heightBOX = 125;
% trial = 1;

%%% Zero impedance
% test_mode = 2;
% test_method = 0;
% mins = .2;
% widthBOX = 300;
% heightBOX = 125;

%%% Sweep
% test_mode = 2;
% test_method = 0;
% mins = 2;
% widthBOX = 300;
% heightBOX = 125;

function [t,vo2_rate,vco2_rate] = getmetabolics_V7(test_mode,test_method,mins,widthBOX,heightBOX,trial)
warning('off','all')
fig=figure;
data_num = 1;
data_collection = false;

% Access webcam and display/import videoframes
% Set up new computer with webcam
% camList = webcamlist
% cam_num = input('Which cam?'); % For s[etting up new comp with webcam
% cam = camList(cam_num);

% Already setup
if test_mode == 2
    camList = webcamlist;
    web_cam = camList(1);
    cam = webcam(web_cam{:});
    cam.Resolution = cam.AvailableResolutions{end}; % Set webcam res to max
    currentFrame = snapshot(cam);
    ax = axes(fig); 
    im = image(ax,zeros(size(currentFrame),'uint8'));
    axis(ax,'image');
    try
        preview(cam,im)
    catch
        
    end
end

% Saved Video
if test_mode == 1
    vo2_text = [];
    vco2_text = [];
    v = VideoReader('testvid1trim.mp4');
    currentFrame = readFrame(v);
    imshow(currentFrame)
    videoPlayer = vision.VideoPlayer('Position',[50,300,680,400]);
end
% Visualization for camera position and function
global objectRegion1 objectRegion2 rows1 cols1 rows2 cols2

%%   TEST METHOD 1
% box no tracking

if test_method == 1
    vo2approved = 1;
    vco2approved = 1;
    sizing_run = 1;
    % current frame
    figure(fig);
    currentFrame = snapshot(cam);
    %         imshow(currentFrame)
    I = rgb2gray(currentFrame);
    preview(cam,im)
    while sizing_run ~= 0
        
        %             VO2
        
        %     Create center for VO2
        while vo2approved == 1
            if exist('vo2modified')
                edits = input(['Please enter a new width and height for VO2 as shown: [width,height]\n']);
                widthBOX1 = edits(1);
                heightBOX1 = edits(2);
                figure;
                fprintf('Please click in the center of the VO2 numbers\n')
                [rows1,cols1] = ginput(1);
            else
                fprintf('Please click in the center of the VO2 numbers\n')
                [rows1,cols1] = ginput(1);
                widthBOX1 = widthBOX;
                heightBOX1 = heightBOX;
            end
            % Compute box dimensions for VO2
            VOX1 = (rows1 - widthBOX1/2);
            VOX2 = (rows1 + widthBOX1/2);
            VOY1 = (cols1 - heightBOX1/2);
            VOY2 = (cols1 + heightBOX1/2);
            vo2approved = 0;
        end
        
        %             VCO2
        
        %     Create center for VCO2
        while vco2approved == 1
            if exist('vco2modified')
                edits = input(['Please enter a new width and height for VCO2 as shown: [width,height]\n']);
                widthBOX2 = edits(1);
                heightBOX2 = edits(2);
                f2 = figure;
                fprintf('Please click in the center of the VCO2 numbers\n')
                [rows2,cols2] = ginput(1);
                delete(f2);
            else
                fprintf('Please click in the center of the VCO2 numbers\n')
                [rows2,cols2] = ginput(1);
                widthBOX2 = widthBOX;
                heightBOX2 = heightBOX;
            end
            % Create box dimensions for VCO2
            VCOX1 = (rows2 - widthBOX2/2);
            VCOX2 = (rows2 + widthBOX2/2);
            VCOY1 = (cols2 - heightBOX2/2);
            VCOY2 = (cols2 + heightBOX2/2);
            vco2approved = 0;
        end
        
        
        tic
        while toc <= 10
            % Show inital results    
            vco2_text_frame = I(VCOY1:VCOY2,VCOX1:VCOX2);
            vo2_text_frame = I(VOY1:VOY2,VOX1:VOX2);
            figure(fig);
            currentFrame = snapshot(cam);
            I = rgb2gray(currentFrame);
            subplot(1,2,1); imshow(vo2_text_frame);
            subplot(1,2,2); imshow(vco2_text_frame);
            text1 = ocr(vo2_text_frame);
            text2 = ocr(vco2_text_frame);
            vo2_text = str2double(text1.Text)
            vco2_text = str2double(text2.Text)
        end
        vo2approved = input(['Are you satisfied with the VO2 textboxes?(0 = Yes; 1 = No)\n']);
        vo2modified = 1;
        vco2approved = input(['Are you satisfied with the VCO2 textboxes?(0 = Yes; 1 = No)\n']);
        vco2modified = 1;
        
        if vo2approved == 1 || vco2approved == 1
            sizing_run = 1;
        else
            sizing_run = 0;
        end
        
    end
    
    data_collection = true;
    
    %% TEST METHOD 2
    %tracking and resizing
    
elseif test_method == 2
    tic
    vo2approved = 1;
    vco2approved = 1;
    sizing_run = 1;
    while toc <= 5
        %         current frame
        figure(fig);
        I = rgb2gray(currentFrame);
        while sizing_run ~= 0
            
            % VO2
            
            while vo2approved == 1
                if exist('vo2modified')
                    edits = input(['Please enter a new width and height for VO2 as shown: [width,height]\n']);
                    widthBOX1 = edits(1);
                    heightBOX1 = edits(2);
                    f2 = figure;
                    imshow(currentFrame);
                    fprintf('Please click in the center of the VO2 numbers\n')
                    [rows1,cols1] = ginput(1);
                    delete(f2);
                else
                    if test_method == 2
                        % Create VO2 text box
                        fprintf('Highlight V02ml for anchoring:\nTop Left of V02ml\nBottom right V02ml\n')
                        % Create marker points for tracking VO2
                        objectRegion1 = round(getPosition(imrect));
                        fprintf('Please click in the center of the VO2 numbers\n')
                        [rows1,cols1] = ginput(1);
                    end
                    points1 = detectMinEigenFeatures(rgb2gray(currentFrame),'ROI',objectRegion1);
                    tracker1 = vision.PointTracker('MaxBidirectionalError',1);
                    initialize(tracker1,points1.Location,currentFrame);
                    % Tracker points for VO2
                    [points1,~] = tracker1(currentFrame);
                    widthBOX1 = widthBOX;
                    heightBOX1 = heightBOX;
                end
                Left1 = min(points1);
                Right1 = max(points1);
                Xoffset1 = rows1 - Left1(1);
                Yoffset1 = cols1 - Right1(2);
                % Compute box dimensions for VO2
                VOX1 = (Left1(1) + Xoffset1 - widthBOX1/2);
                VOX2 = (Left1(1) + Xoffset1 + widthBOX1/2);
                VOY1 = (Right1(2) + Yoffset1 - heightBOX1/2);
                VOY2 = (Right1(2) + Yoffset1 + heightBOX1/2);
                vo2_text_frame = I(VOY1:VOY2,VOX1:VOX2);
                vo2approved = 0;
            end
            
            % VCO2
            
            
            while vco2approved == 1
                if exist('vco2modified')
                    edits = input(['Please enter a new width and height for VCO2 as shown: [width,height]\n']);
                    widthBOX2 = edits(1);
                    heightBOX2 = edits(2);
                    f2 = figure;
                    imshow(currentFrame);
                    fprintf('Please click in the center of the VCO2 numbers\n')
                    [rows2,cols2] = ginput(1);
                    delete(f2);
                else
                    if test_method == 2
                        % Create VCO2 text box
                        fprintf('Highlight VC02ml for anchoring:\nTop Left of VC02ml\nBottom right VC02ml\n')
                        % Create marker points for tracking VCO2
                        objectRegion2 = round(getPosition(imrect));
                        fprintf('Please click in the center of the VCO2 numbers\n')
                        [rows2,cols2] = ginput(1);
                    end
                    points2 = detectMinEigenFeatures(rgb2gray(currentFrame),'ROI',objectRegion2);
                    tracker2 = vision.PointTracker('MaxBidirectionalError',1);
                    initialize(tracker2,points2.Location,currentFrame);
                    % Tracker points for VCO2
                    [points2,~] = tracker2(currentFrame);
                    widthBOX2 = widthBOX;
                    heightBOX2 = heightBOX;
                end
                Left2 = min(points2);
                Right2 = max(points2);
                Xoffset2 = rows2 - Left2(1);
                Yoffset2 = cols2 - Right2(2);
                % Compute box dimensions for VCO2
                VCOX1 = (Left2(1) + Xoffset2 - widthBOX2/2);
                VCOX2 = (Left2(1) + Xoffset2 + widthBOX2/2);
                VCOY1 = (Right2(2) + Yoffset2 - heightBOX2/2);
                VCOY2 = (Right2(2) + Yoffset2 + heightBOX2/2);
                vco2_text_frame = I(VCOY1:VCOY2,VCOX1:VCOX2);
                vco2approved = 0;
            end
            %         Show initial results
            subplot(1,2,1); imshow(vo2_text_frame);
            subplot(1,2,2); imshow(vco2_text_frame);
            text1 = ocr(vo2_text_frame);
            text2 = ocr(vco2_text_frame);
            vo2_text = str2double(text1.Text)
            vco2_text = str2double(text2.Text)
            if test_method == 2
                vo2approved = input(['Are you satisfied with the VO2 textboxes?(0 = Yes; 1 = No)\n']);
            end
            vo2modified = 1;
            if test_method == 2
                vco2approved = input(['Are you satisfied with the VCO2 textboxes?(0 = Yes; 1 = No)\n']);
            end
            vco2modified = 1;
            
            if vo2approved == 1 || vco2approved == 1
                sizing_run = 1;
            else
                sizing_run = 0;
            end
        end
    end
    
    data_collection = true;
    
    %% TEST METHOD 0
    %no tracking or resizing
    
elseif test_method == 0
    widthBOX1 = widthBOX;
    heightBOX1 = heightBOX;
    
    VOX1 = (rows1 - widthBOX1/2);
    VOX2 = (rows1 + widthBOX1/2);
    VOY1 = (cols1 - heightBOX1/2);
    VOY2 = (cols1 + heightBOX1/2);
    
    widthBOX2 = widthBOX;
    heightBOX2 = heightBOX;
    
    VCOX1 = (rows2 - widthBOX2/2);
    VCOX2 = (rows2 + widthBOX2/2);
    VCOY1 = (cols2 - heightBOX2/2);
    VCOY2 = (cols2 + heightBOX2/2);
    
    data_collection = true;
end

if data_collection
    tic
    global pauseFlag continueFlag stopFlag; %flags for pause/stop
    while toc < mins*60
        if ~pauseFlag
        
            if test_mode == 1
                currentFrame = readFrame(v);
            elseif test_mode == 2
                currentFrame = snapshot(cam);
            end

            if test_method == 2
                % Recreate points for each frame
                [points1,validity] = tracker1(currentFrame);
                [points2,validity] = tracker2(currentFrame);

                %     VO2

                Left1 = min(points1);
                Right1 = max(points1);
                Xoffset1 = rows1 - Left1(1);
                Yoffset1 = cols1 - Right1(2);
                % Compute box dimensions for VO2
                VOX1 = (Left1(1) + Xoffset1 - widthBOX1/2);
                VOX2 = (Left1(1) + Xoffset1 + widthBOX1/2);
                VOY1 = (Right1(2) + Yoffset1 - heightBOX1/2);
                VOY2 = (Right1(2) + Yoffset1 + heightBOX1/2);

                %       VCO2

                Left2 = min(points2);
                Right2 = max(points2);
                Xoffset2 = rows2 - Left2(1);
                Yoffset2 = cols2 - Right2(2);
                % Compute box dimensions for VCO2
                VCOX1 = (Left2(1) + Xoffset2 - widthBOX2/2);
                VCOX2 = (Left2(1) + Xoffset2 + widthBOX2/2);
                VCOY1 = (Right2(2) + Yoffset2 - heightBOX2/2);
                VCOY2 = (Right2(2) + Yoffset2 + heightBOX2/2);
            end

            I = rgb2gray(currentFrame);

            vo2_text_frame = I(VOY1:VOY2,VOX1:VOX2);
            vco2_text_frame = I(VCOY1:VCOY2,VCOX1:VCOX2);

            subplot(1,2,1); imshow(vo2_text_frame);
            subplot(1,2,2); imshow(vco2_text_frame);
            % OCR Number reading and output
            text1 = ocr(vo2_text_frame);
            text2 = ocr(vco2_text_frame);
            vo2_text = str2double(text1.Text);
            vco2_text = str2double(text2.Text);
            disp(['time  = ' sprintf('%.2f', toc)]);
            disp(['Vo2      Vco2']);
            spaces = char(' ' .* ones(1, 6 - (length(num2str(vo2_text)) - 3)));
            disp([num2str(vo2_text) spaces num2str(vco2_text)]);
            disp('  ');

            % Continue to next frame if OCR does not work
            if isnan(vo2_text) || isnan(vco2_text)
                continue
            end

            % Record new data only and mark time point (OCR based)
            if data_num == 1
                vo2_rate(data_num) = vo2_text; %#ok<*AGROW>
                vco2_rate(data_num) = vco2_text;
                %         image_check{data_num} = I;
                %         vo2_rate(frame_num)
                %         vco2_rate(frame_num)
                t(data_num) = toc;
                data_num = data_num + 1;
            else
                if vco2_text ~= vco2_rate(data_num-1) && vo2_text ~= vo2_rate(data_num-1)
                    vo2_rate(1,data_num) = vo2_text;
                    vco2_rate(1,data_num) = vco2_text;
                    %         vo2_rate(frame_num)
                    %         vco2_rate(frame_num)
                    t(1,data_num) = toc;
                    data_num = data_num + 1;
                end
            end
        else
            pauseFlag = false;
            paused = true;
            while paused
               pause(0.1);
               if continueFlag
                   paused = false;
                   continueFlag = false;
                   stopFlag = false;
               elseif stopFlag
                   paused = false;
                   continueFlag = false;
                   stopFlag = false;
                   error('Experiment Stopped');
               end
            end
        end
          
    end
    
%     save('test.mat', 'vo2_rate' ,'vco2_rate' ,'t')
%     if test_mode == 1
%         release(videoPlayer);
%     elseif test_mode == 0 || test_mode == 2
%         closePreview(cam)
%         clear('cam')
%     end
%     delete(figure(fig));
end
% if isempty(trial)
%     save('Imp2D_Wendy09132019_Opt.mat', 'vo2_rate' ,'vco2_rate' ,'t')
% else
%     save(sprintf('Imp2D_Wendy09132019_Opt_%d.mat',trial), 'vo2_rate' ,'vco2_rate' ,'t')
% end
% [~,estimated_dot_E,~,~,~] = parvo_meta_rate_est_fxn(t,[],vo2_rate,vco2_rate);
% sprintf('Metabolic Cost: %.2f',estimated_dot_E)
if test_mode == 1
    release(videoPlayer);
elseif test_mode == 0 || test_mode == 2
    closePreview(cam)
    clear('cam')
end
delete(figure(fig));    
end