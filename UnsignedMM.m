% Version for 13.02.2023, Tomasz Ko³odziej
% This program calculates uMM angle and plots major axis with displacement
% direction.

function [uMMangle, MAcoordinates, Fxy, OrientationDirection] = UnsignedMM(File, Settings, k, ResultsMain, ResultsAuxiliary, NameAppendix, FrameList);

% First condition is for most of frames (except the last one). All proper
% things are calculated here. For the last frame, only the properties of
% Major Axes are calculated, because there's no MM angle for the last frame


    if k < size(ResultsMain, 1)

    % Load mask file
    Mask=imread([File.NameBase, num2str(ResultsMain(k,2)), '.', Settings.FileType]);
    MaskBW = im2bw(Mask, 0.8);    
    ShapeProps=regionprops(MaskBW ,'Centroid', 'MinorAxisLength','MajorAxisLength','Perimeter','Area', 'Orientation'); % Mask's properties

% Orientation and length of Major Axis of the ellipse that was fitted to
% the shape:
    ShapeOrientation = ShapeProps.Orientation;
    MAlength = ShapeProps.MajorAxisLength;
     
    %Sprawdziæ czy to dalej wystêpuje:
    %WspKierOsi = tand(ShapeOrientation); %Funkcja "tand" wylicza tangens dla argumentu W STOPNIACH.            
            
% Cell centroid in current frame [Sx1; Sy1] and following frame [Sx2; Sy2]
    Sx1 = ShapeProps.Centroid(1);
    Sx2 = ResultsMain ((k+1), 4);
    Sy1 = ShapeProps.Centroid(2);
    Sy2 = ResultsMain ((k+1), 5);
    
% Calculating the slope and intercept of the line that covers displacement
% vector. For further calculations, I change the sign of the slope, because
% of the difference between coordinates defined by the image ([0;0]placed
% in upper left corner) and "natural" coordinate system. It's just more
% intuitive for me.
    
    DirectionSlope = (Sy2 - Sy1)/(Sx2 - Sx1);
    DirectionIntercept = Sy1 - (DirectionSlope * Sx1);
    AngleDisplXaxis = atan(DirectionSlope);
    WspKierWektora2 = -DirectionSlope;
    AngleDisplXaxis2 = -AngleDisplXaxis;

    %% Calculating the displacement direction vector:
% Similarily to drawing angles, I calculate the line that lays on the
% displacement vector, (because the displacement vector itself might be too
% short to be visible. From that, I'll calculate the uMM
% straightly from the image. I do it just to have a double-chack, that
% resulted image is consistent with actually calculated data.
    
% Calculating the ending point of direction vector [Fx; Fy]
   
    %If the initial length is too small to be visible (<10 px), the length
    %of major axis is taken:
        if MAlength < 10;  
            DispLength = MAlength;
        else
            DispLength = MAlength/2;
        end
    
    % Where the Fx should be:           
    % Correction 1 - If displacement is vertical (cos 90 deg = 0), then 
    % Fx = Sx:
        if Sx1 == Sx2;
           Fx = Sx1;
    % Correction 2 - If Sy1 < Sy2 or (Sy1 = Sy2 and Sx1 < Sx2) then we add
    % the length of direction vector to the point (also taking into account
    % that cosinus can be negative
        elseif Sx1 < Sx2;
           Fx = Sx1 + (DispLength*cos(AngleDisplXaxis2));
    % Correction 3 - If Sy1 > Sy2, then we subtract the length of direction
    % vector, because the point will be more to the left side (also taking
    % into account that cosinus can be negative).
        elseif Sx1 > Sx2;
           Fx = Sx1 - (DispLength*cos(AngleDisplXaxis2));
        end
               
    % Where the Fy should be: 
    % Correction 1 - If Sy1 = Sy2 (horizontal line), then Fy = Sy.
        if Sy1 == Sy2;
           Fy = Sy1;
    % Correction 2 - If Sx1 < Sx2 or (Sx1 = Sx2 ad Sy1 < Sy2), then we add
    % the length of direction vector to the point (also taking into account
    % that cosinus can be negative)
        elseif Sx1 > Sx2 | (Sx1 == Sx2 && Sy1 < Sy2);
           Fy = Sy1 + (DispLength*sin(AngleDisplXaxis2));
    % Correction 3 - If Sy1 > Sy2, then we subtract the length of direction
    % vector, because the point will be more to the left side (also taking
    % into account that cosinus can be negative).
        elseif Sx1 < Sx2 | (Sx1 == Sx2 && Sy1 > Sy2);
           Fy = Sy1 - (DispLength*sin(AngleDisplXaxis2));
        end
    
    % If the resulting [Fx; Fy] is outside of the image, then its length is
    % shorten to the length of Minor Axis. Calculations are similar as
    % above:
   
       if (Fx < 0 | Fx > size(Mask, 2)) | Fy < 0 | Fy > size(Mask, 1);
          DispLength = ShapeProps.MinorAxisLength/2;
          
          % For Fx:
           if Sx1 == Sx2;
           Fx = Sx
    
           elseif Sx1 < Sx2;
           Fx = Sx + (DispLength*cos(AngleDisplXaxis2));
   
           elseif Sx1 > Sx2;
           Fx = Sx - (DispLength*cos(AngleDisplXaxis2));
           end
           
          % For Fy: 
           if Sy1 == Sy2;
           Fy = Sy;
    
           elseif Sx1 > Sx2 | (Sx1 == Sx2 && Sy1 < Sy2);
           Fy = Sy + (DispLength*sin(AngleDisplXaxis2));
    
           elseif Sx1 < Sx2 | (Sx1 == Sx2 && Sy1 > Sy2);
           Fy = Sy - (DispLength*sin(AngleDisplXaxis2));
           end
       else 
           Fx = Fx; 
           Fy = Fy;
       end
    % [Fx; Fy] point is now calcualted.   
    %----------------------------------------------------------------------    
        
%% Drawing the image:

% In the first stage program draws mask, fitted ellipse, centroids and both
% ends of major axis.
     
    %%Create mage canvas:
    Obrazek = figure; set(Obrazek, 'name', 'uMM');
    hold on
 
    %1. Painting the mask in yellow xD
    MaskRGB=cat(3, Mask, Mask, Mask);
    
    for i=1:1:size(MaskRGB,1),
        for j=1:1:size(MaskRGB,2),
            if MaskRGB(i,j,1)==255, %Change R color for 253
               MaskRGB(i,j,1)=253;
            end
            
            if MaskRGB(i,j,2)==255, %Change G color for 253
               MaskRGB(i,j,2)=253;
            else
            end
            
            if MaskRGB(i,j,3)==255, %Change B color for 150
               MaskRGB(i,j,3)=150; 
            else
            end
        end
    end
    
    imshow(uint8(MaskRGB), 'InitialMagnification', 100); % Display the yellow mask                                                                              (jakie rzu³te, ja nie wierze, no normalnie jak papie¿e xD)
    hold on
    xlim([0, size(Mask,2)]);
    ylim([0, size(Mask,1)]);
    
    % 2. Drawing the fitted ellipse. 
    Beta = -pi * ShapeOrientation/180; %Orientation of the shape (in radians)
    X = Sx1;                            %Center of the ellipse = center of the mask
    Y = Sy1;                                                              
    W2 = ShapeProps.MajorAxisLength/2;  
    M2 = ShapeProps.MinorAxisLength/2;  
    t = linspace(0,2*pi,100);          %Matrix with angles between 0 and 2Pi (100 elements)
    
        %Formulations for X and Y coordinates of the ellipse
    Xe = X + W2*cos(t)*cos(Beta) - M2*sin(t)*sin(Beta);                    
    Ye = Y + W2*cos(t)*sin(Beta) + M2*sin(t)*cos(Beta);
    plot(Xe,Ye,'Color', [0.38, 0.78, 0.42],'Linewidth',2);  %Draw the green ellipse
 
            
        % Begin and end of the Major Axis:
    MajorAxisBeginX = Sx1 - (W2*cos(Beta));
    MajorAxisBeginY = Sy1 - (W2*sin(Beta));
    MajorAxisEndX = Sx1 + (W2*cos(Beta)); 
    MajorAxisEndY = Sy1 + (W2*sin(Beta)); 
   
    MAcoordinates = [MajorAxisBeginX, MajorAxisBeginY; MajorAxisEndX, MajorAxisEndY]; %Saving for calculation of sMM angle.
    Fxy = [Fx, Fy];    % Saving for calculation of sMM angle.
    
        % Draw Major Axix and mark its ends. 
    plot([MajorAxisBeginX, MajorAxisEndX], [MajorAxisBeginY, MajorAxisEndY], 'Color', [0.38, 0.78, 0.42], 'LineWidth', 2, 'MarkerSize',6);  
    plot([MajorAxisBeginX, MajorAxisEndX], [MajorAxisBeginY, MajorAxisEndY], 'xy', 'Color', [0.4, 1, 0], 'LineWidth', 2, 'MarkerSize',4);
        
        % 3. Draw direction of displacement in blue
    plot([Sx1 Fx], [Sy1 Fy], '-b', 'LineWidth', 2);
    
        % 4. Draw current S1 = [Sx1; Sy1] and following S2 = [Sx2; Sy2] centroids   
    plot (Sx1, Sy1, 'g+','LineWidth',1, 'MarkerSize',4);
    plot (Sx2, Sy2, 'm*','LineWidth',1, 'MarkerSize',4);
     
        %%% Save the image in .png file
     if isempty(NameAppendix)
     print(Obrazek, '-dpng', '-r200', [File.pathname 'uMM', '\', File.NameBase(1:end-1), 'uMM', num2str(FrameList(k,2)), '.png']);    
     else
     print(Obrazek, '-dpng', '-r200', [File.pathname 'uMM_', NameAppendix, '\', File.NameBase(1:end-1), 'uMM_', NameAppendix, '_', num2str(FrameList(k,2)), '.png']);        
     end 
    %print(Obrazek, '-dpng', '-r200', [File.pathname 'KatMiM_',  File.NameBase, '\',  File.NameBase(1:end-1), 'KatMiM', num2str(ResultsMain(k,2)), '.png']);
   
   
%% Calculating uMM angle
% This program calculates the uMM angle starting from the image done
% before. I did it to be more sure that it really reflects the actual
% image. 
   
    %---------------------------------------------------------------------
    
    %Load the plot to variables again.
    
    %Draw displacement vector again, to be "on the top layer" of the image
    %(and not be covered by any other element) and save it to separate
    %matrix.
    plot([Sx1 Fx], [Sy1 Fy], '-b', 'LineWidth', 2)
    DirectionMask = getframe(Obrazek);
    DirectionMask = DirectionMask.cdata;
    
    %Similarily to direction vector, draw the Major Axis again, but in red
    %(which is easier in further handling), to make it not covered by any
    %other element and save it to separate matrix.

    plot([MajorAxisBeginX, MajorAxisEndX], [MajorAxisBeginY, MajorAxisEndY], '-r', 'LineWidth', 2, 'MarkerSize',6);  
    MAmask = getframe(gcf);
    MAmask = MAmask.cdata;
    
    close(Obrazek)
    
    %---------------------------------------------------------------------
    % Process the displacement vector image:
    
    % Cropp the image to get rid of this annoying white frame around the
    % image. First, I find the first and last zeros in central column and
    % row (identification of the black background). - For the image of
    % direction vector
    
    %Border identification:
    i=round (size(DirectionMask,2)/2);
    j=round (size(DirectionMask,1)/2);
    BorderVertDirection = DirectionMask(j,:,:);
    BorderLeftDirection=min(find(BorderVertDirection(:,:,1)==0)); 
    BorderRightDirection=max(find(BorderVertDirection(:,:,1)==0));
    BorderHorizDirection = DirectionMask(:,i,:);
    BorderTopDirection=min(find(BorderHorizDirection(:,:,1)==0));
    BorderBottomDirection=max(find(BorderHorizDirection(:,:,1)==0));
    
    %Cropping: 
    DirectionMask2=DirectionMask(BorderTopDirection:BorderBottomDirection,BorderLeftDirection:BorderRightDirection,:);

    %Displacement vector is the only element coloured in blue. Thanks to
    %that I take only those blue pixels and make them white, while turning
    %the rest of them to black.     
    for i=1:1:size(DirectionMask2,1),
        for j=1:1:size(DirectionMask2,2),
            if DirectionMask2(i,j,1)== 0 && DirectionMask2(i,j,2)== 0 && DirectionMask2(i,j,3)== 255; %ChanZmiana pierwszej wartoœci piksela na 245
               DirectionMask2(i,j,:)=255;
            else
               DirectionMask2(i,j,:)=0;
            end
        end
    end
    DirectionMaskBW = im2bw(DirectionMask2, 0.9);
        
    % Saving mask of displacement vector
    if isempty(NameAppendix)
    imwrite(DirectionMaskBW, [File.pathname 'DirectMasks', '\',  File.NameBase(1:end-1), 'Direction', num2str(ResultsMain(k,2)), '.png']);  
    else
    imwrite(DirectionMaskBW, [File.pathname 'DirectMasks_', NameAppendix, '\',  File.NameBase(1:end-1), 'Direction_', NameAppendix, '_', num2str(ResultsMain(k,2)), '.png']);         
    end
    clear i j;

    
    %---------------------------------------------------------------------
    % Process the Major Axis image:
    
    
    %Crop Major Axis mask, just the same as previously. I am aware that the
    %dimensions of those images are absolutely the same, but I do prefer to 
    %work with specific images separately.e na obrazku z osi¹.
    
    i=round (size(MAmask,2)/2);
    j=round (size(MAmask,1)/2);
    BorderVertMA = MAmask(j,:,:);
    BorderLeftMA=min(find(BorderVertMA(:,:,1)==0));
    BorderRightMA=max(find(BorderVertMA(:,:,1)==0));
    BorderHorizMA = MAmask(:,i,:);
    BorderTopMA=min(find(BorderHorizMA(:,:,1)==0));
    BorderBottomMA=max(find(BorderHorizMA(:,:,1)==0));
    
    MAmask2=MAmask(BorderTopMA:BorderBottomMA,BorderLeftMA:BorderRightMA,:);
    
    %Major Axis is the only elerment coloured in red. Thanks to that I take
    %only those red pixels and make them whit, while turning the rest of
    %them to black.    
    for i=1:1:size(MAmask2,1),
        for j=1:1:size(MAmask2,2),
            if MAmask2(i,j,1)== 255 && MAmask2(i,j,2)== 0 && MAmask2(i,j,3)== 0; %Zmiana pierwszej wartoœci piksela na 245
               MAmask2(i,j,:)=255;
            else
               MAmask2(i,j,:)=0;
            end
        end
    end
    MAmaskBW = im2bw(MAmask2, 0.9);
    
    % Save major axis mask:    

    if isempty(NameAppendix)
    imwrite(MAmaskBW, [File.pathname 'MAmasks', '\',  File.NameBase(1:end-1), 'MA', num2str(ResultsMain(k,2)), '.png']);  
    else
    imwrite(MAmaskBW, [File.pathname 'MAmasks_', NameAppendix, '\',  File.NameBase(1:end-1), 'MA_', NameAppendix, '_', num2str(ResultsMain(k,2)), '.png']);         
    end
    clear i j;
    %imwrite(MAmaskBW, [File.pathname 'MAmasks_', File.NameBase, '\',  File.NameBase(1:end-1), 'MA', num2str(ResultsMain(k,2)), '.png']);
   
    clear i j;


    %---------------------------------------------------------------------
    % Calculate uMM angle by reading orientations of BW masks of
    % displacement and major axis.
   
    OrientationDirection=regionprops(DirectionMaskBW ,'Orientation');
    % Correction for the fact, that it is possible for centroid to not move
    % between frames (no displacement). Then I assume that direction vector
    % is the same as in previous frame (because there is no rationale
    % behind changing the direction if the cell didn't move).    
    if ~isempty(OrientationDirection)
        OrientationDirection = OrientationDirection.Orientation
    else
        OrientationDirection = ResultsAuxiliary(k-1, 13);
    end
    
    OrientationMA=regionprops(MAmaskBW ,'Orientation');   
    OrientationMA = OrientationMA.Orientation
   
    % First I define angles between Direction or MA and X-axis. I am
    % interested only in the angles on the right-hand site. However, in
    % Matlab it is always defined between -90 and 90 degrees, so I need to
    % make a corretion for it:
    if OrientationMA < 0;
        OrientationMA = 180 + OrientationMA;
    else
        OrientationMA = OrientationMA;
    end
       
    if OrientationDirection < 0;
             OrientationDirection = 180 + OrientationDirection;
    else
             OrientationDirection = OrientationDirection;
    end  
    
    
    % Calculate uMM:
    if OrientationMA < OrientationDirection;
        uMMangle = OrientationDirection - OrientationMA;
    else
        uMMangle = OrientationMA - OrientationDirection;
    end
    
    % Making uMM the acute angle:    
    if uMMangle > 90;
        uMMangle = 180-uMMangle;
    else
        uMMangle = uMMangle;
    end
    
    % Just in case I coded something wrong:
    if uMMangle > 90
        error(['uMM angle cannot be higher than 90 degrees, however in the file', '\n',...
            [File.NameBase, num2str(k), Settings.FileType], ' it is equal to ',  num2str(uMMangle), '. Check the script.'])
    elseif uMMangle < 0
        error(['uMM angle cannot be lower than 0 degrees, however in the file', '\n',...
            [File.NameBase, num2str(k), Settings.FileType], ' it is equal to ',  num2str(uMMangle), '. Check the script.'])
    else
        uMMangle = uMMangle
    end        
    
    
    else
        
    Mask=imread([File.NameBase, num2str(ResultsMain(k,2)), '.', Settings.FileType]);
    MaskBW = im2bw(Mask, 0.8);    
    ShapeProps=regionprops(MaskBW ,'Centroid', 'MinorAxisLength','MajorAxisLength','Perimeter','Area', 'Orientation'); % Mask's properties

% Orientation and length of Major Axis of the ellipse that was fitted to
% the shape:
    ShapeOrientation = ShapeProps.Orientation;
    MAlength = ShapeProps.MajorAxisLength;
     
    %Sprawdziæ czy to dalej wystêpuje:
    %WspKierOsi = tand(ShapeOrientation); %Funkcja "tand" wylicza tangens dla argumentu W STOPNIACH.            
            
% Cell centroid in current frame [Sx1; Sy1]
    Sx1 = ShapeProps.Centroid(1);
    Sy1 = ShapeProps.Centroid(2);

%% Drawing the image:

% In the first stage program draws mask, fitted ellipse, centroids and both
% ends of major axis.
     
    %%Create mage canvas:
    Obrazek = figure; set(Obrazek, 'name', 'uMM');
    hold on
 
    %1. Painting the mask in yellow xD
    MaskRGB=cat(3, Mask, Mask, Mask);
    
    for i=1:1:size(MaskRGB,1),
        for j=1:1:size(MaskRGB,2),
            if MaskRGB(i,j,1)==255, %Change R color for 253
               MaskRGB(i,j,1)=253;
            end
            
            if MaskRGB(i,j,2)==255, %Change G color for 253
               MaskRGB(i,j,2)=253;
            else
            end
            
            if MaskRGB(i,j,3)==255, %Change B color for 150
               MaskRGB(i,j,3)=150; 
            else
            end
        end
    end
    
    imshow(uint8(MaskRGB), 'InitialMagnification', 100); % Display the yellow mask                                                                              (jakie rzu³te, ja nie wierze, no normalnie jak papie¿e xD)
    hold on
    xlim([0, size(Mask,2)]);
    ylim([0, size(Mask,1)]);
    
    % 2. Drawing the fitted ellipse. 
    Beta = -pi * ShapeOrientation/180; %Orientation of the shape (in radians)
    X = Sx1;                            %Center of the ellipse = center of the mask
    Y = Sy1;                                                              
    W2 = ShapeProps.MajorAxisLength/2;  
    M2 = ShapeProps.MinorAxisLength/2;  
    t = linspace(0,2*pi,100);          %Matrix with angles between 0 and 2Pi (100 elements)
    
        %Formulations for X and Y coordinates of the ellipse
    Xe = X + W2*cos(t)*cos(Beta) - M2*sin(t)*sin(Beta);                    
    Ye = Y + W2*cos(t)*sin(Beta) + M2*sin(t)*cos(Beta);
    plot(Xe,Ye,'Color', [0.38, 0.78, 0.42],'Linewidth',2);  %Draw the green ellipse
 
            
        % Begin and end of the Major Axis:
    MajorAxisBeginX = Sx1 - (W2*cos(Beta));
    MajorAxisBeginY = Sy1 - (W2*sin(Beta));
    MajorAxisEndX = Sx1 + (W2*cos(Beta)); 
    MajorAxisEndY = Sy1 + (W2*sin(Beta)); 
   
    MAcoordinates = [MajorAxisBeginX, MajorAxisBeginY; MajorAxisEndX, MajorAxisEndY]; %Saving for calculation of sMM angle.
    Fxy = [NaN, NaN];    % Saving for calculation of sMM angle.
    
        % Draw Major Axix and mark its ends. 
    plot([MajorAxisBeginX, MajorAxisEndX], [MajorAxisBeginY, MajorAxisEndY], 'Color', [0.38, 0.78, 0.42], 'LineWidth', 2, 'MarkerSize',6);  
    plot([MajorAxisBeginX, MajorAxisEndX], [MajorAxisBeginY, MajorAxisEndY], 'xy', 'Color', [0.4, 1, 0], 'LineWidth', 2, 'MarkerSize',4);
     
        %%% Save the image in .png file
     if isempty(NameAppendix)
     print(Obrazek, '-dpng', '-r200', [File.pathname 'uMM', '\', File.NameBase(1:end-1), 'uMM', num2str(FrameList(k,2)), '.png']);    
     else
     print(Obrazek, '-dpng', '-r200', [File.pathname 'uMM_', NameAppendix, '\', File.NameBase(1:end-1), 'uMM_', NameAppendix, '_', num2str(FrameList(k,2)), '.png']);        
     end 
    %print(Obrazek, '-dpng', '-r200', [File.pathname 'KatMiM_',  File.NameBase, '\',  File.NameBase(1:end-1), 'KatMiM', num2str(ResultsMain(k,2)), '.png']);
    
    %---------------------------------------------------------------------
    % Process the Major Axis image:
    
    plot([MajorAxisBeginX, MajorAxisEndX], [MajorAxisBeginY, MajorAxisEndY], '-r', 'LineWidth', 2, 'MarkerSize',6);  
    MAmask = getframe(gcf);
    MAmask = MAmask.cdata;
    
    close(Obrazek)

    %Crop Major Axis mask, just the same as previously.
    
    i=round (size(MAmask,2)/2);
    j=round (size(MAmask,1)/2);
    BorderVertMA = MAmask(j,:,:);
    BorderLeftMA=min(find(BorderVertMA(:,:,1)==0));
    BorderRightMA=max(find(BorderVertMA(:,:,1)==0));
    BorderHorizMA = MAmask(:,i,:);
    BorderTopMA=min(find(BorderHorizMA(:,:,1)==0));
    BorderBottomMA=max(find(BorderHorizMA(:,:,1)==0));
    
    MAmask2=MAmask(BorderTopMA:BorderBottomMA,BorderLeftMA:BorderRightMA,:);
    
    %Major Axis is the only elerment coloured in red. Thanks to that I take
    %only those red pixels and make them whit, while turning the rest of
    %them to black.    
    for i=1:1:size(MAmask2,1),
        for j=1:1:size(MAmask2,2),
            if MAmask2(i,j,1)== 255 && MAmask2(i,j,2)== 0 && MAmask2(i,j,3)== 0; %Zmiana pierwszej wartoœci piksela na 245
               MAmask2(i,j,:)=255;
            else
               MAmask2(i,j,:)=0;
            end
        end
    end
    MAmaskBW = im2bw(MAmask2, 0.9);
    
    % Save major axis mask:    

    if isempty(NameAppendix)
    imwrite(MAmaskBW, [File.pathname 'MAmasks', '\',  File.NameBase(1:end-1), 'MA', num2str(ResultsMain(k,2)), '.png']);  
    else
    imwrite(MAmaskBW, [File.pathname 'MAmasks_', NameAppendix, '\',  File.NameBase(1:end-1), 'MA_', NameAppendix, '_', num2str(ResultsMain(k,2)), '.png']);         
    end
    clear i j;
    %imwrite(MAmaskBW, [File.pathname 'MAmasks_', File.NameBase, '\',  File.NameBase(1:end-1), 'MA', num2str(ResultsMain(k,2)), '.png']);
       
    uMMangle = NaN
    OrientationDirection = NaN
        
    end
end
    
    % Print the frame number
    %uMMFrame = k   