% Version for 13.02.2023, Tomasz Kołodziej
% This program calculates displacements, velocities and turning angles,
% plotting the latter as well.

function [ResultsMain, ResultsAuxiliary] = VelocitiesAngles(File, Settings, ResultsMain, ResultsAuxiliary, NameAppendix, FrameList);
   
AnglesInitial = zeros(Settings.NumFrames, 1);

for k=1:1:(size(ResultsMain,1))-1
    
    dX = ResultsMain(k+1,4) - ResultsMain((k), 4); % Calculate X-displacement in pixels
    dY = ResultsMain(k+1,5) - ResultsMain((k), 5); % Calculate Y-displacement in pixels
    dXum = dX*Settings.PixelSize; % Calculate X-displacements in microns
    dYum = dY*Settings.PixelSize; % Calculate Y-displacements in microns
    r = sqrt((dX)^2 + (dY)^2); % Displacement in pixels
    rum = sqrt((dX)^2 + (dY)^2)*Settings.PixelSize; % Displacement in microns
    v = rum/((ResultsMain(k+1,3) - ResultsMain(k, 3))); %Velocity in microns/min
    
    % Copy those values to ResultsMain.
    ResultsMain(k, 8) = r;
    ResultsMain(k, 9) = v;
    
    
%% Turning angle calculation

% Here I calculate turning angles that is defined as an angle between
% two consecutive displacement vectors: first one between previous and
% current frame (n-1 and n frames) and second one between current and
% following frames (n and n+1 frames). 
% As I am not a programmer at all (or et al. xD), so I do it as I would
% do it on a paper. First I check in which direction the
% displacement vector points on a 2D plane. Then, using atan function,
% I calculate the angle between this vector and X axis, always on the
% right side of the vector (rotating it counterclockwise).
% All those wicked conditions below come from the fact that coordinate
% system of the image starts in upper-left corner which comes from the
% numeration of rows and columns. And I wanted to make it just like
% would it be on the paper.
    
    % Calculated angles are stored in temporary matrix
    
    if dXum > 0 && dYum ==0;
        AngleXaxis = 0;
    elseif dXum > 0 && dYum < 0;
        AngleXaxis = -(atand(dYum/dXum));
    elseif dXum == 0 && dYum < 0;
        AngleXaxis = 90;
    elseif dXum < 0 && dYum < 0;
        AngleXaxis = 180 - (atand(dYum/dXum));
    elseif dXum < 0 && dYum == 0;
        AngleXaxis = 180;
    elseif dXum < 0 && dYum > 0;
        AngleXaxis = 180 - (atand(dYum/dXum));
    elseif dXum == 0 && dYum > 0;
        AngleXaxis = 270;
    elseif dXum > 0 && dYum > 0;
        AngleXaxis = 360 - (atand(dYum/dXum));
    end
    
    AnglesInitial(k, 1) = AngleXaxis;
    ResultsAuxiliary(k,4)= AngleXaxis;
    k;

end

for k=2:1:(size(ResultsMain,1)-1)
    % In case there's no displacement, we assume that the angle between the
    % vector and X-axis is the same as in previous frame, because standard
    % formula gives NaN. 
    if ResultsMain(k,7)==0
        AnglesInitial(k,1) = AnglesInitial(k-1,1);
    else
        AnglesInitial(k,1) = AnglesInitial(k,1);
    end
    
    % Calculating turning angle by substracting consecutive angle between
    % displacement vectors and X-axis.
    AnglesDifference = AnglesInitial(k, 1) - AnglesInitial (k-1, 1);
    
    if AnglesDifference<(-180) 
        AnglesDifference = AnglesDifference + 360;
    elseif AnglesDifference > 180 
        AnglesDifference = AnglesDifference - 360;
    else
        AnglesDifference = AnglesDifference;
    end
    
    % Now I change the sign of the angle, because I want it to positive
    % while turning clockwise and negative for counterclockwise.    
    AnglesDifferenceFinal = -AnglesDifference;
    ResultsMain(k, 10) = AnglesDifferenceFinal;
    
    clear k
end

%% Plotting the turning angles
% The loop below is solely for drawing images of turning angles. Can be
% deactivated without any damage to the data.
for k=2:1:ResultsMain(end, 1)-1


    %Loading the previous, current and following mask
    MaskPrev=imread([File.NameBase, num2str(ResultsMain(k-1,2)), '.', Settings.FileType]);
    MaskPrevBW = im2bw(MaskPrev, 0.8);
    MaskCurr=imread([File.NameBase, num2str(ResultsMain(k,2)), '.',  Settings.FileType]);
    MaskCurrBW = im2bw(MaskCurr, 0.5);
    MaskFoll=imread([File.NameBase, num2str(ResultsMain(k+1,2)), '.',  Settings.FileType]);
    MaskFollBW = im2bw(MaskFoll, 0.5);
    
    ShapePropsPrev=regionprops(MaskPrevBW , 'MajorAxisLength', 'MinorAxisLength'); % Load properties of masks
    ShapePropsCurr=regionprops(MaskCurrBW , 'MajorAxisLength', 'MinorAxisLength');
    
    % Sx1, Sy1 - centroids of previous mask
    Sx1 = ResultsMain(k-1, 4); 
    Sy1 = ResultsMain(k-1, 5);
    DispLenghtPrev = ShapePropsPrev.MajorAxisLength/2; %To be visible on the picture, the vector showing displacement direction has to be longer then just a displacement.  It is set for the half of the lenght of the Major Axis.
    
    % Sx2, Sy2 - centroids of current mask
    Sx2 = ResultsMain(k, 4);
    Sy2 = ResultsMain(k, 5);
    DispLenghtCurr = ShapePropsCurr.MajorAxisLength/2;
    
    % Sx3, Sy3 - centroids of following mask
    Sx3 = ResultsMain(k+1, 4);
    Sy3 = ResultsMain(k+1, 5);
    

% --- Calculating first vector of displacement direction, that passes
% --- through previous centroid [Sx1; Sy1] and current centroid [Sx2; Sy2],
% --- having a lenght of DispLenghtPrev
    
    % Calculate the slope:
    SlopePrev = (Sy2 - Sy1)/(Sx2 - Sx1);
    SlopeAnglePrev = -(atan(SlopePrev));

    % Vector of displacement direction ends in point [Fx1, Fy1] that is
    % calculated here:    
    % Calculating Fx
        if Sx1 == Sx2;
           Fx = Sx1;
        elseif Sx1 < Sx2;
           Fx = Sx1 + (DispLenghtPrev*cos(SlopeAnglePrev));
        elseif Sx1 > Sx2;
           Fx = Sx1 - (DispLenghtPrev*cos(SlopeAnglePrev));
        end
               
    % Calculating Fy    
        if Sy1 == Sy2;
           Fy = Sy1;
        elseif Sx1 > Sx2 | (Sx1 == Sx2 && Sy1 < Sy2);
           Fy = Sy1 + (DispLenghtPrev*sin(SlopeAnglePrev));
        elseif Sx1 < Sx2 | (Sx1 == Sx2 && Sy1 > Sy2);
           Fy = Sy1 - (DispLenghtPrev*sin(SlopeAnglePrev));
        end
    
    % Correction for Fx or Fy that may exceed the image dimensions:    
      if (Fx < 0 | Fx > size(MaskCurr, 2)) | Fy < 0 | Fy > size(MaskCurr, 1);
          DispLenghtPrev = ShapePropsCurr.MinorAxisLength;
            % Just the same calculations as before, but lenght of the vector
            % is chosen as Minor and not the Major axis of the fitted
            % ellipse.
            % Calculating Fx:
        if Sx1 == Sx2;
           Fx = Sx1;
        elseif Sx1 < Sx2;
           Fx = Sx1 + (DispLenghtPrev*cos(SlopeAnglePrev));
        elseif Sx1 > Sx2;
           Fx = Sx1 - (DispLenghtPrev*cos(SlopeAnglePrev));
        end
        
            % Calculating Fy:
        if Sy1 == Sy2;
           Fy = Sy1;
        elseif Sx1 > Sx2 | (Sx1 == Sx2 && Sy1 < Sy2);
           Fy = Sy1 + (DispLenghtPrev*sin(SlopeAnglePrev));
        elseif Sx1 < Sx2 | (Sx1 == Sx2 && Sy1 > Sy2);
           Fy = Sy1 - (DispLenghtPrev*sin(SlopeAnglePrev));
        end
      else % if the previous Fx and Fy are fine, they stay the same.
           Fx = Fx; 
           Fy = Fy;
      end
    
% --- First displacement direction is calculated, now it's time to
% --- calculate the second one, that goes through current centroid [Sx2,
% --- Sy2] and following centroid [Sx3, Sy3]. Calculations are done the
% --- same way

      
      
    % Calculate the slope:
    SlopeCurr = (Sy3 - Sy2)/(Sx3 - Sx2);
    SlopeAngleCurrt = -(atan(SlopeCurr));
    % Vector of displacement direction ends in point [Fx2, Fy2] that is
    % calculated here:    
    % Calculating Fx2
        if Sx2 == Sx3;
           Fx2 = Sx2
        elseif Sx2 < Sx3;
           Fx2 = Sx2 + (DispLenghtCurr*cos(SlopeAngleCurrt));
        elseif Sx2 > Sx3;
           Fx2 = Sx2 - (DispLenghtCurr*cos(SlopeAngleCurrt));
        end
               
    % Calculating Fy2   
        if Sy3 == Sy2;
           Fy2 = Sy2;
        elseif Sx2 > Sx3 | (Sx2 == Sx3 && Sy2 < Sy3);
           Fy2 = Sy2 + (DispLenghtCurr*sin(SlopeAngleCurrt));
        elseif Sx2 < Sx3 | (Sx2 == Sx3 && Sy2 > Sy3);
           Fy2 = Sy2 - (DispLenghtCurr*sin(SlopeAngleCurrt));
        end
    
    % Correction for Fx or Fy that may exceed the image dimensions:     
       if (Fx2 < 0 | Fx2 > size(MaskCurr, 2)) | Fy2 < 0 | Fy2 > size(MaskCurr, 1);
          DispLenghtCurr = DispLenghtCurr/2;
            % Just the same calculations as before, but lenght of the vector
            % is chosen as Minor and not the Major axis of the fitted
            % ellipse.
            % Calculating Fx:
           if Sx2 == Sx3;
            Fx2 = Sx2;
           elseif Sx2 < Sx3;
            Fx2 = Sx2 + (DispLenghtCurr*cos(SlopeAngleCurrt));
           elseif Sx2 > Sx3;
           Fx2 = Sx2 - (DispLenghtCurr*cos(SlopeAngleCurrt));
           end
            % Calculating Fy:
           if Sy3 == Sy2;
            Fy2 = Sy2
           elseif Sx2 > Sx3 | (Sx2 == Sx3 && Sy2 < Sy3);
            Fy2 = Sy2 + (DispLenghtCurr*sin(SlopeAngleCurrt));
           elseif Sx2 < Sx3 | (Sx2 == Sx3 && Sy2 > Sy3);
            Fy2 = Sy2 - (DispLenghtCurr*sin(SlopeAngleCurrt));
           end
       end

% --- Drawing the image:
        % Previous and following masks are represented by their perimeters,
        % current mask is filled with white color. 
        
        MaskPrevPerimeter = bwperim(MaskPrevBW, 4);
        MaskCurrBW = MaskCurrBW;
        im1 = zeros(size(MaskCurrBW,1), size(MaskCurrBW,2), 3);
        im1=uint8(im1);
        
        for i=1:1:size(MaskCurrBW,1),%(255,165,0)
            for j=1:1:size(MaskCurrBW,2),
                if MaskCurrBW(i,j)==0 && MaskPrevPerimeter(i,j)==0; % Leave black pixel
                    im1(i,j,1)=0;
                    im1(i,j,2)=0;
                    im1(i,j,3)=0;
                elseif MaskCurrBW(i,j)==0 && MaskPrevPerimeter(i,j)==1; 
                    im1(i,j,1)=255;
                    im1(i,j,2)=0;
                    im1(i,j,3)=0;
                elseif MaskCurrBW(i,j)==1 && MaskPrevPerimeter(i,j)==1;
                    im1(i,j,1)=255;
                    im1(i,j,2)=165;
                    im1(i,j,3)=0;
                elseif MaskCurrBW(i,j)==1 && MaskPrevPerimeter(i,j)==0;
                    im1(i,j,1)=255;
                    im1(i,j,2)=255;
                    im1(i,j,3)=255;
                end
            end
        end
        
        % Adding the following frame (cell perimeter)
        MaskFollPerimeter = bwperim(MaskFollBW, 4);
        im2 = zeros(size(im1));
        im2=uint8(im2);
       
        for i=1:1:size(im1,1),
            for j=1:1:size(im1,2),
                if (im1(i,j,1)==0 && im1(i,j,2)==0 && im1(i,j,3)==0) && MaskFollPerimeter(i,j)==0; % Leave black pixel
                    im2(i,j,1)=0;
                    im2(i,j,2)=0;
                    im2(i,j,3)=0;
                elseif (im1(i,j,1)==255 && im1(i,j,2)==255 && im1(i,j,3)==255) && MaskFollPerimeter(i,j)==0; % Leave white pixel
                    im2(i,j,1)=255;
                    im2(i,j,2)=255;
                    im2(i,j,3)=255;
                elseif (im1(i,j,1)==255 && im1(i,j,2)==0 && im1(i,j,3)==0) && MaskFollPerimeter(i,j)==0; % Leave red pixel
                    im2(i,j,1)=255;
                    im2(i,j,2)=0;
                    im2(i,j,3)=0;   
                elseif (im1(i,j,1)==255 && im1(i,j,2)==165 && im1(i,j,3)==0) && MaskFollPerimeter(i,j)==0; % Leave orange pixel (at the intersections)
                    im2(i,j,1)=255;
                    im2(i,j,2)=0;
                    im2(i,j,3)=0;                             
                elseif (im1(i,j,1)==0 && im1(i,j,2)==0 && im1(i,j,3)==0) && MaskFollPerimeter(i,j)==1; % Colour the perimeter with dark green
                    im2(i,j,1)=0;
                    im2(i,j,2)=128;
                    im2(i,j,3)=0;
                elseif (im1(i,j,1)==255 && im1(i,j,2)==255 && im1(i,j,3)==255) && MaskFollPerimeter(i,j)==1; % Intersection of following perimeter and current mask
                    im2(i,j,1)=50;
                    im2(i,j,2)=205;
                    im2(i,j,3)=50;
                elseif (im1(i,j,1)==255 && im1(i,j,2)==165 && im1(i,j,3)==0) && MaskFollPerimeter(i,j)==1; % Leave orange pixel again
                    im2(i,j,1)=255;
                    im2(i,j,2)=165;
                    im2(i,j,3)=0;
                elseif (im1(i,j,1)==255 && im1(i,j,2)==165 && im1(i,j,3)==0) && MaskPrevPerimeter(i,j)==1; % Intersection of orange and green -> yellow
                    im2(i,j,1)=255;
                    im2(i,j,2)=255;
                    im2(i,j,3)=0;
                elseif (im1(i,j,1)==255 && im1(i,j,2)==0 && im1(i,j,3)==0) && MaskPrevPerimeter(i,j)==1; % Intersection of previous and following frame is pure gold xD
                    im2(i,j,1)=255;
                    im2(i,j,2)=215;
                    im2(i,j,3)=0;
                end
            end
        end

        Figura = figure; 
        imshow(im2, 'InitialMagnification', 100);
        hold on
        
        
        CurrentAngle = ResultsMain(k, 10)
        if CurrentAngle > 0
            Sign = '+';
            Comment = ' (clockwise)';
        elseif CurrentAngle ==0
            Sign = ' ';
            Comment = ' ';
        else
            Sign = '-';
            Comment = ' (counterclockwise)';
        end
        
        CurrentAngle=abs(CurrentAngle); %Adding the sign just as a comment.
        CurrentAngle=sprintf('%.2f',CurrentAngle);

        FrameNumber=num2str(FrameList(k,2))
        % Positioning the legend to make it not cover the image of the
        % cell.
            ImageCenter = round(size(MaskCurr)/2);
            if Sx2 < ImageCenter(2) && Sy2 < ImageCenter(1);
                LegendPosition = 'southeast';
            elseif Sx2 > ImageCenter(2) && Sy2 < ImageCenter(1);
                LegendPosition = 'southwest';
            elseif (Sx2 > ImageCenter(2)) && Sy2 > ImageCenter(1);
                LegendPosition = 'northwest';
            elseif Sx2 < ImageCenter(2) && Sy2 > ImageCenter(1);
                LegendPosition = 'northeast';
            end
                
        
        title(['Frame ', FrameNumber, '. ', 'angle = ', Sign, CurrentAngle, char(176), Comment]);  

        ax = gca
        ax.TitleFontSizeMultiplier = 0.8;
        v1=plot([Sx1 Fx], [Sy1 Fy], '-b', 'LineWidth', 2,'DisplayName', 'Previos direction');
        v2=plot([Sx2 Fx2], [Sy2 Fy2], '-m', 'LineWidth', 2, 'DisplayName', 'Current direction');
        plot (Sx1, Sy1, 'r*','LineWidth',1, 'MarkerSize',3);
        plot (Sx2, Sy2, 'ko','LineWidth',1, 'MarkerSize',4);
        plot (Sx3, Sy3, 'g*','LineWidth',1, 'MarkerSize',3);
        l=legend({'Prev direction', 'Curr Direction', 'Prev Centroid', 'Curr Centroid', 'Follow Centroid'}, 'FontSize',5, 'Location', LegendPosition); 
        
        if isempty(NameAppendix)
        SequenceName = sprintf(File.NameBase)
        else
        SequenceName = [sprintf(File.NameBase), ' ', NameAppendix]
        end
        xlabel([SequenceName], 'Interpreter', 'none')   
     
        
     if isempty(NameAppendix)
     print(Figura, '-dpng', '-r200', [File.pathname 'TurnAngle', '\', File.NameBase(1:end-1), 'TurnAngle', num2str(FrameList(k,2)), '.png']);    
     else
     print(Figura, '-dpng', '-r200', [File.pathname 'TurnAngle_', NameAppendix, '\', File.NameBase(1:end-1), 'TurnAngle_', NameAppendix, '_', num2str(FrameList(k,2)), '.png']);        
     end
     
   hold off
   close(Figura);

   
end
%    clear AnglesInitial
%    clear k
end