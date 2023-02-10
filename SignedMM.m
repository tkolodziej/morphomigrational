% Wersja na 18 maja 2020 roku
function [sMMangle, MArotation, MAdynamics, ResultsAuxiliary] = SignedMM(File, Settings, k, ResultsMain, ResultsAuxiliary, NameAppendix, FrameList);
    

%Current and following centroid:
    Sx1 = ResultsMain (k, 4);
    Sx2 = ResultsMain ((k+1), 4);
    Sy1 = ResultsMain (k, 5);
    Sy2 = ResultsMain ((k+1), 5);

% Wyliczam k¹t pod którym jest u³o¿ona aktualna oœ. 
% Coordinates of Major Axis and and of current direction vector.
    Ax = ResultsAuxiliary(k,5);
    Ay = ResultsAuxiliary(k,6);
    Bx = ResultsAuxiliary(k,7);
    By = ResultsAuxiliary(k,8);
    Fx = ResultsAuxiliary(k,9);
    Fy = ResultsAuxiliary(k,10);
    
% Since the beginning of coordinate system of the image is in upper left
% corner (because there is the 1st row and column), ArcTan for 90 and 0 deg
% is reversed. Moreover, the negative angles are transformed to the angle
% complementary to 180 degrees. 
  if Ax == Bx
    ArcTanMA = 90;
  elseif Ay == By
    ArcTanMA = 0;
  else
    ArcTanMA=-atand((By-Ay)/(Bx-Ax));
  end 
    
    MAend1 = [0,0];
    MAend2 = [0,0];
    % Order points in the way, to have lower y-coordinate in first row. If
    % Ay = By, then smaller x should be in first row.
    if Ay < By
        MAend1(1,1) = Ax;
        MAend1(1,2) = Ay;
        MAend2(1,1) = Bx;
        MAend2(1,2) = By;
    elseif (Ay == By) && (Ax < Bx)
        MAend1(1,1) = Ax;
        MAend1(1,2) = Ay;
        MAend2(1,1) = Bx;
        MAend2(1,2) = By;
    elseif (Ay == By) && (Ax > Bx)
        MAend1(1,1) = Bx;
        MAend1(1,2) = By;
        MAend2(1,1) = Ax;
        MAend2(1,2) = Ay;
    elseif Ay > By
        MAend1(1,1) = Bx;
        MAend1(1,2) = By;
        MAend2(1,1) = Ax;
        MAend2(1,2) = Ay;
    end    

    
% Find rotation point, around which the cell should be rotated to decide about the sign of sMM angle    

% The purpose of choosing particilar rotation poin is chosen in the SI to
% the paper Kolodziej et al., 2023, "Morphomigrational description - a new 
% approach connecting cell migration with its morphology

% For the first frame - rotation point is the one, that is lower in the
% image (higher Y value). Then, the major axis is rotated clockwise to the
% neares horizontal position.
if k==1   
    if MAend1(1,2) ~= MAend2(1,2)
        RotationPoint = [MAend2(1,1), MAend2(1,2)];    
    elseif MAend1(1,2) == MAend2(1,2)
        RotationPoint = [MAend1(1,1), MAend1(1,2)];
    else
        error(['Error in ', [File.NameBase, num2str(ResultsMain(k, 2)), ...
        Settings.FileType], 'cell - rotation point is not properly calculated in 1st frame.',...
        '\n', 'Check the validity of script and input file'])
    end
    
    if ArcTanMA < 0
        MArotation = 180 + ArcTanMA;
    else
        MArotation = ArcTanMA;
    end
    MArotationTemp = MArotation;
    
% For all next frames
elseif k > 1

%Values of sMM angle and MA orientation in previous frame 
        PreviousMM = ResultsMain(k-1, 13);
        OrientationPrevMA = ResultsAuxiliary(k-1, 11);

    % Badam pod jakim k¹tem jest nachylona oœ. Mo¿na przyj¹æ zawsze dwie
    % opcje ró¿ni¹ce sie o 180 stopni. Np. oœ mo¿e byæ pochylona pod k¹tem
    % 45 lub 135 stopni, zale¿y od którego punktu bêdziemy mierzyæ. Tutaj
    % muszê wyliczyæ oba te k¹ty, potem bêdê je porównywaæ.
    
    % What is the orientation of current MA regarding the X-axis? For each
    % line, there are always two options: higrer or lower than 180 degrees
    % (for example 45 or 325 degrees), depending from which point it is
    % measured. Thus both of them is calculated for further comparison.
                
        Angle1 = (360*floor(OrientationPrevMA/360)) + ArcTanMA;
        if (Angle1 - OrientationPrevMA) >= 180
            Angle1 = Angle1 - 360;
        elseif (Angle1 - OrientationPrevMA) <= -180
            Angle1 = Angle1 + 360;
        else
            Angle1 = Angle1;
        end

        Angle2 = (360*floor(OrientationPrevMA/360)) +180 + ArcTanMA;
        if (Angle2 - OrientationPrevMA) >= 180
            Angle2 = Angle2 - 360;
        elseif (Angle2 - OrientationPrevMA) <= -180
            Angle2 = Angle2 + 360;
        else
            Angle2 = Angle2;
        end
    
                                         
  
    
    AngleDiff1 = abs(Angle1 - OrientationPrevMA);
    AngleDiff2 = abs(Angle2 - OrientationPrevMA);
    
    % Choosing the lower difference that is marked as final MArotation    
    if (AngleDiff1 < AngleDiff2) && (ArcTanMA ~= 0)
        MArotation = Angle1;
    elseif (AngleDiff1 < AngleDiff2) && (ArcTanMA == 0)
        MArotation = Angle1;
    elseif (AngleDiff2 < AngleDiff1) && (ArcTanMA ~= 0)
        MArotation = Angle2;
    elseif (AngleDiff2 < AngleDiff1) && (ArcTanMA == 0)
        MArotation = Angle2;
    end
    
    %Here the rotation point is chosen. If it is between [0;180] or [-180;
    %-360] deg, or their multiple, then the higher Y is assumed. If it is
    %between [180; 360] or p0;-180 deg], the lower Y is assumed.
    
    if MArotation < 0 
        if MArotation <= -360 %Reducing the multiple of 360 deg to only one rotation
            MArotationTemp = MArotation + (floor(MArotation/(-360))*360)
            if MArotationTemp > -360 && MArotationTemp <= 0 % It checks if the rotation angle was properly reduced.
            fprintf(['Reduction of the angle below -360 deg done properly in ', '\n'...
                [File.NameBase, num2str(ResultsMain(k, 2)), Settings.FileType], ' from ', num2str(MArotation), ...
                ' to ' num2str(MArotationTemp), ' deg.']);
            else
            error(['Error while reducing Temporary MA rotation in ', '\n',...
            [File.NameBase, num2str(ResultsMain(k, 2)), Settings.FileType], ...
            ' for MA Rotation = ', num2str(MArotation), '. Check the code.']);
            end
        elseif MArotation > -360
            MArotationTemp = MArotation
        else
            error(['Error while reducing Temporary MA rotation in ', '\n',...
            [File.NameBase, num2str(ResultsMain(k, 2)), Settings.FileType], '. Check the code.'])
        end
        
        
        if 0 > MArotationTemp && MArotationTemp > -180
            RotationPoint = [MAend1(1,1), MAend1(1,2)];
        elseif -180 >= MArotationTemp && MArotationTemp > -360
            RotationPoint = [MAend2(1,1), MAend2(1,2)];
        else 
            error(['Error while calculating the rotation point in ', '\n',...
             [File.NameBase, num2str(ResultsMain(k, 2)), Settings.FileType], '. Check the code.'])
        end
        
    elseif MArotation >= 0
        
        if MArotation >= 360 %Reducing the multiple of 360 deg to only one rotation
            MArotationTemp = MArotation - (floor(MArotation/(360))*360)
            if MArotationTemp < 360 && MArotationTemp >= 0 % It checks if the rotation angle was properly reduced.
            fprintf(['Reduction of the angle below -360 deg done properly in ', '\n'...
                [File.NameBase, num2str(ResultsMain(k, 2)), Settings.FileType], ' from ', num2str(MArotation), ...
                ' to ' num2str(MArotationTemp), ' deg.']);
            else
            error(['Error while reducing Temporary MA rotation in ', '\n',...
            [File.NameBase, num2str(ResultsMain(k, 2)), Settings.FileType], ...
            ' for MA Rotation = ', num2str(MArotation), '. Check the code.']);
            end
        elseif MArotation < 360
            MArotationTemp = MArotation
        else
            error(['Error while reducing Temporary MA rotation in ', '\n',...
            [File.NameBase, num2str(ResultsMain(k, 2)), Settings.FileType], '. Check the code.'])
        end
        
        if MArotationTemp == 0
            RotationPoint = [MAend1(1,1), MAend1(1,2)];
        elseif 0 <MArotationTemp && MArotationTemp <= 180
            RotationPoint = [MAend2(1,1), MAend2(1,2)]
        elseif 180 < MArotationTemp && MArotationTemp < 360
            RotationPoint = [MAend1(1,1), MAend1(1,2)]
        else 
            error(['Error while calculating the rotation point in ', '\n',...
            [File.NameBase, num2str(ResultsMain(k, 2)), Settings.FileType], ' sprawdŸ poprawnoœæ kodu.'])
        end
    end   

    %Ok, the rotation angle and rotation points are calculated ...at this
    %point... xD
    
    
end

    %As mentioned in the paper, the MA is rotated clockwise. However, I do
    %not need to change the sign of an angle, because of the different
    %coordinate system (that starts in upper-left corner of the image).
    
    % Rotating the Major Axis at MArotation. First I'll rotate the second
    % and of the MA to find if it has then the same Y value (thus if MA
    % really lies horizontally after the rotation. Just to double-check it.
    Px = RotationPoint(1,1);
    Py = RotationPoint(1,2);
    
    if isequal(RotationPoint,MAend1); % If the first end is the rotation point, then we rotate the second end.
    
        XaxisMAend2a = Px + (MAend2(1,1)-Px)*cosd(MArotationTemp)-(MAend2(1,2)-Py)*sind(MArotationTemp);
        YaxisMAend2a = Py + (MAend2(1,1)-Px)*sind(MArotationTemp)+(MAend2(1,2)-Py)*cosd(MArotationTemp);
    
        if round(YaxisMAend2a) == round(RotationPoint(1,2)) &&...
                (round(sqrt((Px - XaxisMAend2a)^2 + (Py - YaxisMAend2a)^2)) == round(sqrt((Bx - Ax)^2 + (By - Ay)^2)))
            fprintf(['Rotation of MA done properly in ',...
                [File.NameBase, num2str(ResultsMain(k, 2)), Settings.FileType], ' by  ', num2str(MArotationu), ' deg.']);
        else
            error(['Error while rotating the MA in ', [File.NameBase, num2str(ResultsMain(k, 2)), Settings.FileType], ' check the code.']);        
        end
        MAendAfterRotation = [XaxisMAend2a, YaxisMAend2a];
        
    elseif isequal(RotationPoint,MAend2);
        
        Pocz1aXOsi = Px + (MAend1(1,1)-Px)*cosd(MArotationTemp)-(MAend1(1,2)-Py)*sind(MArotationTemp);
        Pocz1aYOsi = Py + (MAend1(1,1)-Px)*sind(MArotationTemp)+(MAend1(1,2)-Py)*cosd(MArotationTemp);
        
        % W warunku poni¿ej sprawdzam czy obrócona oœ jest pozioma, oraz
        % czy nie zmieni³a siê jej d³ugoœæ - takie sprawdzenie kodu.
        if (round(Pocz1aYOsi) == round(RotationPoint(1,2))) &&...
                (round(sqrt((Px - Pocz1aXOsi)^2 + (Py - Pocz1aYOsi)^2)) == round(sqrt((Bx - Ax)^2 + (By - Ay)^2)))
            fprintf(['Rotation of MA done properly in ',...
                [File.NameBase, num2str(ResultsMain(k, 2)), Settings.FileType], ' by  ', num2str(MArotation), ' deg.']);
                
        else
            error(['Error while rotating the MA in ', [File.NameBase, num2str(ResultsMain(k, 2)), Settings.FileType], ' check the code.']);      
        end
       
        MAendAfterRotation = [Pocz1aXOsi, Pocz1aYOsi];
    else
        error(['Something is WROOOOOOOONG in the code while rotating the MA in ', ...
            [File.NameBase, num2str(ResultsMain(k, 2)), Settings.FileType], ' this code has to be fixed!!!']);
        
    end  
    
    
    %Here, the next centroid is rotated around the Rotation Point. Then we
    %chack if it is above or below the Rotation Point in the picture
    %(higher Y value) or above it (lower Y value). 
    Sx1AfterRotation = Px + (Sx1-Px)*cosd(MArotation)-(Sy1-Py)*sind(MArotation);
    Sy1AfterRotation = Py + (Sx1-Px)*sind(MArotation)+(Sy1-Py)*cosd(MArotation);
    Sx2AfterRotation = Px + (Sx2-Px)*cosd(MArotation)-(Sy2-Py)*sind(MArotation);
    Sy2AfterRotation = Py + (Sx2-Px)*sind(MArotation)+(Sy2-Py)*cosd(MArotation);
    FxAfterRotation = Px + (Fx-Px)*cosd(MArotation)-(Fy-Py)*sind(MArotation);
    FyAfterRotation = Py + (Fx-Px)*sind(MArotation)+(Fy-Py)*cosd(MArotation);
    
    S2AfterRotation = [Sx2AfterRotation, Sy2AfterRotation]
        
    if Sy2AfterRotation < Py
        sMMTemp = ResultsMain(k, 12);
    elseif Sy2AfterRotation > Py
        sMMTemp = -(ResultsMain(k, 12));
    else %correction for the situation if next centroid is on the Major Axis. In such situation the previous sign of sMM angle is assigned. What if the previous one was not defined? Well... it will make program crush, but I'll deal with it if I find such situation.
        if sign(PreviousMM) < 0
            sMMTemp = -(ResultsMain(k, 12));
        elseif sign(PreviousMM) > 0
            sMMTemp = (ResultsMain(k, 12));
        end
    end
   
   ResultsAuxiliary(k,10) = sMMTemp; 
   
   % I want to make first sMM angle in the sequence positive - regardless
   % if it is below or above the rotated Major Axis. It will change the
   % sign only in the situation if the direction will point to another side
   % of MA. The ifs below do this job:
   if k==1
       if sMMTemp < 0
        sMMangle = -sMMTemp;
       elseif sMMTemp >=0
        sMMangle = sMMTemp;
       end
   elseif ResultsAuxiliary(1, 10) < 0
       sMMangle = -sMMTemp;
   elseif ResultsAuxiliary(1,10) > 0 
       sMMangle = sMMTemp;
   end
   
   
   % Now it's drawing time! Let's draw mask, with axes and directions
   % before and after rotation! Including the centroid of the following
   % mask!
   
   Mask=imread([File.NameBase, num2str(ResultsMain(k, 2)),  '.', Settings.FileType]); %Load Mask file
   Mask=im2bw(Mask, 0.5);
   
   
   %-These ifs decides about where the legend should be.
   ImageCenter = round(size(Mask)/2);
   if Sx2 < ImageCenter(2) && Sy2 < ImageCenter(1);
      LegendPlacement = 'southeast';
   elseif Sx2 > ImageCenter(2) && Sy2 < ImageCenter(1);
      LegendPlacement = 'southwest';
   elseif (Sx2 > ImageCenter(2)) && Sy2 > ImageCenter(1);
      LegendPlacement = 'northwest';
   elseif Sx2 < ImageCenter(2) && Sy2 > ImageCenter(1);
      LegendPlacement = 'northeast';
   end
   
   %-Ifs for writing the comment about sMM angle
   if sMMTemp > 0
        comment = ' (above the Major Axis)';
   elseif sMMTemp ==0
        comment = ' (on the Major Axis) ';
   elseif sMMTemp < 0
        comment = ' (below the Major Axis)';
   end 
   sMMangleText=sprintf('%.2f',sMMangle);%Turn the number into the text with two decimal places
   MArotationText = sprintf('%.2f',MArotation);
   FrameNumber = num2str(ResultsMain(k, 2))
         
  % Wybór zapisywania obrazka:
        
      
   %-Rysujê wykres:
   Figura = figure, ;
   imshow(Mask, 'InitialMagnification', 100);

   SequenceName = sprintf(File.NameBase);
   
   hold on
        title(['Frame ', FrameNumber, '. ', 'sMM = ', sMMangleText, char(176), comment, '. ', 'Rotation of ', MArotationText, char(176)], 'FontSize', 10);
        xlabel([SequenceName], 'Interpreter', 'none')

        v1=plot([Ax Bx], [Ay By], '-b', 'LineWidth', 2,'DisplayName', 'MA');
        v2=plot([Sx1 Fx], [Sy1, Fy], '--b', 'LineWidth', 2, 'DisplayName', 'Direction');  
        v3=plot([Px MAendAfterRotation(1)], [Py MAendAfterRotation(2)], '-m', 'LineWidth', 2, 'DisplayName', 'Rotated MA');   
        v4=plot([Sx1AfterRotation, FxAfterRotation], [Sy1AfterRotation, FyAfterRotation], '--m', 'LineWidth', 2, 'DisplayName', 'Rotated direction');    

        p1=plot(Sx2, Sy2, 'c+', 'LineWidth', 2, 'MarkerSize',6, 'DisplayName', 'Next centroid');
        p2=plot(Sx2AfterRotation, Sy2AfterRotation, 'g*', 'MarkerSize',5, 'DisplayName', 'Rotated next centroid');        
        l=legend({'MA', 'Direction', 'Rotated MA', 'Rotated direction','Next centroid', 'Rotated next centroid'}, 'FontSize',5, 'Location', LegendPlacement);
        xlabel([SequenceName], 'Interpreter', 'none')
        

        
        if isempty(NameAppendix)
        print(Figura, '-dpng', '-r200', [File.pathname 'sMM', '\', File.NameBase(1:end-1), 'sMM', num2str(FrameList(k,2)), '.png']);    
        else
        print(Figura, '-dpng', '-r200', [File.pathname 'sMM_', NameAppendix, '\', File.NameBase(1:end-1), 'uMM_', NameAppendix, '_', num2str(FrameList(k,2)), '.png']);        
        end 
        %print(Figura, '-dpng', '-r200', [File.pathname 'sMM_', File.PodstawaNazwy, '\',  File.PodstawaNazwy(1:end-1), 'sMM_', num2str(ResultsMain(k,2)), '.png']);
   
   hold off
   close(Figura)        
          
        
   % Calculating the M.A. dynamics
   MAorientationCurr = ResultsAuxiliary(k, 13);
   MAorientationFoll = ResultsAuxiliary(k+1, 13);
   
   MAdynamics = -(MAorientationFoll-MAorientationCurr); % Put minus sign, because I want the positive sign to mark clockwise turn

   sMMangle
   MArotation
   clear Fx Fy FxAfterRotation FyAfterRotation
end