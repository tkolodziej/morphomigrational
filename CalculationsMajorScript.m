% Version for 13.02.2023, Tomasz Ko³odziej
% Script calculates sMM angle and other parameters; creates images of uMM
% angle, sMM angle and turning angle; masks of Major Axis and Displacement.

function [File] = CalculationsMajorScript(File, Settings)

%Add the information of sampling to the file names and folders
if Settings.Sampling > 1
    NameAppendix = ['Every', num2str(Settings.Sampling), 'Frames'];
    % Create folders for additional images
    if ~isdir(['uMM_', NameAppendix]), mkdir(['uMM_', NameAppendix]); end;
    if ~isdir(['MAmasks_', NameAppendix]), mkdir(['MAmasks_', NameAppendix]); end
    if ~isdir(['DirectMasks_', NameAppendix]), mkdir(['DirectMasks_', NameAppendix]); end
    if ~isdir(['TurnAngle_', NameAppendix]), mkdir(['TurnAngle_', NameAppendix]); end
    if ~isdir(['sMM_', NameAppendix]), mkdir(['sMM_', NameAppendix]); end
    if ~isdir(['Results']), mkdir(['Results']); end
else
    NameAppendix = []
    % Create folders for additional images
    if ~isdir(['uMM']), mkdir(['uMM']); end;
    if ~isdir(['MAmasks']), mkdir(['MAmasks']); end
    if ~isdir(['DirectMasks']), mkdir(['DirectMasks']); end
    if ~isdir(['TurnAngle']), mkdir(['TurnAngle']); end
    if ~isdir(['sMM']), mkdir(['sMM']); end
    if ~isdir(['Results']), mkdir(['Results']); end
end




% Depending on the sampling, create the list of binary masks for analysis. 
Settings.LastFrame = floor((Settings.NumFrames-1)/Settings.Sampling)*Settings.Sampling + 1
Temp.FrameList = [1:Settings.Sampling:Settings.LastFrame]'
FrameList = zeros(size(Temp.FrameList,1),3)
for i=1:1:size(FrameList,1)
    FrameList(i,1) = i;
    FrameList(i,2) = Temp.FrameList(i,1);
    if i == 1
        FrameList(i,3) == 0;
    else
        FrameList(i,3) = ((Temp.FrameList(i,1)-Temp.FrameList(i-1,1))*Settings.TimeIntervalSec)/60+FrameList(i-1,3); %Calculating time in minutes
    end
end

    ResultsMain = nan(size(FrameList, 1), 14); % Empty matrix for final data
    ResultsAuxiliary = nan(size(ResultsMain,1), 13);

% First loop is for extracting cell centroids and basic parameters of the
% shape.
for k=1:1:size(FrameList, 1)
    MaskCurrent=imread([File.NameBase, num2str(FrameList(k,2)), '.', Settings.FileType]); % Load the mask
    MaskCurrentBW = im2bw(MaskCurrent, 0.8); %In case of non-binary, but still Black-White image
  
    ShapeProps=regionprops(MaskCurrentBW ,'Centroid', 'MinorAxisLength','MajorAxisLength','Perimeter','Area', 'Orientation'); %Extract shape parameters
    
    
    t = FrameList(k,3); % Time in MINUTES
    X = ShapeProps.Centroid(1); % X-coordinate of shape in pixels
    Y = ShapeProps.Centroid(2); % Y-coordinate of shape in pixels
    Xum = ShapeProps.Centroid(1)*Settings.PixelSize; % X-coordinate of shape in micrometers
    Yum = ShapeProps.Centroid(2)*Settings.PixelSize; % Y-coordinate of shape in pixels
      
    % Calculate elongation using Major and Minor Axes of the shape   
    Major = ShapeProps.MajorAxisLength; 
    Minor = ShapeProps.MinorAxisLength; 
    Elongation = 1-(Minor/Major);
 
     
    % Add those results into the table   
    ResultsMain(k, 1) = FrameList(k,1); % Order number
    ResultsMain(k, 2) = FrameList(k,2); % Frame number
    ResultsMain(k, 3) = t; % czas
    ResultsMain(k, 4) = X;
    ResultsMain(k, 5) = Y;
    ResultsMain(k, 6) = Xum;
    ResultsMain(k, 7) = Yum;
    ResultsMain(k, 8) = NaN; % Assign NaN to displacements to calculate it later
    ResultsMain(k, 9) = NaN; % Assign NaN to velocities, will be calculated later
    ResultsMain(k, 10) = NaN; % Assign NaN to turning angles, will be calculated later
    ResultsMain(k, 11) = Elongation;
    ResultsMain(k, 12) = NaN; % Assign NaN to uMM, will be calculated later
    ResultsMain(k, 13) = NaN; % Assign NaN to sMM, will be calculated later
    ResultsMain(k,14) = NaN; % Assign NaN to M.A. dynamics, will be calculated later
    k
end


% Copy the first three columns of Main results to Auxiliary results
ResultsAuxiliary(:,1) = ResultsMain(:,1);
ResultsAuxiliary(:,2) = ResultsMain(:,2);
ResultsAuxiliary(:,3) = ResultsMain(:,3);


%% Now it's time to calculate displacements, velocities and turning angles, based on the centroids

[ResultsMain, ResultsAuxiliary] = VelocitiesAngles(File, Settings, ResultsMain, ResultsAuxiliary, NameAppendix, FrameList);
   


    %Uruchamiam program, który liczy mi k¹t miêdzy osi¹ d³ug¹ elipsy, a
    %kierunkiem ruchu komórki (wektor ³¹cz¹cy aktualny œrodek masy z
    %nastêpnym w kolejnoœci).
     
for k=1:1:size(ResultsMain, 1)
    

     
     [uMMangle, MAcoordinates, Fxy, OrientationDirection] = UnsignedMM(File, Settings, k, ResultsMain, ResultsAuxiliary, NameAppendix, FrameList);
     ResultsMain(k, 12) = uMMangle; 
     
     % Save Auxuliary results to calculate the sign of sMM angle:
     %Macierz DaneDodatkowe, s³u¿y nadaniu znaku k¹towi MiM
     ResultsAuxiliary(k, 5) = MAcoordinates(1,1);
     ResultsAuxiliary(k, 6) = MAcoordinates(1,2);
     ResultsAuxiliary(k, 7) = MAcoordinates(2,1);
     ResultsAuxiliary(k, 8) = MAcoordinates(2,2);
     ResultsAuxiliary(k, 9) = Fxy(1);
     ResultsAuxiliary(k, 10) = Fxy(2);
     ResultsAuxiliary(k, 13) = OrientationDirection;
     
     
     uMMangle
     k
end
   

    % Next program assigns proper sign to calculate sMM angle
    
for k=1:1:size(ResultsMain, 1)
   
    
     
       [sMMangle, MArotation, ResultsAuxiliary] = SignedMM(File, Settings, k, ResultsMain, ResultsAuxiliary, NameAppendix, FrameList);
          
    ResultsMain (k, 13) = sMMangle;

    ResultsAuxiliary(k,12) = MArotation;
    
    k
    
end

    % Loop for M.A. dynamics calculation:
    
    for k = 1:1:size(ResultsMain,1)-1 
       % Calculating the M.A. dynamics
    MAorientationCurr = ResultsAuxiliary(k, 12);
    MAorientationFoll = ResultsAuxiliary(k+1, 12);
   
    MAdynamics = -(MAorientationFoll-MAorientationCurr);
    ResultsMain (k, 14) = MAdynamics;
       
    end
    
    % Convert data to cell matrix for export to the text file.
    
    ResultsMainCell = num2cell(ResultsMain);
    
    ResultsMainCell{end, 8} = []; % Remove last nonexistent displacement
    ResultsMainCell{end, 9} = []; % Remove last nonexistent velocity
    ResultsMainCell{1, 10} = []; % Remove first nonexistent turning angle
    ResultsMainCell{end, 10} = []; % Remove last nonexistent turning angle
    ResultsMainCell{end, 12} = []; % Remove last nonexistent uMM angle  
    ResultsMainCell{end, 13} = []; % Remove last nonexistent sMM angle  
    ResultsMainCell{end, 14} = []; % Remove last nonexistent M.A. dynamics 

    
    %%% Create dataset with column headers
    %%% WARNING!!! It might be discontinued in newer Matlab versions 
    
    ResultsHeaders= dataset({ResultsMainCell(:,1), 'Order'}, ...
    {ResultsMainCell(:,2),'Frame'}, ...    
    {ResultsMainCell(:,3),'time [min]'}, ...
    {ResultsMainCell(:,4),'X[px]'}, ... 
    {ResultsMainCell(:,5), 'Y[px]'},...
    {ResultsMainCell(:,6), 'X[um]'}, ... 
    {ResultsMainCell(:,7), 'Y[um]'},...
    {ResultsMainCell(:,8), 'Displacement [px]'},...
    {ResultsMainCell(:,9), 'Velocity [um/min]'},...
    {ResultsMainCell(:,10), 'Turning angle [deg]'},...
    {ResultsMainCell(:,11), 'Elongation []'},...
    {ResultsMainCell(:,12), 'uMM angle [deg]'},...
    {ResultsMainCell(:,13), 'sMM angle [deg]'},...
    {ResultsMainCell(:,14), 'M.A. dynamics []'});
   
% Formatting auxiliary data
    ResultsAuxiliaryCell = num2cell(ResultsAuxiliary) 
    ResultsAuxiliaryCell{end, 4} = [];   %Remove last (nonexistent) orientation of displacement direction
    ResultsAuxiliaryCell{end, 9} = []; %Remove last (nonexistent) Fx coordinate
    ResultsAuxiliaryCell{end, 10} = []; %Remove last (nonexistent) Fy coordinate
    ResultsAuxiliaryCell{end, 11} = []; %Remove last (nonexistent) Temporary sMM angle
    ResultsAuxiliaryCell{end, 13} = []; %Remove last (nonexistent) orientation of direction vector (the one read from the image)

    
    ResultsAuxHeaders= dataset({ResultsAuxiliaryCell(:,1), 'Order'},...
    {ResultsAuxiliaryCell(:,2), 'Frame'},...
    {ResultsAuxiliaryCell(:,3), 'Time [min]'}, ...
    {ResultsAuxiliaryCell(:,4), 'Displacement Orientation (1)'}, ...
    {ResultsAuxiliaryCell(:,5), 'Ax'}, ...
    {ResultsAuxiliaryCell(:,6), 'Ay'}, ...
    {ResultsAuxiliaryCell(:,7), 'Bx'}, ...
    {ResultsAuxiliaryCell(:,8), 'By'}, ...
    {ResultsAuxiliaryCell(:,9), 'Fx'}, ...
    {ResultsAuxiliaryCell(:,10), 'Fy'},...
    {ResultsAuxiliaryCell(:,11), 'sMM Temp'},...
    {ResultsAuxiliaryCell(:,12), 'MAorientation'},...
    {ResultsAuxiliaryCell(:,13), 'Displacement Orietation (2)'});

    % Save data to .mat file
    
    if isempty(NameAppendix)
    save([File.pathname 'Results', '\',  File.NameBase(1:end-1), 'ResultsMain', '.mat'],'ResultsHeaders');
    else
    save([File.pathname 'Results', '\',  File.NameBase(1:end-1), 'ResultsMain_', NameAppendix,  '.mat'],'ResultsAuxHeaders');
    end 
    
    if isempty(NameAppendix)
    save([File.pathname 'Results', '\',  File.NameBase(1:end-1), 'ResultsAuxil', '.mat'],'ResultsAuxHeaders');
    else
    save([File.pathname 'Results', '\',  File.NameBase(1:end-1), 'ResultsAuxil_', NameAppendix, '.mat'],'ResultsAuxHeaders');
    end
    
    
    % Export to the .txt file with headers:
        % Main results
    if isempty(NameAppendix)
    export(ResultsHeaders,'File', [File.pathname 'Results', '\',  File.NameBase(1:end-1), 'ResultsMainHeaders', '.txt'],'WriteVarNames',true);
    else
    export(ResultsHeaders,'File', [File.pathname 'Results', '\',  File.NameBase(1:end-1), 'ResultsMainHeaders_', NameAppendix, '.txt'],'WriteVarNames',true);
    end

        % Auxiliary results
    if isempty(NameAppendix)
    export(ResultsAuxHeaders,'File', [File.pathname 'Results', '\', File.NameBase(1:end-1), 'ResultsAuxHeaders', '.txt'],'WriteVarNames',true);
    else
    export(ResultsAuxHeaders,'File', [File.pathname 'Results', '\', File.NameBase(1:end-1), 'ResultsAuxHeaders_', NameAppendix, '.txt'],'WriteVarNames',true);
    end
    
    
    % Export to the .txt file without headers:
        % Main results
    if isempty(NameAppendix)
    export(ResultsHeaders,'File', [File.pathname 'Results', '\',  File.NameBase(1:end-1), 'ResultsMain', '.txt'],'WriteVarNames',false);
    else
    export(ResultsHeaders,'File', [File.pathname 'Results', '\',  File.NameBase(1:end-1), 'ResultsMain_', NameAppendix, '.txt'],'WriteVarNames',false);
    end

        % Auxiliary results
    if isempty(NameAppendix)
    export(ResultsAuxHeaders,'File', [File.pathname 'Results', '\', File.NameBase(1:end-1), 'ResultsAux', '.txt'],'WriteVarNames',false);
    else
    export(ResultsAuxHeaders,'File', [File.pathname 'Results' , '\', File.NameBase(1:end-1), 'ResultsAux_', NameAppendix, '.txt'],'WriteVarNames',false);
    end
    
    % Short comment for the end:
    EndOfCalculationsFor = [File.NameBase]

end
