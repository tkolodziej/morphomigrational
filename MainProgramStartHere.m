% Main program for sMM angle calculation.
% Version for 07.02.2023
% by Tomasz Ko³odziej, Jagielonian University Collegium Medicum, Kraków,
% Poland

% To run the program please go to directory with binary masks of cells. 

% Each file has to have the same file base and suffix number recognizing
% its order. Suffixes has to be single numbers, without additional zero
% fill (File_t1, File_t2... - OK; File_t01, File_t02 - won't work).

% In each frame it has to be only one cell, because this program does not
% follow or calculate parameters for multiple cells!!!

clear all;
close all;
Settings.NumFrames=181; % Number of frames in the whole image sequence (regardless the chosen sampling).
Settings.PixelSize=0.6331; % Pixel size in micrometers
Settings.Sampling = 1; % How densely the sequence should be sampled (1 - all frames, 2 - each second frame, 3 - each third frame, etc.)
Settings.TimeIntervalSec=10;% Time interval of the original time-lapse experiment IN SECONDS (regardless the chosen sampling)
Settings.FileType='tif' %Extension of files with black-white masks (without the dot before an extension)
File.pathname=1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

File = OpenFiles(File, Settings); %Program otwiera pliki .mat z maskami.

File = CalculationsMajorScript(File, Settings); %Program liczy to co trzeba :D

