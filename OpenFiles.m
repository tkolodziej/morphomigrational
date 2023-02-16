% Version for 13.02.2023, Tomasz Ko³odziej
% Script opens the binary (black-white) masks and recognizes their names.
% Structure of this script (for opening the files) comes from TFM program 
% by Xavier Trepat (IBEC Barcelona).
function [File] = OpenFiles(File, Settings);

% Open first mask file.
[Temp.Mask, File.pathname] = uigetfile(['*.' Settings.FileType], 'Open first binary mask');

% Find the base of files names
File.NameBase = Temp.Mask(1:end-5); %Can be modified, depending on the length of file extension.

%How many frames in Masks sequence?
dirTrac = dir;
imNum = 0;
for i = 3:length(dirTrac)
    if ~isempty(findstr(dirTrac(i).name, File.NameBase))
        imNum = imNum + 1;
    end
end
File.FileNumber = imNum

% Defining name of Masks files:
for i=1:File.FileNumber
    prefix = [];
    File.Nazwa(i).Mask = strcat(File.pathname,File.NameBase, num2str(prefix), num2str(i),'.', Settings.FileType);
end



save File File
save Settings Settings
