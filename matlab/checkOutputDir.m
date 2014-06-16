function [success,message,messageID] = checkOutputDir(outputDirectory)
%
% This script checks to see if the outputDirectory already exists. If
% it exists, the user is prompted to see whether or not they want to
% overwrite the directory.

directoryCreate = 0;

if exist(outputDirectory,'dir')
    userInput = input('Directory already exists - overwrite? [Y/n]:','s');
    if strcmp(userInput,'Y');
        directoryCreate = 1;
    elseif strcmp(userInput,'n');
        disp('Not overwriting directory.');
    else
        userInput = input('Try again - enter [Y/n]');
        if strcmp(userInput,'Y');
            directoryCreate = 1;
        elseif strcmp(userInput,'n');
            disp('Not overwriting directory.');
        else
            disp('Learn how to enter yes or no!')
            return
        end
    end
else
    directoryCreate = 1;
end

if directoryCreate
    [success, message]          = rmdir(outputDirectory,'s'); % recursively remove old directory
    [success,message,messageID] = mkdir(outputDirectory); % make new directory
else
    success   = 0;
    message   = 0;
    messageID = 0;
end

return

