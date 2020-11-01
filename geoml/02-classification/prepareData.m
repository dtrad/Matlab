function [trainData,validateData,testData] = prepareData(nTrain,nValidate,nTest, saveData)
% prepareData Extract relevant data from Human Activity Recognition dataset
% This function downloads the data, reads the total acceleration data in the original Human
% Activity Recognition dataset, splits the original train vs test set into three sets,
% train, validate, and test, based on the counts (#subjects) given as parameters to this function,
% removes the 'Standing' category
% and optionally saves it in a new set of MATLAB data (.mat) files.
% 
% The data in the original dataset is available courtesy of:
% 
% Davide Anguita, Alessandro Ghio, Luca Oneto, Xavier Parra
% and Jorge L. Reyes-Ortiz. A Public Domain Dataset for Human Activity
% Recognition Using Smartphones. 21th European Symposium on Artificial
% Neural Networks, Computational Intelligence and Machine Learning,
% ESANN 2013. Bruges, Belgium 24-26 April 2013. 
% 
% Copyright 2014-2019 The MathWorks, Inc.


datalocation = fullfile(fileparts(which('DataPreparation.mlx')),'Data');
zipfilename = fullfile(datalocation,'UCI HAR Dataset.zip');

mkdir(datalocation);
websave(zipfilename, 'http://archive.ics.uci.edu/ml/machine-learning-databases/00240/UCI HAR Dataset.zip');

% Identify archived data file
[sourcePath, sourceFilename] = fileparts(zipfilename);

% Unzip source archive and create nested folder with same name
fprintf('Unzipping HAR dataset...')
unzip(zipfilename,sourcePath)
fprintf('Done.\n')

% Get hold of data
sourceDataPath = fullfile(sourcePath,sourceFilename);
sourceTrainDataPath = fullfile(sourceDataPath,'train');
sourceTestDatapath = fullfile(sourceDataPath,'test');

% Load acceleration data (Total)
ax_train = importAccelerationComponentFile(fullfile(sourceTrainDataPath,'Inertial Signals','total_acc_x_train.txt'));
ay_train = importAccelerationComponentFile(fullfile(sourceTrainDataPath,'Inertial Signals','total_acc_y_train.txt'));
az_train = importAccelerationComponentFile(fullfile(sourceTrainDataPath,'Inertial Signals','total_acc_z_train.txt'));
ax_test = importAccelerationComponentFile(fullfile(sourceTestDatapath,'Inertial Signals','total_acc_x_test.txt'));
ay_test = importAccelerationComponentFile(fullfile(sourceTestDatapath,'Inertial Signals','total_acc_y_test.txt'));
az_test = importAccelerationComponentFile(fullfile(sourceTestDatapath,'Inertial Signals','total_acc_z_test.txt'));

% Load labels - subject and activity id
s_train = importSubjectBufferList(fullfile(sourceTrainDataPath,'subject_train.txt'));
y_train = importClassBufferList(fullfile(sourceTrainDataPath,'y_train.txt'));
s_test = importSubjectBufferList(fullfile(sourceTestDatapath,'subject_test.txt'));
y_test = importClassBufferList(fullfile(sourceTestDatapath,'y_test.txt'));

% Concatenate acceleration
atx = [ax_train; ax_test];
aty = [ay_train; ay_test]; 
atz = [az_train; az_test]; 

% Concatenate labels - subject and activity id
subid = [s_train; s_test]; 
actid = [y_train; y_test]; 

% Create activity labels
actnames = {'Walking'  'WalkingUpstairs'  'WalkingDownstairs'  'Sitting'...
    'Standing'  'Laying'};

% Define sample frequency (in Hz) and time vector for each buffer
% Note: though the original data set is already buffered. So we aren't
% really using these parameters below.
fs = 50;
t = (1/fs)*(0:size(atx,2)-1); %#ok<NASGU>

% Convert the data into sorted table 
bufferedData = table(subid, atx,aty,atz, categorical(actid));
bufferedData = sortrows(bufferedData,'subid','ascend');
bufferedData.Properties.VariableNames = {'subject', 'X','Y','Z','activity'};
bufferedData.activity = renamecats(bufferedData.activity, actnames);
bufferedData(bufferedData.activity == 'Standing', :) = [];    % remove 'Standing' category
bufferedData.activity = removecats(bufferedData.activity);

% nTrain subjects for training and nValidate subjects for validation 
trainData = bufferedData(bufferedData.subject <= nTrain, :);
validateData = bufferedData(bufferedData.subject > nTrain & ...
    bufferedData.subject <= nTrain+nValidate, :);
testData = bufferedData(bufferedData.subject > nTrain+nValidate & ...
    bufferedData.subject <= nTrain+nValidate+nTest, :);

trainData = removevars(trainData, 'subject');
validateData = removevars(validateData, 'subject');
testData = removevars(testData, 'subject');

if saveData
  % Save relevant arrays to MAT file bufferedData.mat
  fprintf('Saving buffered data...')
  save('bufferedData.mat','trainData','validateData','testData','-v7.3');
  fprintf('Done.\n')
end

if true  % UNCOMMENT this section to re-create unbuffered data
  subjects = createUnbufferedData(subid, actid, atx, aty, atz); 
 
  T = struct2table(subjects);
  ta = cell2mat(T.totalacc);
  ai = cell2mat(T.actid);
 
  unbufferedData = array2table(ta,'VariableNames',{'X','Y','Z'});
  unbufferedData.activity = categorical(ai); 
  unbufferedData.activity = renamecats(unbufferedData.activity, actnames);
 
 unbufferedData.subject = zeros(size(ai));
 unbufferedData = movevars(unbufferedData,5,'Before',1);
 
 n = length(subjects);
 p = 0;
 
 for i = 1:n
     m = size(subjects(i).actid,1);
     unbufferedData{(p+1):p+m,'subject'} = i;
     p = p+m;     
 end
 
 unbufferedData(unbufferedData.activity == 'Standing', :) = [];
 unbufferedData.activity = removecats(unbufferedData.activity);
 
 % separate into train, validation, and test per subject counts given above
 unbufferedTrain = unbufferedData(unbufferedData.subject <= nTrain, :);
 unbufferedValidate = unbufferedData(unbufferedData.subject > nTrain & ...
     unbufferedData.subject <= nTrain+nValidate, :);
 unbufferedTest = unbufferedData(unbufferedData.subject > nTrain+nValidate & ...
     unbufferedData.subject <= nTrain+nValidate+nTest, :);
 
 if saveData
     fprintf('Saving unbuffered data...')
     save('unbufferedData.mat','unbufferedTrain', 'unbufferedValidate','unbufferedTest','-v7.3');
     fprintf('Done.\n')
 end
end


% --- Helper functions

function datastruct = createUnbufferedData(subid, bufactid, ...
    totaldatax, totaldatay, totaldataz)

% Initialize structure array - one component per subject
datastruct = struct('actid',[],'totalacc',[]);

s = unique(subid);
% Loop through subjects
for ks = 1:length(s)

    subject = s(ks);
    
     
    subidx = subid==subject;
    totx = totaldatax(subidx,:);
    toty = totaldatay(subidx,:);
    totz = totaldataz(subidx,:);
    aid = bufactid(subidx,:);
    
    [acc3d, unbufactid] = unbufferSubjectData(totx, toty, totz, aid);
    
    datastruct(subject).actid = unbufactid;
    % (Rescale all acelerations to real-world values before saving)
    g = 9.80665;
    datastruct(subject).totalacc = g * acc3d;

        
end

end

function [acc3d, actid] = unbufferSubjectData(axbuf, aybuf, azbuf, bufactid)

if(isempty(axbuf))
    fprintf('Warning: not enough data for subject %d\n',subjectNumber)
    return
end

L = size(axbuf,2);

% Get rid of overlapping regions and make a single long vector

% Data X
totx = axbuf.';
t1 = totx(1:L/2,1);
totx(1:L/2,:) = [];
tx = [t1;totx(:)];

% Data Y
toty = aybuf.';
t1 = toty(1:L/2,1);
toty(1:L/2,:) = [];
ty = [t1;toty(:)];

% Data Z
totz = azbuf.';
t1 = totz(1:L/2,1);
totz(1:L/2,:) = [];
tz = [t1;totz(:)];

% Create new vector of activities id with same length
tmp = [bufactid(1);bufactid];
aextended = tmp * ones(1,L/2);
actid = reshape(aextended.',[],1);

% Concatenate all acceleration components together
acc3d = [tx, ty, tz];

end
end