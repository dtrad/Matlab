function featureT = extractWaveletFeatures(rawData, sf, N)
% Apply Wavelet scattering to (raw) data provided (3 dimensional signal)
% and convert into a feature table, including class label as last column
%
%   <sf> is the scattering function, <N> the signal length.
%   The raw data is expected to have 3 (signal) columns, plus one
%   indicating the subject and another the label ("activity").
%   With "observation" below we refer to data from a specific activity for a specific subject

% Extract X, Y, Z from raw data as a matrix (columns 2-4)
signalData = table2array(rawData(:,2:4));
[gTrain, ~,activityTrain] = findgroups(rawData.subject, rawData.activity);
   % figure out which rows in the (unbuffered) data belong to the same subject <gTrain>, and
   % remember the corresponding activity labels <activityTrain>

% Apply wavelet scattering to train data 
waveletMatrix = splitapply(@(x) {featureMatrix(sf,(x(1:N,:)))},signalData,gTrain);
   % apply wavelet scattering featureMatrix on all rows from each subject
   % giving us a matrices of 15 features X 9 time intervals X 3 signals

% build up feature table (need to line up the features from each of the 3
% signals into a row, and will get 9 rows for each observation
featureT = table;

% loop over each of the #subjects X #activities (125=25 x 3 for the train set) 
% three dimensional matrices that were created above
for i = 1 : size(waveletMatrix,1)
    
    oneO = waveletMatrix{i};   % process this observation
    thisObservation = [oneO(:,:,1); oneO(:,:,2); oneO(:,:,3)];
    thisObservation = array2table(thisObservation');  % don't forget to convert the features from rows into columns
    
    featureT = [featureT; thisObservation];
end

% get labels by duplicating label for each row of "wavelet" features obtained for  subject x activity
featureT.activity = repelem(activityTrain,size(waveletMatrix{1},2));
   
end

