

% load My4Emotions;
load My4EmotionsMultiTaper;

numTest = 7*40;
numTrain = 0.9*numTest;

% numSamplesperEmotion = 40;
% numFeats = 12*2;
% numEmotions = 7;


Model = Model(1:numTest,:);
speciesEmo = speciesEmo(1:numTest,1);

% Model = zscore(Model);
% for k=1:24
%     tmp = Model(:,k);
%     Model(:,k) = zscore(tmp);
% end

%Shuffle the data
rand('twister',0);
perm = randperm(numTest);
Model = Model(perm,:);
speciesEmo = speciesEmo(perm,:);

featureNames = {'Cent_mu','Cent_med','Spread_mu','Spread_med','Skew_mu','Skew_med',...
    'Kurt_mu','Kurt_med','Slope_mu','Slope_med','Decr_mu','Decr_med',...
    'Rolloff_mu','Rolloff_med','Flux_mu','Flux_med','Erg_mu','Erg_med',...
    'Flat_mu','Flat_med','Crest_mu','Crest_med','Entropy_mu','Entropy_med','class'};

%Prepare test and training sets. 
data = [num2cell(Model),speciesEmo];

train = data(1:numTrain,:);
test  = data(numTrain+1:end,:);

classindex = 25;

%Convert to weka format
train = matlab2weka('4Emo-train',featureNames,train);
test =  matlab2weka('4Emo-test',featureNames,test);

%Train the classifier
nb = trainWekaClassifier(train,'functions.SMO');

%Test the classifier
predicted = wekaClassify(test,nb);

%The actual class labels (i.e. indices thereof)
actual = test.attributeToDoubleArray(classindex-1); %java indexes from 0

errorRate = sum(actual ~= predicted)/length(actual);

Accuracy = (1 - errorRate)*100

