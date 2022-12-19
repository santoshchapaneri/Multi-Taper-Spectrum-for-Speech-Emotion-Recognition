
% load My4Emotions;
load My4EmotionsMultiTaper;

Model = Model(1:280,:);
speciesEmo = speciesEmo(1:280,1);

% % Anger v/s Happy
% angerset = find(strcmp(speciesEmo,'Fear'));
% happyset = find(strcmp(speciesEmo,'Happy'));
% 
% Model2 = [Model(angerset,:);Model(happyset,:)];
% speciesEmo2 = [speciesEmo(angerset,:);speciesEmo(happyset,:)];
% numTest = 2*40;

angerset = find(strcmp(speciesEmo,'Anger'));
boredomset = find(strcmp(speciesEmo,'Boredom'));
fearset = find(strcmp(speciesEmo,'Fear'));
happyset = find(strcmp(speciesEmo,'Happy'));
sadset = find(strcmp(speciesEmo,'Sad'));
neutralset = find(strcmp(speciesEmo,'Neutral'));
disgustset = find(strcmp(speciesEmo,'Disgust'));

% % Positive v/s Negative Activation
% Model2 = [Model(angerset,:);Model(fearset,:);Model(happyset,:);Model(disgustset,:);Model(sadset,:);Model(boredomset,:)];
% for k=1:4*40
% 	speciesEmo2{k,1} = 'PosActive';
% end
% for k=4*40+1:6*40
% 	speciesEmo2{k,1} = 'NegActive';
% end
% numTest = 6*40;

% % Positive v/s Negative Valence
% Model2 = [Model(happyset,:);Model(angerset,:);Model(boredomset,:);Model(fearset,:);Model(sadset,:);Model(disgustset,:)];
% for k=1:1*40
% 	speciesEmo2{k,1} = 'PosValence';
% end
% for k=1*40+1:6*40
% 	speciesEmo2{k,1} = 'NegValence';
% end
% numTest = 6*40;

% Neutral v/s Emotion
Model2 = [Model(neutralset,:);Model(angerset,:);Model(boredomset,:);Model(fearset,:);Model(happyset,:);Model(sadset,:);Model(disgustset,:)];
for k=1:1*40
	speciesEmo2{k,1} = 'Neutral';
end
for k=1*40+1:7*40
	speciesEmo2{k,1} = 'Emotion';
end
numTest = 7*40;

%Shuffle the data
rand('twister',0);
perm = randperm(numTest);
Model2 = Model2(perm,:);
speciesEmo2 = speciesEmo2(perm,:);

featureNames = {'Cent_mu','Cent_med','Spread_mu','Spread_med','Skew_mu','Skew_med',...
    'Kurt_mu','Kurt_med','Slope_mu','Slope_med','Decr_mu','Decr_med',...
    'Rolloff_mu','Rolloff_med','Flux_mu','Flux_med','Erg_mu','Erg_med',...
    'Flat_mu','Flat_med','Crest_mu','Crest_med','Entropy_mu','Entropy_med','class'};

%Prepare test and training sets. 
data = [num2cell(Model2),speciesEmo2];
numTrain = 0.05*numTest;

train = data(1:numTrain,:);
test  = data(numTrain+1:end,:);

classindex = 25;

%Convert to weka format
train = matlab2weka('4Emo-train',featureNames,train);
test =  matlab2weka('4Emo-test',featureNames,test);

%Train the classifier
% options = '-K "weka.classifiers.functions.supportVector.PolyKernel -C 250007 -E1.0"';
nb = trainWekaClassifier(train,'functions.SMO');

%Test the classifier
predicted = wekaClassify(test,nb);

%The actual class labels (i.e. indices thereof)
actual = test.attributeToDoubleArray(classindex-1); %java indexes from 0

errorRate = sum(actual ~= predicted)/length(actual);

Accuracy = (1 - errorRate)*100

