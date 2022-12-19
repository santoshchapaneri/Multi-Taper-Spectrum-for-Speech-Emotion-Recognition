% Toy problem - only 2 emotions - binary SVM classifier

% Codes for Emotion
% 1: Anger; 2: Boredom; 

% Using Timbre Toolbox

clc;

numSamplesperEmotion = 40;
numFeats = 12*2;
numEmotions = 7;
k = 1;
Model = zeros(numSamplesperEmotion*numEmotions, numFeats);

%% === Parameters for Timbre Toolbox
do_s.b_TEE				= 0;    % descriptors from the Temporal Energy Envelope
do_s.b_STFTmag			= 1;    % descriptors from the STFT magnitude
do_s.b_STFTMultiTapermag = 0;
do_s.b_STFTpow			= 0;    % descriptors from the STFT power
do_s.b_Harmonic			= 0;	% descriptors from Harmonic Sinusoidal Modeling representation
do_s.b_ERBfft			= 0;    % descriptors from ERB representation (ERB being computed using FFT)
do_s.b_ERBgam			= 0;    % descriptors from ERB representation (ERB being computed using Gamma Tone Filter)
config_s.xcorr_nb_coeff = 12;	% === defines the number of auto-correlation coefficients that will be sued
config_s.threshold_harmo= 0.3;	% === defines the threshold [0,1] below which harmonic-features are not computed
config_s.nb_harmo		= 20;	% === defines the number of harmonics that will be extracted
% ====================


%% 1: Anger
Feat = zeros(numFeats,1);
for numFile=1:numSamplesperEmotion
    reffile = sprintf('DB2/Anger/1_%d.wav', numFile);
    ALLDESC_s{numFile}	= Gget_desc_onefile(reffile, do_s, config_s);
	ALLTM_s{numFile}	= Gget_temporalmodeling_onefile(ALLDESC_s{numFile});
    Feat = [];
    fieldname_c = fieldnames(ALLTM_s{numFile});
    for f=1:length(fieldname_c)
        value	= ALLTM_s{numFile}.(fieldname_c{f});
        Feat = [Feat value];
    end
    if(~isreal(Feat))
        Feat = abs(Feat);
    end
    Model(k,:) = Feat;
    speciesEmo{k,1} = 'Anger';
    k = k + 1;
end

% 2: Boredom
for numFile=1:numSamplesperEmotion
    reffile = sprintf('DB2/Boredom/2_%d.wav', numFile);
    ALLDESC_s{numFile}	= Gget_desc_onefile(reffile, do_s, config_s);
	ALLTM_s{numFile}	= Gget_temporalmodeling_onefile(ALLDESC_s{numFile});
    Feat = [];
    fieldname_c = fieldnames(ALLTM_s{numFile});
    for f=1:length(fieldname_c)
        value	= ALLTM_s{numFile}.(fieldname_c{f});
        Feat = [Feat value];
    end
    if(~isreal(Feat))
        Feat = abs(Feat);
    end
    Model(k,:) = Feat;
    speciesEmo{k,1} = 'Boredom';
    k = k + 1;
end

% 3: Fear
for numFile=1:numSamplesperEmotion
    reffile = sprintf('DB2/Fear/3_%d.wav', numFile);
    ALLDESC_s{numFile}	= Gget_desc_onefile(reffile, do_s, config_s);
	ALLTM_s{numFile}	= Gget_temporalmodeling_onefile(ALLDESC_s{numFile});
    Feat = [];
    fieldname_c = fieldnames(ALLTM_s{numFile});
    for f=1:length(fieldname_c)
        value	= ALLTM_s{numFile}.(fieldname_c{f});
        Feat = [Feat value];
    end
    if(~isreal(Feat))
        Feat = abs(Feat);
    end
    Model(k,:) = Feat;
    speciesEmo{k,1} = 'Fear';
    k = k + 1;
end

% 4: Happiness
for numFile=1:numSamplesperEmotion
    reffile = sprintf('DB2/Happiness/4_%d.wav', numFile);
    ALLDESC_s{numFile}	= Gget_desc_onefile(reffile, do_s, config_s);
	ALLTM_s{numFile}	= Gget_temporalmodeling_onefile(ALLDESC_s{numFile});
    Feat = [];
    fieldname_c = fieldnames(ALLTM_s{numFile});
    for f=1:length(fieldname_c)
        value	= ALLTM_s{numFile}.(fieldname_c{f});
        Feat = [Feat value];
    end
    if(~isreal(Feat))
        Feat = abs(Feat);
    end
    Model(k,:) = Feat;
    speciesEmo{k,1} = 'Happy';
    k = k + 1;
end

% 5: Sadness
for numFile=1:numSamplesperEmotion
    reffile = sprintf('DB2/Sadness/5_%d.wav', numFile);
    ALLDESC_s{numFile}	= Gget_desc_onefile(reffile, do_s, config_s);
	ALLTM_s{numFile}	= Gget_temporalmodeling_onefile(ALLDESC_s{numFile});
    Feat = [];
    fieldname_c = fieldnames(ALLTM_s{numFile});
    for f=1:length(fieldname_c)
        value	= ALLTM_s{numFile}.(fieldname_c{f});
        Feat = [Feat value];
    end
    if(~isreal(Feat))
        Feat = abs(Feat);
    end
    Model(k,:) = Feat;
    speciesEmo{k,1} = 'Sad';
    k = k + 1;
end

% 6: Neutral
for numFile=1:numSamplesperEmotion
    reffile = sprintf('DB2/Neutral/6_%d.wav', numFile);
    ALLDESC_s{numFile}	= Gget_desc_onefile(reffile, do_s, config_s);
	ALLTM_s{numFile}	= Gget_temporalmodeling_onefile(ALLDESC_s{numFile});
    Feat = [];
    fieldname_c = fieldnames(ALLTM_s{numFile});
    for f=1:length(fieldname_c)
        value	= ALLTM_s{numFile}.(fieldname_c{f});
        Feat = [Feat value];
    end
    if(~isreal(Feat))
        Feat = abs(Feat);
    end
    Model(k,:) = Feat;
    speciesEmo{k,1} = 'Neutral';
    k = k + 1;
end

% 7: Disgust
for numFile=1:numSamplesperEmotion
    reffile = sprintf('DB2/Disgust/7_%d.wav', numFile);
    ALLDESC_s{numFile}	= Gget_desc_onefile(reffile, do_s, config_s);
	ALLTM_s{numFile}	= Gget_temporalmodeling_onefile(ALLDESC_s{numFile});
    Feat = [];
    fieldname_c = fieldnames(ALLTM_s{numFile});
    for f=1:length(fieldname_c)
        value	= ALLTM_s{numFile}.(fieldname_c{f});
        Feat = [Feat value];
    end
    if(~isreal(Feat))
        Feat = abs(Feat);
    end
    Model(k,:) = Feat;
    speciesEmo{k,1} = 'Disgust';
    k = k + 1;
end

save('My4Emotions.mat','Model','speciesEmo');

% 
% %% Prepare Training set
% numTrainSamples = 3;
% numCoeffsTrain = 10;
% 
% % ModelTrain = zeros(numTrainSamples*numEmotions, numCoeffs + (numCoeffs^2));
% ModelTrain = zeros(numTrainSamples*numEmotions, numCoeffsTrain);
% LabelTrain = zeros(numTrainSamples*numEmotions,1);
% j = 1; ind = 1;
% for k = 1:numTrainSamples*numEmotions
%     ModelTrain(k,:) = Model(j,1:numCoeffsTrain);
%     LabelTrain(k,1) = ind;
%     if (mod(k,numTrainSamples)==0)
%         j = j + (numSamplesperEmotion-numTrainSamples+1);
%         ind = ind + 1;
%     else
%         j = j + 1;
%     end
% end
% 
% %% Prepare Testing set
% numTestSamples = 2;
% numCoeffsTest = 10;
% 
% ModelTest = zeros(numTestSamples*numEmotions, numCoeffsTest);
% LabelTest = zeros(numTestSamples*numEmotions,1);
% j = numTrainSamples+1; ind = 1;
% for k = 1:numTestSamples*numEmotions
%     ModelTest(k,:) = Model(j,1:numCoeffsTest);
%     LabelTest(k,1) = ind;
%     if (mod(k,numTestSamples)==0)
%         j = j + (numSamplesperEmotion-numTestSamples+1);
%         ind = ind + 1;
%     else
%         j = j + 1;
%     end
% end
% 
% 
% %% Train SVM on training set
% nbclass=2;
% %   Learning and Learning Parameters
% c = 1000;
% lambda = 1e-7;
% kerneloption = 2;
% kernel = 'gaussian';
% verbose = 1;
% 
% %---------------------One Against One algorithms----------------
% [xsup,w,b,nbsv,classifier,pos]=svmmulticlassoneagainstone(ModelTrain,LabelTrain,nbclass,c,lambda,kernel,kerneloption,verbose);
% 
% %% Test the SVM multi-class classifier
% [ypred,maxi] = svmmultivaloneagainstone(ModelTest,xsup,w,b,nbsv,kernel,kerneloption);
% ypred;