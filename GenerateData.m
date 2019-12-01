clear all

% Paramter definition

TrainingSize = 100000;  % number of points used for training

Elements = 4; % number of antenna elements in linear array
range = 0.9;   % range of AoAs (between 0 and 1; 1 == full -pi/2 to pi/2)
             % field of view = 180*range
             
NumSignals = 4;

SNR = inf;   % SNR in dB

DataType = 'CovMat';    % either directly received samples 'RxSamp'
                        % or covariance matrix 'CovMat'

CovSampSize = 100;	   % Number of rx samples over which to 
				   % estimate the covariance matrix
                        

%--------------------------------------------------------------
%% No changes should be needed beyond this line		%%%%%%%%%
%% Of course you are welcome to do so if you wish	%%%%%%%%%
%--------------------------------------------------------------

d_lambda = 0.5*(0:Elements-1);


if DataType == 'RxSamp'
    Target = (rand(TrainingSize,1)*pi*range - pi/2*range);
    Features = exp(j*2*pi*sin(Target)*d_lambda) + sqrt(1/2/10^(SNR/10))*randn(1,4) + j*sqrt(1/2/10^(SNR/10))*randn(1,4);
    
    for ii=2:NumSignals
        tmp = (rand(TrainingSize,1)*pi*range - pi/2*range);
        Target = [Target,tmp];
        Features = Features + exp(j*2*pi*sin(tmp)*d_lambda);
    end
    FeatureVectors = [real(Features),imag(Features)];
    
    FileName1 = strcat('FeaturesTrainingSNR',num2str(SNR),'Range',num2str(range),'NumSig',num2str(NumSignals),'.csv');
    FileName2 = strcat('TargetTrainingSNR',num2str(SNR),'Range',num2str(range),'NumSig',num2str(NumSignals),'.csv');
    csvwrite(FileName1,FeatureVectors);
    csvwrite(FileName2,Target);


elseif DataType == 'CovMat'
    
    AoA = (rand(TrainingSize,1)*pi*range - pi/2*range);
    Target = AoA;

    for i=1:length(Target)
        tmp = repmat([exp(j*2*pi*sin(Target(i,1))*d_lambda)].',1,CovSampSize).*repmat(sign(randn(1,CovSampSize)),length(d_lambda),1) + sqrt(1/2/10^(SNR/10))*randn(4,CovSampSize) + j*sqrt(1/2/10^(SNR/10))*randn(4,CovSampSize);
            for ii=2:NumSignals
                Target(i,ii) = (rand*pi*range - pi/2*range);
                %Target(i,ii) = AllowedVals(ceil(rand*5));
                tmp = tmp + repmat([exp(j*2*pi*sin(Target(i,ii))*d_lambda)].',1,CovSampSize);
            end
        CovMat = 1/CovSampSize/length(d_lambda)*tmp*tmp';
        FeatureVectors(i,:) = [real(CovMat(1,2)),imag(CovMat(1,2)),real(CovMat(1,3)),imag(CovMat(1,3)),real(CovMat(1,4)),imag(CovMat(1,4)),...
        real(CovMat(2,3)),imag(CovMat(2,3)),real(CovMat(2,4)),imag(CovMat(2,4)),real(CovMat(3,4)),imag(CovMat(3,4))];
    end
    
    Target = sort(Target,2);
    
    FileName1 = strcat('CovFeaturesTrainingSNR',num2str(SNR),'Range',num2str(range),'NumSig',num2str(NumSignals),'.csv');
    FileName2 = strcat('CovTargetTrainingSNR',num2str(SNR),'Range',num2str(range),'NumSig',num2str(NumSignals),'.csv');
    csvwrite(FileName1,FeatureVectors);
    csvwrite(FileName2,Target);

 
else 
    error('DataType misdefined.')
end

