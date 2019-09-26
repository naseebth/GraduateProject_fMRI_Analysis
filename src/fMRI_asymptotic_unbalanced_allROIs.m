clear all;

load data-starplus-04847-v7.mat;

%Function to get the unbalanced label dataset(i.e.only two labels 1 or 2 aretaken
[info,data,meta] = transformIDM_selectROIVoxels(info,data,meta,{'CALC','LFEF','LIPL','LIT','LPPREC','LSPL','LTRIA','RFEF','RIPS','ROPER','RSGA','RT','SMA','LDLPFC','LIPS','LOPER','LSGA','LT','RDLPFC','RIPL','RIT','RPPREC','RSPL','RTRIA'});

%Function to get the unbalanced label dataset(i.e.only two labels 1 or 2 aretaken
[example1,labels,expInfo] = idmToExamples_fixation(info,data,meta,'full');

%Extracting number of rows and columns from the dataset
nrows2 = size(example1,1);
ncols = size(example1,2);

%Compute the Markov-chain state changes in the fMRI signals(such that the spike in value
%or signal should be represented by 1 from one signal to another, and,
%decrease in value is represented by 0)
scmatrix = zeros(nrows2, ncols-1); %initialize the state change matrix. 
%(Note:state changes can only be computed 92609 times for 92610 features)
for row = 1:nrows2 %loop to include all 50 trials
     mat_1=example1(row,1:ncols-1); 
     mat_2=example1(row,2: ncols);
     difference_mat=mat_2-mat_1; %calculate difference betweeen the previous and next feature to compute state changes.
     mat_sign=sign(difference_mat); %convert positive values to 1 and negative values to -1
     scmatrix(row,:)=mat_sign>0; %assings all the 1's values to denote positive state changes in fMRI signals.
end


% Compute 2X2 frequency matrix for all 50 trials to count the state changes
% from 0 to 0, 0 to 1, 1 to 0, and, 1 to 1
Fr_matrix = zeros(2,2); %initiliaze 2 by 2 matrix to store the frequency of state changes
for row = 1:50 %loop to include all 50 trials
    x_trials = scmatrix(row,:); %get the trials from the scmatrix one at a time.
    x_trials1 = x_trials(:,1:ncols-2);
    x_trials2 = x_trials(:,2:ncols-1);
    add_xtrials = x_trials1 + x_trials2; %To figure out the state change from 0 to 0 or 1 to 1. i.e. the value remained the same.
    difference_xtrials = x_trials2-x_trials1; %To figure out the state change from 0 to 1 or 1 to 0. i.e. the value changes from one state to another.
    %It counts state changes 0 to 0, 0 to 1, 1 to 0 and 1 to 1
    sc00 = sum(add_xtrials<1);
    sc01 = sum(difference_xtrials>0);
    sc11 = sum(add_xtrials>1);
    sc10 = sum(difference_xtrials<0);
    Fr_matrix(1,1) = sc00; % gives frequency where state change is same from 0 to 0
    Fr_matrix(1,2) = sc01; % gives frequency where state change is different from 0 to 1
    Fr_matrix(2,1) = sc10; % gives frequency where state change is different from 1 to 0
    Fr_matrix(2,2) = sc11; % % gives frequency where state change is same from 1 to 1
    Fr_matrix_final{row} = Fr_matrix; 
end

% Create 2x2 Transition probabilities matrix for all 50 trials
global pr_matrix;
pr_matrix = zeros(2,2);%initialize 2 by 2 matrix to represent transition probabilities
for col = 1:length( Fr_matrix_final)
    global pr_matrix_final;
    fr_mat = Fr_matrix_final{col}; %Get the counts of state changes for probability calculation
    pr_matrix(1,1) = fr_mat(1,1)/(fr_mat(1,1)+fr_mat(1,2)); %pr that if the current state is 0, then the next state is 0 
    pr_matrix(1,2) = 1-pr_matrix(1,1); %pr that if the current state is 0, then next state is 1 
    pr_matrix(2,1) = fr_mat(2,1)/(fr_mat(2,1)+fr_mat(2,2)); %pr that if the current state is 1, then next state is 0
    pr_matrix(2,2) = 1-pr_matrix(2,1); %pr that if the current state is 1, then next state is 1
    pr_matrix_final{col} = pr_matrix;
end  


% Computing two step transition probability matrix for all 50 trials;
% pr_matrix_2 = pr_matrix^2;
pr_matrix_2 = zeros(2,2);
for col = 1:length(pr_matrix_final)
    pr_mat = pr_matrix_final{col}; %Get the counts of state changes for probability calculation
    pr_matrix_2(1,1) = pr_mat(1,1)^2 +(pr_mat(1,2)*pr_mat(2,1)); 
    pr_matrix_2(1,2) = pr_mat(1,2)*(pr_mat(1,1)+pr_mat(2,2)); 
    pr_matrix_2(2,1) = pr_mat(2,1)*(pr_mat(1,1)+pr_mat(2,2));
    pr_matrix_2(2,2) = pr_mat(1,1)^2 +(pr_mat(1,2)*pr_mat(2,1)); 
    pr_matrix_2_final{col} = pr_matrix_2;
end  




%Computing nstep transition probability matrix for all 50 trails for
%asymptotic analysis
%A function naseeb_TransitionProbability is created so that the value of n
%is taken as input arguments to compute the n-step transition probability.
[pr_matrix_n_final] = naseeb_TransitionProbability(100);

% Create a new matrix to represent all the n-step transition probability
% matrix.
asymp_matrix = zeros(200, ncols-1); %initializing asymp matrix.

col_1 = scmatrix(:,1); %extract the first column of state change matrix to initialize the value of the new matrix.

%Now, generate sequence clones for each of 50 observations in the
%pr_matrix_n_final observation.
for obs = 1:50
    for row = 1:4  %since, 50 should be looped 4 times for 200 values required.
        for column = 1:92608
            new_row = row + (4*(obs-1)); % create new row of all the probability values such that there are total of 200 rows
            asymp_matrix(new_row,1) = col_1(obs,1);
            new_mat= pr_matrix_n_final{obs}; %assign new_mat matrix the value of pr_matrix_n_final for all probabilities
       if (asymp_matrix(new_row,column) == 0)
                r=rand; %generate random numbers whose elements are uniformly distributed in the interval(0,1)     
                if (r <= new_mat(1,2)) % access values of probabilities at place (1,2), and if the value is equal or greater to random number between 0 and 1, then assign the values to 1 or else, assign to 0
                    asymp_matrix(new_row,column+1)=1; 
                else
                    asymp_matrix(new_row,column+1)=0;
                end
            elseif (asymp_matrix(new_row,column) == 1)
                r=rand;
                if (r <= new_mat(2,1)) % access values of probabilities at place (1,2), and if the value is equal or greater to random number between 0 and 1, then assign the values to 1 or else, assign to 0
                    asymp_matrix(new_row,column+1)=0; 
                else
                    asymp_matrix(new_row,column+1)=1;
                end
            end
        end
    end
end
% asymp_matrix is our final compressed sensing matrix for asympotic analysis, as 0 in the matrix will eliminate all the unwanted sequence from the dataset

nrows2 = 200; %Since, we have 50 observations for probabilities of each trials. We need total of 200 observations for the probability sequence.

%Create labels from the original dataset such that it matches the new
%matrix. In this case, we create new label set for all the 200 rows.
new_labels = zeros(nrows2,1);
new_labels(1:(nrows2/2),1) = 1; %assign label 1 for half of the observations
new_labels((nrows2/2)+1:nrows2) = 2; %assign label 2 for remaining observations.

% Now, the original dataset has only 50 rows, but we need 200 rows since
% our probabilities matrix has 200 rows. So, we will create a new matrix to
% represent the original matrix by repeating the rows.
new_example1 = zeros(nrows2,ncols);
for obs = 1:50
    for row = 1:4
        new_row = row + (4*(obs-1));
        new_example1(new_row,:) = example1(obs,:);
    end
end

% Now, let us perform matrix multiplication by using our compressed sensing
% matrix "asymp_matrix" with the new data set new_example1.
new_example2 = new_example1(:,2:ncols);
final_example1 = new_example2 .* asymp_matrix;  %this performs element by element multiplication. We need to perform machine learning on this dataset.


%Assigning new labels to the final dataset
ncolumns2 = size(final_example1,2);
NN = 1; 
MM = ncolumns2;
examples = final_example1(:,NN:MM);
% Add the labels to the data matrix
examples(:,92610) = new_labels;
new_example2(:,92610) = new_labels;

% Testing/training percentages, 80-20
dataA = examples;  
p = .80 ;   % proportion of rows to select for training
N = size(examples,1);  % total number of rows 
tf = false(N,1); % create logical index vector
tf(1:round(p*N)) = true;  
tf = tf(randperm(N));  % randomise order
dataTraining = dataA(tf,:); % Training set using degraded data
%dataTesting = dataA(~tf,:);
dataTesting = new_example2(~tf,:); % Testing using original data


%Create traing and testing dataset
trainExamples = dataTraining(:,1:MM);
trainLabels   = dataTraining(:,92610);
testExamples  = dataTesting(:,1:MM);
testLabels    = dataTesting(:,92610); 


% train a classifier
62/5[classifier] = trainClassifier(trainExamples,trainLabels,'nbayes');
%[classifier] = trainClassifier(trainExamples,trainLabels,'logisticRegression');
% [classifier] = trainClassifier(trainExamples,trainLabels,'SMLR');
% [classifier] = trainClassifier(trainExamples,trainLabels,'neural');

% apply a classifier
[predictions] = applyClassifier(testExamples,classifier);

% summarizePredictions
[result,predictedLabels,trace] = summarizePredictions(predictions,classifier,'accuracy',testLabels);
result{1}