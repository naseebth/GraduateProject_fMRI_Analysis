clear all;

load data-starplus-04847-v7.mat;

% only returns the IDM for the 7 ROIs{'CALC' 'LDLPFC' 'LIPL' 'LIPS' 'LOPER' 'LT' 'LTRIA'}
[info,data,meta] = transformIDM_selectROIVoxels(info,data,meta,{'CALC' 'LDLPFC' 'LIPL' 'LIPS' 'LOPER' 'LT' 'LTRIA'});

%Function to get the Balanced label dataset(i.e.only two labels 1 or 2 aretaken
[example1,labels,expInfo] = idmToExamples_stimulus(info,data,meta,'full');

%Extracting number of rows and columns from the dataset
nrows = size(example1,1);
ncols = size(example1,2);

%Compute the Markov-chain state changes in the fMRI signals(such that the spike in value
%or signal should be represented by 1 from one signal to another, and,
%decrease in value is represented by 0)
scmatrix = zeros(nrows, ncols-1); %initialize the state change matrix. 
%(Note:state changes can only be computed 92609 times for 92610 features)
for row = 1:nrows %loop to include all 50 trials
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
naseeb_TransitionProbability(100);



