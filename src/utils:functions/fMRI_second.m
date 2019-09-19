close all;

load data-starplus-04847-v7.mat;

% returns an IDM where the data contains just the voxels belonging to
% {'CALC' 'LDLPFC' 'LIPL' 'LIPS' 'LOPER' 'LT' 'LTRIA'}
[info,data,meta] = transformIDM_selectROIVoxels(info,data,meta,{'CALC' 'LDLPFC' 'LIPL' 'LIPS' 'LOPER' 'LT' 'LTRIA'});

%Balanced label dataset
[example1,labels,expInfo] = EllisidmToExamples_stimulus(info,data,meta,'full');

%Unbalanced label dataset
% [example1,labels,expInfo] = idmToExamples_fixation(info,data,meta,'full');

%Find number of rows and columns to eliminate hardcoding
ncolumns = size(example1,2);
nrows = size(example1,1);

%It creates MC (state changes) only using the first trial data - must be done for all 50
tmps=zeros(nrows,ncolumns-1);
for row = 1:1:nrows
     tmp1=example1(row,1:ncolumns-1);
     tmp2=example1(row,2: ncolumns);
     tmpd=tmp2-tmp1;
     tmpt=sign(tmpd);
     tmps(row,:)=tmpt>0;
end

% %It duplicates the two-state MC; but we need to generate distinct ones using transition probabilities
csmatrix = tmps;

% %Create 50 2x2 transition frequency (count) matrix
Frmatrix = zeros(2,2);
for row = 1:50
    D = csmatrix(row,:);
    t1 = D(:,1:ncolumns-2);
    t2 = D(:,2:ncolumns-1);
    tt3=t1+t2;
    tt4=t2-t1;
    %It counts state changes 0 to 0, 0 to 1, 1 to 0 and 1 to 1
    s00 = sum(tt3<1);
    s11 = sum(tt3>1);
    s01 = sum(tt4>0);
    s10 = sum(tt4<0);
    Frmatrix(1,1) = s00; % count of 00
    Frmatrix(1,2) = s01; % count of 01
    Frmatrix(2,1) = s10; % count of 10
    Frmatrix(2,2) = s11; % count of 11
    matrix{row} = Frmatrix;
end

% Create 50 2x2 Transition probabilities matrix 
prmatrix = zeros(2,2);
for col = 1:length(matrix)
    mat = matrix{col};
    prmatrix(1,1) = mat(1,1)/(mat(1,1)+mat(1,2)); %p00
    prmatrix(1,2) = 1-prmatrix(1,1); %p01
    prmatrix(2,2) = mat(2,2)/(mat(2,2)+mat(2,1)); %p11
    prmatrix(2,1) = 1-prmatrix(2,2); %p10
    matrix2{col} = prmatrix;
end  