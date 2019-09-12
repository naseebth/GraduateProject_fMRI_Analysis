function [pr_matrix_n] = naseeb_TransitionProbability(n)
%Calculate n-step probability for simple asymtotic analys Summary of this function goes here
% The function inputs n-parameter which is the number of times the
% transition probability matrix is to be multiplied.
% The function outputs value P(n) i.e. pr_matrix_n, which gives the value
% of matrix when transition probability is multiplied n-times.

global pr_matrix_final; %Declare as global variable so that it can be used in the main function
for col = 1:length(pr_matrix_final)
    pr_matn = pr_matrix_final{col};
    pr_matrix_n{col} = pr_matn^n; %multiply by n times. This is a short-cut way to do it for analysis. Will be displayed in an elaborate manner also.
end  

end