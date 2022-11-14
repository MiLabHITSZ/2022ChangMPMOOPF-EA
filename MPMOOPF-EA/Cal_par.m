function [par] = Cal_par(prob_k)
% CEC2021 Real world multi-objective Constrained Optimization Test Suite 
% Abhishek Kumar (email: abhishek.kumar.eee13@iitbhu.ac.in, Indian Institute of Technology (BHU), Varanasi) 

% prob_k -> Index of problem.
% par.n  -> Dimension of the problem.
% par.fn -> Number of objective.
% par.g  -> Number of inequility constraints.
% par.h  -> Number of equality constraints.
% par.xmin -> lower bound of decision variables.
% par.xmax -> upper bound of decision variables.


D        = [4,5,3,4,4,7,4,7,4,2,3,4,7,5,3,2,6,3,10,4,6,9,6,9,2,3,3,7,7,25,25,25,30,30,30,28,28,28,28,34,34,34,34,34,34,34,18,18,18,6];
par.n    = D(prob_k); 
O        = [2,2,2,2,2,2,2,3,2,2,5,2,3,2,2,2,3,2,3,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,3,2,2,3,3,4,2,2,3,2];
par.fn   = O(prob_k);
gn       = [2,5,3,4,4,11,1,9,1,2,7,1,11,8,8,2,9,3,10,7,4,2,1,0,2,1,3,4,9,24,24,24,29,29,29,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
hn       = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,6,0,1,0,4,0,0,0,0,0,0,0,24,24,24,24,26,26,26,26,26,26,26,12,12,12,1];
par.gn    = gn;
par.hn    = hn;
par.g     = gn(prob_k);
par.h     = hn(prob_k);
%% range
% bound constraint definitions for all 18 test functions
xmin46   = -1*ones(1,par.n);
xmin46(27:34) = 0;
xmax46   = +1*ones(1,par.n);
eval(['par.xmin=xmin' int2str(prob_k) ';']);
eval(['par.xmax=xmax' int2str(prob_k) ';' ]);
end