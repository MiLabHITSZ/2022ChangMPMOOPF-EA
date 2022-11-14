classdef RWMOP46 < handle
% <problem> <RWMOPs>

%------------------------------- Reference --------------------------------
% Abhishek Kumar, Guohua Wu, Mostafa Ali, Qizhang Luo, Rammohan Mallipeddi,
% Ponnuthurai Suganthan, and Swagatam Das, A Benchmark-Suite of Real-World 
% Constrained Multi-Objective Optimization Problems and some Baseline Resu-
% -lts, submitted to Swarm and Evolutionary Computation, 2020.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    properties
    %% Initialization
        M;
        D;
        DM;
        min;
        max;
        encoding;
    end
    methods
        %% Initialization
        function obj =RWMOP46()
            par = Cal_par(46);
            obj.M = par.fn;
            obj.D = par.n;
            obj.DM = 2;
            obj.min  = par.xmin;
            obj.max    = par.xmax;
            obj.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalV(obj,X)
            [PopObj,~,~] = CEC2021_func(X,46);    
        end
        %% Calculate constraint violations
        function PopCon = Cal_Con_h(obj,X)
            [~,~,h] = CEC2021_func(X,46);
            PopCon = max([abs(h)-1e-4],0);
        end 
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = load('nadir_46.txt');
        end
    end
end