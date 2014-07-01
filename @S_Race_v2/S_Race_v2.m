classdef S_Race_v2 < handle
% S_Race_v2 class
% S_Race_v2 refers to a combination of the sign test, the discrete Holm's procedure 
% and adaptive alpha scheme 

% COPYRIGHT
%  (C) 2014 Tiantian Zhang, Machine Learning Lab, University of Central Florida
% zhangtt@knights.ucf.edu
% Reference
% T. Zhang, M. Georgiopoulos, G. C. Anagnostopoulos, "Multi-Objective Model
% Selection via Racing", 2014
    
    %% Public Properties
    properties
        M               % number of current models
        currentM        % index of current models, initialized as 1:M        
        batch_size      % number of the test instances every step
        results         % a 3 dimemsional matrix of size (step x batch_size) x M x No_obj matrix, storing the performance vectors of all current models so far  
        Max_step        % assign the maximum number of steps as the stopping criterion
        alpha           % alpha value used in each family of S-Race
        Delta           % overall probability of making any Type I errors
        No_obj          % number of objectives
        fam             % the number of families of current step
    end

    %% Protected Properties
    properties(GetAccess = protected, SetAccess = protected)
        iniM            % number of initial models  
        step            % index of current step, initialized as 1
        dom             % the dominance matrix, dom(i,j) is the no. that i dominates j, dom(i,j) is the no. that j dominates i
    end
    
    %% Public Methods
    methods
        
        % Constructor
        function obj = S_Race_v2(M, Max_step, Delta, batch_size)
            if nargin == 4
                obj.M = M;
                obj.currentM = 1:M;
                % the results are initialized to be empty
                obj.results = [];
                obj.Max_step = Max_step;
                obj.alpha = (1 - Delta)/(M - 1)/Max_step;
                obj.Delta = Delta;
                obj.iniM = M;
                % there are two components in the third dimension of dom
                % the first component corresponds to the number of
                % dominance, and the second one is a flag since we only
                % need to calculate the number of dominance only once in
                % each step. If obj.dom(i,j,2) is 0, it means we haven't
                % calculate obj.dom(i,j) yet; if it is 1, it means we
                % already do the calculation and there is no need to do it
                % anagin
                obj.dom = zeros(obj.iniM, obj.iniM, 2);
                obj.step = 1;
                obj.batch_size = batch_size;
            else
                error('Too few/many arguments');
            end
        end % S_Race
        
        % S-Race
        % the indices of remained models are retuend 
        function retained = Racing(obj, results)
            obj.dom(:,:,2) = zeros(obj.iniM, obj.iniM);
            if obj.step == 1
                obj.No_obj = size(results,3);
            end
            if obj.step > 1
                obj.Delta = obj.Delta + obj.alpha * obj.fam;
                obj.alpha = (1 - obj.Delta)/(obj.M - 1)/(obj.Max_step - obj.step + 1);
            end
            obj.fam = 0;
            % attach new results to previous ones
            obj.results = [obj.results;results];
            % the comparisons start from the first model
            j = 1;
            while(j <= obj.M)
                % initial the p_values for one family
                % the first row stores the p_values for each comparison
                % the second row is a flag, of which 1 means there is a
                % comparison and 0 means there is no comparison
                p_values = zeros(3, obj.M);
                temp_ind = [];
                for t = 1:obj.M
                    if t ~= j
                        % w1 is the number of times that model(j) dominates model(t)
                        % w2 is the number of times that model(t) dominates model(j)
                        if obj.dom(obj.currentM(j),obj.currentM(t),2) == 0 % we haven't calculate the number of dominance at current step
                            if obj.step == 1 
                                [w1, w2] = obj.dominates(obj.results(:,j,:), obj.results(:,t,:));
                                obj.dom(obj.currentM(j),obj.currentM(t),1) = w1;
                                obj.dom(obj.currentM(t),obj.currentM(j),1) = w2;
                                obj.dom(obj.currentM(j),obj.currentM(t),2) = 1;
                                obj.dom(obj.currentM(t),obj.currentM(j),2) = 1;
                            else
                                [w1, w2] = obj.dominates(obj.results(end - obj.batch_size + 1:end,j,:), obj.results(end - obj.batch_size + 1: end,t,:));
                                obj.dom(obj.currentM(j),obj.currentM(t)) = obj.dom(obj.currentM(j),obj.currentM(t)) + w1;
                                obj.dom(obj.currentM(t),obj.currentM(j)) = obj.dom(obj.currentM(t),obj.currentM(j)) + w2;
                                w1 = obj.dom(obj.currentM(j),obj.currentM(t),1);
                                w2 = obj.dom(obj.currentM(t),obj.currentM(j),1); 
                                % change the flag
                                obj.dom(obj.currentM(j),obj.currentM(t),2) = 1;
                                obj.dom(obj.currentM(t),obj.currentM(j),2) = 1;
                            end
                        else % we already calculte the number of dominance at current step
                            w1 = obj.dom(obj.currentM(j),obj.currentM(t));
                            w2 = obj.dom(obj.currentM(t),obj.currentM(j)); 
                        end
                        if(w1 == 0 && w2 == 0) % it means the two models are non-dominated to each other or have same performance
                            if (isequal(obj.results(:,j,:), obj.results(:,t,:)) == 1 && obj.step == obj.Max_step)
                                % same models are removed until the end
                                temp_ind = [temp_ind, t];
                            end
                        elseif(w1 > w2) % we only do comparison if j is better than t
                            % sign test
                            p_values(1,t) = 1 - binocdf(w1-1, w1 + w2, 0.5);          
                            p_values(2,t) = 1;
                            p_values(3,t) = w1 + w2;
                        end
                    end
                end
                % if there is a family, add the family number by 1 
                if(max(max(p_values)) >= 1)
                    obj.fam = obj.fam + 1;
                end
                % Adopt FWER (Holm's procedure 1987) and abandon model(t) where 
                % the null hypothesis is rejected
                [index] = obj.DiscreteHolm(p_values, obj.alpha);
                % If we are at the end of racing, models have excatly the
                % same performance with others are removed
                index = [index, temp_ind];        
                if isempty(index) == 0
                    index = sort(index);
                    obj.currentM(index) = [];
                    for k = 1:length(index)
                        obj.results(:, index(k) - k + 1,:) = [];
                    end
                end
                % the next family is to compare model j with all the other
                % retained models
                j = j + 1 - size(find(index < j), 2);
                obj.M = length(obj.currentM);
            end
            obj.step = obj.step + 1;
            retained = obj.currentM;
        end % Racing              
    end
    
    methods (Access = private)
        % find the number of dominance 
        function [w1, w2] = dominates(~, results1, results2)
            [a,~,~] = size(results1);
            w1 = 0;
            w2 = 0;
            for i = 1:a
                temp = results1(i,1,:) - results2(i,1,:);
                % for minimization problem
                if max(temp) <=0 && min(temp) < 0
                    w1 = w1 + 1;
                elseif min(temp) >=0 && max(temp) > 0
                    w2 = w2 + 1;
                end
            end
        end % dominates
        
        % Discrete Holm's procedure
        function index = DiscreteHolm(obj,p_values, alpha)
            % find out the index that contains hypothesis
            temp_index = find(p_values(2,:) == 1);
            % find out the p_values involved in the family
            temp = p_values(1, temp_index);
            N = p_values(3,temp_index);
            % sor the p_values in ascending order
            [sorted_p,corre_index] = sort(temp);
            temp_in = [];
            % Discrete Holm's procedure
            for i = 1:size(temp,2)
                temp = obj.CalP(sorted_p(i), N(corre_index(i:end)));
                if(temp > alpha)
                    temp_in = corre_index(1:1:i - 1);
                    break;
                end
                if (i == size(temp,2))
                    temp_in = corre_index;
                end
            end
            index = temp_index(temp_in);
        end % DiscreteHolm
        
        % Auxiliary function needed in discrete Holm
        function p = CalP(~,minP,Ns)
            P_s = zeros(1,size(Ns,2));
            for i = 1:size(Ns,2)
                for j = 1:Ns(i) + 1
                    temp_p = 1 - binocdf(j - 2, Ns(i), 0.5);
                    if(temp_p <= minP)
                        P_s(i) = temp_p;
                        break;
                    end
                end
            end
            p = sum(P_s);
        end % CalP
    end
    
end
