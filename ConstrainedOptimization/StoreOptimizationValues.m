classdef StoreOptimizationValues < handle
    %Class to store the optimization processs
    %   Needs to be a handle to be passed to a function
    
    properties
        targetFunctionValues
        variableValues
    end
    
    methods
        function obj = StoreOptimizationValues(maxIterations,numberVariables)
            obj.variableValues = zeros(maxIterations,numberVariables);
            obj.targetFunctionValues = zeros(maxIterations,1);
        end
        
        function addValues(obj,index,fval,variables)
            obj.variableValues(index,:) = variables;
            obj.targetFunctionValues(index) = fval;
        end
    end
end