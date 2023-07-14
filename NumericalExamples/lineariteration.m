classdef lineariteration
    %LINEARITERATION This is a service class for storing the results in
    %terms of linear iterations inside the Riemannian Optimization
    %algorithms. It is *not* an efficient way of doing so and should be
    %used only for investigating convergence.

    properties
        iteration   % Will contain iteration count
        residuals   % Cell variable containing all the convergence histories
    end

    methods
        function obj = lineariteration()
            obj.iteration = [];
            obj.residuals = {};
        end

        function obj = additeration(obj,iter)
            % Stores iteration count
            obj.iteration = [obj.iteration;iter];
        end
        function obj = addresidual(obj,residual)
            % Stores normalized residuals
            obj.residuals{end+1} = residual./residual(1);
        end
        function allresiduals = plotaveragecurve(obj)
            % Many linear solve are performed inside the optimization
            % routine, this function transforms them in a matrix (with Nan
            % for misssing data) so that a statistical analysis of the
            % performance can be done.
            M = max(cellfun(@length, obj.residuals));
            equalized = ...
                cellfun(@(x) [x; NaN(M - numel(x), 1)], ...
                obj.residuals, 'un', 0);
            allresiduals = cell2mat(equalized);
        end

    end
end

