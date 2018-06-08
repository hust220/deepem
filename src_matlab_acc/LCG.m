classdef LCG < handle
    % https://en.wikipedia.org/wiki/Linear_congruential_generator

    properties
        seed_
        a_
        c_
        modules_
    end

    methods
        function obj = LCG(varargin)
            obj.a_ = int64(1103515245);
            obj.modules_ = 2147483647; 
            obj.c_ = 12345;
            if (nargin >= 1)
                obj.seed_ = varargin{1};
            else
                obj.seed_ = 1; 
            end
        end

        function r = rand(obj, varargin)
            n = length(varargin);

            if (n == 0)
                obj.seed_ = mod(obj.a_ * obj.seed_ + obj.c_, obj.modules_);
                r = double(obj.seed_)/obj.modules_;
            elseif (n == 1)
                s = mod(obj.a_ * varargin{1} + obj.c_, obj.modules_);
                r = double(s)/obj.modules_;
            end
        end

        function r = seed(obj, s)
            obj.seed_ = s;
            r = s;
        end
    end
end

