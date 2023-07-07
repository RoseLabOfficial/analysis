classdef Extractor
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        filename
        fhand
        times
    end

    methods
        function load_CED_libs(path)
            if nargin < 2
                path = "C:\CEDMATLAB\CEDS64ML";
            end
            if isempty(getenv('CEDS64ML'))
                setenv('CEDS64ML', path);
            end
            cedpath = getenv('CEDS64ML');
            addpath(cedpath);
            CEDS64LoadLib(cedpath);
        end

        function obj = Extractor(filename)
            if nargin > 1
                obj.filename = filename;
                
            end
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end