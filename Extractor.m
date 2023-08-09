classdef Extractor
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        filename
        saveDirectory
        fhand
        fs
        times
    end

    %% Utility Methods
    methods
        function load_CED_libs(~, path)
            if nargin < 2
                path = "C:\CEDMATLAB\CEDS64ML";
            end
%             if isempty(getenv('CEDS64ML'))
%                 setenv('CEDS64ML', path);
%             end
            setenv('CEDS64ML', path);
            cedpath = getenv('CEDS64ML');
            addpath(cedpath);
            CEDS64LoadLib(cedpath);
        end
        
        function result = removeCloseElements(~, arr, threshold)
            diffArr = diff(arr);  % Calculate the differences between consecutive elements
            closeIndices = abs(diffArr) < threshold;  % Find indices where the difference is less than the threshold
            result = arr([true; ~closeIndices]);  % Include the first element and elements where the difference is not close
        end

    end

    %% Object Methods
    methods
        function obj = Extractor(filename)
            obj.load_CED_libs();
            obj.filename = filename;
            obj.fhand = CEDS64Open(filename, 1);
            if (obj.fhand <= 0); CEDS64ErrorMessage(obj.fhand); unloadlibrary ceds64int; return; end
            obj.fs = 1.0/(CEDS64ChanDiv(obj.fhand, 1)*CEDS64TimeBase(obj.fhand));
        end

        function [output, times] = getAverage(obj, event_channel, wave_channel, timePairs, duration, compDur, minPobs, maxEvents)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            nduration = floor(obj.fs*duration);
            compLength = floor(obj.fs*compDur);
            nPairs = size(timePairs, 1);
            events = [];
            for n = 1: 1: nPairs
                nstart = CEDS64SecsToTicks(obj.fhand, timePairs(n, 1));
                nend = CEDS64SecsToTicks(obj.fhand, timePairs(n, 2));
                [obj, e] = obj.get_events(event_channel, nstart, nend, maxEvents);
                events = cat(1, events, e);
            end
            events = unique(events);
            events = events(1:maxEvents);
            [obj, data] = obj.get_wave(wave_channel, nduration, events);
            if gpuDeviceCount > 0
                subsetIndices = obj.getOptimalAverageGPU(data, compLength, minPobs);
            else
                subsetIndices = obj.getOptimalAverageCPU(data, compLength, minPobs);
            end
            fprintf("Total observations = %d\n", size(data, 1));
            fprintf("Averaged observations = %d\n", length(subsetIndices));
            output = mean(data(subsetIndices, :), 1);
            times = linspace(0, duration, nduration);
        end

        function obj = lock(obj)
            obj.fhand = CEDS64Open(obj.filename, 1);
            if (obj.fhand <= 0); CEDS64ErrorMessage(obj.fhand); unloadlibrary ceds64int; return; end
        end

        function obj = unlock(obj)
            CEDS64Close(obj.fhand);
        end

        function [obj, events] = get_events(obj, channel, nstart, nend, nevents)
            obj = obj.lock();

            [~, events] = CEDS64ReadEvents(obj.fhand, channel, nevents, nstart, nend);
            events = obj.removeCloseElements(events, CEDS64SecsToTicks(2e-3));
            obj = obj.unlock();
        end

        function [obj, data] = get_wave(obj, channel, nduration, events)
            obj = obj.lock();
            nevents = size(events, 1);
            data = zeros(nevents, nduration);
            for i = 1: nevents
                [~, data(i, :), ~] = CEDS64ReadWaveF(obj.fhand, channel, nduration, events(i));
            end
            obj = obj.unlock();
        end

        function subsetIndices = getOptimalAverageCPU(~, data, compLength, minPobs)
            cpuMatrix = data(:, 1:compLength);
            minNobs = floor(minPobs*size(data, 1));
            nobs = size(cpuMatrix, 1);
            binary = dec2bin(0:2^nobs-1, nobs);
            nobsBinary = sum(binary == '1', 2);
            binary = binary(nobsBinary >= minNobs, :);
            varSubsets = zeros(size(binary, 1), 1) + Inf;
            for i = 1: 1: size(binary, 1)
                indices = find(binary(i, :) == '1');
                if ~isempty(indices)
                    varSubsets(i, :) = sum(var(cpuMatrix(indices, :), 0, 1), 2);
                end
            end
            varSubsets = gather(varSubsets);
            [~, minVarIndex] = min(varSubsets);
            subsetIndices = find(binary(minVarIndex, :) == '1');
        end

        function subsetIndices = getOptimalAverageGPU(~, data, compLength, minPobs)
            gpuMatrix = gpuArray(data(:, 1:compLength));
            minNobs = floor(minPobs*size(data, 1));
            nobs = size(gpuMatrix, 1);
            binary = dec2bin(0:2^nobs-1, nobs);
            nobsBinary = sum(binary == '1', 2);
            binary = binary(nobsBinary >= minNobs, :);
            varSubsets = gpuArray(zeros(size(binary, 1), 1) + Inf);
            for i = 1: 1: size(binary, 1)
                indices = find(binary(i, :) == '1');
                if ~isempty(indices)
                    varSubsets(i, :) = sum(var(gpuMatrix(indices, :), 0, 1), 2);
                end
            end
            varSubsets = gather(varSubsets);
            [~, minVarIndex] = min(varSubsets);
            subsetIndices = find(binary(minVarIndex, :) == '1');
        end
    end
end