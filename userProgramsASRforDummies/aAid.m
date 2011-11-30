classdef aAid
    %AAID Summary of this class goes here
    %   Abstract (pure virtual) class defining the necessary interface of all hearing aid
    %   objects.
    
    properties (Abstract)
        input
        %The input stimulus
        output
        %The output stimulus - recommed that is is formatted to precisely
        %match the dimensions of the input
    end
    
    methods (Abstract)
        obj = processStim(obj)
        %Method that is called to process the input and generate output
    end
    
end

