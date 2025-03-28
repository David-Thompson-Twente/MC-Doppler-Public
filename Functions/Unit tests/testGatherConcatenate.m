%%  testGatherConcatenate.m
%   David Thompson, 28-03-2024
%   Unit tests for gatherPhotonsGPU and concatenatePhotonStruct, to see if the gathering and combining of
%   structs works as intended
%
%   Inputs:     struct fail_pass carrying any information about previously tested functions
%
%   Outputs:    PASS_ALL if all subtests succesfull
%               FAIL message indicating which part failed
%%
function fail_pass = testGatherConcatenate(fail_pass,options)
arguments
    fail_pass           struct
    options.numEntries  (1,1)   {mustBeNumeric,mustBePositive}  = 1e6;
end
fail_pass.gatherPhotonsGPU = 'FAIL';
fail_pass.concatenatePhotonStruct = 'FAIL';

testStruct.entry1 = gpuArray(ones(1,options.numEntries));
testStruct.entry2 = ones(1,options.numEntries);
testStruct.entry3 = gpuArray(rand(1,options.numEntries));
testStruct.entry4 = 'StringTest';

gatheredStruct = gatherPhotonsGPU(testStruct);
concatenatedStruct = concatenatePhotonStruct(gatheredStruct,gatheredStruct);
testFieldNames = fieldnames(testStruct);
gatheredFieldNames = fieldnames(gatheredStruct);
concatenatedFieldNames = fieldnames(concatenatedStruct);
for n = 1:length(testFieldNames)
    if ~strcmp(testFieldNames{n},gatheredFieldNames{n})
        fail_pass.gatherPhotonsGPU = [fail_pass.gatherPhotonsGPU '_FIELD_NAME_NOT_SAME'];
    end
    if min(testStruct.(testFieldNames{n}) == gatheredStruct.(gatheredFieldNames{n})) == 0
        fail_pass.gatherPhotonsGPU = [fail_pass.gatherPhotonsGPU '_FIELD_CONTENT_NOT_SAME'];
    end
    if isa(gatheredStruct.(gatheredFieldNames{n}),'gpuArray')
        fail_pass.gatherPhotonsGPU = [fail_pass.gatherPhotonsGPU '_STILL__GPUARRAY'];
    end
    if ~strcmp(concatenatedFieldNames{n},gatheredFieldNames{n})
        fail_pass.concatenatePhotonStruct= [fail_pass.gatherPhotonsGPU '_FIELD_NAME_NOT_SAME'];
    end
    if length(concatenatedStruct.(concatenatedFieldNames{n})) ~= 2*length(gatheredStruct.(gatheredFieldNames{n}))
        fail_pass.concatenatePhotonStruct= [fail_pass.gatherPhotonsGPU '_ENTRIES_OMITTED'];
    end
    if min(concatenatedStruct.(concatenatedFieldNames{n})(1:length(gatheredStruct.(gatheredFieldNames{n}))) == concatenatedStruct.(concatenatedFieldNames{n})(length(gatheredStruct.(gatheredFieldNames{n}))+1:2*length(gatheredStruct.(gatheredFieldNames{n})))) == 0
        fail_pass.concatenatePhotonStruct= [fail_pass.gatherPhotonsGPU '_WRONG_VALUES_CONCATENATED'];
    end

end
if strcmp(fail_pass.gatherPhotonsGPU,'FAIL')
    fail_pass.gatherPhotonsGPU = 'PASS';
end
if strcmp(fail_pass.concatenatePhotonStruct,'FAIL')
    fail_pass.concatenatePhotonStruct = 'PASS';
end
testedFunctions = {'gatherPhotonsGPU','concatenatePhotonStruct'};
passCounter = 0;
fail_pass.gatherConcatenate = '';
for n = 1:length(testedFunctions)
    if strcmp(fail_pass.(testedFunctions{n}),'PASS')
        fail_pass = rmfield(fail_pass,testedFunctions{n});
        passCounter = passCounter + 1;
    else
        fail_pass.gatherConcatenate = [fail_pass.gatherConcatenate 9 testedFunctions{n} '_' fail_pass.(testedFunctions{n})];
        fail_pass = rmfield(fail_pass,testedFunctions{n});
    end
end
if passCounter == length(testedFunctions)
    fail_pass.gatherConcatenate = 'PASS_ALL';
end
end