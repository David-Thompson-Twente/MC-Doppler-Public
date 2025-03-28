function [OutputData, folder, filename] = opendata(directory,varargin)
defaultFileType = 'matrix';
p = inputParser;
defaultSkipLine = 0;

addRequired(p,'directory', @ischar)
addParameter(p,'fileType',defaultFileType, @ischar)
addParameter(p,'skipLine',defaultSkipLine, @isnumeric)

parse(p, directory, varargin{:})
directory = p.Results.directory;
fileType = p.Results.fileType;
skipLine = p.Results.skipLine;
%%
[filename, pathname] = uigetfile('*.*','DefaultName',directory);
if isequal(fileType,'matrix')
    if skipLine == 0
        OutputData = load([pathname filename],'-ascii');
    else
        OutputData = readmatrix([pathname filename],',',skipLine,0);
    end
elseif isequal(fileType,'image')
    OutputData = imread([pathname filename]);
elseif isequal(fileType,'csv')
    OutputData = readmatrix([pathname filename]);
elseif isequal(fileType,'bin')
    Binary = fopen([pathname filename]);
    OutputData = fread(Binary,'*int32');
elseif isequal(fileType,'mat')
    OutputData = load([pathname filename]);
else
    disp('file type not supported')
end
folder = pathname;