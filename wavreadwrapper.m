function [data,sampRate] = wavreadwrapper(filename,endianness)

if nargin < 2
    endianness = 'l';
end

try
    [data,sampRate] = wavload(filename);
catch
    [data,sampRate] = read_NIST_file(filename,endianness);
end

data = normalize(data,'wav');

function [data,sampRate,header] = read_NIST_file(filename,endianness)
% [data,sampRate,header] = read_NIST_file(filename)
% - reads a NIST audio file
% - returns the samples in a 1 x D matrix

% The MIT License (MIT)
%
% Copyright (c) 2013 James Lyons
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of
% this software and associated documentation files (the "Software"), to deal in
% the Software without restriction, including without limitation the rights to
% use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
% the Software, and to permit persons to whom the Software is furnished to do so,
% subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
% FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
% COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
% IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


if nargin < 2
    endianness = 'l';
end

file = fopen(filename,'r',endianness);
header = fread(file,1024,'char');
fseek(file,0,'bof');

line = fgetl(file);
if ~strcmp(line(1:4),'NIST')
    fprintf('read_NIST_file: not a valid nist file, "%s"\n',filename);
    data = [];
    sampRate = [];
    return
end
line = fgetl(file);
while ~strncmp(line,'end_head',8)
    if strncmp(line,'sample_count',12)
        [a b nSamples] = strread(line,'%s %s %f');
    elseif strncmp(line,'sample_rate',11)
        [a b sampRate] = strread(line,'%s %s %f');
    end
    line = fgetl(file);
end
fseek(file,1024,'bof');
data = fread(file,[1,nSamples],'int16');
fclose(file);
