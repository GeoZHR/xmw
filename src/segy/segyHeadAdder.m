addpath('S4M/Geophysics_3.0/')

% input a dataset (without header) and a segy file, 
% and replace segy traces with but keep its header.
% the input dataset should have the same dimensions as the segy file.

nx = 66;   %inline numbers
ny = 46;   %crossline numbers
nt = 230;  %vertical samples per trace

segyFile = '../../../data/seis/test/sd2mshead.sgy'; % input segy
dataFile = '../../../data/seis/test/synsd60amp.dat';% input data
sgyWrite = '../../../data/seis/test/sd2msheadOut.sgy'; % output segy

dataId = fopen(dataFile);
data = fread(dataId,nx*ny*nt,'single');
segy = read_segy_file(segyFile);

dataR = reshape(data,[nt,nx*ny]);  % reshape a column vector to a 2D matrix

segy.traces = dataR; % replace data
write_segy_file(segy,sgyWrite); % write segy file


sgyout = read_segy_file(sgyWrite); % for test only




