addpath('../segy/S4M/Geophysics_3.0/')

% input a dataset (without header) and a segy file, 
% and replace segy traces with but keep its header.
% the input dataset should have the same dimensions as the segy file.

nx = 7000;   %inline numbers
nt = 1500;  %vertical samples per trace

segyFile = '../../../data/seis/tjxd/2d/interp/seis.sgy'; % input segy
dataFile = '../../../data/seis/tjxd/2d/interp/deni.dat';% input data
sgyWrite = '../../../data/seis/tjxd/2d/interp/deni.sgy'; % output segy

dataId = fopen(dataFile);
data = fread(dataId,nx*nt,'single','b');
segy = read_segy_file(segyFile);

dataR = reshape(data,[nt,nx]);  % reshape a column vector to a 2D matrix

segy.traces = dataR; % replace data
write_segy_file(segy,sgyWrite); % write segy file


sgyout = read_segy_file(sgyWrite); % for test only




