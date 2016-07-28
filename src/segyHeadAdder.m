% input a dataset (without header) and a segy file, 
% and replace segy traces with but keep its header.
% the input dataset should have the same dimensions as the segy file.

nx = 66;   %inline numbers
ny = 46;   %crossline numbers
nt = 230;  %vertical samples per trace

segyFile = '../../../data/seis/test/sd2mshead.sgy';     % name of the segy file
dataFile = '../../../data/seis/test/synsd60amp.dat';    % name of the data file

dataId = fopen(dataFile);
data = fread(dataId,nx*ny*nt,'single');
segy = read_segy_file(segyFile);
