dataType = 'humanAirflow';
bmObj = breathmetrics(respiratoryTrace, srate, dataType);
bmObj.estimateAllFeatures();
fig = bmObj.plotCompositions('raw');
fig = bmObj.plotFeatures({'extrema','maxflow'});
fig = bmObj.plotCompositions('normalized');
fig = bmObj.plotCompositions('line');
addpath 'C:\Toolboxes\breathmetrics-master'
dataType = 'humanAirflow';
bmObj = breathmetrics(respiratoryTrace, srate, dataType);
bmObj.estimateAllFeatures();
fig = bmObj.plotCompositions('raw');
fig = bmObj.plotFeatures({'extrema','maxflow'});
fig = bmObj.plotCompositions('normalized');
fig = bmObj.plotCompositions('line');
clc
dataType = 'humanAirflow';
bmObj = breathmetrics(respiratoryTrace, srate, dataType);
bmObj.estimateAllFeatures();
fig = bmObj.plotCompositions('raw');
fig = bmObj.plotFeatures({'extrema','maxflow'});
fig = bmObj.plotCompositions('normalized');
fig = bmObj.plotCompositions('line');
bmObj.inhaleOnsets
bmObj.time(bmObj.inhaleOnsets)
bmObj.time(bmObj.inhaleOffsets)
editedBm=bmGui(bmObj)
bmObj