%%
clear;
filename = './spike2files/2021-12-10_4_550Hz_ICN_NBQX_Intra.smrx';
extractor = Extractor(filename, saveDirectory);

%%
eventChannel = 23;
waveChannel = 5;

tstart = 515.11;
tend = 530.64;

duration = 0.98;
compDur = 0.3;

minPobs = 0.5;
maxEvents = 1000;

%%
[output, times] = extractor.getAverage(eventChannel, waveChannel, tstart, tend, duration, compDur, minPobs, maxEvents);

%%
figure();
plot(times, output);
title(saveFilename);
xlabel('sec');
ylabel('mV');
