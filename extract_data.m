%%
clear;
filename = './spike2files/2021-12-10_4_550Hz_ICN_NBQX_Intra.smrx';
extractor = Extractor(filename);

%%
eventChannel = 11;
waveChannel = 15;

timePairs = [849.10, 867.03; 1010.96, 1034.31; 1734.96, 1748.10; 1937.33, 1941.94];

duration = 0.98;
compDur = 0.67;

minPobs = 0.5;
maxEvents = 20;

%%
[output, times] = extractor.getAverage(eventChannel, waveChannel, timePairs, duration, compDur, minPobs, maxEvents);

%%
figure();
plot(times, output);
xlabel('sec');
ylabel('mV');