% THIS CODE TESTS SYNTHETIC SIGNALS 

% first signal is a sine wave with amplitude 4
% second signal is a higher frequency sine wave
% third signal is the sum of the two

% first "band" added is eq14_1 -- corresponds to higher frequency band
% second "band" added is eq14_2 -- lower frequency band
% sum of the two

% run over each signal to see if variance of individual signals
% recovered in sum of two signals
%%

save_on = false;
n=2; i = 13;



% SIGNAL 1
deltad = pi/50;
t=0:deltad:50*pi;
T = t(end)+deltad;
synth = 4*sin(t);
figure()
plot(t,4*sin(t))
dowave_sat(synth,deltad,n,t,synth,[],save_on,[],i);

% SIGNAL 2
deltad = pi/50;
t=0:deltad:50*pi;
T = t(end)+deltad;
synth = sin(t/pi/2);
figure()
plot(t,sin(t/pi/2))
dowave_sat(synth,deltad,n,t,synth,[],save_on,[],i);

% SIGNAL 3
deltad = pi/50;
t=0:deltad:50*pi;
T = t(end)+deltad;
synth = 4*sin(t) + sin(t/pi/2);
figure()
plot(t,4*sin(t),t,sin(t/pi/2),t,4*sin(t) + sin(t/pi/2))
dowave_sat(synth,deltad,n,t,synth,[],save_on,[],i);



