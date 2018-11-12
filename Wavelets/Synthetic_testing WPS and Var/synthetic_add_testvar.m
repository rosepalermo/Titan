% test synthetic

save_on = false;
n=2; i = 13;



% ADD SINE WAVES
deltad = pi/50;
t=0:deltad:50*pi;
T = t(end)+deltad;
synth = 4*sin(t);
figure()
plot(t,4*sin(t))
dowave_sat(synth,deltad,n,t,synth,[],save_on,[],i);

deltad = pi/50;
t=0:deltad:50*pi;
T = t(end)+deltad;
synth = sin(t/pi/2);
figure()
plot(t,sin(t/pi/2))
dowave_sat(synth,deltad,n,t,synth,[],save_on,[],i);

deltad = pi/50;
t=0:deltad:50*pi;
T = t(end)+deltad;
synth = 4*sin(t) + sin(t/pi/2);
figure()
plot(t,4*sin(t),t,sin(t/pi/2),t,4*sin(t) + sin(t/pi/2))
dowave_sat(synth,deltad,n,t,synth,[],save_on,[],i);



