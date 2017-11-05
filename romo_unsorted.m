load data/rr014/prefront/R14001_001.mat;

% auxiliary parameters
Ntrials = size(result,1)-1;
T      = 12000;
Toff   = 5000;
dt     = 2;
t      = (0:dt:T)-Toff;
Nt     = length(t);

% extract spike trains from one electrode
eno    = 1;                    % electrode number
spk    = zeros(Ntrials,Nt);    % array for spike times
for k=1:Ntrials
    spiketimes = Toff + result{k+1,6}{eno} - result{k+1,9};
    spk(k,round(spiketimes/dt)) = 1/dt;
end

% rasterplot
[x,y] = find(spk>0);
plot( t(y), x, 'k.' );
