load data/rr014/prefront/R14001_001.mat;

% auxiliary parameters
Ntrials = size(result,1)-1;
T      = 12000;
Toff   = 5000;
dt     = 2;
t      = (0:dt:T)-Toff;
Nt     = length(t);

% extract spike trains from one electrode
nno    = 2;                    % electrode number
spk    = zeros(Ntrials,Nt);    % array for spike times
f1     = zeros(Ntrials,1);     % array of stimuli f1
f2     = zeros(Ntrials,1);     % array of stimuli f2
d      = zeros(Ntrials,1);     % array of decisions (-1,1)
for k=1:Ntrials
    f1(k) = result{k+1,4};
    f2(k) = result{k+1,5};
    d(k)  = sign(f1(k)-f2(k)) * (2*result{k+1,3}-1);
    spiketimes = Toff + result{k+1,6}{nno} - result{k+1,9};
    spk(k,round(spiketimes/dt)) = 1/dt;
end

% sort spike trains into groups according to trial type
f1un = unique(f1);
spksort = cell(length(f1un),2); % cell array for each group
for k=1:length(f1un)
    spksort{k,1} = spk(  f1==f1un(k) & d==-1, : );
    spksort{k,2} = spk(  f1==f1un(k) & d== 1, : );
end

% plot sorted spike trains
clf; hold on;
offs = 0;
for k=1:length(f1un)
    
    % decision left in red
    [x,y]  = find(spksort{k,1}>0);
    plot( t(y), x+offs, 'r.' );
    offs = offs + size(spksort{k,1},1);
    
    % decision right in blue
    [x,y] = find(spksort{k,2}>0);
    plot( t(y), x+offs, 'b.' );
    offs = offs + size(spksort{k,2},1) + 5;

end
