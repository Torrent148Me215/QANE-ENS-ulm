
% time parameters (same for all sessions!)
T       = 15000;
Toff    = 5000;
dt      = 2;
t       = (0:dt:T)-Toff;
Nt      = length(t);

% auxiliary: gaussian kernel for smoothing with conv
Nk      = round(500/dt);            % kernel length
knl     = exp( - ( 1.6*(-Nk:2:Nk)/Nk ).^2 );
knl     = knl / sum(knl); 

% loop over all data files
path = 'data/rr014/prefront/';
dd = dir( path );
c = 1;                              % cell-array counter;
for pp=3:length(dd)

    load( [path, dd(pp).name] );    % load one file
    Ntrials = size(result,1)-1;     % # trials

    % loop over seven electrodes
    for nno=1:7;              
        spk = zeros(Ntrials,Nt);    % matrix of spikes
        f1  = zeros(Ntrials,1);     % stimuli f1
        f2  = zeros(Ntrials,1);     % stimuli f2
        d   = zeros(Ntrials,1);     % decisions
        for k=1:Ntrials
            f1(k) = result{k+1,4};
            f2(k) = result{k+1,5};
            d(k)  = sign(f1(k)-f2(k)) * (2*result{k+1,3}-1);
            spiketimes = Toff + result{k+1,6}{nno} - result{k+1,9};
            spk(k,round(spiketimes/dt)) = 1/dt;
        end
        
        % sort spike trains
        f1un = unique(f1);
        Nf1  = length(f1un);
        spksort = cell(Nf1,2);      % groups for f1 and d
        for k=1:Nf1
            spksort{k,1} = spk(  f1==f1un(k) & d==-1, : );
            spksort{k,2} = spk(  f1==f1un(k) & d== 1, : );
        end            

        % compute psth
        psth = zeros(Nf1*2,Nt);     % array for PSTHs
        for k=1:Nf1                 % loop over f1
            for l =1:2              % loop over d
                
                % check that there are one or more trials / group
                if ~isempty( spksort{k,l} )
                    avrate = mean( spksort{k,l}, 1 );
                else
                    avrate = zeros( 1, Nt );
                end
                psth( (l-1)*Nf1+k,:) = conv( avrate, knl, 'same' );

            end
        end
        psth = psth * 1000;         % conversion to Hz
        
        % plot psth
        mp = colormap;
        mp = mp( round(linspace(1,64,Nf1)), : );
        figure(1); clf; hold on;
        for k=1:Nf1                 % loop over f1
            plot( t, psth(k,:), 'Color', mp(k,:) );
            plot( t, psth(Nf1+k,:), '--', 'Color', mp(k,:) );
        end
        pause(0.1);
        
        % store psth
        allpsth{c} = psth;
        allf1{c} = f1un;
        c = c+1;
    end
end

% Generate one big array for all cells with f1=[10 14 18 24 30 34]
f1s = [10 14 18 24 30 34];
c = 1;
for k=1:numel(allpsth)
    cp = (allf1{k}==f1s);
    if (sum(sum(cp)) == 6 )
        thispsth = allpsth{k}( logical([sum(cp') sum(cp')]), : );
        X(c,:,:) = reshape(thispsth,1,12,Nt);
        c=c+1;
    end
end

save romo_allpsth.mat X t f1s
