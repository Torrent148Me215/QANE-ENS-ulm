load romo_allpsth.mat

% reshape three-dim array into matrix: (#neurons) x (#data points)
sz    = size(X);
Z     = reshape( X, sz(1), sz(2)*sz(3) );

% center the data
Z     = Z - mean(Z')';

% compute covariance matrix
C     = cov( Z' );

% compute eigenvalues and eigenvectors
[V,D] = eig( C );

% extract first nine principal components
PCZ   = V(:,end-8:end)'*Z;

% reshape back into three-d array: (# pcs) x (# conditions) x (# time points)
PCX   = reshape( PCZ, 9, sz(2), sz(3) );

% create color map for plotting
mp = colormap;
mp = mp( round(linspace(1,64,6)), : );

% plot first K principal components (in order)
K = 5;
Nf1   = length(f1s);
for i=1:K
    figure(i); clf;
    hold on;
    for k=1:Nf1                 % loop over f1
        plot( t, squeeze(PCX(7-i,k,:)), 'Color', mp(k,:) );
        plot( t, squeeze(PCX(7-i,Nf1+k,:)), '--', 'Color', mp(k,:) );
    end
    axis([ -5000 12000 -130 130] );
    set(gca,'XTick', [0 500 3500 4000] );
    set(gca,'XTickLabel', {'0', '0.5', '3.5', '4'} );
end

