function plotOnlineDlsMetric2(metric,nItersi,nIters,dt,nReps)
% Syntax:   plotOnlineDlsMetric2(metric,nItersi,nIters,dt,nReps);

% Parse inputs
nIters = [nIters(:)', repmat(nIters(end),1,nReps - numel(nIters))];
nItersi = [nItersi(:)', repmat(nItersi(end),1,nReps - numel(nItersi))];

% Re-arrange metric
nt = (numel(metric) - sum(nItersi)) / sum(nIters);
vals = cell(1,nReps);
off = 0;
for j = 1:nReps
    nk = nItersi(j) + nt * nIters(j);
    vals{j} = metric((off + 1):(off + nk));
    off = off + nk;
end

% Plot online metric
phndl = zeros(1,nReps);
pstr = cell(1,nReps);
cm = linspecer(nReps);
for j = 1:nReps
    tt = 1:nIters(j);
    phndl(j) = plot(1:nItersi(j),vals{j}(1:nItersi(j)),'-','Color',cm(j,:));
    pstr{j} = sprintf('pass %d',j);
    hold on;
    for k = 1:nt
        x0 = k * dt;
        y0 = nItersi(j) + (k - 1) * nIters(j);
        plot(x0 + tt,vals{j}(y0 + tt),'-','Color',cm(j,:));
    end
end
%legend(phndl,pstr{:});
