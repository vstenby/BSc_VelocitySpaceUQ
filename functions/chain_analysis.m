function chain_analysis(delta, lambda)

% Output for lambda- and delta-chains using sample_plot.
lamdel_chain = [lambda, delta]';
names        = cell(2,1);
names{1}     = '\lambda';
names{2}     = '\delta';
names{3}     = '\alpha';

[taux,acfun] = sample_plot(lamdel_chain,names);
% Plot the autocorrelations functions together.
[~,nacf] = size(acfun);

plot([1:nacf],acfun(1,:),'r',[1:nacf],acfun(2,:),'k--','LineWidth',2)
axis([1,nacf,0,1])
title('ACFs for \lambda and \delta.')
legend(['\tau_{\rm int}(\lambda)=',num2str(taux(1))],['\delta: \tau_{\rm int}(\delta)=',num2str(taux(2))])

end