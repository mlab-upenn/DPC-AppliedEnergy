function plot_ml_results(yTrue, yPred, argStr, N)

switch argStr
    
    case 'true-predicted'
        % true-predicted plot
        figure; grid on;
        plot(yTrue(:,1), yPred(1,:)', '.');
        
    case 'predictive'
        % quantitative output comparison
        figure; hold on; grid on;
        l = length(yTrue(:,1));
        h1 = plot(1:l, yPred(1,:), 'LineWidth', 2);
        h3 = plot((1:l)+N-1, yPred(N,:), 'LineWidth', 2);
        h4 = plot(1:l, yTrue(:,1), 'LineWidth', 2);
        legend([h1, h3, h4], 'T', 'T+N', 'ground truth')
        
end
