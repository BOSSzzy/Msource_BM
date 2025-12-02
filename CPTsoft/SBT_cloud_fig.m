clc; clear; close all;

%% Parameter Settings
% Graphic parameters
fontName = 'Times New Roman';
fontSize = 10;
titleSize = 14;
labelSize = 12;
lineWidth = 1.5;
markerSize = 30;

% Number of grid points (increase for higher precision)
numPoints = 100;

% Colormap type
colormap_type = 'parula'; % Options: 'jet', 'parula', 'hsv', 'hot', 'cool', 'spring', 'summer', 'autumn', 'winter', 'gray'

%% Read Data
try
    % Read parameter data from CSV file
    data = readtable('SCdata.csv');
    
    % Extract data columns
    % Assuming CSV columns: soiltype, mean_fs, mean_qe, std_fs, std_qe, cor, Prior
    soil_types = data.soiltype; 
    fs_mean = data.mean_fs;
    qE_mean = data.mean_qe;
    fs_std = data.std_fs;
    qE_std = data.std_qe;
    corr_values = data.cor;
    prior_probs = data.Prior;
    
    % Data validity check
    if any(fs_mean <= 0) || any(qE_mean <= 0)
        warning('Zero or negative values found in data, which may cause issues in logarithmic plotting.');
        % Remove invalid data
        valid_idx = fs_mean > 0 & qE_mean > 0;
        soil_types = soil_types(valid_idx);
        fs_mean = fs_mean(valid_idx);
        qE_mean = qE_mean(valid_idx);
        fs_std = fs_std(valid_idx);
        qE_std = qE_std(valid_idx);
        corr_values = corr_values(valid_idx);
        prior_probs = prior_probs(valid_idx);
    end
    
    % Note: Using 'exam_CPT.csv' instead of 'SCdata.csv' for test data 
    try
        test_data = readtable('exam_CPT.csv');
        test_fs = test_data.fs;  % First column: fs
        test_qE = test_data.qe;  % Second column: qE
        
        % Data validity check for test data
        if any(test_fs <= 0) || any(test_qE <= 0)
            warning('Zero or negative values found in test data, which may cause issues in logarithmic scale.');
            % Remove invalid data
            valid_test_idx = test_fs > 0 & test_qE > 0;
            test_fs = test_fs(valid_test_idx);
            test_qE = test_qE(valid_test_idx);
        end
        
        has_test_data = true;
        fprintf('Successfully read test point data, total %d points.\n', length(test_fs));
    catch ME
        warning('Failed to read test point data: %s', ME.message);
        has_test_data = false;
    end
catch ME
    error('Failed to read data: %s', ME.message);
end

%% Calculate Covariance Matrices for Each Soil Type
unique_types = unique(soil_types);
num_types = length(unique_types);

% Store mean vectors and covariance matrices for each soil type
mu_vectors = cell(num_types, 1);
cov_matrices = cell(num_types, 1);
corr_coeffs = zeros(num_types, 1);
prior_probabilities = zeros(num_types, 1);

for i = 1:num_types
    type_idx = strcmp(soil_types, unique_types{i});
    
    % Mean vector
    mu_fs = mean(fs_mean(type_idx));
    mu_qE = mean(qE_mean(type_idx));
    mu_vectors{i} = [mu_fs; mu_qE];
    
    % Standard deviation
    sigma_fs = mean(fs_std(type_idx)); 
    sigma_qE = mean(qE_std(type_idx)); 
    
    % Correlation coefficient
    corr_coeff = mean(corr_values(type_idx));
    corr_coeffs(i) = corr_coeff;
    
    % Prior probability
    prior_probabilities(i) = mean(prior_probs(type_idx));
    
    % Covariance matrix
    cov_fs_qE = corr_coeff * sigma_fs * sigma_qE; 
    cov_matrices{i} = [sigma_fs^2, cov_fs_qE; 
                       cov_fs_qE, sigma_qE^2];
end

%% Define Colorbar Limits and Axis Limits for Subplots
colorbar_limits = [
    0, 0.15;    % Soil 1
    0, 0.01;    % Soil 2
    0, 0.03;    % Soil 3
    0, 0.0015;  % Soil 4
];

posterior_colorbar_limits = [
    0, 1.0;
    0, 1.0;
    0, 1.0;
    0, 1.0;
];

axis_limits = [
    1, 1000, 0.1, 100;
    1, 1000, 0.1, 100;
    1, 1000, 0.1, 100;
    1, 1000, 0.1, 100;
];

%% Create Probability Density Plots
fig1 = figure('Color', 'white', 'Position', [100, 100, 1000, 800]);
set(fig1, 'Name', 'SBT_Probability_Density', 'NumberTitle', 'off');

% Create log-space grid
[X, Y] = meshgrid(logspace(log10(1), log10(1000), numPoints), ...
                  logspace(log10(0.1), log10(100), numPoints));

% Flatten grid for vectorized calculation
grid_points = [X(:), Y(:)]; % N^2 x 2

% Parameters for boundary lines (SBT chart lines)
% Line [12,1]-[1000,3]
m2 = log10(3/1) / log10(1000/12); 
b2 = log10(1) - m2 * log10(12); 

% Line [7,1]-[1000,10]
m3 = log10(10/1) / log10(1000/7); 
b3 = log10(1) - m3 * log10(7); 

% Line [1,1]-[1000,30]
m4 = log10(30/1) / log10(1000/1); 
b4 = log10(1) - m4 * log10(1); 

density_max_values = zeros(num_types, 1);
density_min_values = zeros(num_types, 1);
Z_all = zeros(numPoints, numPoints, num_types);

for i = 1:min(num_types, 4)
    subplot(2, 2, i);
    
    mu = mu_vectors{i};
    Sigma = cov_matrices{i};
    det_Sigma = det(Sigma);
    inv_Sigma = inv(Sigma);
    
    % Vectorized calculation of PDF
    % X_centered = X - mu
    diff = grid_points - mu'; % N^2 x 2
    % Mahalanobis distance: diag(diff * inv_Sigma * diff')
    % Efficiently: sum((diff * inv_Sigma) .* diff, 2)
    mahal_dist = sum((diff * inv_Sigma) .* diff, 2);
    
    % PDF
    Z_vec = (1/(2*pi*sqrt(det_Sigma))) * exp(-0.5 * mahal_dist);
    Z = reshape(Z_vec, numPoints, numPoints);
    
    Z_all(:,:,i) = Z;
    density_max_values(i) = max(Z(:));
    density_min_values(i) = min(Z(:));
    
    % Plot
    [~, h] = contourf(X, Y, Z, 20);
    set(h, 'LineStyle', 'none');
    colormap(colormap_type);
    
    cb = colorbar;
    if i <= size(colorbar_limits, 1)
        caxis(colorbar_limits(i, :));
    end
    
    set(gca, 'XScale', 'log', 'YScale', 'log');
    if i <= size(axis_limits, 1)
        xlim([axis_limits(i, 1), axis_limits(i, 2)]);
        ylim([axis_limits(i, 3), axis_limits(i, 4)]);
    end
    
    % Draw boundaries
    hold on;
    line([1, 12], [1, 1], 'Color', 'k', 'LineWidth', lineWidth); 
    line([12, 12], [1, 0.1], 'Color', 'k', 'LineWidth', lineWidth); 
    
    x_line3 = logspace(log10(12), log10(1000), numPoints); 
    y_line3 = 10.^(m2 * log10(x_line3) + b2); 
    plot(x_line3, y_line3, 'k-', 'LineWidth', lineWidth); 
    
    x_line4 = logspace(log10(7), log10(1000), numPoints); 
    y_line4 = 10.^(m3 * log10(x_line4) + b3); 
    plot(x_line4, y_line4, 'k-', 'LineWidth', lineWidth); 
    
    x_line5 = logspace(log10(1), log10(1000), numPoints); 
    y_line5 = 10.^(m4 * log10(x_line5) + b4); 
    plot(x_line5, y_line5, 'k-', 'LineWidth', lineWidth);
    hold off;
    
    xlabel('Sleeve Friction, f_s (kPa)', 'FontName', fontName, 'FontSize', labelSize);
    ylabel('Effective Corrected Cone Resistance, q_E (MPa)', 'FontName', fontName, 'FontSize', labelSize);
    
    soil_name = unique_types{i};
    title([soil_name, ' (corr=', num2str(corr_coeffs(i), '%.3f'), ', prior=', num2str(prior_probabilities(i), '%.3f'), ')'], ...
          'FontName', fontName, 'FontSize', titleSize);
    
    grid on;
    box on;
    set(gca, 'LineWidth', 1.2, 'FontName', fontName, 'FontSize', fontSize);
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
    set(gca, 'TickDir', 'in');
end

set(fig1, 'PaperPositionMode', 'auto');
sgtitle('Probability Density Cloud Map', 'FontName', fontName, 'FontSize', titleSize+2, 'FontWeight', 'bold');
saveas(fig1, 'SBT_probability_density.png');

%% Posterior Probability Plots
fig2 = figure('Color', 'white', 'Position', [150, 150, 1000, 800]);
set(fig2, 'Name', 'SBT_Posterior_Probability', 'NumberTitle', 'off');

% Calculate Joint and Total Probability
joint_prob = zeros(numPoints, numPoints, num_types);
for i = 1:num_types
    joint_prob(:,:,i) = prior_probabilities(i) * Z_all(:,:,i);
end

total_prob = sum(joint_prob, 3);

posterior_prob = zeros(numPoints, numPoints, num_types);
for i = 1:num_types
    posterior_prob(:,:,i) = joint_prob(:,:,i) ./ (total_prob + eps);
end

posterior_max_values = zeros(num_types, 1);
posterior_min_values = zeros(num_types, 1);

for i = 1:min(num_types, 4)
    subplot(2, 2, i);
    
    Z_posterior = posterior_prob(:,:,i);
    posterior_max_values(i) = max(Z_posterior(:));
    posterior_min_values(i) = min(Z_posterior(:));
    
    [~, h] = contourf(X, Y, Z_posterior, 20);
    set(h, 'LineStyle', 'none');
    colormap(colormap_type);
    
    cb = colorbar;
    if i <= size(posterior_colorbar_limits, 1)
        caxis(posterior_colorbar_limits(i, :));
    end
    
    set(gca, 'XScale', 'log', 'YScale', 'log');
    if i <= size(axis_limits, 1)
        xlim([axis_limits(i, 1), axis_limits(i, 2)]);
        ylim([axis_limits(i, 3), axis_limits(i, 4)]);
    end
    
    % Draw boundaries
    hold on;
    line([1, 12], [1, 1], 'Color', 'k', 'LineWidth', lineWidth); 
    line([12, 12], [1, 0.1], 'Color', 'k', 'LineWidth', lineWidth); 
    plot(x_line3, y_line3, 'k-', 'LineWidth', lineWidth); 
    plot(x_line4, y_line4, 'k-', 'LineWidth', lineWidth); 
    plot(x_line5, y_line5, 'k-', 'LineWidth', lineWidth);
    hold off;
    
    xlabel('Sleeve Friction, f_s (kPa)', 'FontName', fontName, 'FontSize', labelSize);
    ylabel('Effective Corrected Cone Resistance, q_E (MPa)', 'FontName', fontName, 'FontSize', labelSize);
    
    title([unique_types{i}, ' (Posterior Prob)'], 'FontName', fontName, 'FontSize', titleSize);
    
    grid on; box on;
    set(gca, 'LineWidth', 1.2, 'FontName', fontName, 'FontSize', fontSize);
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
    set(gca, 'TickDir', 'in');
end

set(fig2, 'PaperPositionMode', 'auto');
sgtitle('Posterior Probability Cloud Map', 'FontName', fontName, 'FontSize', titleSize+2, 'FontWeight', 'bold');
saveas(fig2, 'SBT_posterior_probability.png');

% Print statistics
fprintf('\nProbability Density Statistics:\n');
fprintf('%-20s %-20s %-20s\n', 'Soil Type', 'Max Density', 'Min Density');
for i = 1:min(num_types, 4)
    fprintf('%-20s %-20.6f %-20.6f\n', unique_types{i}, density_max_values(i), density_min_values(i));
end

fprintf('\nPosterior Probability Statistics:\n');
fprintf('%-20s %-20s %-20s\n', 'Soil Type', 'Max Posterior', 'Min Posterior');
for i = 1:min(num_types, 4)
    fprintf('%-20s %-20.6f %-20.6f\n', unique_types{i}, posterior_max_values(i), posterior_min_values(i));
end

fprintf('\nAnalysis complete. Images saved.\n');

%% Calculate Test Points Posterior
if has_test_data
    fprintf('Calculating posterior probabilities for test points...\n');
    
    num_test_points = length(test_fs);
    posterior_probs = zeros(num_test_points, num_types);
    
    % Vectorized calculation for test points
    test_points = [test_fs, test_qE]; % N_test x 2
    
    likelihoods = zeros(num_test_points, num_types);
    
    for j = 1:num_types
        mu = mu_vectors{j};
        Sigma = cov_matrices{j};
        det_Sigma = det(Sigma);
        inv_Sigma = inv(Sigma);
        
        diff = test_points - mu';
        mahal_dist = sum((diff * inv_Sigma) .* diff, 2);
        
        likelihoods(:, j) = (1/(2*pi*sqrt(det_Sigma))) * exp(-0.5 * mahal_dist);
    end
    
    % Evidence
    weighted_likelihoods = likelihoods .* prior_probabilities';
    evidence = sum(weighted_likelihoods, 2);
    
    % Posterior
    valid_evidence = evidence > 0;
    posterior_probs(valid_evidence, :) = weighted_likelihoods(valid_evidence, :) ./ evidence(valid_evidence);
    posterior_probs(~valid_evidence, :) = 1 / num_types;
    
    if any(~valid_evidence)
        warning('%d test points had zero evidence and were assigned uniform probability.', sum(~valid_evidence));
    end
    
    % Save results
    fid = fopen('soil_classification_results.txt', 'w');
    fprintf(fid, 'Test Points Classification Results\n');
    fprintf(fid, '%-15s %-15s', 'fs (kPa)', 'qE (MPa)');
    for j = 1:num_types
        fprintf(fid, ' %-15s', ['P(' unique_types{j} ')']);
    end
    fprintf(fid, ' %-15s\n', 'Max Prob Soil');
    
    [max_probs, max_indices] = max(posterior_probs, [], 2);
    
    for i = 1:num_test_points
        fprintf(fid, '%-15.2f %-15.2f', test_fs(i), test_qE(i));
        for j = 1:num_types
            fprintf(fid, ' %-15.4f', posterior_probs(i, j));
        end
        fprintf(fid, ' %-15s\n', unique_types{max_indices(i)});
    end
    fclose(fid);
    fprintf('Results saved to soil_classification_results.txt\n');
    
    % Plot Test Points
    figure('Color', 'white', 'Position', [200, 200, 800, 600]);
    set(gcf, 'Name', 'Test Points on SBT Chart', 'NumberTitle', 'off');
    
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim([1, 1000]);
    ylim([0.1, 100]);
    
    hold on;
    line([1, 12], [1, 1], 'Color', 'k', 'LineWidth', lineWidth); 
    line([12, 12], [1, 0.1], 'Color', 'k', 'LineWidth', lineWidth); 
    plot(x_line3, y_line3, 'k-', 'LineWidth', lineWidth); 
    plot(x_line4, y_line4, 'k-', 'LineWidth', lineWidth); 
    plot(x_line5, y_line5, 'k-', 'LineWidth', lineWidth);
    
    colors = lines(num_types);
    scatter(test_fs, test_qE, 50, colors(max_indices, :), 'filled', 'MarkerEdgeColor', 'k');
    
    legend_entries = unique_types;
    % Create dummy plots for legend to ensure correct colors
    h_legend = zeros(num_types, 1);
    for i = 1:num_types
        h_legend(i) = plot(nan, nan, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
    end
    legend(h_legend, legend_entries, 'Location', 'best');
    
    xlabel('Sleeve Friction, f_s (kPa)', 'FontName', fontName, 'FontSize', labelSize);
    ylabel('Effective Corrected Cone Resistance, q_E (MPa)', 'FontName', fontName, 'FontSize', labelSize);
    title('Test Points Classification', 'FontName', fontName, 'FontSize', titleSize);
    
    grid on; box on;
    set(gca, 'LineWidth', 1.2, 'FontName', fontName, 'FontSize', fontSize);
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
    set(gca, 'TickDir', 'in');
    
    saveas(gcf, 'test_points_classification.png');
    fprintf('Classification chart saved to test_points_classification.png\n');
end
