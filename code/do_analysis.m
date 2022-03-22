decoding_done = false;
if decoding_done == false
  % add EEGlab to the path
  addpath('../eeglab2021.1/');

  % Start EEG path
  eeglab

  % Run the decoding analysis for experiment 1
  addpath('Experiment 1 data and analysis code/')
  cd 'Experiment 1 data and analysis code'/
  SVM_ECOC_ERP_Decoding_PreviousTrial_1()

  % Run the decoding analysis
  addpath('../Experiment 2 data and analysis code/')
  cd '../Experiment 2 data and analysis code'/
  SVM_ECOC_ERP_Decoding_PreviousTrial_2()
  cd '../'
end
% Get a list of files for experiment 1
files = dir('Experiment 1 data and analysis code/OR*.mat');
files = arrayfun(@(x) x.name, files, 'UniformOutput', false);
experiment_1 = cell2mat(arrayfun(@(x) calc_means(x), files, 'UniformOutput',false));

load(files{1})
experiment_1_time = svmECOC.time;

files = dir('Experiment 2 data and analysis code/OR*.mat');
files = arrayfun(@(x) x.name, files, 'UniformOutput', false);
experiment_2 = cell2mat(arrayfun(@(x) calc_means(x), files, 'UniformOutput',false));

load(files{1})
experiment_2_time = svmECOC.time;


[exp1_mean, exp1_ci] = sliding_mean_ci(experiment_1);
[exp2_mean, exp2_ci] = sliding_mean_ci(experiment_2);

makeplot(exp1_mean, exp1_ci, experiment_1_time, 'Experiment 1')

makeplot(exp2_mean, exp2_ci, experiment_2_time, 'Experiment 2')

function mean_acc = calc_means(file)
  % load the file
  load(file{1});
  % calcualte the mean
  mean_acc = squeeze(mean(mean(mean(svmECOC.targets == svmECOC.modelPredict,1),3),4));
end

function [m, ci] = sliding_mean_ci(mean_acc)

  window_width = 5; % sliding window width
  ci_value = 0.95; % 95% confidence interval
  ci_width = norminv(1 - ((1 - ci_value) / 2));
  n = size(mean_acc,1); % sample size

  ci = ci_width * (movmean(std(mean_acc), window_width) / sqrt(n));
  m = movmean(mean(mean_acc),window_width);
end

function makeplot(exp_mean, exp_ci, exp_time, plot_title)

  t = tiledlayout(1,1,'Padding','tight');
  t.Units = 'inches';
  t.OuterPosition = [0.25 0.25 6 3];

  nexttile

  % Create axes
  %xes1 = axes('Parent',gcf);
  %hold(axes1,'on');

  upper_bound = exp_mean + exp_ci;
  lower_bound = exp_mean - exp_ci;
  time = [exp_time, fliplr(exp_time)];
  ci_bound = [upper_bound, fliplr(lower_bound)];
  fill(time, ci_bound, 'y', 'FaceAlpha',0.3,...
      'FaceColor',[1 0.411764705882353 0.16078431372549],...
      'EdgeColor','none');
  hold on;
  plot(exp_time, exp_mean, 'y', 'LineWidth', 2);
  set(gca,'XTick',[0 500 1000 1500],'YTick',[0.06 0.08]);
  % Create ylabel
  ylabel('Decoding accuracy');

  % Create xlabel
  xlabel('Time');

  % Create title
  title(plot_title);

  box on;
  grid on;


  exportgraphics(gcf, [strrep(plot_title, ' ','_'), '.png'], "Resolution",300)

  close all

end

