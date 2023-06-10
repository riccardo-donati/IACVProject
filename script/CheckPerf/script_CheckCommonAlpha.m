%%
% In this file we check the performances of the models "Subset", "SubsetOnlyHF","Fundamental","Homography","FundamentalA"
% we take the best alpha for each model type and save a table of the
% results
%%
%% Check the subset constrained motion segmentation results


clear all;
close all;


max_NumHypoPerFrame = 500;
FrameGap = 1;
model_type = 'Subset';
model_types = ["Subset", "SubsetOnlyHF","Fundamental","Homography"];
best_gammas = [1e-2,1e-2,-1,-1];
gamma_range = [1e-2];
Alpha_Range = 5:15;


% Load Seq Information
temp = load('../../Data/SeqList.mat');
SeqList = temp.SeqList;

seq_range = 1:length(SeqList);

error_avg = [];
error_med = [];
valid_gamma = [];
experiments_errors = [];
model_alphas = containers.Map;
%%
for exp = 1:8
    for  type_i = 1:length(model_types)

    gamma = best_gammas(type_i);

    Para = [];
    Para.MaxItr = 8;   % max iteration
    %     Para.gamma = gamma;
    Para.gamma = gamma;
    %find best alpha
    maxAlpha = -1;
    minError = 100;
    %%
    for t=1:8
        result_path = fullfile('../../Results/MoSeg/Archivio/',int2str(t),'/',model_types(type_i));

        for Alpha = Alpha_Range
            if gamma == -1
                result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g.mat',...
                max_NumHypoPerFrame,Alpha));
            else
                result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g_gamma-%g.mat',...
                max_NumHypoPerFrame,Alpha,gamma));
            end

            if ~exist(result_filepath,'file')
                continue;
            end        

            temp = load(result_filepath);

            error = temp.error;
            if mean(error)< minError
                minError = mean(error);
                maxAlpha = Alpha;
            end
        end
        fprintf("current best alpha for model %s experiment %d is %d\n",model_types(type_i),t,maxAlpha)
        model_alphas(model_types(type_i)) = maxAlpha;
    end
    fprintf("----------------------------------------------\n")
    %% select the desired experiment

    result_path = fullfile('../../Results/MoSeg/Archivio/',int2str(exp),'/',model_types(type_i));
    %%
    %take best alpha
    if gamma == -1
        result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g.mat',...
        max_NumHypoPerFrame,maxAlpha));
    else
        result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g_gamma-%g.mat',...
        max_NumHypoPerFrame,maxAlpha,gamma));
    end
    temp = load(result_filepath);
    error = temp.error;

    if gamma ~= -1
      fprintf("Results for model %s  alpha: %d  gamma:%.2f\n",model_types(type_i),maxAlpha, gamma)
    else
      fprintf("Results for model %s  alpha: %d\n",model_types(type_i),maxAlpha)
    end


    for j = 1:length(error)
        fprintf('Miss Classification of seq %d:  %.2f%%  \n',j,error(j)*100)
        errors(type_i,j) = 100*error(j);
    end
    errors(type_i,23) = mean(error)*100;
    fprintf('Overall Miss Classification Rate =   %.2f%%\n',100*mean(error(end,:)));
    fprintf('----------------------------------------\n')
    error_avg(type_i) = mean(error);
    error_med(type_i) = median(error);

    end
    experiments_errors{exp} = errors;
end
%%
valid_gamma = [valid_gamma gamma];
%
format short g;
n1 = 1;
n2 = 25;

n3 = 2;
n4 = 3;
for exp=1:8
    T=array2table(experiments_errors{exp}');
    T.Properties.VariableNames = model_types;
    T.Properties.RowNames = [SeqList "avg"];
    T
    Talpha = array2table(values(model_alphas))
    Talpha.Properties.VariableNames = keys(model_alphas);
    Talpha.Properties.RowNames = "alpha";
    Talpha

    cell = strcat('A',int2str(n1),':E',int2str(n2));
    n1 = n1 + 26;
    n2 = n2 + 26;
    writetable(T,"tmp.xls",'WriteRowNames',true,'Sheet',1,'Range',cell)
    
    cell = strcat('G',int2str(n3),':K',int2str(n4));
    n3 = n3 + 26;
    n4 = n4 + 26;
    writetable(Talpha,"tmp.xls",'WriteRowNames',true,'Sheet',1,'Range',cell)
end
