%%
% In this file we check the performances of the models "Subset", "SubsetOnlyHF","Fundamental","Homography","FundamentalA"
% we take the best alpha for each model type on average on all the experiments and save a table of the
% results
%%
%% Find the best alphas 

clear all;
close all;


max_NumHypoPerFrame = 500;
model_types = ["Subset", "SubsetOnlyHF","Fundamental","Homography","FundamentalA","FundamentalT"];
best_gammas = [1e-2,1e-2,-1,-1,-1,-1];
Alpha_Range = 5:15;


% Load Seq Information
temp = load('../../Data/SeqList.mat');
SeqList = temp.SeqList;

%% Extract the errors for each model, for each clip, for each alpha
for  type_i = 1:length(model_types)

    gamma = best_gammas(type_i);

    Para = [];
    Para.MaxItr = 8;   % max iteration
    Para.gamma = gamma;
    maxAlpha = -1;
    minError = 100;
    %% Iterate over the experiments
    for t=1:8
        result_path = fullfile('../../Results/MoSeg/Archivio/',int2str(t),'/',model_types(type_i));
        %% Iterate over the alphas
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
            %% save the errors
            for j = 1:length(error)
                errors(t,type_i,Alpha - 4,j) = 100*error(j); %Alpha - 4 perchè inizia da 5 -> 1 = 5,..., 11 = 15
            end
        end
    end
end
%% 
avgErrors = mean(errors,1);
avgErrors = squeeze(avgErrors);
avgErrors = mean(avgErrors,3);

T=array2table(avgErrors);
T.Properties.VariableNames = ["5","6","7","8","9","10","11","12","13","14","15"];
T.Properties.RowNames = model_types;
T
[val,ind] = min(avgErrors');

ind = ind + 4; %Alpha start from 5
model_alphas = containers.Map(model_types,ind);

% save the best alphas
save("best_model_alphas.mat","model_alphas");
%% Check the subset constrained motion segmentation results

max_NumHypoPerFrame = 500;
FrameGap = 1;

Alpha_Range = 5:15;


% Load Seq Information
temp = load('../../Data/SeqList.mat');
SeqList = temp.SeqList;

seq_range = 1:length(SeqList);

errors = [];
error_avg = [];
error_med = [];
valid_gamma = [];
experiments_errors = [];
nexp = 8;
%% Find the errors in each experiments (8 in total)
for exp = 1:nexp
    for  type_i = 1:length(model_types)

    gamma = best_gammas(type_i);

    Para = [];
    Para.MaxItr = 8;   % max iteration
    Para.gamma = gamma;
    %% select the desired experiment

    result_path = fullfile('../../Results/MoSeg/Archivio/',int2str(exp),'/',model_types(type_i));
    %%
    %take best alpha
    if gamma == -1
        result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g.mat',...
        max_NumHypoPerFrame,model_alphas(model_types(type_i))));
    else
        result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g_gamma-%g.mat',...
        max_NumHypoPerFrame,model_alphas(model_types(type_i)),gamma));
    end
    temp = load(result_filepath);
    error = temp.error;

    if gamma ~= -1
      fprintf("Results for model %s  alpha: %d  gamma:%.2f\n",model_types(type_i),model_alphas(model_types(type_i)), gamma)
    else
      fprintf("Results for model %s  alpha: %d\n",model_types(type_i),model_alphas(model_types(type_i)))
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
% for each experiment visualize the table with the results and with the
% alphas
for exp=1:nexp
    T=array2table(experiments_errors{exp}');
    T.Properties.VariableNames = model_types;
    T.Properties.RowNames = [SeqList "avg"];
    T
    Talpha = array2table(values(model_alphas));
    Talpha.Properties.VariableNames = keys(model_alphas);
    Talpha.Properties.RowNames = "alpha";
    Talpha

    cell = strcat('A',int2str(n1),':G',int2str(n2));
    n1 = n1 + 26;
    n2 = n2 + 26;
    writetable(T,"tmp.xls",'WriteRowNames',true,'Sheet',1,'Range',cell)
    
    cell = strcat('I',int2str(n3),':O',int2str(n4));
    n3 = n3 + 26;
    n4 = n4 + 26;
    writetable(Talpha,"tmp.xls",'WriteRowNames',true,'Sheet',1,'Range',cell)
end

