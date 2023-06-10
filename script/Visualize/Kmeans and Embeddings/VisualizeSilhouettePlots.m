%%
% In this file we visualize the silhouette plots for every model of interest
% with the possibility of choosing a specific silhouette plot
% Moreover we also visualize a bargraph where we can view the results for
% all the sequences
%%
%% Check the subset constrained motion segmentation results
clear all;
close all;

max_NumHypoPerFrame = 500;
gamma_range = [1e-2];

%% Load Seq Information
temp = load('../../../Data/SeqList.mat');
SeqList = temp.SeqList;
FrameGap = 1;
addpath '..\..\Tools\multigs\model_specific';


%Input the desired model
model_types = ["Subset", "SubsetOnlyHF","Fundamental","Homography","FundamentalA","FundamentalT"];

%Load the best alpha for every model
model_alphas = load("../../CheckPerf/best_model_alphas.mat").model_alphas;


%% kmeans silhouette (compute the silhouette wrt the obtained clustering and not the Ground Truth)
%{
for i = 1:length(model_types)
    
    Alpha = model_alphas(model_types(i));
    title_figure ="kmeans "+ model_types(i)+" best_alpha="+Alpha;
    figure('Position', [50+(i*15), 50+(i*15), 1024, 512], 'Name', title_figure , 'NumberTitle','off');
    result_path = fullfile('../../../Results/MoSeg/',model_types(i));

    if model_types(i) == "Fundamental" || model_types(i) == "Homography" || model_types(i) == "FundamentalA" || model_types(i) == "FundamentalT"
        gamma = -1;
    else
        gamma = gamma_range(1);
    end

    %load the desired segmentation file
    if gamma == -1
       result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g.mat',...
       max_NumHypoPerFrame,Alpha));
    else
       result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g_gamma-%g.mat',...
       max_NumHypoPerFrame,Alpha,gamma));
    end
    
    temp = load(result_filepath);
    
    %extract observations and clustering resulting from kmeans
    ClusterIdx = temp.ClusterIdx; %output of kmeans
    observations = temp.observations;

    temps(i) = temp;
    
    for seq = 1:length(SeqList)   
        subplot(3,8, seq);

        [s,h] = silhouette(observations{seq}, ClusterIdx{seq}); %silhouette of kmeans
        xline(mean(s),'--r');
        title(SeqList{seq},'Interpreter','none');
    end

end
%}

%% Ground Truth Silhouette
for i = 1:length(model_types)
    
    Alpha = model_alphas(model_types(i));
    title_figure ="Ground Truth "+ model_types(i)+" best_alpha="+Alpha;
    figure('Position', [50+(i*15), 50+(i*15), 1024, 512], 'Name', title_figure , 'NumberTitle','off');
    result_path = fullfile('../../../Results/MoSeg/',model_types(i));

    % gamma is a parameter only for Subset and SubsetHF
    if model_types(i) == "Fundamental" || model_types(i) == "Homography" || model_types(i) == "FundamentalA" || model_types(i) == "FundamentalT"
        gamma = -1;
    else
        gamma = gamma_range(1);
    end

    %load the desired segmentation file
    if gamma == -1
       result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g.mat',...
       max_NumHypoPerFrame,Alpha));
    else
       result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g_gamma-%g.mat',...
       max_NumHypoPerFrame,Alpha,gamma));
    end
    
    temp = load(result_filepath);
    
    %extract observations  from kmeans
    observations = temp.observations;
    temps(i) = temp;
    
    for seq = 1:length(SeqList)
        %load the ground truth
        gt =  load(strcat("../../../Data/",SeqList{seq},'_Tracks.mat'));
        Gt = gt.Data.GtLabel; %ground truth
        subplot(3,8, seq);
        
        [s,h] = silhouette(observations{seq}, Gt); %silhouette of ground truth
        xline(mean(s),'--r');
        title(SeqList{seq},'Interpreter','none');
    end
    
end


%%
%Input the desired clip to analyze
map = containers.Map([1,2,3,4,5,6],model_types);
prompt = "What is desired the model to analyze [1 = Subset, 2 = SubsetHF, 3 = Fundamental, 4 = Homography, FundamentalA = 5, FundamentalT = 6] ?  ";
x = input(prompt);
model_type = map(x);
fprintf("%s selected\n",model_type);

%WSS
wss= [];
for i=1:22
    tmp = temps(x).allWSS{i}';
    while size(tmp) ~= 5
        tmp = [tmp 0];
    end
    wss = [wss; tmp];
end
figure(13);

bar(categorical(SeqList),wss,1);
set(gca,'TickLabelInterpreter','none')

title(strcat("WSS of the Sequences for model ",model_type));

%% plot single silhouette
prompt = "What is desired sequence to analyze (from 1 to 22) ? ";
for i=1:22
    fprintf("%d: %s\n",i,SeqList{i});
end
clip = input(prompt);
fprintf("%s selected\n",SeqList{clip});

%plot the silhouette
figure(14);
observations = temps(x).observations;
%ClusterIdx = temps(x).ClusterIdx; for the kmeans silhouette

gt =  load(strcat("../../../Data/",SeqList{clip},'_Tracks.mat'));
Gt = gt.Data.GtLabel; %for ground truth silhouette

[s,h] = silhouette(observations{clip}, Gt); 
xline(mean(s),'--r','LineWidth',2.0);

title(strcat(SeqList{clip}," ",model_type),'Interpreter','none');

set(gca,'FontSize',17)


