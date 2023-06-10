%%
% In this file we find the best lambdas for GRIC with 
% hyperparameter tuning of the parameters lambda1, lambda2 and lambda3
%%
%% Obtain silhouette graphs and tables for all clips, for all models
clear all;
close all;

max_NumHypoPerFrame = 500;
motions = ["motion1", "motion2","motion3","motion4", "motion5"];
n1 = 1;
n2 = 7;

%% Load Seq Information
temp = load('../../../Data/SeqList.mat');
SeqList = temp.SeqList;
FrameGap = 1;
addpath '../..\..\Tools\multigs\model_specific';
%% models silhouette
%Input theg desired geometric model to analyze
models = ["fundamental", "fundamentala", "fundamentalt", "homography","affine"];
scoresMotion = []; 
silhouetteMotion = [];
silhouetteMotionAvg = [];
silhouetteAvg = [];

ps = [7, 4, 2, 4, 3]; % also called p: minimum number of correspondences needed to estimate the constraint
ds = [3, 3, 3, 2, 2]; % dimension of the constraint:3 for fundamentals 2 for homography
ks = [7, 4, 2, 8, 6]; % number of parameters of the model:7 4 2 for fundamentals 8 for homography

% model selection algorithms
names = ["GRIC"];

gric = [];

lambdas = [];
if ~exist("BestSilhouette.mat", 'file')
    bestAvgSilhouette = -1;
else
    bestAvgSilhouette = load("BestSilhouette").bestAvgSilhouette; %is the avg silhouette of GRIC with default lambdas
end
if ~exist("BestLambdas.mat", 'file')
    bestLambdas = [];
else
    bestLambdas = load("BestLambdas").lambdas;
end
tic
for j=1:1000
    if ~isempty(bestLambdas)
        fprintf("-------------------------------------\n");
        fprintf("Current best silhouette="+num2str(bestAvgSilhouette)+"\n");
        fprintf("with lambda1 ="+num2str(bestLambdas(1))+" lambda2 = "+num2str(bestLambdas(2))+" lambda3 = "+num2str(bestLambdas(3))+"\n");
        fprintf("-------------------------------------\n");
    end
    fprintf("Iteration n."+int2str(j)+"\n");
    % hyperparameter tuning for GRIC
    lambda1 = rand*10;
    lambda2 = rand*10;
    lambda3 = rand*10;
    fprintf("Testing lambda1="+num2str(lambda1)+" lambda2="+num2str(lambda2)+" lambda3="+num2str(lambda3)+"\n");
    lambdas = [lambda1, lambda2, lambda3];
    for seq = 1:length(SeqList)
        seqData =  load(strcat("../../../Data/",SeqList{seq},'_Tracks.mat'));
        Data = seqData.Data;
        gt = seqData.Data.GtLabel;
        nMotions = max(gt);
      
        for f_i = 1:Data.nFrames-FrameGap
            for model=1:size(models,2)
                model_type = models(model);
                [fitModel, resModel, ~, ~, ~] = getModelParam(model_type);
                %% Select points for every motion visible on both frames
                visible_pts_ind = Data.visibleSparse(:,f_i) & Data.visibleSparse(:,f_i+1);
                Models = [];
                matchedPoints1 = [];
                matchedPoints2 = [];
                % obtain a binary array containing indices of point divided per
                % motions and for every motion estimate the F
                for m = 1:nMotions
                    pointsInMotion = visible_pts_ind & (gt==m);
                    % matched points of motion m visible in frames f_i and f_i+1
                    matchedPoints1{m} = Data.ySparse(:,pointsInMotion,f_i);
                    matchedPoints2{m} = Data.ySparse(:,pointsInMotion,f_i+FrameGap);
                    if size(matchedPoints1{m}, 2) > 1
                        matchedPoints1{m} = normalise2dpts(Data.ySparse(:,pointsInMotion,f_i));
                        matchedPoints2{m} = normalise2dpts(Data.ySparse(:,pointsInMotion,f_i+FrameGap));
                    end
                    try
                        Models{m} = feval(fitModel, [matchedPoints1{m}; matchedPoints2{m}]);
                    catch 
                        Models{m} = NaN;
                    end
                end     
    
                ReswrtMotion = [];
                variances = [];
                ReswrtMotionNormalized = [];
                robustRes = [];
                r = 4;
                d = ds(model); 
                k = ks(model);
                p = ps(model);

                
                for m=1:nMotions
                    if isnan(Models{m})
                        grictmp{f_i,model,m} = Inf;
                    else
                        ReswrtMotion{m} = feval(resModel, Models{m},[matchedPoints1{m}; matchedPoints2{m}]);
                        variances{m} = var(ReswrtMotion{m});
                        ReswrtMotionNormalized{m} = ReswrtMotion{m}.^2./variances{m};
                        n = size(ReswrtMotionNormalized{m},1);


                        tmp = lambda2*(r-d);
                        robustRes{m} = min(ReswrtMotionNormalized{m}, tmp);
                        grictmp{f_i,model,m} = sum(robustRes{m}) + lambda1*d*log(n) + lambda3*k;
        
                    end
                end
            end
        end
        gric{seq} = grictmp;

    end 
    % --------------------------------------------------------------------------------

    for seq=1:length(SeqList)
        seqData =  load(strcat("../../../Data/",SeqList{seq},'_Tracks.mat'));
        Data = seqData.Data;
        nMotions = max(Data.GtLabel);
        results = gric{seq};
        best_models = [];
        fr_list = [];
        for fr=1:Data.nFrames-FrameGap
            fr_list = [fr_list; "frame "+int2str(fr)+"-"+int2str(fr+1)];
            m = cell2mat(squeeze(results(fr,:,:)));
            [val,idx] = min(m);
            curr_best=[];
            for motion=1:nMotions
                best_model = models(idx(motion));
                curr_best=[curr_best best_model];
%                        fprintf(names(i)+": seq "+int2str(seq)+" frame "+int2str(fr)+" best model for motion "+int2str(motion)+" is "+best_model+"\n");
            end
            best_models=[best_models;curr_best];
        end

        seq_best_models{seq} = best_models;
    end

    
    %% Silhouette
    for seq = 1:length(SeqList)
        seqData =  load(strcat("../../../Data/",SeqList{seq},'_Tracks.mat'));
        Data = seqData.Data;
        gt = seqData.Data.GtLabel;
        nMotions = max(gt);
        silhouettes = [];
        curr_best_models = seq_best_models{ seq};
        sizes = [];
        for f_i = 1:Data.nFrames-FrameGap
            %% Select points for every motion visible on both frames
            visible_pts_ind = Data.visibleSparse(:,f_i) & Data.visibleSparse(:,f_i+1);
            Models = [];
            matchedPoints1 = [];
            matchedPoints2 = [];
            % obtain a binary array containing indices of point divided per
            % motions and for every motion estimate the F
            for m = 1:nMotions
                pointsInMotion = visible_pts_ind & (gt==m);
                % matched points of motion m visible in frames f_i and f_i+1
                matchedPoints1{m} = Data.ySparse(:,pointsInMotion,f_i);
                matchedPoints2{m} = Data.ySparse(:,pointsInMotion,f_i+FrameGap);
                if size(matchedPoints1{m}, 2) > 1
                    matchedPoints1{m} = normalise2dpts(Data.ySparse(:,pointsInMotion,f_i));
                    matchedPoints2{m} = normalise2dpts(Data.ySparse(:,pointsInMotion,f_i+FrameGap));
                end
                sizes = [sizes size(matchedPoints1{m},2)];

%                 NON NORMALIZED VERSION
%                 model_type = curr_best_models(f_i,m);
%                 [fitModel, resModel, ~, ~, ~] = getModelParam(model_type);
%                 try
% 
%                     Models{m} = feval(fitModel, [matchedPoints1{m}; matchedPoints2{m}]);
%                 catch 
%                     Models{m} = NaN;
%                 end
            end   
            
            % NORMALIZED VERSION
            sizeMin =max(8,min(sizes));
            for m=1:nMotions
                model_type = curr_best_models(f_i,m);
                [fitModel, resModel, ~, ~, ~] = getModelParam(model_type);
                try
                    indexes = randsample(size(matchedPoints1{m},2), sizeMin, true);
                    Models{m} = feval(fitModel, [matchedPoints1{m}(:,indexes); matchedPoints2{m}(:,indexes)]); %stochastic
                    %Models{m} = feval(fitModel, [matchedPoints1{m}(:,1:sizeMin); matchedPoints2{m}(:,1:sizeMin)]);
                catch 
                    Models{m} = NaN;
                end
            end

            % now compute residuals of every point wrt the model for every motion.
            % This is needed to compute a and b for every point
            % ReswrtMotion: matrix m x n: every row corresponds to the residuals of
            % a motion wrt all the m motions. e.g. 2 motions => 2x2:
            % [res(points_motion_1,model_1) res(points_motion_1,model_2);
            % res(points_motion_2,model_1) res(points_motion_2,model_2) ]
            ReswrtMotion = [];
            for m=1:nMotions
                for n=1:nMotions
                    model_type = curr_best_models(f_i,n);
                    [~, resModel, ~, ~, ~] = getModelParam(model_type);
                    if isnan(Models{n})
                         ReswrtMotion{m}{n} = Inf(size(matchedPoints1{m},2),1);
                    else
                        ReswrtMotion{m}{n} = feval(resModel, Models{n},[matchedPoints1{m}; matchedPoints2{m}]);
                    end
                end
            end

            %now compute a and b for every point and its silhouette score     
            for m=1:nMotions
                b = [];
                for n=1:nMotions
                    if m==n % vector of a's for motion m
                        a = ReswrtMotion{m}{n};                    
                    else %residuals of points of motion m wrt motion n
                        if isempty(b)
                            b = ReswrtMotion{m}{n};
                        else
                            b = min(b, ReswrtMotion{m}{n});
                        end
                    end
                end
                % we have a and b for every point of motion m. Now silhouette
                scoresMotion{seq}{f_i}{m} = (b - a) ./ max(a, b);
                silhouetteMotion{seq}{f_i}{m} = mean(scoresMotion{seq}{f_i}{m});
                silhouettes(m, f_i) = mean(scoresMotion{seq}{f_i}{m});
                frames(f_i) = "frames "+f_i+" "+(f_i+1);
            end
        end
        silhouetteMotionAvg{seq} = nanmean(silhouettes,2); 
        silhouetteMethodAvg{seq} = mean(silhouetteMotionAvg{seq});
    end
    
    meanSilhouettes = mean(cell2mat(silhouetteMethodAvg));
    fprintf("Silhouette = "+num2str(meanSilhouettes(1))+"\n");
    if meanSilhouettes(1) > bestAvgSilhouette 
        fprintf("NEW BEST GRIC AVGSILHOUETTE: %f\n",meanSilhouettes(1))
        bestAvgSilhouette = meanSilhouettes(1);
        save("BestSilhouette","bestAvgSilhouette");
        bestLambdas = lambdas;
        save("BestLambdas","lambdas");

    end
end
toc










