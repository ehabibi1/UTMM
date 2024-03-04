project = '20220706beads'; rounds = [5];
%project = '20220615'; rounds = [1 2 3 4 5];
home_dir = '\\sodium\broad_thechenlab\ehsan\analysis\InSitu_preImpEmbryo\';


Names_batch=["WED1-1","Wed1-2",...
    "WED1-3","WED1-4","WED1-5","WED1-6","WED1-7",...
    "THU1-1","THU1-2","THU1-3",...
    "THU1-4","THU1-5","THU1-6","THU1-7", ...
    "FRI1-1","FRI1-2","FRI1-3",...
    "ICM","TE"];

Corresponding_stage=["2","2",...
    "4","4","4","4","4",...
    "PreComp8","PreComp8","PreComp8",...
    "Compact8","Compact8","Compact8","Compact8",...
    "Late","Late","Late",...
    "ICM","TE"];
% Create a container map
map = containers.Map(Names_batch, Corresponding_stage);
for round = rounds
%results_dir = sprintf('%s/projects/mouse/%s/mechanics/experiments/mechanical_measurements/round%d/results/',home_dir,project,round);
results_dir = sprintf('%s/projects/mouse/%s/round%d/',home_dir,project,round);

files = dir(fullfile(results_dir, 'XY*.mat'));
for i = 1:length(files) 
    file = files(i);
    load(sprintf('%s/%s', results_dir, file.name));
    
    
    % Extract the name batch from the filename
    [~, name_batch, ~] = fileparts(file.name);
    name_batch = erase(name_batch, 'XY_MSD ');
    
    % Check if the name batch exists in the map
    if isKey(map, name_batch)
        stage = map(name_batch);
    else
        stage = 'Unknown';  % Use this for name batches that are not in the map
    end
    
    
    % Find unique numbers in the first column.
    uniqueSeries = unique(S_c);
    % Find out how many times each number occurs.
    mmXc=[];
    mmYc=[];
    mmZc=[];
    mmMsd=[];
    seriesName={};
    for i=1:length(uniqueSeries),
    seriesName{end+1} = sprintf('S%03d',uniqueSeries(i));
    mmXc(end+1)=mean(X_c(S_c==uniqueSeries(i)));
    mmYc(end+1)=mean(Y_c(S_c==uniqueSeries(i)));
    mmZc(end+1)=mean(Z_c(S_c==uniqueSeries(i)));
    mmMsd(end+1)=mean(compliance(S_c==uniqueSeries(i)));
    end
    data = table(seriesName', mmXc', mmYc', mmZc', mmMsd', repmat({name_batch}, length(seriesName), 1), repmat({stage}, length(seriesName), 1));
    data.Properties.VariableNames = {'File_name'    'x' 'y' 'z' 'MSD' 'Name_batch' 'Stage'};
    data
    outputName = regexprep(file.name,{'XY_MSD\s+','.mat'},{'',''})

    writetable(data,sprintf('%s/%s.xlsx', results_dir, outputName),'Sheet',1,'WriteVariableNames',true);
end
end
