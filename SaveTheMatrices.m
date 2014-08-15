
% THIS CODE SAVES THE NECESSARY MATRICES TO THE CORRECT NAME
% MATRICES ARE NAMED ACCORDING TO THE SCENARIO

% save the networkperformancematrix
if strcmp(char(generationdemandscenario),'scenario_baseline')
    networkperformancematrixB = networkperformancematrix;
elseif strcmp(char(generationdemandscenario),'scenario_centralized_0')
    networkperformancematrixC00 = networkperformancematrix;
elseif strcmp(char(generationdemandscenario),'scenario_centralized_1.5')
    networkperformancematrixC15 = networkperformancematrix;
elseif strcmp(char(generationdemandscenario),'scenario_centralized_3')
    networkperformancematrixC30 = networkperformancematrix;
elseif strcmp(char(generationdemandscenario),'scenario_distributed_0')
    networkperformancematrixD00 = networkperformancematrix;
elseif strcmp(char(generationdemandscenario),'scenario_distributed_1.5')
    networkperformancematrixD15 = networkperformancematrix;
elseif strcmp(char(generationdemandscenario),'scenario_distributed_3')
    networkperformancematrixD30 = networkperformancematrix;
elseif strcmp(char(generationdemandscenario),'scenario_offshorewind_0')
    networkperformancematrixW00 = networkperformancematrix;
elseif strcmp(char(generationdemandscenario),'scenario_offshorewind_1.5')
    networkperformancematrixW15 = networkperformancematrix;
elseif strcmp(char(generationdemandscenario),'scenario_offshorewind_3')
    networkperformancematrixW30 = networkperformancematrix;
elseif strcmp(char(generationdemandscenario),'scenario_import_0')
    networkperformancematrixI00 = networkperformancematrix;
elseif strcmp(char(generationdemandscenario),'scenario_import_1.5')
    networkperformancematrixI15 = networkperformancematrix;
elseif strcmp(char(generationdemandscenario),'scenario_import_3')
    networkperformancematrixI30 = networkperformancematrix;
elseif strcmp(char(generationdemandscenario),'scenario_export_0')
    networkperformancematrixE00 = networkperformancematrix;
elseif strcmp(char(generationdemandscenario),'scenario_export_1.5')
    networkperformancematrixE15 = networkperformancematrix;
elseif strcmp(char(generationdemandscenario),'scenario_export_3')
    networkperformancematrixE30 = networkperformancematrix;
else
    disp('ERROR: CANNOT SAVE NETWORKPERFORMANCEMATRIX')
end