
%% TODO

% See the following links for data to include climate change effects:
% Overstromingsrisico's in Nederland in een veranderend klimaat: http://www.rivm.nl/bibliotheek/digitaaldepot/WL_rapport_Overstromingsrisicos_Nederland.pdf
% The effects of climate change in the Netherlands 2012 (PBL): http://www.pbl.nl/sites/default/files/cms/publicaties/PBL_2013_The%20effects%20of%20climate%20change%20in%20the%20Netherlands_957.pdf
% The effects of climate change in the Netherlands: http://www.rivm.nl/bibliotheek/rapporten/773001037.pdf
% New cooling water discharge guidelines in the Netherlands: http://www.aquator.nl/files/CIW%20nwe%20systematiek.pdf
% Effecten van klimaatverandering op de waterkwaliteit in de Rijn en Maas: http://www.delftcluster.nl/website/files/Waterkwaliteit_/2008-U-R0629_Michelle_van_Vliet_definitieve_versie.pdf
% Aard, ernst en omvang van watertekorten in Nederland (pg. 91), by RIZA

% Add comments to generation dispatch code

% OPTIONAL:

% Assumption: We are ignoring the flood risk of substations in the
% buitendijksgebied I think.

% Incorporate code for calculating movement through the state space

% Calculate resilience values, also max and min values, for each case, both
% with the risk factors and without the risk factors.

% Validate the model by comparing a default case with line capacities and
% maybe with load flow data.  If there's not time for this, you can leave
% this blank and do it later.

% Include power plant maintenance scenarios based on planned unavailability
% data?

% Fill in missing values in distribution grid demand data based on the
% known total peak consumption in 2010?

% Include sudden failures due to line sag/soil movement?

%% SET GLOBAL VARIABLES

clear

OS = 1; % Linux
%OS = 2; % Windows

%add paths
if OS == 1
    addpath('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersionCh8');
    addpath(genpath('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersionCh8/matpower4.0b4/'));
    addpath('/home/andrewbollinger/AndrewBollinger/projects/NetlogoModels/NetworkEvolution/modelFiles/outputdata/v10/n2');
    outputpath = '/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersionCh8/output/1000iterations/';
end

if OS == 2
    addpath('C:\\Users\\Andrew\\Documents\\TUD\\AndrewBollinger\\projects\\MatlabResilienceModel\\ThesisVersionCh8');
    addpath(genpath('C:\\Users\\Andrew\\Documents\\TUD\\AndrewBollinger\\projects\\MatlabResilienceModel\\ThesisVersionCh8\\matpower4.0b4'));
    addpath('C:\\Users\\Andrew\\Documents\\TUD\\AndrewBollinger\\projects\\NetlogoModels\\NetworkEvolution\\modelFiles\\outputdata\\v10\\n2')
    outputpath = 'C:\\Users\\Andrew\\Documents\\TUD\\AndrewBollinger\\projects\\MatlabResilienceModel\\ThesisVersionCh8\\output\\';
end

%set options
printstuff = false;
allowcascadingfailures = true;
allowgeneratorredispatching = true;
redispatchattemptslimit = 100;
mpopt = mpoption('PF_DC', 1, 'OUT_ALL', 0, 'VERBOSE', 0);
iterations = 1000;
plotstuff = 1;

%% EVALUATE FLOOD RESILIENCE

% create a matrix to hold the results
meanperformancematrix = [];
generationdemandscenarios = {'scenario_baseline',...
    'scenario_centralized_0','scenario_centralized_1.5','scenario_centralized_3',...
    'scenario_distributed_0','scenario_distributed_1.5','scenario_distributed_3',...
    'scenario_offshorewind_0','scenario_offshorewind_1.5','scenario_offshorewind_3',...
    'scenario_import_0','scenario_import_1.5','scenario_import_3',...
    'scenario_export_0','scenario_export_1.5','scenario_export_3'};
    
for generationdemandscenario = generationdemandscenarios

    % create a matrix to hold the results
    LoadNetwork
    networkperformancematrix = [];

    i = 1;
    while i <= iterations

        % load the network data
        LoadNetwork

        % create some variables to hold the results
        lostloadfraction = 0;
        networkperformance = [];
        disp(strcat('paramvalues=',char(generationdemandscenario),' iteration=',num2str(i)));

        % set the bus demand values
        mindemandvector = mpc.bus(:,16);
        maxdemandvector = mpc.bus(:,15);
        demandvector = mindemandvector + (maxdemandvector - mindemandvector) * rand(1);
        mpc.bus(:,3) = demandvector;

        % create a random generation dispatch distribution
        generatoroutputs = random('uniform', 0, 1, length(mpc.gen(:,1)), 1);
        generatoroutputs = generatoroutputs .* mpc.gen(:,9);
        generatoroutputs = generatoroutputs * sum(mpc.bus(:,3)) / sum(generatoroutputs);
        mpc.gen(:,2) = generatoroutputs;

        % create the weight matrix, which will be used to determine the removal of substations
        % The imported vulnerability values vary between 0 and 5, which relate
        % to different maximum flood height levels, as follows: 0 = 0m, 1 =
        % 0-0.2m, 2 = 0.2-0.5m, 3 = 0.5 - 0.8m, 4 = 0.8 - 2m, 5 = 2 - 5m. We
        % assume that substations are protected up to 2m (the supposed actual
        % value is 2.5m), so all values less than 3 are considered to constitute
        % a risk of 0.
        busfloodvulnerabilities = mpc.bus(:,14);
        busfloodvulnerabilities(busfloodvulnerabilities == 1) = 0;
        busfloodvulnerabilities(busfloodvulnerabilities == 2) = 0;
        busfloodvulnerabilities(busfloodvulnerabilities == 3) = 0;
        busfloodvulnerabilities(busfloodvulnerabilities == 4) = 0;
        busfloodvulnerabilities(busfloodvulnerabilities == 5) = 1;
        busfloodvulnerabilities_normalized = busfloodvulnerabilities / (5 - 4);

        mpc.bus(:,14) = busfloodvulnerabilities_normalized;

        % successively remove buses and calculate network performance
        numberofbusestoremove = length(find(mpc.bus(:,14)));
        %numberofbusestoremove = length(mpc.bus(:,1));
        busestoremove = mpc.bus;
        for numberofbusesremoved = 1:numberofbusestoremove

            % get the number of branches remaining
            numrows = length(mpc.branch(:,1));

            redispatchattempts = 0;

            % calculate the power flows
            CalculatePowerFlows;

            % calculate the network performance
            lostload = sum(mpc.bus(:,3)) - sum(busresults2(:,2));
            lostloadfraction = lostload ./ sum(mpc.bus(:,3));
            networkperformance = vertcat(networkperformance, 1 - lostloadfraction);

            % determine the bus to be removed based on the vulnerability and remove the bus
            if sum(busestoremove(:,14)) > 0
                x = randsample(1:length(busestoremove(:,1)),1,true,busestoremove(:,14));
                y = busestoremove(x,1);
                mpc.branch((mpc.branch(:,1) == y),:) = [];
                mpc.branch((mpc.branch(:,2) == y),:) = [];
                mpc.gen((mpc.gen(:,1) == y),9) = 0;
                mpc.gen((mpc.gen(:,1) == y),2) = 0;
                busestoremove(x,:) = [];
            end
        end

        % plot the network performance (optional)
        if plotstuff == 1
            plot(networkperformance,'LineWidth',1,'Color',[0.7 0.7 0.7]);
            hold on
        end

        % update the network performance matrix
        networkperformancematrix = horzcat(networkperformancematrix, networkperformance);

        % update the iterator
        i = i + 1;

    end

    SaveTheMatrices;

    % calculate and plot the mean network performance
    networkperformancematrix = networkperformancematrix';
    meanperformance = mean(networkperformancematrix);

    if plotstuff == 1
        plot(meanperformance, 'LineWidth',3,'Color','black');
        xlabel('Event magnitude (no .substations failed)')
        ylabel('Network performance')
        ylim([0 1.05])
        hold off
    end

    meanperformancematrix = vertcat(meanperformancematrix,meanperformance);

end

% calculate the resilience
resiliencevector1 = mean(meanperformancematrix');
resiliencevector1 = resiliencevector1';

if plotstuff == 1
    
    % plot mean network performance across the adaptation scenarios
    meanperformancematrix = meanperformancematrix';
    xvalues = 1:numberofbusestoremove;
    plot(xvalues, meanperformancematrix(:,1),xvalues, meanperformancematrix(:,2),...
        xvalues, meanperformancematrix(:,3),xvalues, meanperformancematrix(:,4),...
        xvalues, meanperformancematrix(:,5),xvalues, meanperformancematrix(:,6),...
        xvalues, meanperformancematrix(:,7),xvalues, meanperformancematrix(:,8),...
        xvalues, meanperformancematrix(:,9),xvalues, meanperformancematrix(:,10),...
        xvalues, meanperformancematrix(:,11),xvalues, meanperformancematrix(:,12),...
        xvalues, meanperformancematrix(:,13),xvalues, meanperformancematrix(:,14),...
        xvalues, meanperformancematrix(:,15),xvalues, meanperformancematrix(:,16));
    legend('scenariobaseline',...
        'scenariocentralized0','scenariocentralized1.5','scenariocentralized3',...
        'scenariodistributed0','scenariodistributed1.5','scenariodistributed3',...
        'scenariooffshorewind0','scenariooffshorewind1.5','scenariooffshorewind3',...
        'scenarioimport0','scenarioimport1.5','scenarioimport3',...
        'scenarioexport0','scenarioexport1.5','scenarioexport3','Location','NorthEast')
    xlabel('Event magnitude (# substations failed)')
    ylabel('Network performance')
    saveas(gcf,strcat(outputpath,'floodperformance.fig'));
    saveas(gcf,strcat(outputpath,'floodperformance.png'));
    hold off

end

save(strcat(outputpath,'floodperformance'));


%% EVALUATE HEAT WAVE RESILIENCE

capacitytosubtracteachiteration = 100;

% create a matrix to hold the results
meanperformancematrix = [];

generationdemandscenarios = {...
    'scenario_baseline',...
    'scenario_centralized_0','scenario_centralized_1.5','scenario_centralized_3',...
    'scenario_distributed_0','scenario_distributed_1.5','scenario_distributed_3',...
    'scenario_offshorewind_0','scenario_offshorewind_1.5','scenario_offshorewind_3',...
    'scenario_import_0','scenario_import_1.5','scenario_import_3',...
    'scenario_export_0','scenario_export_1.5','scenario_export_3'};

%find the maximum generatorsensitivecapacities in the tested scenarios
for generationdemandscenario = generationdemandscenarios
    
    LoadNetwork
    maxsum = 0;
    generatorssensitivecapacities = mpc.gen(:,24);
    if sum(generatorssensitivecapacities) > maxsum
        maxsum = sum(generatorssensitivecapacities);
    end
end
sizeofnetworkperformancematrix = ceil(maxsum / capacitytosubtracteachiteration);

for generationdemandscenario = generationdemandscenarios
    
    % create a matrix to hold the results
    LoadNetwork
    networkperformancematrix = zeros(sizeofnetworkperformancematrix, iterations);
    %networkperformancematrix = [];

    i = 1;
    while i <= iterations

        % load the network data
        LoadNetwork

        % create some variables to hold the results
        lostloadfraction = 0;
        networkperformance = [];
        disp(strcat('paramvalue=',char(generationdemandscenario),' iteration=',num2str(i)));

        % set the bus demand values
        mindemandvector = mpc.bus(:,16);
        maxdemandvector = mpc.bus(:,15);
        %adjust demand vector to account for lower demand during summer
        %14483 is the maximum May-Sep 2013 demand
        %16101 is the maximum Apr-Dec 2013 demand
        maxdemandvector = maxdemandvector * 14483 / 16101;
        demandvector = mindemandvector + (maxdemandvector - mindemandvector) * rand(1);
        %demandvector = mindemandvector + (maxdemandvector - mindemandvector); %TEST
        mpc.bus(:,3) = demandvector;
        
        % modify line capacities depending on the temperature
        Ta = 33.2; %ambient temp
        Taref = 25; %reference temp (rating temp)
        Tcmax = 100; %maximum conductor temp
        capacitylist(:,3) = capacitylist(:,3) * sqrt((Tcmax - Ta) / (Tcmax - Taref));

        % modify demand magnitudes depending on the temperature
        Ta = 33.2; %ambient temp
        Taref = 20.15; %reference temp = mean maximum daily temp at De Bilt between 1/5/2010 and 31/09/2010
        mpc.bus(:,3) = mpc.bus(:,3) * (1 + 0.005 * (Ta - Taref));

        % create a random generation dispatch distribution
        generatoroutputs = random('uniform', 0, 1, length(mpc.gen(:,1)), 1);
        generatoroutputs = generatoroutputs .* mpc.gen(:,9);
        %generatoroutputs = mpc.gen(:,9); %TEST
        generatoroutputs = generatoroutputs * sum(mpc.bus(:,3)) / sum(generatoroutputs);
        mpc.gen(:,2) = generatoroutputs;

        % successively reduce generation capacity
        generatorssensitivecapacities = mpc.gen(:,24);
        %while sum(generatorssensitivecapacities) > 0
        for removesensitivecapacity = 1:sizeofnetworkperformancematrix
        
            % get the number of branches
            numrows = length(mpc.branch(:,1));
            
            redispatchattempts = 0;

            % calculate the power flows
            CalculatePowerFlows;

            % calculate the network performance
            lostload = sum(mpc.bus(:,3)) - sum(busresults2(:,2));
            lostloadfraction = lostload ./ sum(mpc.bus(:,3));
            networkperformance = vertcat(networkperformance, 1 - lostloadfraction);
            
            if sum(isnan(networkperformance)) > 0
                disp('NaNs found')
            end

            % determine the amount of capacity to subtract this iteration
            if sum(generatorssensitivecapacities) >= capacitytosubtracteachiteration
                amounttosubtract = capacitytosubtracteachiteration;
            else 
                amounttosubtract = sum(generatorssensitivecapacities); 
            end

            % get a random list of generator indices
            generatorindices = 1:length(mpc.gen(:,9));
            generatorindices = generatorindices';
            generatorsensitivitycapacityvector = horzcat(generatorindices, generatorssensitivecapacities);
            generatorsensitivitycapacityvector = horzcat(generatorsensitivitycapacityvector, mpc.gen(:,9));
            generatorsensitivitycapacityvector = generatorsensitivitycapacityvector(randperm(size(generatorsensitivitycapacityvector,1)),:);
            generatorssensitivecapacities = [generatorsensitivitycapacityvector(:,1) generatorsensitivitycapacityvector(:,2)];
            generatorstotalcapacities = [generatorsensitivitycapacityvector(:,1) generatorsensitivitycapacityvector(:,3)];
            
            % starting at the top of the random list, subtract capacity from
            % each generator until the desired amount has been subtracted
            % or until no more can be subtracted
            j = 1;
            while amounttosubtract > 0 && sum(generatorssensitivecapacities(:,2)) > 0 && j <= length(generatorstotalcapacities(:,1))
                if generatorssensitivecapacities(j,2) >= amounttosubtract
                    amounttosubtractfromthisgenerator = amounttosubtract;
                else
                    amounttosubtractfromthisgenerator = generatorssensitivecapacities(j,2);
                end
                generatorstotalcapacities(j,2) = generatorstotalcapacities(j,2) - amounttosubtractfromthisgenerator;
                generatorssensitivecapacities(j,2) = generatorssensitivecapacities(j,2) - amounttosubtractfromthisgenerator;
                %mpc.gen(generatorindices(j),9) = mpc.gen(generatorindices(j),9) - amounttosubtractfromthisgenerator;
                %generatorssensitivecapacities(j) = generatorssensitivecapacities(j) - amounttosubtractfromthisgenerator;
                amounttosubtract = amounttosubtract - amounttosubtractfromthisgenerator;
                j = j + 1;
            end
            generatorssensitivecapacities = sortrows(generatorssensitivecapacities,1);
            generatorssensitivecapacities = generatorssensitivecapacities(:,2);
            generatorstotalcapacities = sortrows(generatorstotalcapacities,1);
            mpc.gen(:,9) = generatorstotalcapacities(:,2);
            if min(mpc.gen(:,9)) < 0
                disp(min(mpc.gen(:,9)));
            end

            % reset the output of generators whose output exceeds their capacity
            generatorswithoutputexceedingcapacity_indices = find(mpc.gen(:,9) < mpc.gen(:,2),2);
            mpc.gen(generatorswithoutputexceedingcapacity_indices,2) = mpc.gen(generatorswithoutputexceedingcapacity_indices,9);

        end

        % plot the network performance (optional)
        if plotstuff == 1
            xincrement = capacitytosubtracteachiteration / 1000;
            xvalues = xincrement:xincrement:length(networkperformance)*xincrement;
            plot(xvalues,networkperformance,'LineWidth',1,'Color',[0.7 0.7 0.7]);
            hold on
        end

        % update the network performance matrix
        for z = 1:length(networkperformance(:,1))
            networkperformancematrix(z,i) = networkperformance(z,1);
        end
        %networkperformancematrix = horzcat(networkperformancematrix, networkperformance);

        % update the iterator
        i = i + 1;

    end
    
    SaveTheMatrices;

    % plot mean network performance
    networkperformancematrix = networkperformancematrix';
    meanperformance = mean(networkperformancematrix);
    xincrement = capacitytosubtracteachiteration / 1000;
    xvalues = xincrement:xincrement:length(meanperformance)*xincrement;

    plot(xvalues, meanperformance, 'LineWidth',3,'Color','black');
    xlabel('Event magnitude (GW generation capacity disabled)')
    ylabel('Network performance')
    ylim([0 1.05])

    hold off
    
    meanperformancematrix = vertcat(meanperformancematrix,meanperformance);

end

resiliencevector = mean(meanperformancematrix');
resiliencevector = resiliencevector';

% plot mean network performance across the adaptation scenarios
xincrement = capacitytosubtracteachiteration / 1000;
xvalues = xincrement:xincrement:length(meanperformance)*xincrement;
meanperformancematrix = meanperformancematrix';
plot(xvalues, meanperformancematrix(:,1),xvalues, meanperformancematrix(:,2),...
        xvalues, meanperformancematrix(:,3),xvalues, meanperformancematrix(:,4),...
        xvalues, meanperformancematrix(:,5),xvalues, meanperformancematrix(:,6),...
        xvalues, meanperformancematrix(:,7),xvalues, meanperformancematrix(:,8),...
        xvalues, meanperformancematrix(:,9),xvalues, meanperformancematrix(:,10),...
        xvalues, meanperformancematrix(:,11),xvalues, meanperformancematrix(:,12),...
        xvalues, meanperformancematrix(:,13),xvalues, meanperformancematrix(:,14),...
        xvalues, meanperformancematrix(:,15),xvalues, meanperformancematrix(:,16));
    legend('scenariobaseline',...
        'scenariocentralized0','scenariocentralized1.5','scenariocentralized3',...
        'scenariodistributed0','scenariodistributed1.5','scenariodistributed3',...
        'scenariooffshorewind0','scenariooffshorewind1.5','scenariooffshorewind3',...
        'scenarioimport0','scenarioimport1.5','scenarioimport3',...
        'scenarioexport0','scenarioexport1.5','scenarioexport3','Location','NorthEast')
xlabel('Event magnitude (GW)')
ylabel('Network performance')
saveas(gcf,strcat(outputpath,'heatwaveperformance.fig'));
saveas(gcf,strcat(outputpath,'heatwaveperformance.png'));

save(strcat(outputpath,'heatwaveperformance'))

hold off


%% VALIDATE POWER FLOWS

clear

addpath('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion');
addpath(genpath('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/matpower4.0b4/'));
addpath('/home/andrewbollinger/AndrewBollinger/projects/NetlogoModels/NetworkEvolution/modelFiles/outputdata/');
outputpath = '/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/';

%set options
printstuff = false;
allowcascadingfailures = false;
mpopt = mpoption('PF_DC', 1, 'OUT_ALL', 0, 'VERBOSE', 0);
allowgeneratorredispatching = true;
redispatchattemptslimit = 100;


overloadedbranches = [];

for i = 1:100
    
    LoadNetwork

    % set the bus demand values
    mindemandvector = mpc.bus(:,17);
    maxdemandvector = mpc.bus(:,16);
    demandvector = mindemandvector + (maxdemandvector - mindemandvector) * rand(1);
    %demandvector = mindemandvector;
    mpc.bus(:,3) = demandvector;

    % create a random generation dispatch distribution
    generatoroutputs = random('uniform', 0, 1, length(mpc.gen(:,1)), 1);
    generatoroutputs = generatoroutputs .* mpc.gen(:,9);
    generatoroutputs = generatoroutputs * sum(mpc.bus(:,3)) / sum(generatoroutputs);
    mpc.gen(:,2) = generatoroutputs;

    % get the number of branches
    numrows = length(mpc.branch(:,1));
    
    redispatchattempts = 0;

    CalculatePowerFlows

    for s = 1:length(branchflows2(:,1))
        for t = 1:length(capacitylist(:,1))
            if branchflows2(s,1) == capacitylist(t,1) && branchflows2(s,2) == capacitylist(t,2) 
                branchflows2(s,4) = capacitylist(t,3);
            end
        end
    end

    capacityflowcomparison = 1.2 * branchflows2(:,4) - branchflows2(:,3);
    brancheswithnonzerocapacity_indices = find(branchflows2(:,4) > 0);
    overloadedbranches_indices = find(capacityflowcomparison < 0);
    overloadedbranches_adjusted = intersect(brancheswithnonzerocapacity_indices, overloadedbranches_indices);
    overloadedbranches_buses = [branchflows2(overloadedbranches_adjusted,1), branchflows2(overloadedbranches_adjusted,2), ...
        branchflows2(overloadedbranches_adjusted,3), branchflows2(overloadedbranches_adjusted,4), ...
        branchflows2(overloadedbranches_adjusted,3)./branchflows2(overloadedbranches_adjusted,4)];
    
    overloadedbranches = vertcat(overloadedbranches, overloadedbranches_buses);

end

overloadedbranches = unique(overloadedbranches, 'rows');
disp(overloadedbranches)

%bar(capacityflowcomparison)


%% CALCULATE CAPACITY UTILIZATION AND THE RATIO OF DEMAND TO SUPPLY UNDER NO-FAULT CONDITIONS

clear

OS = 1; % Linux
%OS = 2; % Windows

%add paths
if OS == 1
    addpath('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersionCh8');
    addpath(genpath('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersionCh8/matpower4.0b4/'));
    addpath('/home/andrewbollinger/AndrewBollinger/projects/NetlogoModels/NetworkEvolution/modelFiles/outputdata/v10/n2');
    outputpath = '/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersionCh8/output/';
end

if OS == 2
    addpath('C:\\Users\\Andrew\\Documents\\TUD\\AndrewBollinger\\projects\\MatlabResilienceModel\\ThesisVersionCh8');
    addpath(genpath('C:\\Users\\Andrew\\Documents\\TUD\\AndrewBollinger\\projects\\MatlabResilienceModel\\ThesisVersionCh8\\matpower4.0b4'));
    addpath('C:\\Users\\Andrew\\Documents\\TUD\\AndrewBollinger\\projects\\NetlogoModels\\NetworkEvolution\\modelFiles\\outputdata\\v10\\n2')
    outputpath = 'C:\\Users\\Andrew\\Documents\\TUD\\AndrewBollinger\\projects\\MatlabResilienceModel\\ThesisVersionCh8\\output\\';
end

%set options
printstuff = false;
allowcascadingfailures = true;
allowgeneratorredispatching = true;
redispatchattemptslimit = 100;
mpopt = mpoption('PF_DC', 1, 'OUT_ALL', 0, 'VERBOSE', 0);

generationdemandscenarios = {'scenario_baseline',...
    'scenario_centralized_0','scenario_centralized_1.5','scenario_centralized_3',...
    'scenario_distributed_0','scenario_distributed_1.5','scenario_distributed_3',...
    'scenario_offshorewind_0','scenario_offshorewind_1.5','scenario_offshorewind_3',...
    'scenario_import_0','scenario_import_1.5','scenario_import_3',...
    'scenario_export_0','scenario_export_1.5','scenario_export_3'};
    
for generationdemandscenario = generationdemandscenarios

    % create a matrix to hold the results
    LoadNetwork
    meancapacityutilization_total = [];
    demandsupplyratio_total = [];

    i = 1;
    while i <= 50
        
        % load the network data
        LoadNetwork

        % create some variables to hold the results
        lostloadfraction = 0;
        networkperformance = [];
        %disp(strcat('paramvalues=',char(generationdemandscenario),' iteration=',num2str(i)));

        % set the bus demand values
        mindemandvector = mpc.bus(:,16);
        maxdemandvector = mpc.bus(:,15);
        demandvector = mindemandvector + (maxdemandvector - mindemandvector) * rand(1);
        mpc.bus(:,3) = demandvector;

        % create a random generation dispatch distribution
        generatoroutputs = random('uniform', 0, 1, length(mpc.gen(:,1)), 1);
        generatoroutputs = generatoroutputs .* mpc.gen(:,9);
        generatoroutputs = generatoroutputs * sum(mpc.bus(:,3)) / sum(generatoroutputs);
        mpc.gen(:,2) = generatoroutputs;
        
        demandsupplyratio = sum(mpc.bus(:,3)) / sum(mpc.gen(:,9));
        
        numrows = length(mpc.branch(:,1));
        redispatchattempts = 0;
        CalculatePowerFlows;
        
        capacityutilization = branchflows2(:,3) ./ branchflows2(:,4);
        capacityutilization = capacityutilization(branchflows2(:,4) > 0,:); %ignore the transformers
        meancapacityutilization = mean(capacityutilization);
        meancapacityutilization_total = vertcat(meancapacityutilization_total, meancapacityutilization);
        demandsupplyratio_total = vertcat(demandsupplyratio_total, demandsupplyratio);
        
        i = i + 1;
    end
    
    meancapacityutilization_mean = mean(meancapacityutilization_total);
    demandsupplyratio_mean = mean(demandsupplyratio_total);
    
    disp(' ');
    disp(char(generationdemandscenario));
    disp(strcat('mean capacity utilization = ', num2str(meancapacityutilization_mean)))
    disp(strcat('mean demand supply ratio = ', num2str(demandsupplyratio_mean)))
end
        

%% CHECK CONSISTENCY OF LOADED NETWORK

LoadNetwork;

linkvoltages1 = mpc.branch(:,15);
linkcapacities1 = capacitylist(:,3);
linkvalues1 = horzcat(linkvoltages1, linkcapacities1);
linkvalues1 = sortrows(linkvalues1,[1 2]);

generatorcapacities1 = mpc.gen(:,9);
generatorcapacities1 = sortrows(generatorcapacities1);

busdemandvalues1 = mpc.bus(:,3);
busdemandvalues1 = sortrows(busdemandvalues1);

%%

LoadNetwork;

linkvoltages2 = mpc.branch(:,15);
linkcapacities2 = capacitylist(:,3);
linkvalues2 = horzcat(linkvoltages2, linkcapacities2);
linkvalues2 = sortrows(linkvalues2,[1 2]);

generatorcapacities2 = mpc.gen(:,9);
generatorcapacities2 = sortrows(generatorcapacities2);

busdemandvalues2 = mpc.bus(:,3);
busdemandvalues2 = sortrows(busdemandvalues2);

%%

linkdiff = linkvalues1 - linkvalues2;
gendiff = generatorcapacities1 - generatorcapacities2;
busdiff = busdemandvalues1 - busdemandvalues2;

maxlinkdiff = max(linkdiff);
maxgendiff = max(gendiff);
maxbusdiff = max(busdiff);

disp(maxlinkdiff)
disp(maxgendiff)
disp(maxbusdiff)

