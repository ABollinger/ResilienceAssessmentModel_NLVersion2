%% LOAD THE FLOOD RESULTS

clear

% load the results
load('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersionCh8/output/1000iterations/floodperformance.mat')
outputpath = ('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersionCh8/output/finalplots_1000iterations/');

%load('C:\Users\Andrew\Documents\TUD\AndrewBollinger\projects\MatlabResilienceModel\ThesisVersionCh8\output\floodperformance.mat')
%outputpath = ('C:\Users\Andrew\Documents\TUD\AndrewBollinger\projects\MatlabResilienceModel\ThesisVersionCh8\output\finalplots\');

%% PLOT THE FLOOD RESULTS - LINE PLOTS

figure(1)
z = figure(1);

% plot the default case results
xvalues = 1:numberofbusestoremove;
xvalues = xvalues - 1;
subplot(1,2,1);
plot(xvalues, networkperformancematrix00,'LineWidth',1,'Color',[0.7 0.7 0.7])
hold on
plot(xvalues, mean(networkperformancematrix00'),'LineWidth',3,'Color','black')
hold off
title({'Infrastructure performance vs. Event magnitude';'(no adaptation measures)'},'FontSize',16)
xlabel('Event magnitude (# substations failed)','FontSize',15)
ylabel('Performance (fraction demand served)','FontSize',15)
xlim([0 max(xvalues)])
ylim([min(min(networkperformancematrix00)) - 0.05 1.01])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',14)
set(gca,'XTickMode','manual','YTickMode','manual')
saveas(gcf,strcat(outputpath,'floodperformance_noadaptations.fig'));
saveas(gcf,strcat(outputpath,'floodperformance_noadaptations.png'));

% plot the combined results
subplot(1,2,2);
plot(xvalues, meanperformancematrix(:,1),...
    xvalues, meanperformancematrix(:,2),...
    xvalues, meanperformancematrix(:,3),...
    xvalues, meanperformancematrix(:,4),...
    xvalues, meanperformancematrix(:,5),...
    xvalues, meanperformancematrix(:,6),'LineWidth',3);
title({'Infrastructure performance vs. Event magnitude';'(comparison of measures)'},'FontSize',16)
xlabel('Event magnitude (# substations failed)','FontSize',15)
ylabel('Performance (fraction demand served)','FontSize',15)
leg = legend('0% protected','20% protected','40% protected','60% protected','80% protected','100% protected','Location','SouthWest');
set(leg,'FontSize',12)
v = get(leg,'title');
set(v,'string','Adaptation measure','FontSize',14);
xlim([0 max(xvalues)])
ylim([min(mean(networkperformancematrix00)) - 0.01 1.01])
colormap('Lines')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',14)
set(gca,'XTickMode','manual','YTickMode','manual')
saveas(gcf,strcat(outputpath,'floodperformance_adaptationscomparison.fig'));
saveas(gcf,strcat(outputpath,'floodperformance_adaptationscomparison.png'));

set(z,'PaperUnits','inches','PaperPosition',[0 0 14 6])
print -dpng floodresults.png;

min(min(networkperformancematrix00))
min(mean(networkperformancematrix00))
max(mean(networkperformancematrix00))
std(mean(networkperformancematrix00))

%% PLOT THE FLOOD RESULTS - SUMMARY BOX PLOT

% box plot of combined results
networkperformancematrix_total = [mean(networkperformancematrixB);...
    mean(networkperformancematrixC00);...
    mean(networkperformancematrixC15);...
    mean(networkperformancematrixC30);...
    mean(networkperformancematrixD00);...
    mean(networkperformancematrixD15);...
    mean(networkperformancematrixD30);...
    mean(networkperformancematrixW00);...
    mean(networkperformancematrixW15);...
    mean(networkperformancematrixW30);...
    mean(networkperformancematrixI00);...
    mean(networkperformancematrixI15);...
    mean(networkperformancematrixI30);...
    mean(networkperformancematrixE00);...
    mean(networkperformancematrixE15);...
    mean(networkperformancematrixE30)];
h = boxplot(networkperformancematrix_total','boxstyle','filled','symbol','k+',...
    'labels',{'Baseline ','Centralized 0%','Centralized 1.5% ','Centralized 3% ',...
    'Distributed 0% ','Distributed 1.5% ','Distributed 3% ',...
    'Offsh. wind 0%   ','Offsh. wind 1.5%   ','Offsh. wind 3%   ',...
    'Import 0% ','Import 1.5% ','Import 3% ',...
    'Export 0% ','Export 1.5% ','Export 3% '},'orientation','horizontal');
set(h,'color','black')
title({'Flood resilience of the infrastructure';'(comparison of scenarios)'},'FontSize',16)
ylabel({'Generation / demand development scenario'},'FontSize',15)
xlabel('Resilience','FontSize',15)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',14)
b = get(gca,'YTickLabel');
set(gca,'YTickLabel',b,'FontSize',14)
set(gca,'Ydir','reverse')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 7 9]);
saveas(gcf,strcat(outputpath,'floodperformance_scenarioscomparison_box.fig'));
saveas(gcf,strcat(outputpath,'floodperformance_scenarioscomparison_box.png'));


%% UNUSED FLOOD RESULTS PLOTS 

% create a histogram and kernel density plot
hist(mean(networkperformancematrix00),20)
set(get(gca,'child'),'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
hold on
ksdensity(mean(networkperformancematrix00))
set(get(gca,'child'),'Color','black','LineWidth',2);
hold off

% create a bar plot with error bars
addpath('C:\Users\Andrew\Documents\TUD\AndrewBollinger\projects\MatlabResilienceModel\ThesisVersion')

%   Symmetric Example:
y = [mean(mean(networkperformancematrix00));...
    mean(mean(networkperformancematrix20));...
    mean(mean(networkperformancematrix40));...
    mean(mean(networkperformancematrix60));...
    mean(mean(networkperformancematrix80));...
    mean(mean(networkperformancematrix100))];
errY = zeros(6,1,2);
errY(:,:,1) = [mean(mean(networkperformancematrix00)) - min(mean(networkperformancematrix00));...
    mean(mean(networkperformancematrix20)) - min(mean(networkperformancematrix20));...
    mean(mean(networkperformancematrix40)) - min(mean(networkperformancematrix40));...
    mean(mean(networkperformancematrix60)) - min(mean(networkperformancematrix60));...
    mean(mean(networkperformancematrix80)) - min(mean(networkperformancematrix80));...
    mean(mean(networkperformancematrix100)) - min(mean(networkperformancematrix100))];
errY(:,:,2) = [max(mean(networkperformancematrix00)) - mean(mean(networkperformancematrix00));...
    max(mean(networkperformancematrix20)) - mean(mean(networkperformancematrix20));...
    max(mean(networkperformancematrix40)) - mean(mean(networkperformancematrix40));...
    max(mean(networkperformancematrix60)) - mean(mean(networkperformancematrix60));...
    max(mean(networkperformancematrix80)) - mean(mean(networkperformancematrix80));...
    max(mean(networkperformancematrix100)) - mean(mean(networkperformancematrix100))];
h = barwitherr(errY, y);% Plot with errorbars
set(gca,'XTickLabel',{'0%','20%','40%','60%','80%','100%'})
ylim([min(mean(networkperformancematrix00)) - 0.1 1.01])
set(h(1),'FaceColor',[0.7 0.7 0.7]);
title('Flood resilience of the infrastructure (comparison of measures)','FontSize',17)
xlabel('Adaptation measure (% vulnerable substations protected)','FontSize',14)
ylabel('Resilience','FontSize',14)

%% LOAD THE HEAT WAVE RESULTS

clear

% load the results
load('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersionCh8/output/1000iterations/heatwaveperformance.mat')
outputpath = ('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersionCh8/output/finalplots_1000iterations/');

%load('C:\Users\Andrew\Documents\TUD\AndrewBollinger\projects\MatlabResilienceModel\ThesisVersionCh8\output\heatwaveperformance_repaired.mat')
%outputpath = ('C:\Users\Andrew\Documents\TUD\AndrewBollinger\projects\MatlabResilienceModel\ThesisVersionCh8\output\finalplots\');

%% PLOT THE HEAT WAVE RESULTS - LINE PLOTS

figure(1)
z = figure(1);
axislabelfontsize = 10;
titlelabelfontsize = 12;

% plot the default case results
xincrement = capacitytosubtracteachiteration / 1000;
xvalues = xincrement:xincrement:length(meanperformance)*xincrement;
subplot(1,2,1);
plot(xvalues, networkperformancematrix00,'LineWidth',1,'Color',[0.7 0.7 0.7])
hold on
plot(xvalues, mean(networkperformancematrix00'),'LineWidth',3,'Color','black')
hold off
title({'Infrastructure performance vs. Event magnitude';'(no adaptation measures)'},'FontSize',16)
xlabel('Event magnitude (GW generation capacity disabled)','FontSize',15)
ylabel('Performance (fraction demand served)','FontSize',15)
xlim([0 max(xvalues)])
ylim([min(min(networkperformancematrix00)) - 0.01 1.01])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',14)
set(gca,'XTickMode','manual','YTickMode','manual')
saveas(gcf,strcat(outputpath,'heatwaveperformance_noadaptations.fig'));
saveas(gcf,strcat(outputpath,'heatwaveperformance_noadaptations.png'));

% plot the combined results
subplot(1,2,2);
plot(xvalues, meanperformancematrix(:,1),...
    xvalues, meanperformancematrix(:,2),...
    xvalues, meanperformancematrix(:,3),...
    xvalues, meanperformancematrix(:,4),...
    xvalues, meanperformancematrix(:,5),...
    xvalues, meanperformancematrix(:,6),'LineWidth',3);
title({'Infrastructure performance vs. Event magnitude';'(comparison of measures)'},'FontSize',16)
xlabel('Event magnitude (GW generation capacity disabled)','FontSize',15)
ylabel('Performance (fraction demand served)','FontSize',15)
leg = legend('0% reduction','5% reduction','10% reduction','15% reduction','20% reduction','25% reduction','Location','SouthWest');
set(leg,'FontSize',12)
v = get(leg,'title');
set(v,'string','Adaptation measure','FontSize',14);
xlim([0 max(xvalues)])
ylim([mean(min(networkperformancematrix00)) - 0.001 1.0003])
colormap('Lines')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',14)
set(gca,'XTickMode','manual','YTickMode','manual')
saveas(gcf,strcat(outputpath,'heatwaveperformance_adaptationscomparison.fig'));
saveas(gcf,strcat(outputpath,'heatwaveperformance_adaptationscomparison.png'));

set(z,'PaperUnits','inches','PaperPosition',[0 0 14 6])
print -dpng heatwaveresults.png;

min(min(networkperformancematrix00))
min(mean(networkperformancematrix00))
max(mean(networkperformancematrix00))
std(mean(networkperformancematrix00))

%% PLOT THE HEAT WAVE RESULTS - SUMMARY BOX PLOT

% box plot of combined results
networkperformancematrix_total = [mean(networkperformancematrixB);...
    mean(networkperformancematrixC00);...
    mean(networkperformancematrixC15);...
    mean(networkperformancematrixC30);...
    mean(networkperformancematrixD00);...
    mean(networkperformancematrixD15);...
    mean(networkperformancematrixD30);...
    mean(networkperformancematrixW00);...
    mean(networkperformancematrixW15);...
    mean(networkperformancematrixW30);...
    mean(networkperformancematrixI00);...
    mean(networkperformancematrixI15);...
    mean(networkperformancematrixI30);...
    mean(networkperformancematrixE00);...
    mean(networkperformancematrixE15);...
    mean(networkperformancematrixE30)];
h = boxplot(networkperformancematrix_total','boxstyle','filled','symbol','k+',...
    'labels',{'Baseline ','Centralized 0% ','Centralized 1.5% ','Centralized 3% ',...
    'Distributed 0% ','Distributed 1.5% ','Distributed 3% ',...
    'Offsh. wind 0%   ','Offsh. wind 1.5%   ','Offsh. wind 3%   ',...
    'Import 0% ','Import 1.5% ','Import 3% ',...
    'Export 0% ','Export 1.5% ','Export 3% '},'orientation','horizontal');
set(h,'color','black')
title({'Heat wave resilience of the infrastructure';'(comparison of scenarios)'},'FontSize',16)
ylabel({'Generation / demand development scenario'},'FontSize',15)
xlabel('Resilience','FontSize',15)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',14)
b = get(gca,'YTickLabel');
set(gca,'YTickLabel',b,'FontSize',14)
set(gca,'Ydir','reverse')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 7 9]);
saveas(gcf,strcat(outputpath,'heatwaveperformance_scenarioscomparison_box.fig'));
saveas(gcf,strcat(outputpath,'heatwaveperformance_scenarioscomparison_box.png'));

%% PLOT THE FLOOD RESULTS -- FREQUENCY VS MAGNITUDE (CURRENT SITUATION)

clear

% load the results
load('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/100iterations_flooddefensescenarios/floodperformance.mat')
outputpath = ('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/finalplots/');

% line plot
xvalues = 1:numberofbusestoremove;
xvalues = xvalues - 1;
semilogy(xvalues,mean(eventprobabilitymatrix00a'),xvalues,mean(eventprobabilitymatrix00b'),xvalues,mean(eventprobabilitymatrix00c'),xvalues,mean(eventprobabilitymatrix00d'))
xlabel('Event frequency (# substations failed)')
ylabel('Event probability')

%% PLOT THE FLOOD RESULTS -- FRQUENCY-ADJUSTED PERFORMANCE VS MAGNITUDE (CURRENT SITUATION)

clear

% load the results
load('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/100iterations_flooddefensescenarios/floodperformance.mat')
outputpath = ('/home/andrewbollinger/AndrewBollinger/projects/MatlabResilienceModel/ThesisVersion/output/finalplots/');

xvalues = 1:numberofbusestoremove;
xvalues = xvalues - 1;
%xvalues(1) = [];
meanfrequencyadjustedperformancematrix_cropped = meanfrequencyadjustedperformancematrix;
meanfrequencyadjustedperformancematrix_cropped(1,:) = [];

plot(xvalues, meanfrequencyadjustedperformancematrix(:,1))

xlabel('Event magnitude (# substations failed)')
ylabel('Frequency-adjusted performance')


%% PLOT THE FLOOD RESULTS -- FRQUENCY VS MAGNITUDE (FUTURE SCENARIOS)



%% PLOT THE FLOOD RESULTS -- FREQUENCY-ADJUSTED PERFORMANCE VS. MAGNITUDE (FUTURE SCENARIOS)



%%