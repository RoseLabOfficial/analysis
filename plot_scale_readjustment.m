%% Y limits
ylims = [-10e-10, 10e-10];
yylims = [-4, 4];

%% Get Figures
gh1 = figure(3);
gh2 = figure(2);

%% set limits of y axis:
figure(gh1);
yyaxis left;
ylim(ylims);
yyaxis right;
ylim(yylims);

figure(gh2);
yyaxis left;
ylim(ylims);
yyaxis right;
ylim(yylims);
