% Discrete Tool v.2
% Author: Jasmin Rizko
% Written: 4/17/2016

% This code aids in the study of the 1-D discrete dynamical system
% x_n+1=f(x_n; a), where a is a parameter. Time series, cobweb, and
% bifurcation diagrams appear simultaneously to allow easy exploration of
% parameter space.

% This code was adapted from Professor Darryl Yong's code (Discrete Tool v.1
% 11/20/2012), which was in turn adapted from "Three Views of the Logistic
% Map" by Hiroki Sayama.


% Make a dialog box to prompt user for inputs (for more information on the
% function inputdlg and its arguments, see MATLAB's documentation website).
function NewDiscreteTool2016
prompt = {'Enter discrete dynamical system, f(p)=', 'Enter population range start:', 'Enter population range end:',...
    'Enter number of steps:', 'Enter parameter range start:', 'Enter parameter range end:', 'Enter parameter starting value: '};
dlg_title = 'Inputs';
num_lines = 1;
defaultans = {'a*p*(1-p)', '0', '1', '20', '-2', '4', '1.0'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans, 'on');

% Extract user inputs from the variable called answer, in which the user
% inputs were stored
S = vectorize(char(answer{1})); % first, convert the input function from a string to a vector equation
C.f = str2func(['@(a, p) ' S]); % then, assign the vector equation to a function f capable of taking values of a and p
C.g = @(x) x;

% Then, extract the ...
C.startR = str2num(answer{2}); % interval start value
C.endR = str2num(answer{3}); % interval end value
C.seed = (C.startR + C.endR)/2;
C.steps = str2num(answer{4}); % number of steps to take in cobweb diagram
paramStart = str2num(answer{5}); % parameter range start value
paramEnd = str2num(answer{6}); % parameter range end value
C.a = str2num(answer{7}); % parameter start value

% And finally, set and calculate the ...
numReps = 200; % the incrementation of the parameter value for the bifurcation diagram
C.pInc = (C.endR-C.startR)/numReps; % increments of population values
paramInc = (paramEnd - paramStart)/numReps; % increment for bifurcation diagram
bmaxit = 400; % max iterates for bifurcation
bminit = 200; % min iterates for bifurcation


                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     % BIFURCATION DIAGRAM CODE %
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We will now make a 2-D array that will hold our f values as the parameter a changes, as
% well as the f values with each iteration. In this 2-D array, the columns
% will represent an iteration while the rows represent f values where the
% same value of a was used for each element in the row.
acounter = 1; % the counter will serve as our index when creating the bifurcation data array
p = (C.startR + C.endR)/2; % PICK A P-VALUE
bifData = zeros(length(paramStart:paramInc:paramEnd), bmaxit); % 2-D bifurcation data array
for iterate=1:bmaxit; % iterate through the values of a and the max number of iterates
    for aval=paramStart:paramInc:paramEnd;
         if iterate == 1; % if we're on the first iterate, we use the value for p defined above
             bifData(acounter, iterate) = C.f(aval, p); % store the value f(a, p) in the array
         else % if we're not on the first iterate, we use the f(a, p) value from the previous iterate in place of p
             bifData(acounter, iterate) = C.f(aval, bifData(acounter, iterate - 1));
         end
    acounter = acounter + 1; % increment acounter
    end
    acounter = 1; % reset the counter since we need to go to the top of the next column
end

% Before we plot the points, we need to take into account the settling time
% of the points. As a result of this settling time,
% the first however many iterates (chosen by the user as bminit) will be
% removed from bifData.
bifData = bifData(:, bminit:bmaxit);

% Now we can plot the points
figure(); % opens plot in a new window
hold on % will plot over what is currently on the graph; this enables us to plot all of the column vectors
aval=paramStart:paramInc:paramEnd; % x-axis
for iterate=1:length(bminit:bmaxit); % plots all columns vector in the array against x-axis
    plot(aval, bifData(:, iterate));
end

% aesthetics of the graph
xlabel('Parameter a');
ylabel('Population Function f(a, p)');
title('Bifurcation Diagram of f(a, p)');


                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     % COBWEB AND TIME SERIES CODE %
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
% Create an object C with the following properties
C.g1 = figure();
C.ax1 = axes('Parent',C.g1,'position',[.1 .4 .8 .5]);
C.g2 = figure();
C.x_axis = 1:C.steps; % on the x-axis are discrete time steps

% Aesthetics
plotall(C);
xlabel (C.ax1,'population at time "t"');
ylabel (C.ax1,'population at time "t+1"');
title(C.ax1,'Cobweb Diagram of f(a, p)');
set(groot,'CurrentFigure',C.g1);

% Make the sliders along with their labels
b = uicontrol('Parent', C.g1, 'Style', 'slider', 'Position', [81, 84, 419, 23], ...
    'value', C.a, 'min', paramStart, 'max', paramEnd, 'Callback', {@callback,C});
bgcolor = C.g1.Color;
bl1 = uicontrol('Parent', C.g1, 'Style', 'text', 'Position', [50, 84, 23, 23], ...
    'String', paramStart, 'BackgroundColor', bgcolor);
bl2 = uicontrol('Parent', C.g1, 'Style', 'text', 'Position', [500, 84, 23, 23], ...
    'String', paramEnd, 'BackgroundColor', bgcolor);
bl3 = uicontrol('Parent', C.g1, 'Style', 'text', 'Position', [240, 65, 100, 23], ...
    'String', 'Parameter a', 'BackgroundColor', bgcolor);

s = uicontrol('Parent', C.g1, 'Style', 'slider', 'Position', [81, 34, 419, 23], ...
    'value', C.seed, 'min', 0, 'max', 1, 'Callback', {@seed_callback,C});
sgcolor = C.g1.Color;
sl1 = uicontrol('Parent', C.g1, 'Style', 'text', 'Position', [50, 34, 23, 23], ...
    'String', 0, 'BackgroundColor', bgcolor);
sl2 = uicontrol('Parent', C.g1, 'Style', 'text', 'Position', [500, 34, 23, 23], ...
    'String', 1, 'BackgroundColor', bgcolor);
sl3 = uicontrol('Parent', C.g1, 'Style', 'text', 'Position', [240, 15, 100, 23], ...
    'String', 'Initial Population', 'BackgroundColor', bgcolor);
end

% Make a callback function for each slider
function [] = callback(varargin)
[h,C] = varargin{[1,3]};
C.a = get(h,'value');
plotall(C);
end

function [] = seed_callback(varargin)
[h,C] = varargin{[1,3]};
C.seed = get(h,'value');
plotall(C);
end

% Make a function to plot everything, which the callback functions will
% call in order to change what is viewed on the plots
function plotall(C)
set(groot,'CurrentFigure',C.g1); % first change the population plot
cla % clear what's currently on the plot

% Make population vector y
x = C.startR:C.pInc:C.endR; % population
counter = 1; % for indexing
y = zeros(1, length(C.startR:C.pInc:C.endR)); % initialize population vector
for i=C.startR:C.pInc:C.endR; % loop through x-values
    y(counter) = C.f(C.a, i); % compute
    counter = counter + 1; % increment counter
end
h1 =plot(C.ax1,x, y, x, x); % plot the population and y=x line

% Now, we want to display a cobweb from a value of p and see how it changes
% as a changes.
hold on % allows us to overlay plots on the same set of axes
time_series_x = zeros(1, C.steps); % hold the time series data x-values
time_series_y = zeros(1, C.steps); % hold the time series data y-values
x0 = C.seed; % initial x-value
y0 = 0; % initial y-value
for step=1:C.steps;
    % vertical half-step
    x1 = x0; % next point has same x-value
    y1 = C.f(C.a, x0); % but different y-value
    time_series_x(step) = x0;
    time_series_y(step) = y1;
    plot([x0, x0], [y0, y1], 'k'); % plot a line between the initial point and the next point
    % horizontal half-step
    x2 = y1; % next point lies on the line y=x, so both x and y values must be equal
    y2 = y1;
    plot([x1, x2], [y1, y2], 'k'); % plot a line between the second point and the third point
    % reset values of x0 and y0 for next iteration
    x0 = x2;
    y0 = y2;
end

% Now change the time series plot
set(groot,'CurrentFigure',C.g2);
h2 = scatter(C.x_axis, time_series_x, 'blue', 'filled'); % scatter plot
xlabel('Steps');
ylabel('f(a, p)');
title('Time Series');

end