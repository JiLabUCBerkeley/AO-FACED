n = 2^8;
% green_colormap = [linspace(0, 0, n)', linspace(0, 1, n)', linspace(0, 0, n)'];
% red_colormap = [linspace(0, 1, n)' zeros(n, 1) zeros(n, 1)];

% Define the red-hot colormap
redHotMap = [
    0, 0, 0;        % Black
    0.5, 0, 0;      % Dark Red
    1, 0, 0;        % Red
    1, 0.5, 0;      % Orange
    1, 1, 0;        % Yellow
    1, 1, 1         % White (optional, for brightness)
];
% Interpolate to create a smooth colormap
redHotMap = interp1(linspace(0, 1, size(redHotMap, 1)), redHotMap, linspace(0, 1, n));

%%
greenHotMap = [
    0, 0, 0;        % Black
    0, 0.5, 0;      % Dark Red
    0, 1, 0;        % Green
    % 0, 0.5, 0;      % Orange
    1, 1, 0;        % Yellow
    1, 1, 1         % White (optional, for brightness)
];
% Interpolate to create a smooth colormap
greenHotMap = interp1(linspace(0, 1, size(greenHotMap, 1)), greenHotMap, linspace(0, 1, n));

%% Define the yellow-hot colormap
yellowHotMap = [
    0, 0, 0;          % Black
    0.4, 0.4, 0;      % Dark Yellow (olive-like)
    0.7, 0.7, 0;      % Lime Yellow (greenish)
    1, 1, 0;          % Bright Yellow
    1, 0.85, 0;       % Golden Yellow (deeper)
    1, 1, 0.8;        % Light Yellow (softer)
    1, 1, 1;          % White (optional for brightness)
];

% Interpolate to create a smooth colormap
yellowHotMap = interp1(linspace(0, 1, size(yellowHotMap, 1)), yellowHotMap, linspace(0, 1, n),'linear');

%% cyan hot
cyanHotMap = [
    0, 0, 0;        % Black
    0, 0.5, 0.5;    % Dark Cyan
    0, 1, 1;        % Cyan
    1, 1, 1         % White
];
cyanHotMap = interp1(linspace(0, 1, size(cyanHotMap, 1)), cyanHotMap, linspace(0, 1, n),'linear');

%%
redmap = [linspace(0,1,n)', zeros(n,1), zeros(n,1)]; colormap(redmap);