% Datacube animation example

% Initialise a random cube of data, which is made up of 100 lots of 60 wide
% by 40 tall images.

raw_cube = rand(40,60,100); % "height, width, depth" convention

% Import cube data into AODataCube class format
AOcube = AODataCube(raw_cube);

% Create gif which plays at 10fps, and is titled 'example.gif'
AOcube.cube2gif('example',0.1); % (0.1Hz) = 1/(10fps)

% Check to see that the file was created. Note that the warning about 1st
% and 2nd dimensions is only relevant in specific AO cases that the
% AODataCube class is also used for.

% Open it in your default application for viewing gifs:
winopen 'example.gif'

% Any issues contact jesse.cranney@uon.edu.au