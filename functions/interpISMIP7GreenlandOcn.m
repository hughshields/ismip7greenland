function frontalforcing = interpISMIP7GreenlandOcn(md,model_name,scenario)
%interpISMIP7GreenlandOcn - interpolate ISMIP7 climate model forcing 
%									 to prepare frontal forcings via recommended parameterization
%
%   Parameterization from Rignot et al. 2016:
%							https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016GL068784	
%					  see also Slater et al. 2020:
%							https://tc.copernicus.org/articles/14/985/2020/
%
%   Thermal forcing (tf) and subglacial discharge (sgd) are each stored as
%   one NetCDF file per year; this function globs every file in the
%   relevant directory, loads each year, and concatenates them into one
%   multi-year time series before interpolating onto the model mesh.
%
%   Input:
%     - md (model object)
%     - model_name (string): name of the climate model
%     - scenario (string): climate scenario
%     - supported combinations:
%
%             Model              Scenarios
%             -----------------------------------------------
%             CESM2-WACCM        historical, ssp370, ssp126, ssp585, calving-calibration
%             MRI-ESM2-0         historical, ssp370, ssp126, ssp585, calving-calibration
%
%     - 'calving-calibration' is not a raw model scenario: for the given
%       model_name it concatenates every available year of 'historical'
%       forcing with 'ssp370' forcing restricted to 2015-2024, producing one
%       continuous, observationally-grounded time series (e.g. for
%       calibrating a calving law over 2007-2022). See SPLICE_SCENARIO,
%       SPLICE_START_YEAR and SPLICE_END_YEAR below to change the splice
%       scenario or window.
%
%   Units:
%     - thermal forcing (tf):        deg C
%     - subglacial discharge (sgd):  m^3 s^-1 in the source files, converted
%                                    to m^3 d^-1 in the returned frontalforcing
%                                    (frontalforcingsrignot expects m^3/day)
%     - x, y (grid and mesh coords): meters
%     - time (last row of each matrix returned): decimal year
%     - basin_id:                    unitless integer basin label
%
%   Note: the sgd source data is a basin-integrated total discharge --
%   every grid cell within a given basin carries the same value (that
%   basin's total), not a spatially-varying local discharge. Interpolating
%   it onto the mesh (as done below) therefore does not add spatial detail
%   within a basin, it only resamples the same repeated value.
%
%   Output:
%     - frontalforcing: prepared to be input directly into md.frontalforcings
%
%   Examples:
%      md.frontalforcings = interpISMIP7GreenlandOcn(md,'CESM2-WACCM','ssp370');
%      md.frontalforcings = interpISMIP7GreenlandOcn(md,'MRI-ESM2-0','historical');
%      md.frontalforcings = interpISMIP7GreenlandOcn(md,'CESM2-WACCM','calving-calibration');

% ---------------------------------------------------------------------
% Resolve data directory and check inputs
% ---------------------------------------------------------------------
switch oshostname(),
	case {'totten'}
		path='/totten_1/ModelData/ISMIP7/';
	case {'epica'}
		path='/data2/issm/shields/ismip7greenland/ModelData/ISMIP7/';
	otherwise
		error('machine not supported yet, please provide your own path');
end

valid_models = {'CESM2-WACCM', 'MRI-ESM2-0'};
valid_scenarios = {'historical', 'ssp370', 'ssp126', 'ssp585', 'calving-calibration'};

if ~ismember(model_name, valid_models)
	error('Model ''%s'' not supported. Valid models: %s', model_name, strjoin(valid_models, ', '));
end

if ~ismember(scenario, valid_scenarios)
	error('Scenario ''%s'' not supported. Valid scenarios: %s', scenario, strjoin(valid_scenarios, ', '));
end

% Directory structure: model/scenario/variable-type/version/*.nc, one file per year.
%
% 'calving-calibration' is not itself a scenario directory on disk: it
% concatenates the given model's 'historical' files (every year available)
% with its SPLICE_SCENARIO files restricted to SPLICE_START_YEAR-SPLICE_END_YEAR,
% giving one historical+near-term time series for calving calibration work.
SPLICE_SCENARIO   = 'ssp370';
SPLICE_START_YEAR = 2015;
SPLICE_END_YEAR   = 2024;

switch model_name
	case 'CESM2-WACCM'
		tf_subpath  = 'ocean-1000m/tf/v2';
		sgd_subpath = 'SDBN1-1000m/sgd/v2';
	case 'MRI-ESM2-0'
		tf_subpath  = 'ocean-1000m/tf/v1';
		sgd_subpath = 'GEMB-SDBN1-1000m/sgd/v1';
end

if strcmp(scenario, 'calving-calibration')

	% -- historical portion: use every year that is available --
	hist_root        = [path 'GrIS/' model_name '/historical'];
	hist_tf_pattern  = [hist_root '/' tf_subpath  '/*.nc'];
	hist_sgd_pattern = [hist_root '/' sgd_subpath '/*.nc'];
	hist_tffiles  = dir(hist_tf_pattern);
	hist_sgdfiles = dir(hist_sgd_pattern);

	if isempty(hist_tffiles)
		error('No historical thermal forcing files found. Pattern: %s', hist_tf_pattern);
	end
	if isempty(hist_sgdfiles)
		error('No historical subglacial discharge files found. Pattern: %s', hist_sgd_pattern);
	end
	if length(hist_tffiles) ~= length(hist_sgdfiles)
		error(['Number of historical thermal forcing files (%d) does not match number of ' ...
			'historical subglacial discharge files (%d). Check that both directories cover ' ...
			'the same years.'], length(hist_tffiles), length(hist_sgdfiles));
	end

	% -- future portion: splice in SPLICE_SCENARIO, restricted to the splice window --
	splice_root        = [path 'GrIS/' model_name '/' SPLICE_SCENARIO];
	splice_tf_pattern  = [splice_root '/' tf_subpath  '/*.nc'];
	splice_sgd_pattern = [splice_root '/' sgd_subpath '/*.nc'];
	splice_tffiles  = dir(splice_tf_pattern);
	splice_sgdfiles = dir(splice_sgd_pattern);

	if isempty(splice_tffiles)
		error('No %s thermal forcing files found. Pattern: %s', SPLICE_SCENARIO, splice_tf_pattern);
	end
	if isempty(splice_sgdfiles)
		error('No %s subglacial discharge files found. Pattern: %s', SPLICE_SCENARIO, splice_sgd_pattern);
	end

	% restrict the splice scenario to the splice window years only
	splice_tf_years  = cellfun(@(n) str2double(regexp(n, '\d{4}(?=\.nc$)', 'match', 'once')), {splice_tffiles.name});
	splice_sgd_years = cellfun(@(n) str2double(regexp(n, '\d{4}(?=\.nc$)', 'match', 'once')), {splice_sgdfiles.name});
	splice_tffiles   = splice_tffiles(splice_tf_years  >= SPLICE_START_YEAR & splice_tf_years  <= SPLICE_END_YEAR);
	splice_sgdfiles  = splice_sgdfiles(splice_sgd_years >= SPLICE_START_YEAR & splice_sgd_years <= SPLICE_END_YEAR);

	if isempty(splice_tffiles)
		error('No %s thermal forcing files found in the splice window %d-%d. Pattern: %s', ...
			SPLICE_SCENARIO, SPLICE_START_YEAR, SPLICE_END_YEAR, splice_tf_pattern);
	end
	if isempty(splice_sgdfiles)
		error('No %s subglacial discharge files found in the splice window %d-%d. Pattern: %s', ...
			SPLICE_SCENARIO, SPLICE_START_YEAR, SPLICE_END_YEAR, splice_sgd_pattern);
	end
	if length(splice_tffiles) ~= length(splice_sgdfiles)
		error(['Number of %s splice-window thermal forcing files (%d) does not match number of ' ...
			'%s splice-window subglacial discharge files (%d). Check that both directories cover ' ...
			'the same years.'], SPLICE_SCENARIO, length(splice_tffiles), SPLICE_SCENARIO, length(splice_sgdfiles));
	end

	disp(['   == ''calving-calibration'': splicing ' model_name ' historical with ' ...
		SPLICE_SCENARIO ' (' num2str(SPLICE_START_YEAR) '-' num2str(SPLICE_END_YEAR) ')']);

	tffiles  = [hist_tffiles;  splice_tffiles];
	sgdfiles = [hist_sgdfiles; splice_sgdfiles];

else
	rootname    = [path 'GrIS/' model_name '/' scenario];
	tf_pattern  = [rootname '/' tf_subpath  '/*.nc'];
	sgd_pattern = [rootname '/' sgd_subpath '/*.nc'];

	tffiles  = dir(tf_pattern);
	sgdfiles = dir(sgd_pattern);

	if isempty(tffiles)
		error('No thermal forcing files found. Pattern: %s', tf_pattern);
	end
	if isempty(sgdfiles)
		error('No subglacial discharge files found. Pattern: %s', sgd_pattern);
	end
end

% filenames end in the 4-digit year (e.g. ..._v2_2007.nc). Sort on the
% extracted year itself (rather than the raw filename string) so that the
% historical+splice file set built above -- whose filenames mix two
% different scenario names -- still ends up in true chronological order.
tf_years  = cellfun(@(n) str2double(regexp(n, '\d{4}(?=\.nc$)', 'match', 'once')), {tffiles.name});
[~, idx]  = sort(tf_years);
tffiles   = tffiles(idx);

sgd_years = cellfun(@(n) str2double(regexp(n, '\d{4}(?=\.nc$)', 'match', 'once')), {sgdfiles.name});
[~, idx]  = sort(sgd_years);
sgdfiles  = sgdfiles(idx);

if length(tffiles) ~= length(sgdfiles)
	error(['Number of thermal forcing files (%d) does not match number of ' ...
		'subglacial discharge files (%d). Check that both directories cover ' ...
		'the same years.'], length(tffiles), length(sgdfiles));
end

% ---------------------------------------------------------------------
% Load all years and concatenate into one multi-year time series
% ---------------------------------------------------------------------
nfiles = length(tffiles);

% filenames are chronologically sorted, so the first/last give the year range
yr_first = regexp(tffiles(1).name,   '\d{4}(?=\.nc$)', 'match', 'once');
yr_last  = regexp(tffiles(end).name, '\d{4}(?=\.nc$)', 'match', 'once');
disp(['   == Loading ISMIP7 ocean forcing (TF/SGD) for ' model_name ' ' scenario ', ' ...
	yr_first '-' yr_last]);

x_n = [];
y_n = [];
time_days = zeros(0,1);            % days since 1900-1-1, common axis across years
tf_list  = cell(1, nfiles);        % thermal forcing,       [x, y, time], deg C
sgd_list = cell(1, nfiles);        % subglacial discharge,  [x, y, time], m^3 s^-1
progress_msg = '';

for f = 1:nfiles
	% print the year being read, erasing the previous one in place
	yr_str = regexp(tffiles(f).name, '\d{4}(?=\.nc$)', 'match', 'once');
	fprintf(repmat('\b', 1, length(progress_msg)));
	progress_msg = sprintf('   -- year %s (%d/%d)', yr_str, f, nfiles);
	fprintf('%s', progress_msg);

	tfnc  = fullfile(tffiles(f).folder,  tffiles(f).name);
	sgdnc = fullfile(sgdfiles(f).folder, sgdfiles(f).name);

	% grid is assumed identical across years, so only read it once
	if f == 1
		x_n = double(ncread(tfnc, 'x')); % meters
		y_n = double(ncread(tfnc, 'y')); % meters
	end

	% each file's time units give its own reference date, e.g.
	% "days since 2007-01-16 00:00:00" -- convert to a common axis
	time_f     = double(ncread(tfnc, 'time'));
	time_units = ncreadatt(tfnc, 'time', 'units');
	tok = regexp(time_units, 'since\s+(\d+)-(\d+)-(\d+)', 'tokens', 'once');
	if isempty(tok)
		error('Could not parse time reference date from units string: %s', time_units);
	end
	refdate  = [str2double(tok{1}), str2double(tok{2}), str2double(tok{3})];
	time_days(end+1:end+length(time_f), 1) = time_f + datenum(refdate) - datenum(1900,1,1);

	tf_list{f}  = double(ncread(tfnc,  'tf'));
	sgd_list{f} = double(ncread(sgdnc, 'sgd'));

	for i = 1:3
		assert(size(sgd_list{f}, i) == size(tf_list{f}, i), ...
			'Dimension mismatch between TF and SGD at dimension %d in year file %s', ...
			i, tffiles(f).name);
	end
end
fprintf('\n');

tf_data  = cat(3, tf_list{:});
sgd_data = cat(3, sgd_list{:});

% guard against the yearly files not already being in chronological order
[time_days, sortidx] = sort(time_days);
tf_data  = tf_data(:, :, sortidx);
sgd_data = sgd_data(:, :, sortidx);

time_year = str2num(datestr(time_days + datenum(1900,1,1), 'yyyy')); %#ok<ST2NM>
nt        = length(time_year);

% ---------------------------------------------------------------------
% Interpolate onto the model mesh
% ---------------------------------------------------------------------
disp(['   == Interpolating on model mesh (' num2str(md.mesh.numberofvertices) ' vertices)']);
TF_matrix  = zeros(md.mesh.numberofvertices, nt); % deg C, per vertex per month
SGD_matrix = zeros(md.mesh.numberofvertices, nt); % m^3 s^-1, per vertex per month (converted to m^3/day below)

for i = 1:nt
	TF_matrix(:, i) = InterpFromGrid(x_n, y_n, tf_data(:, :, i)', ...
		md.mesh.x, md.mesh.y, 'linear');
	SGD_matrix(:, i) = InterpFromGrid(x_n, y_n, sgd_data(:, :, i)', ...
		md.mesh.x, md.mesh.y, 'linear');
end

TF_matrix  = max(0, TF_matrix);
SGD_matrix = max(0, SGD_matrix);

% ---------------------------------------------------------------------
% Assemble the frontalforcingsrignot object to return.
%
% md is passed by value in MATLAB, so writing into md.frontalforcings here
% would not persist back to the caller -- build into 'frontalforcing' (the
% function's actual output variable) instead, matching the documented call
% pattern: md.frontalforcings = interpISMIP7GreenlandOcn(md, model, scenario);
% ---------------------------------------------------------------------
time_decyear = date2decyear(time_days + datenum(1900,1,1))'; % decimal year, row vector

frontalforcing = frontalforcingsrignot();

% one row per mesh vertex, plus a final row holding time (decimal years) --
% matches ISSM's marshall() for this class ('timeserieslength', numberofvertices + 1)
frontalforcing.thermalforcing       = zeros(md.mesh.numberofvertices + 1, nt);
frontalforcing.subglacial_discharge = zeros(md.mesh.numberofvertices + 1, nt);

frontalforcing.thermalforcing(1:end-1, :)       = TF_matrix;
frontalforcing.subglacial_discharge(1:end-1, :) = SGD_matrix * 86400; % m^3 s^-1 -> m^3 d^-1 (frontalforcingsrignot expects m^3/day)

frontalforcing.thermalforcing(end, :)       = time_decyear;
frontalforcing.subglacial_discharge(end, :) = time_decyear;

% ---------------------------------------------------------------------
% Basin IDs (one per mesh element), for grouping calving-front vertices
% into discharge/melt basins
% ---------------------------------------------------------------------
disp('   == Reading basin data from NetCDF file');
basin_file = [path 'tools/ismip7-gris-ocean-forcing/subglacial_discharge_basins_ismip.nc'];

x_basin    = double(ncread(basin_file, 'x'));      % meters
y_basin    = double(ncread(basin_file, 'y'));      % meters
basin_data = double(ncread(basin_file, 'basin'));  % unitless integer basin label

if ~isequal(size(basin_data), [length(x_basin), length(y_basin)])
	basin_data = basin_data';
end

unique_basins = unique(basin_data(:));
unique_basins(unique_basins == 0) = [];  % 0 = background / no basin
num_basins = length(unique_basins);
disp(['   -- found ' num2str(num_basins) ' basins']);

% remap to contiguous 1..N labels in case the raster's raw ids aren't already
basin_data_remapped = zeros(size(basin_data));
for b = 1:num_basins
	basin_data_remapped(basin_data == unique_basins(b)) = b;
end

xc = mean(md.mesh.x(md.mesh.elements), 2); % element centroid x, meters
yc = mean(md.mesh.y(md.mesh.elements), 2); % element centroid y, meters

% 'nearest' avoids blending basin labels into meaningless intermediate values
frontalforcing.basin_id   = InterpFromGrid(x_basin, y_basin, basin_data_remapped', ...
	xc, yc, 'nearest');
frontalforcing.num_basins = num_basins;
end
