function smb = interpISMIP7GreenlandSMB(md, modelname, scenario, start_end)
	%{
	%interpISMIP7GreenlandSMB - interpolate chosen ISMIP7 SMB forcing to model.
	%
	%Example
	%-------
	%.. code-block:: matlab
	%	md.smb = interpISMIP7GreenlandSMB(md,'CESM2-WACCM','ssp126');
	%	md.smb = interpISMIP7GreenlandSMB(md,'CESM2-WACCM','ssp585',[2015, 2100]);
	%	md.smb = interpISMIP7GreenlandSMB(md,'OCX');
	%	md.smb = interpISMIP7GreenlandSMB(md,'OCX',[2007, 2020]);
	%
	%Inputs
	%------
	%md: ISSM model
	%
	%modelname: string
	%	CMIP7 model, or 'OCX'. Supported: CESM2-WACCM, MRI-ESM2-0, OCX
	%
	%scenario: string
	%	Scenario in CMIP7 model. Supported: historical, ssp126, ssp370, ssp585, ctrl
	%	Only the files under this single scenario's directory are loaded --
	%	e.g. 'historical' loads just the historical years, and an SSP does
	%	not automatically get the historical record prepended.
	%	Not used for 'OCX' (it has no scenario level) -- when modelname is
	%	'OCX' this argument may be omitted or left empty, and start_end may
	%	be passed here instead (see below).
	%
	%start_end (optional int array): two entry array of [start_year end_year].
	%	If omitted, every yearly file found for the requested model/scenario
	%	is used. If given, the found files are subset to that range (e.g.
	%	[2015, 2100] out of files spanning 2015-2300 for a given SSP).
	%	For 'OCX', this defaults to [2007, 2025] instead, since OCX has no
	%	scenario-level range of its own. Because 'OCX' has no scenario
	%	argument, start_end may be passed as the 3rd argument instead of
	%	the 4th when modelname is 'OCX'.
	%
	%Output
	%------
	%smb: SMBgradients class in ISSM. smb.smbref holds the reference SMB --
	%	for CESM2-WACCM/MRI-ESM2-0, RACMO23p2 climatology + ISMIP7 SMB
	%	anomaly (monthly); for OCX, the RACMO23p2-ERA SMB directly (it is
	%	not an anomaly) -- and smb.b_pos/b_neg hold the SMB-vs-elevation
	%	gradient (dacabfdz, annual), each stacked with its own trailing
	%	time row. smb.href is set to md.geometry.surface.
	%}

	if nargin < 4
		start_end = [];
	end
	if nargin < 3
		scenario = '';
	end

	% OCX has no scenario level, so allow start_end to be passed
	% positionally as the 3rd argument instead of the 4th, e.g.
	% interpISMIP7GreenlandSMB(md,'OCX',[2007, 2020]).
	if strcmpi(modelname, 'OCX') && isnumeric(scenario) && ~isempty(scenario)
		start_end = scenario;
		scenario  = '';
	end

	% Find appropriate directory
	switch oshostname(),
		case {'totten'}
			datadir='/totten_1/ModelData/ISMIP7/ISMIP7/GrIS/';
		case {'epica'}
			datadir='/data2/issm/shields/ismip7greenland/ModelData/ISMIP7/GrIS/';
		otherwise
			error('machine not supported yet, please provide your own path');
	end

	valid_models    = {'CESM2-WACCM', 'MRI-ESM2-0', 'OCX'};
	valid_scenarios = {'historical', 'ssp126', 'ssp370', 'ssp585', 'ctrl'};

	if ~ismember(modelname, valid_models)
		error('Model ''%s'' not supported. Valid models: %s', modelname, strjoin(valid_models, ', '));
	end
	if ~strcmpi(modelname, 'OCX') && ~ismember(scenario, valid_scenarios)
		error('Scenario ''%s'' not supported. Valid scenarios: %s', scenario, strjoin(valid_scenarios, ', '));
	end

	% -----------------------------------------------------------------
	% OCX: no scenario level, and no separate climatology + anomaly --
	% interpRACMO23p2MonthlySMB.m already returns the actual SMB (not
	% an anomaly), monthly, stacked with its own trailing time row,
	% which is exactly the format smb.smbref needs, so it is used
	% as-is. The SMB-vs-elevation gradient (dacabfdz) is pulled from
	% its own OCX/RACMO2.3p2-ERA product, the same way as the other
	% models' dacabfdz.
	% -----------------------------------------------------------------
	if strcmpi(modelname, 'OCX')
		if ~isempty(scenario)
			warning('ISMIP7:ocxScenarioIgnored', ...
				'OCX has no scenario level; ignoring scenario argument ''%s''.', scenario);
		end

		ocx_acabf_dir        = fullfile(datadir, 'OCX', 'RACMO2.3p2-ERA', 'SDBN1-1000m', 'acabf',    'v1');
		ocx_dacabfdz_pattern = fullfile(datadir, 'OCX', 'RACMO2.3p2-ERA', 'SDBN1-1000m', 'dacabfdz', 'v1', 'dacabfdz*.nc');

		if isempty(start_end)
			ocx_start_year = 2007;
			ocx_end_year   = 2025;
		else
			ocx_start_year = start_end(1);
			ocx_end_year   = start_end(2);
		end

		disp(['   == Loading RACMO23p2-ERA SMB (OCX), ' num2str(ocx_start_year) '-' num2str(ocx_end_year)]);
		% NOTE: unlike CESM2-WACCM/MRI-ESM2-0, this is the actual SMB, not
		% an anomaly relative to a climatology.
		smb_ocx = interpRACMO23p2MonthlySMB(md.mesh.x, md.mesh.y, ocx_start_year, ocx_end_year, ocx_acabf_dir);

		[smb_dz, time_dz] = loadISMIP7YearlyField(ocx_dacabfdz_pattern, 'dacabfdz', ...
			[ocx_start_year ocx_end_year], md, 'SMB-vs-elevation gradient (dSMB/dz) for OCX RACMO2.3p2-ERA');

		smb        = SMBgradients();
		smb.href   = md.geometry.surface;
		smb.smbref = smb_ocx; % already (nVertices+1) x nTime, with trailing time row
		smb.b_pos  = [smb_dz; time_dz'];
		smb.b_neg  = [smb_dz; time_dz'];
		return
	end

	% Directory structure: model/scenario/variable-type/version/*.nc, one file per year
	switch modelname
		case 'CESM2-WACCM'
			anomaly_pattern = 'SDBN1-1000m/acabf-anomaly/v3/acabf*.nc';
			dz_pattern      = 'SDBN1-1000m/dacabfdz/v3/dacabfdz*.nc';
		case 'MRI-ESM2-0'
			anomaly_pattern = 'GEMB-SDBN1-1000m/acabf-anomaly/v2/acabf*.nc';
			dz_pattern      = 'GEMB-SDBN1-1000m/dacabfdz/v2/dacabfdz*.nc';
	end

	% Load RACMO24p1_ERA5 dataset and compute climatological mean value of SMB 
	% (Jan. 1960 - Dec. 1989, see Nowicki et al. 2020: https://tc.copernicus.org/articles/14/2331/2020/)
	clim_start_year = 1960;
	clim_end_year   = 1989;
	disp(['   == Loading RACMO23p2 climatology (' num2str(clim_start_year) '-' num2str(clim_end_year) ')']);
	smb_clim = interpRACMO23p2MonthlySMB(md.mesh.x, md.mesh.y, clim_start_year, clim_end_year);
	smb_clim(isnan(smb_clim))=0;

	pos0 = find(smb_clim(1:end-1,1)==0);
	pos  = find(smb_clim(1:end-1,1)~=0);
	if ~isempty(pos0)
		nearest_idx = dsearchn([md.mesh.x(pos), md.mesh.y(pos)], [md.mesh.x(pos0), md.mesh.y(pos0)]);
		for ii = 1:size(smb_clim,2)
			smb_clim(pos0,ii) = smb_clim(pos(nearest_idx),ii);
		end
	end

	disp('   -- computing climatological mean');
	smb_clim = mean(smb_clim(1:end-1,:),2); % climatological mean value (exclude timestamps)

	% Load the SMB anomaly (acabf-anomaly) and the SMB-vs-elevation gradient
	% (dacabfdz), both stored as one yearly file per variable.
	[smb_anom, time_anom] = loadISMIP7YearlyField(fullfile(datadir, modelname, scenario, anomaly_pattern), ...
		'acabf-anomaly', start_end, md, sprintf('SMB anomaly for %s %s', modelname, scenario));
	[smb_dz, time_dz] = loadISMIP7YearlyField(fullfile(datadir, modelname, scenario, dz_pattern), ...
		'dacabfdz', start_end, md, sprintf('SMB-vs-elevation gradient (dSMB/dz) for %s %s', modelname, scenario));

	% NOTE: acabf-anomaly is stored monthly (12 time records per yearly file)
	% while dacabfdz is an annual regression coefficient (1 time record per
	% yearly file), so time_anom and time_dz will have different lengths.
	% This is expected: smbref and b_pos/b_neg each carry their own trailing
	% time row and are interpolated independently in time by ISSM's
	% SMBgradients class, so they do not need to share a common time base.

	% Now, SMB = SMB_ref + SMB_anomaly
	temp_matrix_smb = repmat(smb_clim,1,size(smb_anom,2)) + smb_anom;

	% Save data
	smb = SMBgradients();
	smb.href   = md.geometry.surface;
	smb.smbref = [temp_matrix_smb; time_anom'];
	smb.b_pos  = [smb_dz; time_dz'];
	smb.b_neg  = [smb_dz; time_dz'];
end

function [data_matrix, time_vec] = loadISMIP7YearlyField(field_pattern, varname, start_end, md, label)
	%loadISMIP7YearlyField - find, subset, load and interpolate one ISMIP7
	%yearly-file variable (e.g. acabf-anomaly or dacabfdz) onto the model mesh.
	%
	%field_pattern is the full dir()-style glob for the yearly files, e.g.
	%'<datadir>/CESM2-WACCM/ssp126/SDBN1-1000m/dacabfdz/v3/dacabfdz*.nc', or
	%for OCX (no model/scenario subdirectories),
	%'<datadir>/OCX/RACMO2.3p2-ERA/SDBN1-1000m/dacabfdz/v1/dacabfdz*.nc'.
	%
	%Returns data_matrix (numberofvertices x ntime), already converted from
	%kg m-2 s-1 to ice m yr-1, and time_vec (ntime x 1) in decimal years.

	field_file = dir(field_pattern);

	if isempty(field_file)
		error('No files found for %s. Pattern: %s', label, field_pattern);
	end

	% filenames end in the 4-digit year (e.g. ..._v3_2007.nc), so a plain
	% string sort puts the file list in chronological order
	[~, idx] = sort({field_file.name});
	field_file = field_file(idx);

	% By default, use every yearly file found for this forcing. If start_end
	% is given, subset to that range instead (e.g. 2015-2100 out of the
	% 2015-2300 files available for a given SSP).
	%NOTE: File format: acabf_AIS_CESM2-WACCM_historical_SDBN1_v3_2014.nc
	if ~isempty(start_end)
		years = start_end(1):start_end(2);
		keep  = false(numel(field_file), 1);
		for i = 1:numel(field_file)
			yr_str = regexp(field_file(i).name, '\d{4}(?=\.nc$)', 'match', 'once');
			keep(i) = any(years == str2double(yr_str));
		end
		field_file = field_file(keep);
	end

	if isempty(field_file)
		error('No files found for %s within the requested start_end range. Pattern: %s', label, field_pattern);
	end

	% Recover file names in cell.
	field_file = fullfile({field_file.folder}, {field_file.name});

	% field_file is already chronologically sorted, so the first/last give the year range
	yr_first = regexp(field_file{1},   '\d{4}(?=\.nc$)', 'match', 'once');
	yr_last  = regexp(field_file{end}, '\d{4}(?=\.nc$)', 'match', 'once');
	disp(['   == Loading ISMIP7 ' label ', ' yr_first '-' yr_last]);
	x_n = double(ncread(field_file{1},'x'));
	y_n = double(ncread(field_file{1},'y'));

	% Crop to the mesh's bounding box (+ a small buffer) instead of reading
	% the full grid for every year/time-slice, same approach as
	% interpISMIP6GreenlandSMB.m. All yearly files for a given model/scenario
	% share the same x/y grid, so the crop indices only need computing once.
	offset = 2;

	xmin = min(md.mesh.x(:)); xmax = max(md.mesh.x(:));
	posx = find(x_n <= xmax);
	id1x = max(1, find(x_n >= xmin,1) - offset);
	id2x = min(numel(x_n), posx(end) + offset);

	ymin = min(md.mesh.y(:)); ymax = max(md.mesh.y(:));
	posy = find(y_n <= ymax);
	id1y = max(1, find(y_n >= ymin,1) - offset);
	id2y = min(numel(y_n), posy(end) + offset);

	x_n = x_n(id1x:id2x);
	y_n = y_n(id1y:id2y);

	data_matrix = [];
	time_vec    = [];
	progress_msg = '';
	for i = 1:length(field_file)
		% configure out starting year of current file.
		[~, fname] = fileparts(field_file{i});      % strip path and .nc
		tok = strsplit(fname, '_');
		temp_time_start = str2double(tok{end});

		% print the year being read, erasing the previous one in place
		fprintf(repmat('\b', 1, length(progress_msg)));
		progress_msg = sprintf('   -- year %d (%d/%d)', temp_time_start, i, length(field_file));
		fprintf('%s', progress_msg);

		%NOTE: unit for acabf-based fields in netcdf file: kg m-2 s-1 (or per m elevation for dacabfdz)
		% Only read the cropped x/y window computed above; Inf reads every
		% time step present in this file (12 for monthly acabf-anomaly, 1
		% for annual dacabfdz).
		field_data = double(ncread(field_file{i}, varname, ...
			[id1x id1y 1], [id2x-id1x+1, id2y-id1y+1, Inf])); % dimension = (x,y,time)
		field_data = field_data/md.materials.rho_ice*md.constants.yts; % kg m-2 s-1 -> ice m yr-1

		%Load time data
		%NOTE: the 'time' variable's reference date is NOT the same across
		%all ISMIP7 GrIS files: historical/ssp files use a fixed epoch of
		%1850-01-01, but 'ctrl' files instead encode time relative to each
		%file's own year (e.g. dacabfdz: 'days since 2206-12-31 ...',
		%acabf-anomaly: 'days since 2206-01-16'). Hardcoding the 1850-01-01
		%epoch decodes 'ctrl' files to ~1850 instead of their real year, so
		%the reference date is read from each file's own time:units
		%attribute instead of being assumed.
		%date2decyear expects a MATLAB datenum, so convert the parsed
		%reference date + raw day offset into one before calling it (same
		%pattern used for the ISMIP7 ocean forcing in
		%interpISMIP7GreenlandOcn.m).
		time_units   = ncreadatt(field_file{i}, 'time', 'units'); % e.g. 'days since 1850-01-01' or 'days since 2206-01-16'
		ref_date_str = regexp(time_units, '\d{4}-\d{2}-\d{2}', 'match', 'once');
		if isempty(ref_date_str)
			error('Could not parse a reference date out of time:units ''%s'' for %s', time_units, field_file{i});
		end
		ref_datenum   = datenum(ref_date_str, 'yyyy-mm-dd');
		temp_time_raw = double(ncread(field_file{i},'time')); % days since ref_datenum
		temp_time     = date2decyear(ref_datenum + temp_time_raw);

		% sanity check: decoded year should be within ~1 year of the
		% filename's year (catches any future change in the time encoding)
		if any(abs(temp_time - temp_time_start) > 1.5)
			warning(['ISMIP7:timeDecodeMismatch: decoded time (' num2str(temp_time(1)) ...
				') is far from the year implied by the filename (' num2str(temp_time_start) ...
				') for ' field_file{i} '. Check the time units/encoding.']);
		end
		time_vec = cat(1,time_vec, temp_time); % concatenate time series

		% Now, interpolate onto the mesh
		for j = 1:size(field_data,3)
			temp_interp = InterpFromGrid(x_n,y_n,field_data(:,:,j)',double(md.mesh.x),double(md.mesh.y));

			% Concatenate dataset
			data_matrix = [data_matrix, temp_interp];
			clear temp_interp;
		end
	end
	fprintf('\n');

	clear field_data x_n y_n;
end

