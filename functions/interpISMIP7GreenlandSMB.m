function smb = interpISMIP7GreenlandSMB(md, modelname, scenario, start_end)
	%{
	%interpISMIP7GreenlandSMB - interpolate chosen ISMIP7 SMB forcing to model.
	%
	%Example
	%-------
	%.. code-block:: matlab
	%	md.smb = interpISMIP7GreenlandSMB(md,'CESM2-WACCM','ssp126');
	%	md.smb = interpISMIP7GreenlandSMB(md,'CESM2-WACCM','ssp585',[2015, 2100]);
	%
	%Inputs
	%------
	%md: ISSM model
	%
	%modelname: string
	%	CMIP7 model. Supported: CESM2-WACCM, MRI-ESM2-0
	%
	%scenario: string
	%	Scenario in CMIP7 model. Supported: historical, ssp126, ssp370, ssp585
	%	Only the files under this single scenario's directory are loaded --
	%	e.g. 'historical' loads just the historical years, and an SSP does
	%	not automatically get the historical record prepended.
	%
	%start_end (optional int array): two entry array of [start_year end_year].
	%	If omitted, every yearly file found for the requested model/scenario
	%	is used. If given, the found files are subset to that range (e.g.
	%	[2015, 2100] out of files spanning 2015-2300 for a given SSP).
	%
	%Output
	%------
	%smb: SMBgradients class in ISSM. smb.smbref holds the reference SMB
	%	(RACMO23p2 climatology + ISMIP7 SMB anomaly, monthly), and
	%	smb.b_pos/b_neg hold the ISMIP7 SMB-vs-elevation gradient
	%	(dacabfdz, annual), each stacked with its own trailing time row.
	%	smb.href is set to md.geometry.surface.
	%}

	if nargin < 4
		start_end = [];
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

	valid_models    = {'CESM2-WACCM', 'MRI-ESM2-0'};
	valid_scenarios = {'historical', 'ssp370', 'ssp126', 'ssp585'};

	if ~ismember(modelname, valid_models)
		error('Model ''%s'' not supported. Valid models: %s', modelname, strjoin(valid_models, ', '));
	end
	if ~ismember(scenario, valid_scenarios)
		error('Scenario ''%s'' not supported. Valid scenarios: %s', scenario, strjoin(valid_scenarios, ', '));
	end

	% Directory structure: model/scenario/variable-type/version/*.nc, one file per year
	switch modelname
		case 'CESM2-WACCM'
			anomaly_pattern = 'SDBN1-1000m/acabf-anomaly/v3/acabf*.nc';
			dz_pattern      = 'SDBN1-1000m/dacabfdz/v3/acabf*.nc';
		case 'MRI-ESM2-0'
			anomaly_pattern = 'GEMB-SDBN1-1000m/acabf-anomaly/v2/acabf*.nc';
			dz_pattern      = 'GEMB-SDBN1-1000m/dacabfdz/v2/acabf*.nc';
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
	[smb_anom, time_anom] = loadISMIP7YearlyField(datadir, modelname, scenario, ...
		anomaly_pattern, 'acabf-anomaly', start_end, md, 'SMB anomaly');
	[smb_dz, time_dz] = loadISMIP7YearlyField(datadir, modelname, scenario, ...
		dz_pattern, 'dacabfdz', start_end, md, 'SMB-vs-elevation gradient (dSMB/dz)');

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

function [data_matrix, time_vec] = loadISMIP7YearlyField(datadir, modelname, scenario, pattern, varname, start_end, md, label)
	%loadISMIP7YearlyField - find, subset, load and interpolate one ISMIP7
	%yearly-file variable (e.g. acabf-anomaly or dacabfdz) onto the model mesh.
	%
	%Returns data_matrix (numberofvertices x ntime), already converted from
	%kg m-2 s-1 to ice m yr-1, and time_vec (ntime x 1) in decimal years.

	field_pattern = fullfile(datadir, modelname, scenario, pattern);
	field_file    = dir(field_pattern);

	if isempty(field_file)
		error('No %s %s files found. Pattern: %s', scenario, label, field_pattern);
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

	% Recover file names in cell.
	field_file = fullfile({field_file.folder}, {field_file.name});

	% field_file is already chronologically sorted, so the first/last give the year range
	yr_first = regexp(field_file{1},   '\d{4}(?=\.nc$)', 'match', 'once');
	yr_last  = regexp(field_file{end}, '\d{4}(?=\.nc$)', 'match', 'once');
	disp(['   == Loading ISMIP7 ' label ' for ' modelname ' ' scenario ', ' ...
		yr_first '-' yr_last]);
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
		%NOTE: confirmed via ncdump that 'time' is stored as absolute days
		%since 1850-01-01 (standard/Gregorian calendar) for both
		%acabf-anomaly and dacabfdz files -- NOT days-of-year -- so it must
		%be decoded against that epoch rather than added to the filename year.
		temp_time_raw = double(ncread(field_file{i},'time')); % days since 1850-01-01
		temp_time = ismip7_days1850_to_decyear(temp_time_raw);

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
			temp_interp = InterpFromGridToMesh(x_n,y_n,field_data(:,:,j)',md.mesh.x,md.mesh.y,NaN);

			% Concatenate dataset
			data_matrix = [data_matrix, temp_interp];
			clear temp_interp;
		end
	end
	fprintf('\n');

	clear field_data x_n y_n;
end

function decyear = ismip7_days1850_to_decyear(days_since_1850)
	%ismip7_days1850_to_decyear - convert ISMIP7 'time' values (days since
	%1850-01-01, standard/Gregorian calendar -- confirmed via ncdump) into a
	%decimal year, accounting for leap years exactly rather than assuming a
	%fixed 365 or 365.25-day year.

	ref_date = datetime(1850,1,1);
	dates    = ref_date + days(days_since_1850);
	yrs      = year(dates);
	doy      = day(dates,'dayofyear');

	is_leap = (mod(yrs,4)==0 & (mod(yrs,100)~=0 | mod(yrs,400)==0));
	days_in_year = 365 + double(is_leap);

	decyear = yrs + (doy-1)./days_in_year;
end
