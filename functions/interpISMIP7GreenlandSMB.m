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
	%smb: SMBforcing class in ISSM.
	%}

	% Find appropriate directory
	switch oshostname(),
		case {'totten'}
			datadir='/totten_1/ModelData/ISMIP7/GrIS/';
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
		case 'MRI-ESM2-0'
			anomaly_pattern = 'GEMB-SDBN1-1000m/acabf-anomaly/v2/acabf*.nc';
	end

	smb_pattern = fullfile(datadir, modelname, scenario, anomaly_pattern);
	smb_file    = dir(smb_pattern);

	if isempty(smb_file)
		error('No %s SMB anomaly files found. Pattern: %s', scenario, smb_pattern);
	end

	% filenames end in the 4-digit year (e.g. ..._v3_2007.nc), so a plain
	% string sort puts the file list in chronological order
	[~, idx] = sort({smb_file.name});
	smb_file = smb_file(idx);

	% By default, use every yearly file found for this forcing. If start_end
	% is given, subset to that range instead (e.g. 2015-2100 out of the
	% 2015-2300 files available for a given SSP).
	%NOTE: File format: acabf_AIS_CESM2-WACCM_historical_SDBN1_v3_2014.nc
	if nargin == 4
		years = start_end(1):start_end(2);
		keep  = false(numel(smb_file), 1);
		for i = 1:numel(smb_file)
			yr_str = regexp(smb_file(i).name, '(\d{4})\.nc$', 'match', 'once');
			keep(i) = any(years == str2double(yr_str));
		end
		smb_file = smb_file(keep);
	end

	% Recover file names in cell.
	smb_file = fullfile({smb_file.folder}, {smb_file.name});

	% Load RACMO24p1_ERA5 dataset
	% Compute climatological mean value of SMB (Jan. 1995 - Dec. 2014 in Nowicki et al. (2020@TC))
	smb_clim = interpRACMO23p2MonthlySMB(md.mesh.x, md.mesh.y, 1995, 2014);
	smb_clim(smb_clim==-9999)=0;
	for ii = 1:size(smb_clim,2)
		pos0 = find(smb_clim(1:end-1,ii)==0);
		pos = find(smb_clim(1:end-1,ii)~=0);
		smb_clim(pos0,ii) = griddata(md.mesh.x(pos),md.mesh.y(pos),smb_clim(pos,ii),md.mesh.x(pos0),md.mesh.y(pos0),'nearest');
	end
	disp('Computing climatological mean')
	smb_clim = mean(smb_clim(1:end-1,:),2); % climatological mean value (exclude timestamps)

	% Load data from files
	disp('   == loading SMB anomaly data');
	x_n = double(ncread(smb_file{1},'x'));
	y_n = double(ncread(smb_file{1},'y'));	

	temp_matrix_smb_anon = [];
	temp_matrix_time= [];
	for i = 1:length(smb_file)
		fprintf('    processing file %d/%d \r',i,length(smb_file));

		%NOTE: unit for acabf in netcdf file: kg m-2 s-1
		smb_data = double(ncread(smb_file{i},'acabf-anomaly')); % dimension = (x,y,time)
		smb_data = smb_data/md.materials.rho_ice*md.constants.yts; % kg m-2 s-1 -> ice m yr-1

		%Load time data
		temp_time = double(ncread(smb_file{i},'time')); % time since year from current file...

		% configure out starting year of current file.
		[~, fname] = fileparts(smb_file{i});      % strip path and .nc
		tok = strsplit(fname, '_');
		temp_time_start = str2double(tok{end});
		
		% convert days in year decimal
		%FIXME: standard calendar for time is 365 days in year (with noleap)?
		temp_time = temp_time/365 + temp_time_start;
		temp_matrix_time = cat(1,temp_matrix_time, temp_time); % concatenate time series

		% Now, interpolate SMB 
		for j = 1:size(smb_data,3)
			temp_smb_anon = InterpFromGridToMesh(x_n,y_n,smb_data(:,:,j)',md.mesh.x,md.mesh.y,NaN);

			% Concatenate dataset
			temp_matrix_smb_anon = [temp_matrix_smb_anon, temp_smb_anon];
			clear temp_smb_anon;
		end
	end

	clear smb_data x_n y_n;

	% Now, SMB = SMB_ref + SMB_anomaly
	temp_matrix_smb = repmat(smb_clim,1,size(temp_matrix_smb_anon,2)) + temp_matrix_smb_anon;	

	% Save data
	smb = SMBforcing();
	smb.mass_balance = [temp_matrix_smb; temp_matrix_time'];
end
