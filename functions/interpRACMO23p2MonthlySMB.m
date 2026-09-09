function [output] = interpRACMO23p2MonthlySMB(X,Y,t_start,t_end,ncpath),
	%   Input:
	%     - X:       model x values
	%     - Y:       model y values
	%     - t_start: first decimal year time for which you want SMB data (float)
	%     - t_end:   last decimal year time for which you want SMB data (float)
	%     - ncpath:  (optional) path to netCDF files directory
	%                default: '/totten_1/ModelData/ISMIP7/GrIS/OCX/RACMO2.3p2-ERA/SDBN1-1000m/acabf/v1'
	%
	%   Output:
	%     - output:  matrix of size number_of_vertices+1 x number of time steps
	%                it is a P1 timeseries, values in m ice/yr
	%
	%   Examples:
	%      smb = interpRACMO23p2MonthlySMB(md.mesh.x,md.mesh.y,2007.0,2022.5);
	%      smb = interpRACMO23p2MonthlySMB(md.mesh.x,md.mesh.y,2007.0,2022.5,'/custom/path/to/data');

	% Set default path if not provided
	if nargin < 5 || isempty(ncpath)
		switch oshostname(),
			case {'totten'}
				ncpath = '/totten_1/ModelData/ISMIP7/GrIS/OCX/RACMO2.3p2-ERA/SDBN1-1000m/acabf/v1';
			case {'epica'}
				ncpath = '/data2/issm/shields/ismip7greenland/ModelData/ISMIP7/GrIS/OCX/RACMO2.3p2-ERA/SDBN1-1000m/acabf/v1';
			otherwise
				error('machine not supported yet, please provide your own path');
		end
	end
	filestruct = dir([ncpath '/*.nc']);
	filenames = {filestruct(:).name};
	directories = {filestruct(:).folder};

	% Identify up front which yearly files fall within the requested time
	% range, so the progress line below can show "year (n/total)" the same
	% way interpISMIP7GreenlandOcn does.
	% Files look like: acabf_GrIS_RACMO2.3p2-ERA_OCX_SDBN1-1000m_v1_1995.nc
	match_idx  = [];
	match_year = {};
	for ii = 1:length(filenames)
		tok = regexp(filenames{ii}, '^acabf_.*_(\d{4})\.nc$', 'once', 'tokens');
		if isempty(tok)
			continue
		end
		yr_str = tok{1};
		if (str2double(yr_str) < floor(t_start)) || (str2double(yr_str) > t_end)
			continue
		end
		match_idx(end+1)  = ii; %#ok<AGROW>
		match_year{end+1} = yr_str; %#ok<AGROW>
	end
	nfiles = length(match_idx);

	% Initialize output matrix
	output = NaN(length(X)+1, 1);

	count = 1;
	progress_msg = '';
	for kk = 1:nfiles
		ii       = match_idx(kk);
		filename = filenames{ii};
		yr_str   = match_year{kk};
		yr_date  = datetime([yr_str '-01-01']);

		% print the year being read, erasing the previous one in place
		fprintf(repmat('\b', 1, length(progress_msg)));
		progress_msg = sprintf('   -- year %s (%d/%d)', yr_str, kk, nfiles);
		fprintf('%s', progress_msg);

		% Load x and y coordinate vectors (already in meters)
		x_full = double(ncread([directories{ii} '/' filename],'x'))';
		y_full = double(ncread([directories{ii} '/' filename],'y'))';

		% Find subset of netCDF file to use
		offset = 2;
		xmin = min(X(:)); xmax = max(X(:));
		ymin = min(Y(:)); ymax = max(Y(:));

		inX = find(x_full >= xmin & x_full <= xmax);
		id1x = max(1, min(inX) - offset);
		id2x = min(numel(x_full), max(inX) + offset);

		inY = find(y_full >= ymin & y_full <= ymax);
		id1y = max(1, min(inY) - offset);
		id2y = min(numel(y_full), max(inY) + offset);

		x = x_full(id1x:id2x);
		y = y_full(id1y:id2y);

		% Load time (units are 'days since 1850-01-01 00:00:00')
		time_temp = double(ncread([directories{ii} '/' filename],'time'));
		dateyear = datenum(1850,1,1) - 1 + time_temp;
		time = date2decyear(dateyear);

		% Days in the calendar year, used below to annualize the per-second flux
		leapyear = @(year) mod(year, 4) == 0 & (mod(year, 100) ~= 0 | mod(year, 400) == 0);
		daysinyear = 365 + leapyear(year(yr_date));
		seconds_per_year = daysinyear * 86400;
		unit_transformation = seconds_per_year / 917; % ice density, kg/m^3

		% Load acabf (land ice surface specific mass balance flux, kg m-2 s-1)
		data = double(ncread([directories{ii} '/' filename],'acabf',[id1x id1y 1],[id2x-id1x+1 id2y-id1y+1 length(time)],[1 1 1]));
		data(abs(data) >= 1.e20) = NaN;

		% Loop through months, transform units, regrid, and put into output matrix
		% Transform units from kg m-2 s-1 to m ice/yr, assuming ice density 917 kg/m3:
		% (kg/m2/s * s/yr) / (917 kg/m3) = m ice/yr.
		for jj = 1:size(data,3)
			datamat = data(:,:,jj)' * unit_transformation;

			% Put data and times into output matrix
			if count == 1
				output(1:end-1,1) = InterpFromGrid(x, y, datamat, double(X), double(Y));
			else
				output(1:end-1,end+1) = InterpFromGrid(x, y, datamat, double(X), double(Y));
			end
			output(end,end) = time(jj);

			count = count + 1;
		end
	end
	fprintf('\n');
end
