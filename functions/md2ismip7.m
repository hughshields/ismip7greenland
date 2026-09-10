function results=md2ismip7(md,directoryname,icesheetname,source_id,ism_id,ism_member_id,esm_id,forcing_member_id,experiment_id,set_counter,time_range,resolution_km,output_interval_yr);
	%Create netcdf files for an experiment following ISMIP7 conventions
	%
	%directoryname:      base output directory (ISMIP7 subdirectories are created below it)
	%icesheetname:       domain_id, 'GrIS' or 'AIS'
	%source_id:          modelling group name (no '_', '.' or special characters)
	%ism_id:             ice sheet model name and version (no '_', '.' or special characters)
	%ism_member_id:      ISM choice variant, e.g. 'm001'
	%esm_id:             CMIP/ESM model name, e.g. 'CESM2-WACCM'
	%forcing_member_id:  forcing choice variant, e.g. 'f001'
	%experiment_id:      'historical', 'ctrl', 'ssp126', ...
	%set_counter:        unique index within a set, e.g. 'C001', 'E001', 'P001'
	%time_range:         start-end year of the experiment, e.g. '2015-2300' (optional;
	%                    auto-derived from the model output years if empty/omitted)
	%resolution_km:      output grid resolution in km (optional, default 8).
	%                    Allowed: GrIS 1/2/4/8/16 km, AIS 2/4/8/16 km.
	%output_interval_yr: interval in years between 2-D (gridded) output snapshots
	%                    (optional, default 1). Scalars are always annual.

	experiment = experiment_id; %kept for the init single-snapshot branching below

	if nargin<12 || isempty(resolution_km),
		resolution_km=8; %default resolution valid for both GrIS and AIS
	end
	if nargin<13 || isempty(output_interval_yr),
		output_interval_yr=1;
	end

	if exist(directoryname)~=7,
		error(['directory ' directoryname ' does not exist']);
	end

	if ~strcmp(icesheetname,'GrIS') & ~strcmp(icesheetname,'AIS'),
		error('icesheetname (domain_id) should be GrIS or AIS');
	end

	if isempty(experiment_id),
		error('experiment_id must be provided (e.g. historical, ctrl, ssp126)');
	end

	%Allowed ISMIP7 resolutions (km): GrIS 1/2/4/8/16, AIS 2/4/8/16
	if strcmp(icesheetname,'GrIS'),
		allowed_res=[1 2 4 8 16];
	else
		allowed_res=[2 4 8 16];
	end
	if ~any(resolution_km==allowed_res),
		error(['resolution_km=' num2str(resolution_km) ' km is not allowed for ' icesheetname ' (allowed: ' num2str(allowed_res) ' km)']);
	end


	%ISMIP7 output directory structure: <domain_id>/<source_id>/<ism_id>/<set_id>/<set_counter>/
	switch upper(set_counter(1))
		case 'C', set_id='CORE';
		case 'E', set_id='ESM';
		case 'P', set_id='PPE';
		otherwise, error('set_counter must start with C, E or P (e.g. C001)');
	end
	outdir = fullfile(directoryname,icesheetname,source_id,ism_id,set_id,set_counter);
	if exist(outdir)~=7,
		mkdir(outdir);
	end

	%Derive <time_range> from the model output years if it was not supplied. Uses
	%floor(TransientSolution.time), so it assumes .time is in (absolute) calendar years;
	%if your times are relative to the run start, pass time_range explicitly instead.
	if nargin<11 || isempty(time_range),
		ytmp = floor(arrayfun(@(k) md.results.TransientSolution(k).time, 1:numel(md.results.TransientSolution)));
		if strcmp(experiment,'init'),
			y0=ytmp(end); y1=ytmp(end);            %single snapshot: final stored year
		else
			y0=min(ytmp); y1=max(ytmp);
		end
		if y0==y1,
			time_range=sprintf('%04d',y0);          %single year, e.g. 2014
		else
			time_range=sprintf('%04d-%04d',y0,y1);  %year range, e.g. 1850-2014
		end
	end

	%ISMIP7 filename builder:
	%<variable_id>_<domain_id>_<source_id>_<ism_id>_<ISM_member_id>_<ESM_id>_<forcing_member_id>_<experiment_id>_<set_counter>_<time_range>.nc
	mkfname = @(variable_id) fullfile(outdir,[variable_id '_' icesheetname '_' source_id '_' ism_id '_' ism_member_id '_' esm_id '_' forcing_member_id '_' experiment_id '_' set_counter '_' time_range '.nc']);

	%Scalar variables %{{{
	%ISMIP7 wants annual output, with two different time-registration rules
	%(see https://github.com/ismip/ismip7-time-encoding):
	%  - state variables: instantaneous snapshot registered at the END of the year
	%  - flux variables:  yearly AVERAGE registered at the MIDDLE of the year, with
	%                      time bounds spanning the start to the end of that year
	%NOTE: grouping is by floor(TransientSolution.time), so .time is assumed to be
	%      in (absolute) decimal years - see the time-encoding note further below.
	Nsol    = numel(md.results.TransientSolution);
	alltime = arrayfun(@(k) md.results.TransientSolution(k).time, 1:Nsol);
	if strcmp(experiment,'init'),
		uyears  = floor(alltime(end));   %single snapshot: treat it as its own "year"
		yearidx = {Nsol};
	else
		years   = floor(alltime);
		uyears  = unique(years);
		yearidx = cell(1,numel(uyears));
		for j=1:numel(uyears),
			yearidx{j} = find(years==uyears(j));   %all solutions in this calendar year
		end
	end
	nyears = numel(uyears);

	%STATE annual index: last solution stored in each year (end-of-year snapshot)
	state_idx = zeros(1,nyears);
	for j=1:nyears,
		kk = yearidx{j};
		state_idx(j) = kk(end);
	end

	%--- state scalar diagnostics: instantaneous value at the end-of-year snapshot ---
	results.mass         = zeros(1,nyears);
	results.massaf        = zeros(1,nyears);
	results.groundedarea  = zeros(1,nyears);
	results.floatingarea  = zeros(1,nyears);
	for i=1:nyears,
		k=state_idx(i);
		results.mass(i)        = md.results.TransientSolution(k).IceVolume*md.materials.rho_ice;
		results.massaf(i)       = md.results.TransientSolution(k).IceVolumeAboveFloatation*md.materials.rho_ice;
		results.groundedarea(i) = md.results.TransientSolution(k).GroundedArea;
		results.floatingarea(i) = md.results.TransientSolution(k).FloatingArea;
	end

	%--- flux scalar diagnostics: yearly average over every solution stored within
	%    that calendar year (Gt/yr -> kg/s where relevant) ---
	raw_smbtot       = arrayfun(@(k) md.results.TransientSolution(k).TotalSmb*10^12/md.constants.yts, 1:Nsol);
	raw_bmbgr        = zeros(1,Nsol);     %grounded basal mass balance not in this model setup - set to 0
	raw_bmbfl        = zeros(1,Nsol);     %floating basal mass balance not in this model setup - set to 0
	%TODO: not provided by this model setup - populate before submission (kg/s)
	raw_calvtot      = NaN*ones(1,Nsol);  %tendlicalvf   (total calving flux)
	raw_frontmelttot = NaN*ones(1,Nsol);  %tendlifmassbf (total ice-front melt flux)
	raw_gltot        = NaN*ones(1,Nsol);  %tendligroundf (total grounding-line flux)

	results.smbtot       = zeros(1,nyears);
	results.bmbgr         = zeros(1,nyears);
	results.bmbfl         = zeros(1,nyears);
	results.calvtot       = NaN*ones(1,nyears);
	results.frontmelttot  = NaN*ones(1,nyears);
	results.gltot         = NaN*ones(1,nyears);
	for j=1:nyears,
		kk = yearidx{j};
		results.smbtot(j)       = mean(raw_smbtot(kk));
		results.bmbgr(j)        = mean(raw_bmbgr(kk));
		results.bmbfl(j)        = mean(raw_bmbfl(kk));
		results.calvtot(j)      = mean(raw_calvtot(kk));
		results.frontmelttot(j) = mean(raw_frontmelttot(kk));
		results.gltot(j)        = mean(raw_gltot(kk));
	end

	%ISMIP7 default netCDF4 f4 fill value
	fillval = single(9.969209968386869e+36);

	%Time coordinate in "days since 1850-01-01" (standard/Gregorian calendar).
	%NOTE: this assumes md.results.TransientSolution(k).time holds ABSOLUTE calendar
	%      years. If your model stores time relative to the run start, add the start
	%      year (e.g. parsed from time_range) before conversion. The reference tool at
	%      https://github.com/ismip/ismip7-time-encoding should be used to generate the
	%      day values and the accompanying time bounds for the official submission.
	days_since_1850 = @(yr) datenum(floor(yr),1,1) - datenum(1850,1,1) + (yr-floor(yr)).*(datenum(floor(yr)+1,1,1)-datenum(floor(yr),1,1));

	%STATE time: actual model time of the end-of-year snapshot
	tvec_state = arrayfun(@(k) md.results.TransientSolution(k).time, state_idx);
	time_state = days_since_1850(tvec_state);

	%FLUX time: middle of the calendar year, with bounds spanning the whole year.
	%(for the 'init' single-snapshot case there is no averaging period, so the
	%flux time collapses to the same instant used for the state variables.)
	if strcmp(experiment,'init'),
		time_flux      = time_state;
		time_flux_bnds = [time_state; time_state];
	else
		time_flux      = zeros(1,nyears);
		time_flux_bnds = zeros(2,nyears);
		for j=1:nyears,
			yr = uyears(j);
			b0 = datenum(yr,1,1)  -datenum(1850,1,1);   %start of the year
			b1 = datenum(yr+1,1,1)-datenum(1850,1,1);   %end of the year (=start of next year)
			time_flux_bnds(:,j) = [b0;b1];
			time_flux(j)        = (b0+b1)/2;             %middle of the year
		end
	end

	%One file per scalar variable (ISMIP7 requires a single main variable per file).
	%Columns: variable_id, standard_name, units, data, is_flux
	%  is_flux=false (ST): end-of-year snapshot, no time bounds
	%  is_flux=true  (FL): yearly average at mid-year, with time bounds
	scalar_vars = {
		'lim',             'land_ice_mass',                          'kg',    results.mass,         false;
		'limnsw',          'land_ice_mass_not_displacing_sea_water', 'kg',    results.massaf,       false;
		'iareagr',         'grounded_ice_sheet_area',                'm^2',   results.groundedarea, false;
		'iareafl',         'floating_ice_shelf_area',                'm^2',   results.floatingarea, false;
		'tendacabf',       'tendency_of_land_ice_mass_due_to_surface_mass_balance', 'kg s-1', results.smbtot,       true;
		'tendlibmassbfgr', 'tendency_of_land_ice_mass_due_to_basal_mass_balance',   'kg s-1', results.bmbgr,        true;
		'tendlibmassbffl', 'tendency_of_land_ice_mass_due_to_basal_mass_balance',   'kg s-1', results.bmbfl,        true;
		'tendlicalvf',     'tendency_of_land_ice_mass_due_to_calving',              'kg s-1', results.calvtot,      true;
		'tendlifmassbf',   '',                                                      'kg s-1', results.frontmelttot, true;
		'tendligroundf',   '',                                                      'kg s-1', results.gltot,        true;
	};

	for v=1:size(scalar_vars,1),
		if scalar_vars{v,5},
			t=time_flux; b=time_flux_bnds;
		else
			t=time_state; b=[];
		end
		write_scalar_var(mkfname(scalar_vars{v,1}),scalar_vars{v,1},scalar_vars{v,2},scalar_vars{v,3},scalar_vars{v,4},t,b,fillval,source_id,ism_id,icesheetname);
	end
	%}}}

	%Field variables {{{
	%Output periods: subsample the annual STATE index list every output_interval_yr
	%years; each output period covers the block of calendar years between two
	%consecutive selected years (a single year when output_interval_yr==1).
	sel      = 1:output_interval_yr:nyears;
	noutput  = numel(sel);
	results.timegrid = state_idx(sel);   %STATE snapshot index for each output period (end-of-year)
	blockidx = cell(1,noutput);          %raw solution indices covered by each output period (for flux averaging)
	block_y0 = zeros(1,noutput);         %first calendar year in each output period
	block_y1 = zeros(1,noutput);         %last calendar year in each output period
	for j=1:noutput,
		if j==1,
			yfirst = 1;
		else
			yfirst = sel(j-1)+1;
		end
		ylast = sel(j);
		block_y0(j) = uyears(yfirst);
		block_y1(j) = uyears(ylast);
		blockidx{j} = [yearidx{yfirst:ylast}];
	end
	%Some setups do not save the spatial sub-shelf melt rate. If absent, treat floating
	%basal melt (libmassbffl) as 0 rather than erroring.
	has_bmbfl = isfield(md.results.TransientSolution,'BasalforcingsFloatingiceMeltingRate');
	if ~has_bmbfl,
		warning('ISMIP7:nobmbfl','BasalforcingsFloatingiceMeltingRate not found in the solutions - gridded floating basal melt (libmassbffl) treated as 0.');
	end
	%ISMIP7 standard grids (ISMIP6 domains, EPSG:3413 for GrIS, EPSG:3031 for AIS).
	%Allowed resolutions (multiples of 2 km): GrIS 1/2/4/8/16 km, AIS 2/4/8/16 km.
	posting = resolution_km*1000;
	if strcmp(icesheetname,'GrIS'),
		%Cell centres from (-720000,-3450000) to (960000,-570000); 1681 x 2881 at 1 km
		results.gridx  = -720000 :posting: 960000;
		results.gridy  = -3450000:posting:-570000;
	elseif strcmp(icesheetname,'AIS'),
		%Domain corners (-3040000,-3040000) to (3040000,3040000); 761 x 761 at 8 km
		results.gridx  = -3040000:posting:3040000;
		results.gridy  = -3040000:posting:3040000;
	else
		error('ice sheet not supported yet');
	end

	%Parameters for InterpFromMeshToGrid (index,x,y,data,xgrid,ygrid,default_value).
	%The interpolated grid starts at (xmin,ymax) and steps down in y, so the stored
	%arrays are ordered (x, y-descending); the y coordinate variable is written to match.
	xmin   = results.gridx(1);
	ymax   = results.gridy(end);
	ncols  = numel(results.gridx);
	nlines = numel(results.gridy);
	xgrid  = results.gridx;                      %x grid values (ascending)
	ygrid  = results.gridy(end:-1:1);            %y grid values (descending from ymax)
	results.xcoord = results.gridx;              %x coordinate variable (ascending)
	results.ycoord = results.gridy(end:-1:1);    %y coordinate variable (descending, matches data)

	if strcmp(icesheetname,'GrIS'),
		results.thickness= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.surface= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.bed= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.geoflux= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.smb= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.bmb= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.dhdt= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.vxsurf= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.vysurf= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		%results.vzsurf= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.vxbase= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.vybase= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		%results.vzbase= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.vxmean= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.vymean= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.surftemp= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.basetemp= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.drag= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.calving= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.groundedice= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.floatingice= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.mask= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
	elseif strcmp(icesheetname,'AIS'),
		results.thickness= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.surface= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.bed= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.geoflux= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.smb= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.bmb= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.dhdt= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.vxsurf= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.vysurf= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.vzsurf= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.vxbase= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.vybase= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.vzbase= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.vxmean= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.vymean= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.surftemp= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.basetemp= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.drag= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.calving= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.groundedice= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.floatingice= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
		results.mask= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid)); %y,x,time
	else
		error('ice sheet not supported yet');
	end

	%Extra mandatory ISMIP7 fields handled outside the per-domain loops
	results.base       = NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid));
	results.sftgif     = NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid));
	%Mandatory, but not produced by this converter - written as fill values (TODO: populate before submission)
	results.libmassbfgr= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid));
	results.ligroundf  = NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid));
	results.lifmassbf  = NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid));

	if strcmp(icesheetname,'GrIS'),
		for i=1:length(results.timegrid),
			thickness=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Thickness(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
			results.thickness(:,:,i)=transpose(thickness);
			base=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Base(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
			results.base(:,:,i)=transpose(base);
			surface=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Surface(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
			results.surface(:,:,i)=transpose(surface);
			bed=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.bed(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
			results.bed(:,:,i)=transpose(bed);
			geoflux=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.basalforcings.geothermalflux(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
			results.geoflux(:,:,i)=transpose(geoflux);
			if i==1,
				dhdt=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,(md.results.TransientSolution(results.timegrid(i)).Thickness(1:md.mesh.numberofvertices)-md.geometry.thickness(1:md.mesh.numberofvertices))/(md.results.TransientSolution(results.timegrid(i)).time-0),xgrid,ygrid,NaN);
			else
				dhdt=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,(md.results.TransientSolution(results.timegrid(i)).Thickness(1:md.mesh.numberofvertices)-md.results.TransientSolution(results.timegrid(i-1)).Thickness(1:md.mesh.numberofvertices))/(md.results.TransientSolution(results.timegrid(i)).time-md.results.TransientSolution(results.timegrid(i-1)).time),xgrid,ygrid,NaN);
			end
			results.dhdt(:,:,i)=transpose(dhdt)/(md.constants.yts);
			vxsurf=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vx(end-md.mesh.numberofvertices+1:end),xgrid,ygrid,NaN);
			results.vxsurf(:,:,i)=transpose(vxsurf)/md.constants.yts;
			vysurf=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vy(end-md.mesh.numberofvertices+1:end),xgrid,ygrid,NaN);
			results.vysurf(:,:,i)=transpose(vysurf)/md.constants.yts;
			%vzsurf=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vz(end-md.mesh.numberofvertices+1:end),xgrid,ygrid,NaN);
			%results.vzsurf(:,:,i)=transpose(vzsurf)/md.constants.yts;
			vxbase=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vx(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
			results.vxbase(:,:,i)=transpose(vxbase)/md.constants.yts;
			vybase=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vy(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
			results.vybase(:,:,i)=transpose(vybase)/md.constants.yts;
			%vzbase=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vz(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
			%results.vzbase(:,:,i)=transpose(vzbase)/md.constants.yts;
			vxmean=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vx,xgrid,ygrid,NaN);
         results.vxmean(:,:,i)=transpose(vxmean)/md.constants.yts;
			vymean=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vy,xgrid,ygrid,NaN);
			results.vymean(:,:,i)=transpose(vymean)/md.constants.yts;
			surftemp=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.initialization.temperature(end-md.mesh.numberofvertices+1:end),xgrid,ygrid,NaN);
			results.surftemp(:,:,i)=transpose(surftemp);
			basetemp=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.initialization.temperature(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
			results.basetemp(:,:,i)=transpose(basetemp);
			drag=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.friction.coefficient(1:md.mesh.numberofvertices).^2*md.constants.g.*(md.materials.rho_ice*md.results.TransientSolution(results.timegrid(i)).Thickness(1:md.mesh.numberofvertices)+md.materials.rho_water*md.results.TransientSolution(results.timegrid(i)).Base(1:md.mesh.numberofvertices)).*sqrt(md.results.TransientSolution(results.timegrid(i)).Vx(1:md.mesh.numberofvertices).^2+md.results.TransientSolution(results.timegrid(i)).Vy(1:md.mesh.numberofvertices).^2)/md.constants.yts,xgrid,ygrid,NaN);
			results.drag(:,:,i)=transpose(drag);
			calving=NaN*ones(nlines,ncols);
			results.calving(:,:,i)=transpose(calving);
			mask=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,-md.mask.ice_levelset(1:md.mesh.numberofvertices),xgrid,ygrid,-1);
			mask(find(mask>0))=1;
			mask(find(mask<0))=0;
			results.mask(:,:,i)=transpose(mask);
			groundedice=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).MaskOceanLevelset(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
			groundedice(find(groundedice>0))=1;
			groundedice(find(groundedice<0))=0;
			floatingice=1-groundedice;
			groundedice(find(isnan(groundedice)))=0;
			floatingice(find(isnan(floatingice)))=0;
			results.groundedice(:,:,i)=transpose(groundedice.*mask);
			results.floatingice(:,:,i)=transpose(floatingice.*mask);
			%--- flux fields for this output period: averaged over every raw solution
			%    stored within it, rather than a single end-of-period snapshot ---
			bidx = blockidx{i};
			nb   = numel(bidx);
			meltfl_sum = zeros(md.mesh.numberofvertices,1);
			smb_sum    = zeros(md.mesh.numberofvertices,1);
			for ib=1:nb,
				kb = bidx(ib);
				if has_bmbfl,
					meltfl_sum = meltfl_sum + md.results.TransientSolution(kb).BasalforcingsFloatingiceMeltingRate(1:md.mesh.numberofvertices);
				end
				smb_kb = md.results.TransientSolution(kb).SmbMassBalance(end-md.mesh.numberofvertices+1:end);
				smb_kb(find(smb_kb<-1000))=0;   %clip spurious values (as in the original snapshot approach)
				smb_sum = smb_sum + smb_kb;
			end
			if has_bmbfl,
				meltfl = meltfl_sum/nb;
			else
				meltfl = zeros(md.mesh.numberofvertices,1);   %floating basal melt not saved - treated as 0
			end
			bmb=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,meltfl,xgrid,ygrid,NaN);
			results.bmb(:,:,i)=-(transpose(bmb)*md.materials.rho_ice/md.constants.yts).*transpose(1-groundedice);
			smb=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,smb_sum/nb,xgrid,ygrid,NaN);
			results.smb(:,:,i)=transpose(smb)*md.materials.rho_ice/md.constants.yts;
		end
	elseif strcmp(icesheetname,'AIS'),
		for i=1:length(results.timegrid),
			thickness=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Thickness,xgrid,ygrid,NaN);
			results.thickness(:,:,i)=transpose(thickness);
			base=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Base,xgrid,ygrid,NaN);
			results.base(:,:,i)=transpose(base);
			surface=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Surface,xgrid,ygrid,NaN);
			results.surface(:,:,i)=transpose(surface);
			bed=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.bed,xgrid,ygrid,NaN);
			results.bed(:,:,i)=transpose(bed);
			geoflux=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.basalforcings.geothermalflux,xgrid,ygrid,NaN);
			results.geoflux(:,:,i)=transpose(geoflux);
			if i==1,
				dhdt=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,(md.results.TransientSolution(results.timegrid(i)).Thickness-md.geometry.thickness)/(md.results.TransientSolution(results.timegrid(i)).time-0),xgrid,ygrid,NaN);
			else
				dhdt=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,(md.results.TransientSolution(results.timegrid(i)).Thickness-md.results.TransientSolution(results.timegrid(i-1)).Thickness)/(md.results.TransientSolution(results.timegrid(i)).time-md.results.TransientSolution(results.timegrid(i-1)).time),xgrid,ygrid,NaN);
			end
			results.dhdt(:,:,i)=transpose(dhdt)/(md.constants.yts);
			vxsurf=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vx,xgrid,ygrid,NaN);
			results.vxsurf(:,:,i)=transpose(vxsurf)/md.constants.yts;
			vysurf=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vy,xgrid,ygrid,NaN);
			results.vysurf(:,:,i)=transpose(vysurf)/md.constants.yts;
			%vzsurf=zeros(nlines,ncols);
			%vzsurf=NaN*ones(nlines,ncols);
			%results.vzsurf(:,:,i)=transpose(vzsurf)/md.constants.yts;
			vxbase=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vx,xgrid,ygrid,NaN);
			results.vxbase(:,:,i)=transpose(vxbase)/md.constants.yts;
			vybase=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vy,xgrid,ygrid,NaN);
			results.vybase(:,:,i)=transpose(vybase)/md.constants.yts;
			%vzbase=NaN*ones(nlines,ncols);
			%results.vzbase(:,:,i)=transpose(vzbase)/md.constants.yts;
			vxmean=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vx,xgrid,ygrid,NaN);
			results.vxmean(:,:,i)=transpose(vxmean)/md.constants.yts;
			vymean=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vy,xgrid,ygrid,NaN);
			results.vymean(:,:,i)=transpose(vymean)/md.constants.yts;
			surftemp=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.initialization.temperature,xgrid,ygrid,NaN);
			results.surftemp(:,:,i)=transpose(surftemp);
			basetemp=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.initialization.temperature,xgrid,ygrid,NaN);
			basetemp=NaN*ones(nlines,ncols);
			results.basetemp(:,:,i)=transpose(basetemp);
			drag=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.friction.coefficient.^2*md.constants.g.*(md.materials.rho_ice*md.results.TransientSolution(results.timegrid(i)).Thickness+md.materials.rho_water*md.results.TransientSolution(results.timegrid(i)).Base).*md.results.TransientSolution(results.timegrid(i)).Vel/md.constants.yts,xgrid,ygrid,NaN);
			results.drag(:,:,i)=transpose(drag);
			calving=NaN*ones(nlines,ncols);
			results.calving(:,:,i)=transpose(calving);
			groundedice=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).MaskGroundediceLevelset,xgrid,ygrid,NaN);
			groundedice(find(groundedice>0))=1;
			groundedice(find(groundedice<0))=0;
			results.groundedice(:,:,i)=transpose(groundedice);
			results.floatingice(:,:,i)=transpose(1-groundedice);
			mask=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,-md.mask.ice_levelset,xgrid,ygrid,-1);
			mask(find(mask>0))=1;
			mask(find(mask<0))=-1;
			results.mask(:,:,i)=transpose(mask);
			%--- flux fields for this output period: averaged over every raw solution
			%    stored within it, rather than a single end-of-period snapshot ---
			bidx = blockidx{i};
			nb   = numel(bidx);
			meltfl_sum = zeros(md.mesh.numberofvertices,1);
			smb_sum    = zeros(md.mesh.numberofvertices,1);
			for ib=1:nb,
				kb = bidx(ib);
				if has_bmbfl,
					meltfl_sum = meltfl_sum + md.results.TransientSolution(kb).BasalforcingsFloatingiceMeltingRate(1:md.mesh.numberofvertices);
				end
				smb_sum = smb_sum + md.results.TransientSolution(kb).SmbMassBalance;
			end
			if has_bmbfl,
				meltfl = meltfl_sum/nb;
			else
				meltfl = zeros(md.mesh.numberofvertices,1);   %floating basal melt not saved - treated as 0
			end
			bmb=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,meltfl,xgrid,ygrid,NaN);
			results.bmb(:,:,i)=-(transpose(bmb)*md.materials.rho_ice/md.constants.yts).*transpose(1-groundedice);
			smb=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,smb_sum/nb,xgrid,ygrid,NaN);
			results.smb(:,:,i)=transpose(smb)*md.materials.rho_ice/md.constants.yts;
		end
	else
		error('ice sheet not supported yet');
	end
	%}}}

	%Gridded output: derive mask fractions, build the gridded time axes, then write
	%one file per mandatory variable (ISMIP7 requires a single main variable per file).
	icepresent = double(results.mask>0);                %1 where ice present, 0 elsewhere (both domains)
	results.sftgif = icepresent;                         %land ice area fraction
	results.sftgrf = results.groundedice .* icepresent;  %grounded fraction, 0 outside ice
	results.sftflf = results.floatingice .* icepresent;  %floating fraction, 0 outside ice

	%STATE grid time: end-of-year snapshot (actual model time)
	tvec_grid_state = arrayfun(@(k) md.results.TransientSolution(k).time, results.timegrid);
	time_grid_state = days_since_1850(tvec_grid_state);

	%FLUX grid time: middle of each output period, with bounds spanning its start
	%to its end (a single calendar year when output_interval_yr==1).
	if strcmp(experiment,'init'),
		time_grid_flux      = time_grid_state;
		time_grid_flux_bnds = [time_grid_state;time_grid_state];
	else
		time_grid_flux      = zeros(1,noutput);
		time_grid_flux_bnds = zeros(2,noutput);
		for j=1:noutput,
			b0 = datenum(block_y0(j),1,1)  -datenum(1850,1,1);
			b1 = datenum(block_y1(j)+1,1,1)-datenum(1850,1,1);
			time_grid_flux_bnds(:,j) = [b0;b1];
			time_grid_flux(j)        = (b0+b1)/2;
		end
	end

	warning('ISMIP7:placeholder',['These mandatory variables are not produced by this ' ...
		'converter and are written as fill/NaN - populate before submission: libmassbfgr, ' ...
		'licalvf, ligroundf, lifmassbf (gridded); tendlicalvf, tendlifmassbf, tendligroundf (scalar).']);

	%Columns: variable_id, standard_name, units, data, is_mask, is_flux
	%  is_flux=false (ST): end-of-year snapshot, no time bounds
	%  is_flux=true  (FL): yearly average at mid-year, with time bounds
	grid_vars = {
		'lithk',       'land_ice_thickness',                         'm',          results.thickness,   false, false;
		'orog',        'surface_altitude',                           'm',          results.surface,     false, false;
		'topg',        'bedrock_altitude',                           'm',          results.bed,         false, false;
		'base',        '',                                           'm',          results.base,        false, false;
		'acabf',       'land_ice_surface_specific_mass_balance_flux','kg m-2 s-1', results.smb,         false, true;
		'libmassbfgr', 'land_ice_basal_specific_mass_balance_flux',  'kg m-2 s-1', results.libmassbfgr, false, true;
		'libmassbffl', 'land_ice_basal_specific_mass_balance_flux',  'kg m-2 s-1', results.bmb,         false, true;
		'dlithkdt',    'tendency_of_land_ice_thickness',             'm s-1',      results.dhdt,        false, true;
		'xvelmean',    'land_ice_vertical_mean_x_velocity',          'm s-1',      results.vxmean,      false, false;
		'yvelmean',    'land_ice_vertical_mean_y_velocity',          'm s-1',      results.vymean,      false, false;
		'strbasemag',  'land_ice_basal_drag',                        'Pa',         results.drag,        false, false;
		'licalvf',     'land_ice_specific_mass_flux_due_to_calving', 'kg m-2 s-1', results.calving,     false, true;
		'ligroundf',   '',                                           'kg m-2 s-1', results.ligroundf,   false, true;
		'lifmassbf',   '',                                           'kg m-2 s-1', results.lifmassbf,   false, true;
		'sftgif',      'land_ice_area_fraction',                     '1',          results.sftgif,      true,  false;
		'sftgrf',      'grounded_ice_sheet_area_fraction',           '1',          results.sftgrf,      true,  false;
		'sftflf',      'floating_ice_shelf_area_fraction',           '1',          results.sftflf,      true,  false;
	};

	for v=1:size(grid_vars,1),
		if grid_vars{v,6},
			t=time_grid_flux; b=time_grid_flux_bnds;
		else
			t=time_grid_state; b=[];
		end
		write_gridded_var(mkfname(grid_vars{v,1}),grid_vars{v,1},grid_vars{v,2},grid_vars{v,3},grid_vars{v,4},results.xcoord,results.ycoord,t,b,fillval,grid_vars{v,5},source_id,ism_id,icesheetname);
	end

function ismip7_global_attributes(ncid,source_id,ism_id,domain_id)
	%Mandatory ISMIP7 global attributes (see conventions section 5)
	if strcmp(domain_id,'GrIS'),
		crs='epsg:3413';
	else
		crs='epsg:3031';
	end
	g=netcdf.getConstant('NC_GLOBAL');
	netcdf.putAtt(ncid,g,'group',source_id);
	netcdf.putAtt(ncid,g,'model',ism_id);
	netcdf.putAtt(ncid,g,'contact_name','Hugh Shields');                 %EDIT: set the actual contact name(s)
	netcdf.putAtt(ncid,g,'contact_email','nikolaus.h.shields.gr@dartmouth.edu');   %EDIT: set the actual contact email(s)
	netcdf.putAtt(ncid,g,'crs',crs);

function write_scalar_var(fname,variable_id,standard_name,units,data,time_days,time_bnds,fillval,source_id,ism_id,domain_id)
	%Write one scalar (time-only) ISMIP7 variable to its own netCDF4 file.
	%time_bnds: [] for state variables (instantaneous snapshot), or a 2 x ntime
	%array of [start;end] bounds (days since 1850) for flux variables (yearly
	%average registered at the middle of the year).
	mode = netcdf.getConstant('NETCDF4');
	mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
	ncid = netcdf.create(fname,mode);
	ismip7_global_attributes(ncid,source_id,ism_id,domain_id);
	ntime_id = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
	time_var_id = netcdf.defVar(ncid,'time','NC_FLOAT',ntime_id);
	netcdf.putAtt(ncid,time_var_id,'standard_name','time');
	netcdf.putAtt(ncid,time_var_id,'units','days since 1850-01-01');
	netcdf.putAtt(ncid,time_var_id,'calendar','standard');
	netcdf.putAtt(ncid,time_var_id,'axis','T');
	has_bnds = ~isempty(time_bnds);
	if has_bnds,
		netcdf.putAtt(ncid,time_var_id,'bounds','time_bnds');
		nv_id = netcdf.defDim(ncid,'nv',2);
		bnds_var_id = netcdf.defVar(ncid,'time_bnds','NC_FLOAT',[nv_id ntime_id]);
	end
	var_id = netcdf.defVar(ncid,variable_id,'NC_FLOAT',ntime_id);
	if ~isempty(standard_name),
		netcdf.putAtt(ncid,var_id,'standard_name',standard_name);
	end
	netcdf.putAtt(ncid,var_id,'units',units);
	netcdf.putAtt(ncid,var_id,'_FillValue',fillval);
	netcdf.endDef(ncid);
	d = single(data);
	d(isnan(d)) = fillval;
	netcdf.putVar(ncid,time_var_id,0,numel(time_days),single(time_days));
	netcdf.putVar(ncid,var_id,0,numel(d),d);
	if has_bnds,
		netcdf.putVar(ncid,bnds_var_id,[0 0],[2 numel(time_days)],single(time_bnds));
	end
	netcdf.close(ncid);

function write_gridded_var(fname,variable_id,standard_name,units,data,xcoord,ycoord,time_days,time_bnds,fillval,is_mask,source_id,ism_id,domain_id)
	%Write one gridded (x,y,t) ISMIP7 variable to its own netCDF4 file.
	%Data is expected ordered (x, y) with y matching ycoord (descending).
	%time_bnds: [] for state variables (instantaneous snapshot), or a 2 x ntime
	%array of [start;end] bounds (days since 1850) for flux variables (yearly
	%average registered at the middle of the year).
	mode = netcdf.getConstant('NETCDF4');
	mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
	ncid = netcdf.create(fname,mode);
	ismip7_global_attributes(ncid,source_id,ism_id,domain_id);
	nx_id    = netcdf.defDim(ncid,'x',numel(xcoord));
	ny_id    = netcdf.defDim(ncid,'y',numel(ycoord));
	ntime_id = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
	x_var_id = netcdf.defVar(ncid,'x','NC_FLOAT',nx_id);
	netcdf.putAtt(ncid,x_var_id,'standard_name','projection_x_coordinate');
	netcdf.putAtt(ncid,x_var_id,'units','m');
	netcdf.putAtt(ncid,x_var_id,'axis','X');
	y_var_id = netcdf.defVar(ncid,'y','NC_FLOAT',ny_id);
	netcdf.putAtt(ncid,y_var_id,'standard_name','projection_y_coordinate');
	netcdf.putAtt(ncid,y_var_id,'units','m');
	netcdf.putAtt(ncid,y_var_id,'axis','Y');
	time_var_id = netcdf.defVar(ncid,'time','NC_FLOAT',ntime_id);
	netcdf.putAtt(ncid,time_var_id,'standard_name','time');
	netcdf.putAtt(ncid,time_var_id,'units','days since 1850-01-01');
	netcdf.putAtt(ncid,time_var_id,'calendar','standard');
	netcdf.putAtt(ncid,time_var_id,'axis','T');
	has_bnds = ~isempty(time_bnds);
	if has_bnds,
		netcdf.putAtt(ncid,time_var_id,'bounds','time_bnds');
		nv_id = netcdf.defDim(ncid,'nv',2);
		bnds_var_id = netcdf.defVar(ncid,'time_bnds','NC_FLOAT',[nv_id ntime_id]);
	end
	var_id = netcdf.defVar(ncid,variable_id,'NC_FLOAT',[nx_id ny_id ntime_id]);
	if ~isempty(standard_name),
		netcdf.putAtt(ncid,var_id,'standard_name',standard_name);
	end
	netcdf.putAtt(ncid,var_id,'units',units);
	netcdf.putAtt(ncid,var_id,'_FillValue',fillval);
	netcdf.endDef(ncid);
	d = single(data);
	if is_mask,
		d(isnan(d)) = 0;                %ice masks: set missing to 0, not fill (ISMIP7 convention)
	else
		d(isnan(d)) = fillval;
	end
	netcdf.putVar(ncid,x_var_id,single(xcoord));
	netcdf.putVar(ncid,y_var_id,single(ycoord));
	netcdf.putVar(ncid,time_var_id,0,numel(time_days),single(time_days));
	netcdf.putVar(ncid,var_id,[0 0 0],[numel(xcoord) numel(ycoord) numel(time_days)],d);
	if has_bnds,
		netcdf.putVar(ncid,bnds_var_id,[0 0],[2 numel(time_days)],single(time_bnds));
	end
	netcdf.close(ncid);
