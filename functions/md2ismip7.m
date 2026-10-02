function results=md2ismip7(md,directoryname,source_id,ism_id,ism_member_id,forcing_member_id,set_counter_override,time_range,resolution_km,output_interval_yr);
	%Create netcdf files for an experiment following ISMIP7 conventions, for Greenland (GrIS) only.
	%
	%directoryname:        base output directory (ISMIP7 subdirectories are created below it)
	%source_id:            modelling group name (no '_', '.' or special characters)
	%                      (optional, default 'Dartmouth')
	%ism_id:               ice sheet model name and version (no '_', '.' or special characters)
	%                      (optional, default 'ISSM')
	%ism_member_id:        ISM choice variant, e.g. 'm001' (optional, default 'm001')
	%forcing_member_id:    forcing choice variant, e.g. 'f001' (optional, default 'f001')
	%set_counter_override: ISMIP7 set_counter (e.g. 'E001', 'P001') to use INSTEAD of the Core
	%                      experiment automatically matched below (optional; leave empty/omitted
	%                      to default to the matching Core run, e.g. 'C003')
	%time_range:           start-end year of the experiment, e.g. '2015-2300' (optional;
	%                      auto-derived from the model output years if empty/omitted)
	%resolution_km:        output grid resolution in km (optional, default 1).
	%                      Allowed: 1/2/4/8/16 km.
	%output_interval_yr:   interval in years between 2-D (gridded) output snapshots
	%                      (optional, default 1). Scalars are always annual.
	%
	%esm_id and experiment_id are no longer passed in directly - they are parsed from
	%md.miscellaneous.name, which every model input to this function must set, following:
	%   Greenland_ISMIP7Run_<esm_id>_<experiment_id>[_<anything else, e.g. a year range>]
	%   Greenland_ISMIP7Prep_OCX                       (the ocean-forcing prep run; no experiment_id)
	%e.g. Greenland_ISMIP7Run_CESM2-WACCM_SSP370_2015-2101
	%     Greenland_ISMIP7Run_MRI-ESM2-0_Historical
	%     Greenland_ISMIP7Prep_OCX
	%
	%set_counter is likewise no longer passed in: by default this function assumes a Core
	%experiment and looks up the matching C0xx code from esm_id/experiment_id using the table
	%below; pass set_counter_override to use an ESM-sensitivity ('Exxx') or PPE ('Pxxx') code
	%instead.
	%  Core  Experiment   Start year  End year  ESM
	%  C001  Historical   >=1850      2014      CESM2-WACCM
	%  C002  Historical   >=1850      2014      MRI-ESM2-0
	%  C003  SSP370       2015        2100      CESM2-WACCM
	%  C004  SSP370       2015        2100      MRI-ESM2-0
	%  C005  SSP126       2015        2300      CESM2-WACCM
	%  C006  SSP126       2015        2300      MRI-ESM2-0
	%  C007  SSP585       2015        2300      CESM2-WACCM
	%  C008  SSP585       2015        2300      MRI-ESM2-0
	%  C009  CTRL2015     >=1850      2300      CESM2-WACCM
	%  C010  CTRL2015     >=1850      2300      MRI-ESM2-0
	%  C011  OCX          1990-2015   2025      -

	if nargin<3 || isempty(source_id),
		source_id='Dartmouth';
	end
	if nargin<4 || isempty(ism_id),
		ism_id='ISSM';
	end
	if nargin<5 || isempty(ism_member_id),
		ism_member_id='m001';
	end
	if nargin<6 || isempty(forcing_member_id),
		forcing_member_id='f001';
	end
	if nargin<7 || isempty(set_counter_override),
		set_counter_override='';
	end
	if nargin<9 || isempty(resolution_km),
		resolution_km=1; %default GrIS resolution
	end
	if nargin<10 || isempty(output_interval_yr),
		output_interval_yr=1;
	end

	if exist(directoryname)~=7,
		error(['directory ' directoryname ' does not exist']);
	end

	%This function now only supports Greenland - the domain_id is fixed rather than passed in.
	icesheetname='GrIS';

	%Parse esm_id and experiment_id out of md.miscellaneous.name (see the header comment above
	%for the expected naming convention). md is an ISSM model object (not a plain struct), so
	%isfield() can't be used to check for the property - try/catch works for both a struct and
	%a classdef object.
	try
		runname = md.miscellaneous.name;
	catch
		error('md.miscellaneous.name must be set, e.g. Greenland_ISMIP7Run_CESM2-WACCM_SSP370');
	end
	if isempty(runname),
		error('md.miscellaneous.name must be set, e.g. Greenland_ISMIP7Run_CESM2-WACCM_SSP370');
	end
	nameparts = strsplit(runname,'_');
	if numel(nameparts)<3,
		error(['md.miscellaneous.name (''' runname ''') does not follow the expected ' ...
			'Greenland_ISMIP7Run_<esm_id>_<experiment_id> or Greenland_ISMIP7Prep_OCX pattern']);
	end
	esm_id = nameparts{3};
	if strcmpi(esm_id,'OCX'),
		experiment_id=''; %the OCX ocean-forcing prep run has no experiment_id
	else
		if numel(nameparts)<4,
			error(['md.miscellaneous.name (''' runname ''') is missing the experiment_id ' ...
				'token (expected Greenland_ISMIP7Run_<esm_id>_<experiment_id>...)']);
		end
		experiment_id = lower(nameparts{4}); %ISMIP7 naming wants this lowercase, e.g. 'historical', 'ssp370'
	end

	%Calving/ice-front-melt diagnostics are only physically meaningful when the
	%front position is dynamically simulated. Historical runs in this workflow
	%prescribe the front (observed extent imposed each year, no calving law
	%active), so calving-derived outputs are left at fill value throughout.
	is_historical = strcmpi(experiment_id,'historical');

	%Allowed ISMIP7 resolutions (km) for GrIS
	allowed_res=[1 2 4 8 16];
	if ~any(resolution_km==allowed_res),
		error(['resolution_km=' num2str(resolution_km) ' km is not allowed for ' icesheetname ' (allowed: ' num2str(allowed_res) ' km)']);
	end

	%set_counter: default to the Core experiment matching esm_id/experiment_id (see table in the
	%header comment above), unless set_counter_override was provided.
	if ~isempty(set_counter_override),
		set_counter = set_counter_override;
	elseif strcmpi(esm_id,'OCX'),
		set_counter = 'C011';
	else
		if ~isempty(regexpi(esm_id,'CESM')),
			esmgroup='CESM';
		elseif ~isempty(regexpi(esm_id,'MRI')),
			esmgroup='MRI';
		else
			error(['esm_id ''' esm_id ''' (parsed from md.miscellaneous.name) is not recognized as a Core ' ...
				'ESM (expected a CESM2-WACCM or MRI-ESM2-0 variant, or OCX) - pass set_counter_override explicitly.']);
		end
		core_table = {
			'historical', 'CESM', 'C001';
			'historical', 'MRI',  'C002';
			'ssp370',     'CESM', 'C003';
			'ssp370',     'MRI',  'C004';
			'ssp126',     'CESM', 'C005';
			'ssp126',     'MRI',  'C006';
			'ssp585',     'CESM', 'C007';
			'ssp585',     'MRI',  'C008';
			'ctrl2015',   'CESM', 'C009';
			'ctrl2015',   'MRI',  'C010';
		};
		match = find(strcmpi(core_table(:,1),experiment_id) & strcmp(core_table(:,2),esmgroup));
		if isempty(match),
			error(['No Core ISMIP7 experiment found for experiment_id=''' experiment_id ''' and esm_id=''' esm_id '''' ...
				' - pass set_counter_override explicitly for a non-Core (ESM/PPE) run.']);
		end
		set_counter = core_table{match,3};
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
	if nargin<8 || isempty(time_range),
		%Drop the last stored solution before deriving the year range: it is the
		%branch/handoff instant shared with whatever comes next (e.g. a historical
		%run's array ending at exactly 2015.0 really just marks "as of Jan 1 2015",
		%i.e. through end of 2014) rather than a genuine extra year of its own.
		Nuse_tmp = max(numel(md.results.TransientSolution)-1,1);
		ytmp = floor(arrayfun(@(k) md.results.TransientSolution(k).time, 1:Nuse_tmp));
		y0=min(ytmp); y1=max(ytmp);
		if y0==y1,
			time_range=sprintf('%04d',y0);          %single year, e.g. 2014
		else
			time_range=sprintf('%04d-%04d',y0,y1);  %year range, e.g. 1850-2014
		end
	end

	%ISMIP7 filename builder:
	%<variable_id>_<domain_id>_<source_id>_<ism_id>_<ISM_member_id>_<ESM_id>_<forcing_member_id>_<experiment_id>_<set_counter>_<time_range>.nc
	%(the experiment_id token is omitted for the OCX run, which has none)
	fname_tokens = {icesheetname,source_id,ism_id,ism_member_id,esm_id,forcing_member_id};
	if ~isempty(experiment_id),
		fname_tokens{end+1} = experiment_id;
	end
	fname_tokens(end+1:end+2) = {set_counter,time_range};
	fname_base = fname_tokens{1};
	for k=2:numel(fname_tokens),
		fname_base = [fname_base '_' fname_tokens{k}];
	end
	mkfname = @(variable_id) fullfile(outdir,[variable_id '_' fname_base '.nc']);

	%Scalar variables %{{{
	%ISMIP7 wants annual output, with two different time-registration rules
	%(see https://github.com/ismip/ismip7-time-encoding):
	%  - state variables: instantaneous snapshot registered at the END of the year
	%  - flux variables:  yearly AVERAGE registered at the MIDDLE of the year, with
	%                      time bounds spanning the start to the end of that year
	%NOTE: grouping is by floor(TransientSolution.time), so .time is assumed to be
	%      in (absolute) decimal years - see the time-encoding note further below.
	%Drop the last stored solution before doing anything else - see the matching
	%note by the time_range derivation above. Everything downstream (state_idx,
	%yearidx, blockidx) derives from Nsol/alltime, so dropping it here keeps it
	%out of every computation, not just out of the filename.
	Nsol    = max(numel(md.results.TransientSolution)-1,1);
	alltime = arrayfun(@(k) md.results.TransientSolution(k).time, 1:Nsol);
	years   = floor(alltime);
	uyears  = unique(years);
	yearidx = cell(1,numel(uyears));
	for j=1:numel(uyears),
		yearidx{j} = find(years==uyears(j));   %all solutions in this calendar year
	end
	nyears = numel(uyears);

	%STATE annual index: last solution stored in each year (end-of-year snapshot)
	state_idx = zeros(1,nyears);
	for j=1:nyears,
		kk = yearidx{j};
		state_idx(j) = kk(end);
	end

	%--- state scalar diagnostics: instantaneous value at the end-of-year snapshot ---
	%Using the *Scaled ISSM diagnostics (already in md.transient.requested_outputs): the
	%unscaled versions count an element as fully in/out of the ice/grounded/floating
	%domain (binary, whole-element), while Scaled weights each element by its actual
	%fractional coverage - more accurate right at the ice margin/grounding line, where
	%elements are often only partially covered.
	results.mass         = zeros(1,nyears);
	results.massaf        = zeros(1,nyears);
	results.groundedarea  = zeros(1,nyears);
	results.floatingarea  = zeros(1,nyears);
	for i=1:nyears,
		k=state_idx(i);
		results.mass(i)        = md.results.TransientSolution(k).IceVolumeScaled*md.materials.rho_ice;
		results.massaf(i)       = md.results.TransientSolution(k).IceVolumeAboveFloatationScaled*md.materials.rho_ice;
		results.groundedarea(i) = md.results.TransientSolution(k).GroundedAreaScaled;
		results.floatingarea(i) = md.results.TransientSolution(k).FloatingAreaScaled;
	end

	%--- flux scalar diagnostics: yearly average over every solution stored within
	%    that calendar year (Gt/yr -> kg/s where relevant) ---
	%Scalar grounded/floating basal mass balance totals - fall back to 0 if these
	%weren't requested (matches the physical 0 already used for grounded bmb, and
	%avoids a hard crash on older/different result sets).
	has_totalbmbgr = isfield(md.results.TransientSolution,'TotalGroundedBmbScaled');
	if ~has_totalbmbgr,
		warning('ISMIP7:nototalbmbgr','TotalGroundedBmbScaled not found in the solutions - tendlibmassbfgr treated as 0.');
	end
	has_totalbmbfl = isfield(md.results.TransientSolution,'TotalFloatingBmbScaled');
	if ~has_totalbmbfl,
		warning('ISMIP7:nototalbmbfl','TotalFloatingBmbScaled not found in the solutions - tendlibmassbffl treated as 0.');
	end
	raw_smbtot       = arrayfun(@(k) md.results.TransientSolution(k).TotalSmbScaled*10^12/md.constants.yts, 1:Nsol);
	if has_totalbmbgr,
		raw_bmbgr = arrayfun(@(k) md.results.TransientSolution(k).TotalGroundedBmbScaled*10^12/md.constants.yts, 1:Nsol);  %expect ~0 - no grounded-ice basal melt in this model setup
	else
		raw_bmbgr = zeros(1,Nsol);
	end
	if has_totalbmbfl,
		raw_bmbfl = arrayfun(@(k) md.results.TransientSolution(k).TotalFloatingBmbScaled*10^12/md.constants.yts, 1:Nsol);
	else
		raw_bmbfl = zeros(1,Nsol);
	end
	%Calving/front-melt totals: TotalCalvingFluxLevelset and TotalCalvingMeltingFluxLevelset
	%follow the same Gt/yr -> kg/s convention as TotalSmb above (same underlying rho_ice-scaled
	%ISSM diagnostic family). TotalCalvingMeltingFluxLevelset is the COMBINED calving+melt flux,
	%so the melt-only term is isolated by subtracting the calving-only total from it.
	%Left at fill value for historical runs, where the front is prescribed (see is_historical above).
	%NOTE: both Total* diagnostics use ISSM's 'positive = mass flux OUT through the front'
	%convention (the same one compute_calvingflux mirrored before being negated above) - negate
	%here too so tendlicalvf/tendlifmassbf (negative = mass lost) stay consistent with the sign
	%of licalvf/lifmassbf, the gridded fields they are the spatial integral of.
	if ~is_historical,
		raw_calvtot      = arrayfun(@(k) -md.results.TransientSolution(k).TotalCalvingFluxLevelset*10^12/md.constants.yts, 1:Nsol);
		raw_frontmelttot = arrayfun(@(k) -(md.results.TransientSolution(k).TotalCalvingMeltingFluxLevelset-md.results.TransientSolution(k).TotalCalvingFluxLevelset)*10^12/md.constants.yts, 1:Nsol);
	else
		raw_calvtot      = NaN*ones(1,Nsol);  %tendlicalvf   - fill value (front prescribed)
		raw_frontmelttot = NaN*ones(1,Nsol);  %tendlifmassbf - fill value (front prescribed)
	end
	%tendligroundf: GroundinglineMassFlux follows the same Gt/yr -> kg/s convention as
	%TotalSmb/TotalCalvingFluxLevelset above (same underlying rho_ice-scaled ISSM
	%diagnostic family). Not tied to the prescribed-front issue, so filled in every
	%experiment including historical.
	raw_gltot = arrayfun(@(k) md.results.TransientSolution(k).GroundinglineMassFlux*10^12/md.constants.yts, 1:Nsol);

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

	%Time coordinate is "days since 1850-01-01" (standard/Gregorian calendar),
	%computed below via datenum(...)-datenum(1850,1,1). The reference tool at
	%https://github.com/ismip/ismip7-time-encoding should be used to generate
	%the day values and the accompanying time bounds for the official submission.

	%STATE time: EXACT end-of-year instant (Jan 1 of the year after the nominal
	%year), per the ISMIP7 ST convention. NOTE: this is deliberately NOT the raw
	%model save time - ISSM's last solution of a year is typically saved a few
	%days before the exact year boundary, which the compliance checker treats as
	%a different (wrong) nominal year rather than a rounding error.
	time_state = zeros(1,nyears);
	for j=1:nyears,
		time_state(j) = datenum(uyears(j)+1,1,1) - datenum(1850,1,1);   %Jan 1 of year+1
	end

	%FLUX time: EXACTLY July 1 of the calendar year (the ISMIP7 FL convention),
	%with bounds spanning the whole year. NOTE: this is deliberately not a
	%computed (start+end)/2 midpoint - that lands on Jul 2 in a non-leap year
	%because 365 days split in half rounds past Jul 1, which the checker flags.
	time_flux      = zeros(1,nyears);
	time_flux_bnds = zeros(2,nyears);
	for j=1:nyears,
		yr = uyears(j);
		b0 = datenum(yr,1,1)  -datenum(1850,1,1);   %start of the year
		b1 = datenum(yr+1,1,1)-datenum(1850,1,1);   %end of the year (=start of next year)
		time_flux_bnds(:,j) = [b0;b1];
		time_flux(j)        = datenum(yr,7,1)-datenum(1850,1,1);   %exactly Jul 1
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

	%tendlicalvf/tendlifmassbf are not computed for historical runs (the front is prescribed -
	%see is_historical above), so rather than writing a file full of fill value, skip it entirely.
	historical_skip_vars = {'tendlicalvf','tendlifmassbf'};
	for v=1:size(scalar_vars,1),
		if is_historical && any(strcmp(scalar_vars{v,1},historical_skip_vars)),
			continue;
		end
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
	%Grounding-line flux needs the depth-averaged velocity on a 3D/layered mesh (this matches
	%what ISSM's own GroundinglineMassFlux diagnostic uses internally - see movingfront_core.cpp).
	%On a genuinely 2D model Vx/Vy already ARE the depth average, so this only matters for
	%GrIS-style 3D setups.
	has_vxavg = isfield(md.results.TransientSolution,'VxAverage');
	if ~has_vxavg,
		warning('ISMIP7:novxavg','VxAverage/VyAverage not found in the solutions - using Vx/Vy directly for the grounding-line flux calculation (only correct for a depth-uniform/2D model).');
	end
	%Some setups leave md.basalforcings.geothermalflux as its uninitialized scalar NaN
	%(or an all-NaN array) rather than a real per-vertex field - indexing that with
	%(1:numberofvertices) errors outright. Treat a missing/all-NaN field as 0 instead.
	geothermalflux_missing = all(isnan(md.basalforcings.geothermalflux(:)));
	if geothermalflux_missing,
		warning('ISMIP7:nogeoflux','md.basalforcings.geothermalflux is NaN - gridded geothermal flux (hfgeoubed) treated as 0.');
	end
	%ISMIP7 standard grid (ISMIP6 domain, EPSG:3413 for GrIS).
	%Allowed resolutions (multiples of 1 km): 1/2/4/8/16 km.
	posting = resolution_km*1000;
	%Cell centres from (-720000,-3450000) to (960000,-570000); 1681 x 2881 at 1 km
	results.gridx  = -720000 :posting: 960000;
	results.gridy  = -3450000:posting:-570000;

	%Parameters for InterpFromMeshToGrid (index,x,y,data,xgrid,ygrid,default_value).
	%Both x and y are ascending: the compliance checker computes grid resolution
	%as coord(2)-coord(1) and only accepts positive values, so a descending axis
	%(negative resolution) fails that check even though the data itself is fine.
	ncols  = numel(results.gridx);
	nlines = numel(results.gridy);
	xgrid  = results.gridx;                      %x grid values (ascending)
	ygrid  = results.gridy;                      %y grid values (ascending)
	results.xcoord = results.gridx;              %x coordinate variable (ascending)
	results.ycoord = results.gridy;              %y coordinate variable (ascending, matches data)

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

	%Extra mandatory ISMIP7 fields handled outside the main per-timestep loop
	results.base       = NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid));
	results.sftgif     = NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid));
	%ligroundf: filled below in the per-timestep loops via compute_ligroundf, in every experiment.
	%lifmassbf: filled below in the per-timestep loops from CalvingMeltingFluxLevelset-CalvingFluxLevelset
	%           (non-historical only - stays fill value for historical/prescribed-front runs).
	results.ligroundf  = NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid));
	results.lifmassbf  = NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid));
	%libmassbfgr: set below (0 under grounded ice - this model has no grounded-ice basal
	%melt/refreeze - fill value elsewhere) once sftgrf is available.
	results.libmassbfgr= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid));

	%Accumulator for the "exceeds 1e6 Pa checker cap" drag warning: rather than firing a
	%warning every time this happens (it can trigger on many vertices in many timesteps for a
	%single run), the per-timestep detail below is logged silently and a single summary
	%warning is raised once, after the loop, with the worst offender across the whole run.
	drag_cap_summary = struct('noccurrences',0,'nvertices_total',0,'worst_value',-Inf, ...
		'worst_vertex',[],'worst_time_index',[],'worst_x',[],'worst_y',[],'worst_C',[],'worst_N',[],'worst_vel',[]);

	for i=1:length(results.timegrid),
		%One informative line per output period instead of letting InterpFromMeshToGrid
		%print its own uninformative "interpolation progress: XX.XX%" for each of the ~20
		%fields interpolated below (that call is wrapped in interp_quiet, which swallows it).
		if block_y0(i)==block_y1(i),
			yearlabel = sprintf('%d',block_y0(i));
		else
			yearlabel = sprintf('%d-%d',block_y0(i),block_y1(i));
		end
		fprintf('md2ismip7: gridding GrIS fields for %s (period %d/%d)\n',yearlabel,i,noutput);
		thickness=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Thickness(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
		results.thickness(:,:,i)=transpose(thickness);
		base=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Base(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
		results.base(:,:,i)=transpose(base);
		surface=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Surface(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
		results.surface(:,:,i)=transpose(surface);
		bed=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.bed(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
		results.bed(:,:,i)=transpose(bed);
		if geothermalflux_missing,
			geoflux_mesh = zeros(md.mesh.numberofvertices,1);
		else
			geoflux_mesh = md.basalforcings.geothermalflux(1:md.mesh.numberofvertices);
		end
		geoflux=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,geoflux_mesh,xgrid,ygrid,NaN);
		results.geoflux(:,:,i)=transpose(geoflux);
		if i==1,
			dhdt=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,(md.results.TransientSolution(results.timegrid(i)).Thickness(1:md.mesh.numberofvertices)-md.geometry.thickness(1:md.mesh.numberofvertices))/(md.results.TransientSolution(results.timegrid(i)).time-0),xgrid,ygrid,NaN);
		else
			dhdt=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,(md.results.TransientSolution(results.timegrid(i)).Thickness(1:md.mesh.numberofvertices)-md.results.TransientSolution(results.timegrid(i-1)).Thickness(1:md.mesh.numberofvertices))/(md.results.TransientSolution(results.timegrid(i)).time-md.results.TransientSolution(results.timegrid(i-1)).time),xgrid,ygrid,NaN);
		end
		results.dhdt(:,:,i)=transpose(dhdt)/(md.constants.yts);
		vxsurf=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vx(end-md.mesh.numberofvertices+1:end),xgrid,ygrid,NaN);
		results.vxsurf(:,:,i)=transpose(vxsurf)/md.constants.yts;
		vysurf=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vy(end-md.mesh.numberofvertices+1:end),xgrid,ygrid,NaN);
		results.vysurf(:,:,i)=transpose(vysurf)/md.constants.yts;
		%vzsurf=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vz(end-md.mesh.numberofvertices+1:end),xgrid,ygrid,NaN);
		%results.vzsurf(:,:,i)=transpose(vzsurf)/md.constants.yts;
		vxbase=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vx(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
		results.vxbase(:,:,i)=transpose(vxbase)/md.constants.yts;
		vybase=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vy(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
		results.vybase(:,:,i)=transpose(vybase)/md.constants.yts;
		%vzbase=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vz(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
		%results.vzbase(:,:,i)=transpose(vzbase)/md.constants.yts;
		vxmean=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vx,xgrid,ygrid,NaN);
         results.vxmean(:,:,i)=transpose(vxmean)/md.constants.yts;
		vymean=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(results.timegrid(i)).Vy,xgrid,ygrid,NaN);
		results.vymean(:,:,i)=transpose(vymean)/md.constants.yts;
		surftemp=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,md.initialization.temperature(end-md.mesh.numberofvertices+1:end),xgrid,ygrid,NaN);
		results.surftemp(:,:,i)=transpose(surftemp);
		basetemp=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,md.initialization.temperature(1:md.mesh.numberofvertices),xgrid,ygrid,NaN);
		results.basetemp(:,:,i)=transpose(basetemp);
		%Effective pressure clamped at 0: right at flotation, floating-point roundoff
		%produced tiny negative drag values that failed the ISMIP7 checker.
		Neff_drag=max(md.constants.g*(md.materials.rho_ice*md.results.TransientSolution(results.timegrid(i)).Thickness(1:md.mesh.numberofvertices)+md.materials.rho_water*md.results.TransientSolution(results.timegrid(i)).Base(1:md.mesh.numberofvertices)),0);
		vel_drag=sqrt(md.results.TransientSolution(results.timegrid(i)).Vx(1:md.mesh.numberofvertices).^2+md.results.TransientSolution(results.timegrid(i)).Vy(1:md.mesh.numberofvertices).^2)/md.constants.yts;
		drag_mesh=md.friction.coefficient(1:md.mesh.numberofvertices).^2.*Neff_drag.*vel_drag;
		%With linear Budd friction confirmed (p=q=1), the formula above is the right
		%physics, so a vertex this far above the checker's 1e6 Pa cap is a data/
		%inversion artifact (a runaway friction coefficient or a locally noisy
		%velocity right at a shear margin) rather than a formula error. The ISMIP7
		%checker tolerates missing values, so these vertices are masked to NaN rather
		%than clamped to an arbitrary, still-fictitious number.
		baddrag=find(drag_mesh>1e6);
		if ~isempty(baddrag),
			[maxdrag,maxi]=max(drag_mesh(baddrag)); worst=baddrag(maxi);
			%Logged into drag_cap_summary rather than warned here - see the summary
			%warning issued once, after the loop, below.
			drag_cap_summary.noccurrences    = drag_cap_summary.noccurrences + 1;
			drag_cap_summary.nvertices_total = drag_cap_summary.nvertices_total + numel(baddrag);
			if maxdrag > drag_cap_summary.worst_value,
				drag_cap_summary.worst_value      = maxdrag;
				drag_cap_summary.worst_vertex      = worst;
				drag_cap_summary.worst_time_index  = i;
				drag_cap_summary.worst_x           = md.mesh.x(worst);
				drag_cap_summary.worst_y           = md.mesh.y(worst);
				drag_cap_summary.worst_C           = md.friction.coefficient(worst);
				drag_cap_summary.worst_N           = Neff_drag(worst);
				drag_cap_summary.worst_vel         = vel_drag(worst)*md.constants.yts;
			end
			drag_mesh(baddrag)=NaN;
		end
		drag=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,drag_mesh,xgrid,ygrid,NaN);
		drag(drag<0)=0; %clamp any residual floating-point noise from interpolation (leaves NaN untouched)
		results.drag(:,:,i)=transpose(drag);
		%--- calving / ice-front-melt: only meaningful when the front is dynamically
		%    simulated (see is_historical note above). Built from the raw calving-law
		%    rate fields (compute_calvingflux, mirroring compute_ligroundf) rather than
		%    ISSM's own CalvingFluxLevelset/CalvingMeltingFluxLevelset, which are
		%    normalized to the vertical ice-front FACE area, not the horizontal
		%    grid-cell area ISMIP7 wants - see compute_calvingflux's header for detail.
		%    The SCALAR totals below (tendlicalvf/tendlifmassbf) never had this problem
		%    - Total* never divides by area - and still come from TotalCalvingFluxLevelset/
		%    TotalCalvingMeltingFluxLevelset as before.
		if ~is_historical,
			icels_calv = md.results.TransientSolution(results.timegrid(i)).MaskIceLevelset(1:md.mesh.numberofvertices);
			crx_calv   = md.results.TransientSolution(results.timegrid(i)).Calvingratex(1:md.mesh.numberofvertices);
			cry_calv   = md.results.TransientSolution(results.timegrid(i)).Calvingratey(1:md.mesh.numberofvertices);
			mr_calv    = md.results.TransientSolution(results.timegrid(i)).CalvingMeltingrate(1:md.mesh.numberofvertices);
			vx_calv    = md.results.TransientSolution(results.timegrid(i)).Vx(1:md.mesh.numberofvertices);
			vy_calv    = md.results.TransientSolution(results.timegrid(i)).Vy(1:md.mesh.numberofvertices);
			thickness_calv = md.results.TransientSolution(results.timegrid(i)).Thickness(1:md.mesh.numberofvertices);
			[calv_elem,melt_elem] = compute_calvingflux(md.mesh.elements,md.mesh.x,md.mesh.y,icels_calv,thickness_calv,crx_calv,cry_calv,mr_calv,vx_calv,vy_calv,md.materials.rho_ice);
			%compute_calvingflux mirrors ISSM's Tria::CalvingFluxLevelset sign convention
			%(positive = mass flux OUT through the ice front). licalvf/lifmassbf instead use
			%the CF mass-balance convention their -10..0 kg m-2 s-1 allowed range implies
			%(negative = mass lost), so the sign is flipped here before interpolation.
			%licalvf/lifmassbf's fill_policy is 'forbidden' (defined over the whole grid, 0
			%where not applicable) - like ligroundf below, default interp_quiet's output to 0
			%rather than NaN for grid points outside the mesh footprint, so the checker's
			%missing-value test is satisfied.
			calving=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,-calv_elem/md.constants.yts,xgrid,ygrid,0);  %compute_calvingflux uses md.materials.rho_ice directly (no Gt scaling) - only /yts needed
			results.calving(:,:,i)=transpose(calving);
			frontmelt=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,-melt_elem/md.constants.yts,xgrid,ygrid,0);
			results.lifmassbf(:,:,i)=transpose(frontmelt);
		end
		mask=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,-md.mask.ice_levelset(1:md.mesh.numberofvertices),xgrid,ygrid,-1);
		mask(find(mask>0))=1;
		mask(find(mask<0))=0;
		results.mask(:,:,i)=transpose(mask);
		gl_raw=md.results.TransientSolution(results.timegrid(i)).MaskOceanLevelset(1:md.mesh.numberofvertices);
		groundedice=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,gl_raw,xgrid,ygrid,NaN);
		groundedice(find(groundedice>0))=1;
		groundedice(find(groundedice<0))=0;
		floatingice=1-groundedice;
		groundedice(find(isnan(groundedice)))=0;
		floatingice(find(isnan(floatingice)))=0;
		results.groundedice(:,:,i)=transpose(groundedice.*mask);
		results.floatingice(:,:,i)=transpose(floatingice.*mask);
		%--- grounding-line flux: physically meaningful in every experiment (unlike
		%    calving, it is not tied to a prescribed front) ---
		if has_vxavg,
			vxgl=md.results.TransientSolution(results.timegrid(i)).VxAverage(1:md.mesh.numberofvertices);
			vygl=md.results.TransientSolution(results.timegrid(i)).VyAverage(1:md.mesh.numberofvertices);
		else
			vxgl=md.results.TransientSolution(results.timegrid(i)).Vx(1:md.mesh.numberofvertices);
			vygl=md.results.TransientSolution(results.timegrid(i)).Vy(1:md.mesh.numberofvertices);
		end
		thickness_gl=md.results.TransientSolution(results.timegrid(i)).Thickness(1:md.mesh.numberofvertices);
		lig_elem=compute_ligroundf(md.mesh.elements,md.mesh.x,md.mesh.y,gl_raw,thickness_gl,vxgl,vygl,md.materials.rho_ice);
		%ligroundf's fill_policy is 'forbidden' (never missing, anywhere in the domain) - lig_elem
		%is already 0 for every element the grounding line doesn't cross, so the only source of
		%NaN here would be interp_quiet's default for grid points outside the mesh footprint
		%(most of the output grid, since it spans a bounding box well beyond the ice extent).
		%Default those to 0 flux instead of NaN so the checker's missing-value test is satisfied.
		ligroundf=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,lig_elem/md.constants.yts,xgrid,ygrid,0);  %compute_ligroundf uses md.materials.rho_ice directly (no Gt scaling) - only /yts needed
		results.ligroundf(:,:,i)=transpose(ligroundf);
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
		bmb=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,meltfl,xgrid,ygrid,NaN);
		%libmassbffl's fill policy (no_floating_ice) requires missing (not 0) everywhere
		%that is NOT floating ice - that's wherever sftflf will end up 0, i.e. grounded ice
		%OR no ice at all. Masking only on groundedice==1 missed the ice-free case (cells the
		%raw ocean levelset happens to sign "floating" even though there is no ice there),
		%which is what the checker's "holds a value where sftflf is 0" error was catching.
		bmbval=-(transpose(bmb)*md.materials.rho_ice/md.constants.yts);
		bmbval(transpose(floatingice.*mask)==0)=NaN;
		results.bmb(:,:,i)=bmbval;
		smb=interp_quiet(md.mesh.elements,md.mesh.x,md.mesh.y,smb_sum/nb,xgrid,ygrid,NaN);
		results.smb(:,:,i)=transpose(smb)*md.materials.rho_ice/md.constants.yts;
	end

	%One summary warning for every "exceeds the 1e6 Pa checker cap" drag event accumulated
	%above, instead of one warning per occurrence (see drag_cap_summary init above the loop).
	if drag_cap_summary.noccurrences>0,
		warning('md2ismip7:drag',['Basal drag exceeded the 1e6 Pa checker cap at ' num2str(drag_cap_summary.nvertices_total) ...
			' vertex-timestep(s) across ' num2str(drag_cap_summary.noccurrences) ' of ' num2str(length(results.timegrid)) ...
			' output period(s); all such vertices were masked to NaN rather than clamped. Worst case: vertex ' ...
			num2str(drag_cap_summary.worst_vertex) ' (x=' num2str(drag_cap_summary.worst_x) ', y=' num2str(drag_cap_summary.worst_y) ...
			') at time index ' num2str(drag_cap_summary.worst_time_index) ': ' num2str(drag_cap_summary.worst_value) ' Pa. C=' ...
			num2str(drag_cap_summary.worst_C) ', N=' num2str(drag_cap_summary.worst_N) ' Pa, vel=' num2str(drag_cap_summary.worst_vel) ...
			' m/yr - inspect the friction inversion/velocity solve there.']);
	end
	%}}}

	%Gridded output: derive mask fractions, build the gridded time axes, then write
	%one file per mandatory variable (ISMIP7 requires a single main variable per file).
	icepresent = double(results.mask>0);                %1 where ice present, 0 elsewhere
	results.sftgif = icepresent;                         %land ice area fraction
	results.sftgrf = results.groundedice .* icepresent;  %grounded fraction, 0 outside ice
	results.sftflf = results.floatingice .* icepresent;  %floating fraction, 0 outside ice

	%Physical consistency: grounded ice rests directly on the bed, so wherever a cell
	%is classified wholly grounded (sftgrf==1), base must equal topg. Left alone,
	%'base' and 'topg' are interpolated independently from the mesh, and at grid cells
	%straddling the grounding line the two interpolations don't always land on exactly
	%the same value even though the (separately thresholded) mask rounds that cell to
	%"grounded" - the same kind of interpolate-then-threshold mismatch as the
	%grounded/floating mask fixes above. Force it directly rather than leave a small
	%mismatch for the checker.
	fullygrounded = results.sftgrf==1;
	results.base(fullygrounded) = results.bed(fullygrounded);
	%orog (surface) was interpolated independently from the mesh's Surface field, so forcing
	%base=bed above - without touching surface - breaks the orog=base+lithk identity the
	%checker requires, by however much bed differs from the original interpolated base at
	%that cell. Re-sync surface here so the identity holds exactly wherever base was forced.
	results.surface(fullygrounded) = results.base(fullygrounded) + results.thickness(fullygrounded);

	%libmassbfgr: this model does not simulate basal melt/refreeze under grounded ice, so
	%the value is a genuine physical 0 (not a missing computation) everywhere grounded ice
	%is present - consistent with tendlibmassbfgr above, which is computed from
	%TotalGroundedBmbScaled and expected to come out ~0 for the same reason. Left at fill
	%outside grounded ice (floating ice or ice-free), consistent with libmassbffl covering
	%the floating regime instead and the fill_policy for cells where neither applies.
	results.libmassbfgr(results.sftgrf>0) = 0;

	%Data-request fill policies the generic NaN handling in write_gridded_var
	%doesn't know about:
	%  - 'no_ice' (xvelmean, yvelmean, strbasemag): must be FILL outside the ice
	%    mask. The ISSM mesh can extend slightly beyond the ice extent (bare rock
	%    margins), so InterpFromMeshToGrid returns real, non-NaN values there;
	%    force those back to NaN so write_gridded_var converts them to fill.
	results.vxmean(icepresent==0) = NaN;
	results.vymean(icepresent==0) = NaN;
	results.drag(icepresent==0)   = NaN;

	%STATE grid time: EXACT end-of-year instant (Jan 1 of the year after the last
	%calendar year in each output period), per the ST convention - see the note
	%by the scalar STATE time above for why this isn't the raw model save time.
	time_grid_state = zeros(1,noutput);
	for j=1:noutput,
		time_grid_state(j) = datenum(block_y1(j)+1,1,1) - datenum(1850,1,1);
	end

	%FLUX grid time: EXACTLY July 1 for a single-calendar-year output period (the
	%FL convention - see the scalar FLUX time note above); for a multi-year period
	%(output_interval_yr>1) there is no single "Jul 1" to snap to, so the true
	%midpoint of the period is used instead. Bounds always span the period's start
	%to its end.
	time_grid_flux      = zeros(1,noutput);
	time_grid_flux_bnds = zeros(2,noutput);
	for j=1:noutput,
		b0 = datenum(block_y0(j),1,1)  -datenum(1850,1,1);
		b1 = datenum(block_y1(j)+1,1,1)-datenum(1850,1,1);
		time_grid_flux_bnds(:,j) = [b0;b1];
		if block_y0(j)==block_y1(j),
			time_grid_flux(j) = datenum(block_y0(j),7,1)-datenum(1850,1,1);   %exactly Jul 1
		else
			time_grid_flux(j) = (b0+b1)/2;                                    %multi-year: true midpoint
		end
	end

	if is_historical,
		warning('ISMIP7:placeholder',['licalvf, lifmassbf (gridded) and tendlicalvf, tendlifmassbf (scalar) ' ...
			'are not produced by this converter for historical runs, since the front is prescribed ' ...
			'(no calving law active) and these diagnostics are not physically meaningful - their ' ...
			'.nc files are skipped entirely rather than written out at fill value.']);
	end

	%Columns: variable_id, standard_name, units, data, zero_outside_ice, is_flux
	%  zero_outside_ice=true: data-request fill_policy is 'forbidden' - the
	%    variable must be 0 (not fill) wherever there is no ice, so NaN there is
	%    written as 0. Applies to the mask fractions themselves and to lithk /
	%    dlithkdt (fill_policy=forbidden per the data request).
	%  is_flux=false (ST): end-of-year snapshot, no time bounds
	%  is_flux=true  (FL): yearly average at mid-year, with time bounds
	grid_vars = {
		'lithk',       'land_ice_thickness',                         'm',          results.thickness,   true,  false;
		'orog',        'surface_altitude',                           'm',          results.surface,     false, false;
		'topg',        'bedrock_altitude',                           'm',          results.bed,         false, false;
		'base',        '',                                           'm',          results.base,        false, false;
		'acabf',       'land_ice_surface_specific_mass_balance_flux','kg m-2 s-1', results.smb,         false, true;
		'libmassbfgr', 'land_ice_basal_specific_mass_balance_flux',  'kg m-2 s-1', results.libmassbfgr, false, true;
		'libmassbffl', 'land_ice_basal_specific_mass_balance_flux',  'kg m-2 s-1', results.bmb,         false, true;
		'dlithkdt',    'tendency_of_land_ice_thickness',             'm s-1',      results.dhdt,        true,  true;
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
		if is_historical && any(strcmp(grid_vars{v,1},{'licalvf','lifmassbf'})),
			continue;
		end
		if grid_vars{v,6},
			t=time_grid_flux; b=time_grid_flux_bnds;
		else
			t=time_grid_state; b=[];
		end
		write_gridded_var(mkfname(grid_vars{v,1}),grid_vars{v,1},grid_vars{v,2},grid_vars{v,3},grid_vars{v,4},results.xcoord,results.ycoord,t,b,fillval,grid_vars{v,5},source_id,ism_id,icesheetname);
	end

function out=interp_quiet(elements,x,y,data,xgrid,ygrid,defaultvalue)
	%interp_quiet - thin wrapper around InterpFromMeshToGrid that swallows the
	%"interpolation progress: XX.XX%" line it prints on every call. The main loop above
	%calls InterpFromMeshToGrid roughly 20 times per output period, so letting it print
	%directly floods the console with uninformative lines; md2ismip7 prints its own
	%one-line, per-period progress message instead (see the main loop above).
	evalc('out=InterpFromMeshToGrid(elements,x,y,data,xgrid,ygrid,defaultvalue);');

function ismip7_global_attributes(ncid,source_id,ism_id,domain_id)
	%Mandatory ISMIP7 global attributes (see conventions section 5)
	crs='epsg:3413'; %GrIS only
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

function ligroundf=compute_ligroundf(elements,x,y,gl,thickness,vx,vy,rho_ice)
	%compute_ligroundf - per-element grounding-line mass flux, per unit HORIZONTAL area
	%
	%   elements          : Ne x 3 vertex indices (1-based)
	%   x, y              : Nv x 1 vertex coordinates
	%   gl                : Nv x 1 grounding mask (>0 grounded, <0 floating) - e.g.
	%                       MaskOceanLevelset or MaskGroundediceLevelset
	%   thickness, vx, vy : Nv x 1, at the SAME timestep as gl (vx,vy should be the
	%                       depth-averaged velocity for a 3D/layered mesh)
	%   rho_ice           : scalar, kg/m^3
	%
	%Returns ligroundf: Ne x 1, flux per unit horizontal area, in the same time units
	%as vx/vy (e.g. kg/m^2/yr if vx,vy are m/yr - divide by yts outside this function
	%for kg/m^2/s). Zero for elements the grounding line does not pass through.
	%
	%Geometry and the front-point ordering exactly mirror ISSM's
	%Tria::GroundinglineMassFlux (classes/Elements/Tria.cpp); the normal formula
	%exactly mirrors shared/Numerics/Normals.cpp:LineSectionNormal, so - unlike
	%Tria::GroundinglineMassFlux itself - this version divides by the element's
	%HORIZONTAL area instead of skipping the area division, giving a true
	%per-unit-area flux rather than a raw per-element mass-flux contribution.

	Ne = size(elements,1);
	ligroundf = zeros(Ne,1);
	eps_ = 1e-15;

	for e = 1:Ne,
		idx = elements(e,:);
		xyz = [x(idx), y(idx)];
		g   = gl(idx);
		g(g==0) = eps_;   %avoid exact-zero degenerate case

		if all(g>0) || all(g<0), continue; end   %grounding line does not cross this element

		h = thickness(idx);
		u = vx(idx);
		v = vy(idx);

		xyz_front = zeros(2,2);
		h_front   = zeros(2,1);
		u_front   = zeros(2,1);
		v_front   = zeros(2,1);
		pt1 = 1; pt2 = 2;   %MATLAB 1-based version of the C++ pt1=0,pt2=1

		if g(1)*g(2) > 0,        %nodes 1,2 same sign -> cross edges (3-2) and (3-1)
			s1 = g(3)/(g(3)-g(2));
			s2 = g(3)/(g(3)-g(1));
			if g(3) < 0, pt1 = 2; pt2 = 1; end
			xyz_front(pt2,:) = xyz(3,:) + s1*(xyz(2,:)-xyz(3,:));
			xyz_front(pt1,:) = xyz(3,:) + s2*(xyz(1,:)-xyz(3,:));
			h_front(pt2) = h(3)+s1*(h(2)-h(3)); h_front(pt1) = h(3)+s2*(h(1)-h(3));
			u_front(pt2) = u(3)+s1*(u(2)-u(3)); u_front(pt1) = u(3)+s2*(u(1)-u(3));
			v_front(pt2) = v(3)+s1*(v(2)-v(3)); v_front(pt1) = v(3)+s2*(v(1)-v(3));
		elseif g(2)*g(3) > 0,    %nodes 2,3 same sign -> cross edges (1-2) and (1-3)
			s1 = g(1)/(g(1)-g(2));
			s2 = g(1)/(g(1)-g(3));
			if g(1) < 0, pt1 = 2; pt2 = 1; end
			xyz_front(pt1,:) = xyz(1,:) + s1*(xyz(2,:)-xyz(1,:));
			xyz_front(pt2,:) = xyz(1,:) + s2*(xyz(3,:)-xyz(1,:));
			h_front(pt1) = h(1)+s1*(h(2)-h(1)); h_front(pt2) = h(1)+s2*(h(3)-h(1));
			u_front(pt1) = u(1)+s1*(u(2)-u(1)); u_front(pt2) = u(1)+s2*(u(3)-u(1));
			v_front(pt1) = v(1)+s1*(v(2)-v(1)); v_front(pt2) = v(1)+s2*(v(3)-v(1));
		else                      %nodes 1,3 same sign -> cross edges (2-1) and (2-3)
			s1 = g(2)/(g(2)-g(1));
			s2 = g(2)/(g(2)-g(3));
			if g(2) < 0, pt1 = 2; pt2 = 1; end
			xyz_front(pt2,:) = xyz(2,:) + s1*(xyz(1,:)-xyz(2,:));
			xyz_front(pt1,:) = xyz(2,:) + s2*(xyz(3,:)-xyz(2,:));
			h_front(pt2) = h(2)+s1*(h(1)-h(2)); h_front(pt1) = h(2)+s2*(h(3)-h(2));
			u_front(pt2) = u(2)+s1*(u(1)-u(2)); u_front(pt1) = u(2)+s2*(u(3)-u(2));
			v_front(pt2) = v(2)+s1*(v(1)-v(2)); v_front(pt1) = v(2)+s2*(v(3)-v(2));
		end

		d = xyz_front(2,:) - xyz_front(1,:);
		L = norm(d);
		if L < eps_, continue; end
		n = [d(2), -d(1)]/L;   %LineSectionNormal: normal=[dy,-dx], normalized

		%Simpson's rule (exact for this quadratic integrand): t=0, 0.5, 1
		hm = 0.5*(h_front(1)+h_front(2));
		um = 0.5*(u_front(1)+u_front(2));
		vm = 0.5*(v_front(1)+v_front(2));

		f0 = rho_ice*h_front(1)*(u_front(1)*n(1)+v_front(1)*n(2));
		f1 = rho_ice*h_front(2)*(u_front(2)*n(1)+v_front(2)*n(2));
		fm = rho_ice*hm       *(um*n(1)+vm*n(2));
		flux = L/6*(f0 + 4*fm + f1);

		%Horizontal (map-view) element area - NOT the front cross-section area
		area = 0.5*abs((xyz(1,1)-xyz(3,1))*(xyz(2,2)-xyz(1,2)) - ...
		               (xyz(1,1)-xyz(2,1))*(xyz(3,2)-xyz(1,2)));

		ligroundf(e) = flux/area;
	end

function [calvflux,meltflux]=compute_calvingflux(elements,x,y,icels,thickness,calvingratex,calvingratey,meltingrate,vx,vy,rho_ice)
	%compute_calvingflux - per-element calving-only and melt-only mass flux, per unit
	%HORIZONTAL area (fixes the front-face-area normalization ISSM's own
	%CalvingFluxLevelset/CalvingMeltingFluxLevelset use - see md2ismip7's calving
	%comments above).
	%
	%   elements                   : Ne x 3 vertex indices (1-based)
	%   x, y                       : Nv x 1 vertex coordinates
	%   icels                      : Nv x 1 ice-front mask, MaskIceLevelset (<0 ice present)
	%   thickness                  : Nv x 1
	%   calvingratex, calvingratey : Nv x 1, calving-law rate vector (Calvingratex/Calvingratey)
	%   meltingrate                : Nv x 1, frontal/undercutting melt rate MAGNITUDE (CalvingMeltingrate)
	%   vx, vy                     : Nv x 1, ice velocity - gives the melt rate a direction
	%                                (plain Vx/Vy, matching what Tria::CalvingMeltingFluxLevelset
	%                                itself uses - not the depth-averaged velocity grounding
	%                                line needs)
	%   rho_ice                    : scalar, kg/m^3
	%
	%Returns calvflux, meltflux: Ne x 1 each, flux per unit horizontal area, in the same
	%time units as calvingratex/vx (e.g. kg/m^2/yr if those are m/yr - divide by yts
	%outside this function for kg/m^2/s). Zero for elements the ice front does not cross.
	%
	%Geometry, front-point ordering and rate formulas exactly mirror ISSM's
	%Tria::CalvingFluxLevelset / Tria::CalvingMeltingFluxLevelset (classes/Elements/Tria.cpp):
	%calving-only uses (calvingratex,calvingratey); melt-only uses meltingrate projected onto
	%the local ice-velocity direction. The normal is LineSectionNormal's [dy,-dx], NEGATED
	%(Tria::CalvingFluxLevelset flips the sign after calling NormalSection - unlike
	%Tria::GroundinglineMassFlux, which does not). This function divides by the element's
	%HORIZONTAL area instead of skipping the division (Total*) or dividing by the vertical
	%ice-front face area (ISSM's own CalvingFluxLevelset/CalvingMeltingFluxLevelset), and
	%returns calving-only and melt-only separately instead of calving-only and a
	%calving+melt-combined field.
	%
	%ACCURACY NOTE: the calving-only integrand is quadratic in the segment parameter (like
	%compute_ligroundf), so the 3-point Simpson's rule below is exact for it. The melt-only
	%integrand involves meltingrate*v/|v|, a ratio of linear functions, so it is NOT
	%exactly polynomial - Simpson's rule is only an approximation there, evaluated pointwise
	%at the same 3 points ISSM's own 3-point Gauss quadrature would use, so it should be
	%comparably accurate to ISSM's own internal computation, not a downgrade from it.

	eps_ = 1e-15;
	Ne = size(elements,1);
	calvflux = zeros(Ne,1);
	meltflux = zeros(Ne,1);

	for e = 1:Ne,
		idx = elements(e,:);
		xyz = [x(idx), y(idx)];
		g   = icels(idx);
		g(g==0) = eps_;

		if all(g>0) || all(g<0), continue; end   %ice front does not cross this element

		h   = thickness(idx);
		crx = calvingratex(idx);
		cry = calvingratey(idx);
		mr  = meltingrate(idx);
		u   = vx(idx);
		v   = vy(idx);

		xyz_front = zeros(2,2);
		h_front=zeros(2,1); crx_front=zeros(2,1); cry_front=zeros(2,1);
		mr_front=zeros(2,1); u_front=zeros(2,1); v_front=zeros(2,1);
		pt1 = 1; pt2 = 2;

		if g(1)*g(2) > 0,        %nodes 1,2 same sign -> cross edges (3-2) and (3-1)
			s1 = g(3)/(g(3)-g(2));
			s2 = g(3)/(g(3)-g(1));
			if g(3) < 0, pt1 = 2; pt2 = 1; end
			xyz_front(pt2,:) = xyz(3,:) + s1*(xyz(2,:)-xyz(3,:));
			xyz_front(pt1,:) = xyz(3,:) + s2*(xyz(1,:)-xyz(3,:));
			h_front(pt2)  = h(3)+s1*(h(2)-h(3));     h_front(pt1)  = h(3)+s2*(h(1)-h(3));
			crx_front(pt2)= crx(3)+s1*(crx(2)-crx(3)); crx_front(pt1)= crx(3)+s2*(crx(1)-crx(3));
			cry_front(pt2)= cry(3)+s1*(cry(2)-cry(3)); cry_front(pt1)= cry(3)+s2*(cry(1)-cry(3));
			mr_front(pt2) = mr(3)+s1*(mr(2)-mr(3));   mr_front(pt1) = mr(3)+s2*(mr(1)-mr(3));
			u_front(pt2)  = u(3)+s1*(u(2)-u(3));       u_front(pt1)  = u(3)+s2*(u(1)-u(3));
			v_front(pt2)  = v(3)+s1*(v(2)-v(3));       v_front(pt1)  = v(3)+s2*(v(1)-v(3));
		elseif g(2)*g(3) > 0,    %nodes 2,3 same sign -> cross edges (1-2) and (1-3)
			s1 = g(1)/(g(1)-g(2));
			s2 = g(1)/(g(1)-g(3));
			if g(1) < 0, pt1 = 2; pt2 = 1; end
			xyz_front(pt1,:) = xyz(1,:) + s1*(xyz(2,:)-xyz(1,:));
			xyz_front(pt2,:) = xyz(1,:) + s2*(xyz(3,:)-xyz(1,:));
			h_front(pt1)  = h(1)+s1*(h(2)-h(1));     h_front(pt2)  = h(1)+s2*(h(3)-h(1));
			crx_front(pt1)= crx(1)+s1*(crx(2)-crx(1)); crx_front(pt2)= crx(1)+s2*(crx(3)-crx(1));
			cry_front(pt1)= cry(1)+s1*(cry(2)-cry(1)); cry_front(pt2)= cry(1)+s2*(cry(3)-cry(1));
			mr_front(pt1) = mr(1)+s1*(mr(2)-mr(1));   mr_front(pt2) = mr(1)+s2*(mr(3)-mr(1));
			u_front(pt1)  = u(1)+s1*(u(2)-u(1));       u_front(pt2)  = u(1)+s2*(u(3)-u(1));
			v_front(pt1)  = v(1)+s1*(v(2)-v(1));       v_front(pt2)  = v(1)+s2*(v(3)-v(1));
		else                      %nodes 1,3 same sign -> cross edges (2-1) and (2-3)
			s1 = g(2)/(g(2)-g(1));
			s2 = g(2)/(g(2)-g(3));
			if g(2) < 0, pt1 = 2; pt2 = 1; end
			xyz_front(pt2,:) = xyz(2,:) + s1*(xyz(1,:)-xyz(2,:));
			xyz_front(pt1,:) = xyz(2,:) + s2*(xyz(3,:)-xyz(2,:));
			h_front(pt2)  = h(2)+s1*(h(1)-h(2));     h_front(pt1)  = h(2)+s2*(h(3)-h(2));
			crx_front(pt2)= crx(2)+s1*(crx(1)-crx(2)); crx_front(pt1)= crx(2)+s2*(crx(3)-crx(2));
			cry_front(pt2)= cry(2)+s1*(cry(1)-cry(2)); cry_front(pt1)= cry(2)+s2*(cry(3)-cry(2));
			mr_front(pt2) = mr(2)+s1*(mr(1)-mr(2));   mr_front(pt1) = mr(2)+s2*(mr(3)-mr(2));
			u_front(pt2)  = u(2)+s1*(u(1)-u(2));       u_front(pt1)  = u(2)+s2*(u(3)-u(2));
			v_front(pt2)  = v(2)+s1*(v(1)-v(2));       v_front(pt1)  = v(2)+s2*(v(3)-v(2));
		end

		d = xyz_front(2,:) - xyz_front(1,:);
		L = norm(d);
		if L < eps_, continue; end
		n = [-d(2), d(1)]/L;   %LineSectionNormal=[dy,-dx], then NEGATED (matches Tria::CalvingFluxLevelset)

		%midpoint values (linear interpolation along the front - exact since h/crx/cry/mr/u/v are all P1)
		hm   = 0.5*(h_front(1)+h_front(2));
		crxm = 0.5*(crx_front(1)+crx_front(2));
		crym = 0.5*(cry_front(1)+cry_front(2));
		um   = 0.5*(u_front(1)+u_front(2));
		vm   = 0.5*(v_front(1)+v_front(2));
		mrm  = 0.5*(mr_front(1)+mr_front(2));

		%melt rate direction = local ice-velocity direction, magnitude = meltingrate,
		%evaluated pointwise at t=0, 0.5, 1 (matches Tria::CalvingMeltingFluxLevelset)
		vel0 = sqrt(u_front(1)^2+v_front(1)^2)+1e-14;
		vel1 = sqrt(u_front(2)^2+v_front(2)^2)+1e-14;
		velm = sqrt(um^2+vm^2)+1e-14;
		mrx0 = mr_front(1)*u_front(1)/vel0; mry0 = mr_front(1)*v_front(1)/vel0;
		mrx1 = mr_front(2)*u_front(2)/vel1; mry1 = mr_front(2)*v_front(2)/vel1;
		mrxm = mrm*um/velm;                 mrym = mrm*vm/velm;

		%calving-only flux (Simpson's rule - exact, quadratic integrand)
		f0c = rho_ice*h_front(1)*(crx_front(1)*n(1)+cry_front(1)*n(2));
		f1c = rho_ice*h_front(2)*(crx_front(2)*n(1)+cry_front(2)*n(2));
		fmc = rho_ice*hm       *(crxm*n(1)+crym*n(2));
		fluxc = L/6*(f0c + 4*fmc + f1c);

		%melt-only flux (Simpson's rule - approximation, see ACCURACY NOTE above)
		f0m = rho_ice*h_front(1)*(mrx0*n(1)+mry0*n(2));
		f1m = rho_ice*h_front(2)*(mrx1*n(1)+mry1*n(2));
		fmm = rho_ice*hm       *(mrxm*n(1)+mrym*n(2));
		fluxm = L/6*(f0m + 4*fmm + f1m);

		%Horizontal (map-view) element area - NOT the front cross-section area
		area = 0.5*abs((xyz(1,1)-xyz(3,1))*(xyz(2,2)-xyz(1,2)) - ...
		               (xyz(1,1)-xyz(2,1))*(xyz(3,2)-xyz(1,2)));

		calvflux(e) = fluxc/area;
		meltflux(e) = fluxm/area;
	end
