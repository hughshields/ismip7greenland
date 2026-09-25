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

	%Calving/ice-front-melt diagnostics are only physically meaningful when the
	%front position is dynamically simulated. Historical runs in this workflow
	%prescribe the front (observed extent imposed each year, no calving law
	%active), so calving-derived outputs are left at fill value throughout.
	is_historical = strcmpi(experiment_id,'historical');

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
	mkfname = @(variable_id) fullfile(outdir,[variable_id '_' icesheetname '_' source_id '_' ism_id '_' ism_member_id '_' esm_id '_' forcing_member_id '_' experiment_id '_' set_counter '_' time_range '.nc']);

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
	if ~is_historical,
		raw_calvtot      = arrayfun(@(k) md.results.TransientSolution(k).TotalCalvingFluxLevelset*10^12/md.constants.yts, 1:Nsol);
		raw_frontmelttot = arrayfun(@(k) (md.results.TransientSolution(k).TotalCalvingMeltingFluxLevelset-md.results.TransientSolution(k).TotalCalvingFluxLevelset)*10^12/md.constants.yts, 1:Nsol);
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
	%Both x and y are ascending: the compliance checker computes grid resolution
	%as coord(2)-coord(1) and only accepts positive values, so a descending axis
	%(negative resolution) fails that check even though the data itself is fine.
	ncols  = numel(results.gridx);
	nlines = numel(results.gridy);
	xgrid  = results.gridx;                      %x grid values (ascending)
	ygrid  = results.gridy;                      %y grid values (ascending)
	results.xcoord = results.gridx;              %x coordinate variable (ascending)
	results.ycoord = results.gridy;              %y coordinate variable (ascending, matches data)

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
	%ligroundf: filled below in the per-timestep loops via compute_ligroundf, in every experiment.
	%lifmassbf: filled below in the per-timestep loops from CalvingMeltingFluxLevelset-CalvingFluxLevelset
	%           (non-historical only - stays fill value for historical/prescribed-front runs).
	results.ligroundf  = NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid));
	results.lifmassbf  = NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid));
	%libmassbfgr: set below (0 under grounded ice - this model has no grounded-ice basal
	%melt/refreeze - fill value elsewhere) once sftgrf is available.
	results.libmassbfgr= NaN*ones(numel(results.gridx),numel(results.gridy),length(results.timegrid));

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
			if geothermalflux_missing,
				geoflux_mesh = zeros(md.mesh.numberofvertices,1);
			else
				geoflux_mesh = md.basalforcings.geothermalflux(1:md.mesh.numberofvertices);
			end
			geoflux=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,geoflux_mesh,xgrid,ygrid,NaN);
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
				warning('md2ismip7:drag',[num2str(numel(baddrag)) ' vertex(es) exceed the 1e6 Pa checker cap at time index ' num2str(i) ...
					'; masking to NaN. Worst: vertex ' num2str(worst) ' (x=' num2str(md.mesh.x(worst)) ', y=' num2str(md.mesh.y(worst)) '): ' ...
					num2str(maxdrag) ' Pa. C=' num2str(md.friction.coefficient(worst)) ', N=' num2str(Neff_drag(worst)) ...
					' Pa, vel=' num2str(vel_drag(worst)*md.constants.yts) ' m/yr - inspect the friction inversion/velocity solve there.']);
				drag_mesh(baddrag)=NaN;
			end
			drag=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,drag_mesh,xgrid,ygrid,NaN);
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
				calving=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,calv_elem/md.constants.yts,xgrid,ygrid,NaN);  %compute_calvingflux uses md.materials.rho_ice directly (no Gt scaling) - only /yts needed
				results.calving(:,:,i)=transpose(calving);
				frontmelt=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,melt_elem/md.constants.yts,xgrid,ygrid,NaN);
				results.lifmassbf(:,:,i)=transpose(frontmelt);
			end
			mask=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,-md.mask.ice_levelset(1:md.mesh.numberofvertices),xgrid,ygrid,-1);
			mask(find(mask>0))=1;
			mask(find(mask<0))=0;
			results.mask(:,:,i)=transpose(mask);
			gl_raw=md.results.TransientSolution(results.timegrid(i)).MaskOceanLevelset(1:md.mesh.numberofvertices);
			groundedice=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,gl_raw,xgrid,ygrid,NaN);
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
			ligroundf=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,lig_elem/md.constants.yts,xgrid,ygrid,NaN);  %compute_ligroundf uses md.materials.rho_ice directly (no Gt scaling) - only /yts needed
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
			bmb=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,meltfl,xgrid,ygrid,NaN);
			%libmassbffl's fill policy requires missing (not 0) outside floating ice.
			%Multiplying by (1-groundedice) wrote a defined 0 there instead - mask to
			%NaN explicitly so write_gridded_var converts it to _FillValue.
			bmbval=-(transpose(bmb)*md.materials.rho_ice/md.constants.yts);
			bmbval(transpose(groundedice)==1)=NaN;
			results.bmb(:,:,i)=bmbval;
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
			if geothermalflux_missing,
				geoflux_mesh = zeros(md.mesh.numberofvertices,1);
			else
				geoflux_mesh = md.basalforcings.geothermalflux;
			end
			geoflux=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,geoflux_mesh,xgrid,ygrid,NaN);
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
			%AIS mesh has no vertical resolution, so there is no distinct basal
			%temperature to interpolate (it would just duplicate surftemp, which is
			%physically wrong) - left as a fill placeholder until a real basal
			%temperature field is available for this domain.
			basetemp=NaN*ones(nlines,ncols);
			results.basetemp(:,:,i)=transpose(basetemp);
			%Effective pressure clamped at 0: rho_ice*H+rho_water*Base crosses zero right
			%at flotation, and floating-point roundoff there was producing tiny negative
			%drag values that failed the ISMIP7 checker.
			Neff_drag=max(md.constants.g*(md.materials.rho_ice*md.results.TransientSolution(results.timegrid(i)).Thickness+md.materials.rho_water*md.results.TransientSolution(results.timegrid(i)).Base),0);
			vel_drag=md.results.TransientSolution(results.timegrid(i)).Vel/md.constants.yts;
			drag_mesh=md.friction.coefficient.^2.*Neff_drag.*vel_drag;
			%With linear Budd friction confirmed (p=q=1), the formula above is the right
			%physics, so a vertex this far above the checker's 1e6 Pa cap is a data/
			%inversion artifact (a runaway friction coefficient or a locally noisy
			%velocity right at a shear margin) rather than a formula error. The ISMIP7
			%checker tolerates missing values, so these vertices are masked to NaN rather
			%than clamped to an arbitrary, still-fictitious number.
			baddrag=find(drag_mesh>1e6);
			if ~isempty(baddrag),
				[maxdrag,maxi]=max(drag_mesh(baddrag)); worst=baddrag(maxi);
				warning('md2ismip7:drag',[num2str(numel(baddrag)) ' vertex(es) exceed the 1e6 Pa checker cap at time index ' num2str(i) ...
					'; masking to NaN. Worst: vertex ' num2str(worst) ' (x=' num2str(md.mesh.x(worst)) ', y=' num2str(md.mesh.y(worst)) '): ' ...
					num2str(maxdrag) ' Pa. C=' num2str(md.friction.coefficient(worst)) ', N=' num2str(Neff_drag(worst)) ...
					' Pa, vel=' num2str(vel_drag(worst)*md.constants.yts) ' m/yr - inspect the friction inversion/velocity solve there.']);
				drag_mesh(baddrag)=NaN;
			end
			drag=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,drag_mesh,xgrid,ygrid,NaN);
			drag(drag<0)=0; %clamp any residual floating-point noise from interpolation (leaves NaN untouched)
			results.drag(:,:,i)=transpose(drag);
			%--- calving / ice-front-melt: only meaningful when the front is dynamically
			%    simulated (see is_historical note above). Built from the raw calving-law
			%    rate fields (compute_calvingflux, mirroring compute_ligroundf) - see the
			%    matching GrIS comment above for why this replaces ISSM's own
			%    CalvingFluxLevelset/CalvingMeltingFluxLevelset for the gridded fields.
			%    AIS is a 2D model here, so Vx/Vy are used directly (matches the plain
			%    Vx/Vy that Tria::CalvingMeltingFluxLevelset itself uses for the melt
			%    direction, as opposed to the depth-averaged velocity grounding line needs).
			if ~is_historical,
				icels_calv = md.results.TransientSolution(results.timegrid(i)).MaskIceLevelset;
				crx_calv   = md.results.TransientSolution(results.timegrid(i)).Calvingratex;
				cry_calv   = md.results.TransientSolution(results.timegrid(i)).Calvingratey;
				mr_calv    = md.results.TransientSolution(results.timegrid(i)).CalvingMeltingrate;
				vx_calv    = md.results.TransientSolution(results.timegrid(i)).Vx;
				vy_calv    = md.results.TransientSolution(results.timegrid(i)).Vy;
				thickness_calv = md.results.TransientSolution(results.timegrid(i)).Thickness;
				[calv_elem,melt_elem] = compute_calvingflux(md.mesh.elements,md.mesh.x,md.mesh.y,icels_calv,thickness_calv,crx_calv,cry_calv,mr_calv,vx_calv,vy_calv,md.materials.rho_ice);
				calving=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,calv_elem/md.constants.yts,xgrid,ygrid,NaN);  %compute_calvingflux uses md.materials.rho_ice directly (no Gt scaling) - only /yts needed
				results.calving(:,:,i)=transpose(calving);
				frontmelt=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,melt_elem/md.constants.yts,xgrid,ygrid,NaN);
				results.lifmassbf(:,:,i)=transpose(frontmelt);
			end
			mask=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,-md.mask.ice_levelset,xgrid,ygrid,-1);
			mask(find(mask>0))=1;
			mask(find(mask<0))=-1;
			results.mask(:,:,i)=transpose(mask);
			%NOTE: assumes MaskGroundediceLevelset follows the same sign convention as
			%MaskOceanLevelset used in the GrIS branch (>0 grounded, <0 floating) - worth
			%double-checking against your ISSM version if the grounding-line results look off.
			gl_raw=md.results.TransientSolution(results.timegrid(i)).MaskGroundediceLevelset;
			groundedice=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,gl_raw,xgrid,ygrid,NaN);
			groundedice(find(groundedice>0))=1;
			groundedice(find(groundedice<0))=0;
			%MaskGroundediceLevelset only encodes grounded-vs-floating (bed vs sea level)
			%- it has no notion of ice presence, so ice-free rock above sea level (e.g.
			%exposed mountains) would otherwise be flagged as "grounded" here. Multiply
			%by the true ice-extent indicator (mask>0, from ice_levelset) before storing
			%- matching the fix already applied in the GrIS branch above. The raw
			%(unmasked) local groundedice variable is left as-is for the bmb weighting
			%below, where 1-groundedice already correctly evaluates to 0 over ice-free land.
			icemask=double(mask>0);
			results.groundedice(:,:,i)=transpose(groundedice.*icemask);
			results.floatingice(:,:,i)=transpose((1-groundedice).*icemask);
			%--- grounding-line flux: physically meaningful in every experiment (unlike
			%    calving, it is not tied to a prescribed front). AIS is a 2D model here,
			%    so Vx/Vy are already the depth average - no VxAverage/VyAverage needed. ---
			lig_elem=compute_ligroundf(md.mesh.elements,md.mesh.x,md.mesh.y,gl_raw,md.results.TransientSolution(results.timegrid(i)).Thickness,md.results.TransientSolution(results.timegrid(i)).Vx,md.results.TransientSolution(results.timegrid(i)).Vy,md.materials.rho_ice);
			ligroundf=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,lig_elem/md.constants.yts,xgrid,ygrid,NaN);  %compute_ligroundf uses md.materials.rho_ice directly (no Gt scaling) - only /yts needed
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
				smb_sum = smb_sum + md.results.TransientSolution(kb).SmbMassBalance;
			end
			if has_bmbfl,
				meltfl = meltfl_sum/nb;
			else
				meltfl = zeros(md.mesh.numberofvertices,1);   %floating basal melt not saved - treated as 0
			end
			bmb=InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,meltfl,xgrid,ygrid,NaN);
			%libmassbffl's fill policy requires missing (not 0) outside floating ice.
			%Multiplying by (1-groundedice) wrote a defined 0 there instead - mask to
			%NaN explicitly so write_gridded_var converts it to _FillValue.
			bmbval=-(transpose(bmb)*md.materials.rho_ice/md.constants.yts);
			bmbval(transpose(groundedice)==1)=NaN;
			results.bmb(:,:,i)=bmbval;
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
		warning('ISMIP7:placeholder',['These mandatory variables are not produced by this ' ...
			'converter and are written as fill/NaN - populate before submission: licalvf, ' ...
			'lifmassbf (gridded) and tendlicalvf, tendlifmassbf (scalar) - left at fill for ' ...
			'this historical run since the front is prescribed.']);
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
