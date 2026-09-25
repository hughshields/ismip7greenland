steps=[3];
%Run Options{{{

% Cluster Options (for forward transients and inversions)
cluster = 1; % 0 for 'totten', 1 for 'andes', 2 for epica
clusternames = {'totten', 'andes', 'local'};
clustername = clusternames{cluster+1};

% Run parameters
max_run_time = 36; % hours
mem_alloc = 500; % in Gb
num_cpus = 64; % number of cores
interactive = 1;
loadonly = 1;

% Run start and end times
historical_start_time = 2007;
historical_end_time = 2015;
projection_start_time = 2015;
projection_end_time = 2301;

% Climate Models
climate_models = {'CESM2-WACCM', 'MRI-ESM2-0'};
climate_model = climate_models{1}
% Forecast Parameters
scenarios = {'ssp585', 'ssp370', 'ssp126', 'ctrl'};
scenario = scenarios{1};
% Folder to store models
folder = 'models';

addpath functions

%}}}
%Cluster parameters{{{
if strcmpi(clustername,'totten'),
	cluster=generic('name',oshostname(),'np',num_cpus);
elseif strcmpi(clustername,'andes'),
	cluster=andes('numnodes',1,'cpuspernode',num_cpus,'time',max_run_time,'memory', mem_alloc);
	interactive = 0;
else
	cluster=generic('name',oshostname(),'np',num_cpus);
end
%}}}

% Set up organizer to save the model after each step
org=organizer('repository',['./' folder],'prefix', '', 'steps',steps); clear steps;

% Steps 1-2: Historical runs (1 for each climate model)
if perform(org,['Greenland_ISMIP7Prep_' climate_model '_Historical']),% {{{

	md=loadmodel(org,['Greenland_TransientCalibration']);

	% Set the transient parameters
	md.transient.isthermal=0;
	md.transient.isstressbalance=1;
	md.transient.ismasstransport=1;
	md.transient.isgroundingline=1;
	md.groundingline.migration = 'SubelementMigration';
	md.transient.ismovingfront=1;
	
	% SMB
   md.smb = interpISMIP7GreenlandSMB(md, climate_model,'historical');

	% Frontal forcings set to 0, as front is prescribed with observed fronts in levelset
	md.frontalforcings.meltingrate = 0*ones(md.mesh.numberofvertices,1);
	md.calving.calvingrate=0*ones(md.mesh.numberofvertices,1);
	md.frontalforcings.ablationrate=md.calving.calvingrate+md.frontalforcings.meltingrate;

	% Prescribe evolving ice front through the ice levelset
	disp('   Initializing ice mask spclevelset');
	disp('Interpolating Greene ice fronts--->')
	% get observed calving front for all months
	ice_masks = interpMonthlyIceMaskGreene(md.mesh.x, md.mesh.y, [historical_start_time, historical_end_time]);

	disp('Processing Greene ice fronts--->')
	% Ensure no artificial advance
	ref_levelset = md.mask.ice_levelset;
	ref_levelset((md.geometry.bed<-200)&(md.mask.ice_levelset>0))=-1; % allow some advance for large glaciers but prevent slow/grounded from errorneously readvancing

	%for i = 2:size(ice_masks,2)
	for i = 1:size(ice_masks,2)
		fprintf('Processing front %d of %d\r', i, size(ice_masks, 2));
		ice_mask = ice_masks(1:end-1,i);
		pos = find(ice_mask<=0 & ref_levelset>0);
		ice_mask(pos)=1;

		%Reinitialize
		ice_masks(1:end-1,i) = reinitializelevelset(md,ice_mask);
	end
	% transient spc
	md.levelset.spclevelset = ice_masks;
	% fix the mask.ice_levelset to the first ice_mask 
	%(this doesn't cause ice cliff issues when the new front is upstream of the DEM front)
	md.mask.ice_levelset = ice_masks(1:end-1,1);

	% this avoids a bug (in theory, changing this should have no effect,
	% since the levelset does not advect and is entirely prescribed from 
	% the interpolation between the spclevelsets)
	md.levelset.reinit_frequency = 0;
	
	disp('   Prescribing basal melt rate');
	%Basal melt rate from https://tc.copernicus.org/articles/19/2695/2025/
	md.basalforcings=linearbasalforcings();
	md.basalforcings.deepwater_melting_rate=40.; % m/yr ice equivalent
	md.basalforcings.deepwater_elevation=-300;
	%	md.basalforcings.deepwater_elevation=-800;
	md.basalforcings.upperwater_melting_rate=0; % no melting for zb>=0
	md.basalforcings.upperwater_elevation=-100; % sea level
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1); % no melting on grounded ice
	
	% Set parameters
	md.inversion.iscontrol=0;

	savemodel(org,md);
end % }}}
if perform(org,['Greenland_ISMIP7Run_' climate_model '_Historical']),% {{{

	md=loadmodel(org,['Greenland_ISMIP7Prep_' climate_model '_Historical']);

	md.cluster=cluster;
	md.settings.stagingpath = [issmdir() '/execution'];
	md.toolkits.DefaultAnalysis=bcgslbjacobioptions();% biconjugate gradient with block Jacobi preconditioner
	md.settings.solver_residue_threshold = 1.e-4;
	md.verbose.solution=1;
	md.cluster = cluster;
	md.cluster.interactive = interactive;
	md.settings.waitonlock = 0;

	md.transient.requested_outputs={'default','IceVolume','IceVolumeScaled','IceVolumeAboveFloatation','IceVolumeAboveFloatationScaled','MaskIceLevelset','MaskOceanLevelset','SmbMassBalance','GroundedArea','GroundedAreaScaled','FloatingArea','FloatingAreaScaled','TotalSmb','TotalSmbScaled', 'BasalforcingsGroundediceMeltingRate', 'TotalGroundedBmb', 'TotalGroundedBmbScaled', 'BasalforcingsFloatingiceMeltingRate','TotalFloatingBmb', 'TotalFloatingBmbScaled', 'GroundinglineMassFlux'};
	
	md.timestepping.time_step=0.01;
	md.timestepping.start_time=historical_start_time;
	md.timestepping.final_time=historical_end_time;
	md.settings.output_frequency = 10;

	md.masstransport.stabilization = 2;

	md.miscellaneous.name = org.steps(org.currentstep).string;
	md=solve(md,'tr','runtimename',0,'loadonly',loadonly);

	if md.cluster.interactive == 1
		md=loadresultsfromcluster(md);
	end
	if (loadonly == 1) | (md.cluster.interactive == 1)
		savemodel(org,md);
	end
end%}}}

% Step 3-4: Projections (ctrl, ssp126 ssp370, ssp585 for each climate model)
if perform(org,['Greenland_ISMIP7Prep_' climate_model '_' upper(scenario)]),% {{{

	md=loadmodel(org,['Greenland_TransientCalibration']);
	
	% Set the transient parameters
	md.transient.isthermal=0;
	md.transient.isstressbalance=1;
	md.transient.ismasstransport=1;
	md.transient.ismovingfront=1;
	md.transient.isgroundingline=1;
	md.groundingline.migration = 'SubelementMigration';

	% Prescribe frontal calving and melt rate to be 0 for now
	md.frontalforcings.meltingrate = 0*ones(md.mesh.numberofvertices,1);
	md.calving.calvingrate=0*ones(md.mesh.numberofvertices,1);
	md.frontalforcings.ablationrate=md.calving.calvingrate+md.frontalforcings.meltingrate;

	% Prescribe evolving ice front through the ice levelset
	md.transient.ismovingfront=1;

	disp(['	Loading ' climate_model ' historical run']);
	tmp = loadmodel(['models/Greenland_ISMIP7Run_' climate_model '_Historical']);
	end_time = tmp.results.TransientSolution(end).time;
	end_step = size(tmp.results.TransientSolution, 2);
	if end_time ~= projection_start_time
		error(['Projection start time, ' num2str(projection_start_time) ' does not match historical end time, ' num2str(end_time) '.'])
	end
	disp(['    Resetting values and starting projection at ' num2str(projection_start_time)]);
	tmp = transientrestart(tmp,end_step);
	md.initialization.vx = tmp.initialization.vx;
	md.initialization.vy = tmp.initialization.vy;
	md.initialization.vz = tmp.initialization.vz;
	md.initialization.vel = tmp.initialization.vel;
	md.initialization.pressure = tmp.initialization.pressure;
	md.initialization.temperature = tmp.initialization.temperature;

	md.mask.ice_levelset = tmp.mask.ice_levelset;
	md.mask.ocean_levelset = tmp.mask.ocean_levelset;

	md.geometry.base = tmp.geometry.base;
	md.geometry.thickness = tmp.geometry.thickness;
	md.geometry.surface = tmp.geometry.surface;

	disp('   Prescribing SMB');
   md.smb = interpISMIP7GreenlandSMB(md, climate_model,scenario);
	
	disp('   Prescribing basal melt rate');
	%Basal melt rate from https://tc.copernicus.org/articles/19/2695/2025/
	md.basalforcings=linearbasalforcings();
	md.basalforcings.deepwater_melting_rate=40.; % m/yr ice equivalent
	md.basalforcings.deepwater_elevation=-300;
	%	md.basalforcings.deepwater_elevation=-800;
	md.basalforcings.upperwater_melting_rate=0; % no melting for zb>=0
	md.basalforcings.upperwater_elevation=-100; % sea level
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1); % no melting on grounded ice
		
	
	disp('   Prescribing frontal melt');
	% melting rate at the front is prescribed by ISMIP7 parameterization 
   md.frontalforcings = interpISMIP7GreenlandOcn(md, climate_model, scenario);

	disp('   Prescribing von mises calving');
	% Set the calving law to von Mises
	md.calving = calvingvonmises(); 		

	%Default calving threshold
	md.calving.stress_threshold_groundedice=1000*ones(md.mesh.numberofelements,1);
	%md.calving.stress_threshold_groundedice=3500*ones(md.mesh.numberofelements,1);
	md.calving.stress_threshold_floatingice=200*ones(md.mesh.numberofelements,1);

	% === CW Glaciers (DONE) ==== %{{{
	% CW_NONAME1
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==25), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% CW_NONAME2
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==26), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %3000; %1500; %DONE
	% CW_NONAME3
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==27), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% EQIP_SERMIA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==36), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 600; %700; %875; %DONE
	% JAKOBSHAVN_ISBRAE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==86), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1300; %1100; %950; %1200;
	md.calving.stress_threshold_floatingice(pos) = 500; %475; %400; %500; %300; %400;
	% KANGERLUARSUUP_SERMIA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==89), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1500; %1000; %DONE
	% KANGERLUSSUUP_SERMERSUA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==91), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% KANGILERNGATA_SERMIA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==93), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 800; %1500; %2000; %3000; %DONE
	% LILLE_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==111), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE
	% NORDENSKIOLD_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==139), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, grounding line above sea level for now
	% RINK_ISBRAE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==164), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3000; %DONE
	md.calving.stress_threshold_floatingice(pos) = 3000; %DONE
	% SAQQARLIUP_ALANGORLIUP
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==168), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 600; %800; %1000; %1500;%2000; %DONE
	% SERMEQ_AVANNARLEQ
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==181), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 500; %600; %700; %780; %DONE
	% SERMEQ_AVANNARLEQ2
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==182), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3000; %DONE
	% SERMEQ_KUJALLEQ
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==183), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2500; %DONE
	% SERMEQ_SILARLEQ
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==184), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2000; %1000;%3100; %DONE
	% STORE_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==207), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %2500; %1300; %DONE
	%}}}
	% === NW Glaciers (DONE) === %{{{

	% ALISON_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==6), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 650; %750; %DONE
	% BAMSE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==13), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% BOWDOIN
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==16), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% CARLOS
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==21), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% CORNELL_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==23), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 300; %650 %700; %800; %DONE
	% DIEBITSCH
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==30), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% DOCKER_SMITH_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==31), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2500; %3100; %DONE
	% DOCKER_SMITH_GLETSCHER_W
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==32), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% FARQUHAR_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==37), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% GABLE_MIRROR
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==43), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% GADE-MORELL
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==44), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3000; %DONE
	% HARALD_MOLTKE_BRAE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==60), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 200; %500; %600; %700; %DONE
	md.calving.stress_threshold_floatingice(pos) = 100;
	% HART
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==62), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% HART_W
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==63), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% HAYES_GLETSCHER_M_SS
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==64), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2000; %3000; %DONE
	% HAYES_GLETSCHER_N_NN
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==65), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 700; %800; %1200; %1500; %2000; %DONE reduce
	% HEILPRIN_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==66), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% HELLAND
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==70), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% HUBBARD
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==73), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% ILLULLIP_SERMIA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==79), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2000; %DONE
	% INNGIA_ISBRAE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==81), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2800; %2100; %1900; %1750; %DONE
	% ISSUUARSUIT_SERMIA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==83), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 800; %1000; %2000; %2500; %3100; %DONE
	% KAKIVFAAT_SERMIAT
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==88), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1200; %DONE
	% KJER_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==98), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 325; %350; %300; %450; %570; %630; %DONE
	% KONG_OSCAR_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==106), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2000; %DONE
	% LEIDY-MARIE-SERMIARSUPALUK
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==110), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% MEEHAN
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==115), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% MEEHAN_W
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==116), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% MELVILLE_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==117), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% MOHN_GLETSJER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==123), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 500; %1000; %1500; %2000; %DONE
	% MORRIS_JESUP
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==124), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2500; %DONE
	% MORRIS_JESUP_W
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==125), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NANSEN_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==127), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 550; %600; %800; %1000; %1300; %DONE
	% NONAME_NORTH_OSCAR
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==137), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 800; %DONE
	md.calving.stress_threshold_floatingice(pos) = 100;
	% NORDENSKIOLD_GLESCHER_NW
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==138), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 515; %500; %575; %DONE
	% NUNATAKASSAAP_SERMIA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==145), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2000; %DONE
	% NW_NONAME1
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==146), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NW_NONAME2
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==147), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NW_NONAME3
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==148), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2500; %DONE
	% NW_NONAME4
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==149), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% PITUGFIK
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==154), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 200; %400; %700; %900; %1050; %DONE
	% QEQERTARSUUP_SERMIA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==161), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 800; %DONE
	% RINK_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==163), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2000; %2500; %3100; %DONE
	% SAVISSUAQ
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==169), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 400; %DONE
	% SAVISSUAQ_UNNAMED1
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==171), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %1100; %DONE
	% SAVISSUAQ_UNNAMED2
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==172), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SAVISSUAQ_UNNAMED3
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==173), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SAVISSUAQ_UNNAMED4
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==174), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3000; %DONE
	% SAVISSUAQ_W
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==175), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2000; %DONE
	% SAVISSUAQ_WW
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==176), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %2000; %3000; %DONE
	% SAVISSUAQ_WWW
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==177), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 300; %350; %DONE
	% SAVISSUAQ_WWWW
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==178), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 115; %DONE
	% SAVISSUAQ_WWWWW
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==179), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 600; %900; %1100; %1150; %DONE
	% SHARP
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==196), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SHARP_W
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==197), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SIORARSUAQ
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==199), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% STEENSTRUP-DIETRICHSON
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==206), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 650; %600; %800; %DONE
	% SUN
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==211), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SUN_W
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==212), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SVERDRUP_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==213), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 525; %550; %575; %500; %700; %725; %DONE
	% TRACY_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==217), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% TUGTO
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==218), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% UMIAMMAKKU_ISBRAE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==220), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% UPERNAVIK_ISSTROM_C
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==238), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2700; %2500; %2200; %1750; %1500; %1250; %DONE
	% UPERNAVIK_ISSTROM_N
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==239), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 990; %DONE
	% UPERNAVIK_ISSTROM_S
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==240), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2000; %DONE
	% UPERNAVIK_ISSTROM_SS
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==241), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1500; %DONE
	% USSING_BRAEER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==242), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2500; %DONE
	% USSING_BRAEER_N
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==243), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% VERHOEFF
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==245), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% VERHOEFF_W
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==246), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% YNGVAR_NIELSEN_BRAE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==251), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 530; %575; %500; %400; %600; %DONE
	md.calving.stress_threshold_floatingice(pos) = 100; %300;
	% YNGVAR_NIELSEN_BRAE_W
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==252), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1500; %DONE
	%}}}
	% === NO Glaciers (DONE) === %{{{

	% ACADEMY
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==2), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3000; %DONE
	% BRIKKERNE_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==19), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3000; %DONE
	% DODGE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==33), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% HAGEN_BRAE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==59), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE
	% HARDER_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==61), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% HUMBOLDT_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==74), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 300; %350; %450; %600; %700; %800;
	md.calving.stress_threshold_floatingice(pos) = 50; %100;
	% JUNGERSEN_HENSON_NARAVANA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==87), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% MARIE_SOPHIE_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==114), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NEWMAN_BUGT
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==131), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NO_NONAME1
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==142), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NO_NONAME2
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==143), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NO_NONAME3
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==144), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% OSTENFELD_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==150), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE
	% PETERMANN_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==152), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %Domain does not cover shelf
	shelf_pos = find(ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y,'exp/iceshelves/Petermann_Shelf.exp','element',0));
	md.calving.stress_threshold_floatingice(union(pos,shelf_pos)) = 200; %100; %50; %200; %300; %500;
	% PETERMANN_GLETSCHER_N
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==153), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %1000;
	% RYDER_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==166), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	shelf_pos = find(ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y,'exp/iceshelves/Ryder_Shelf.exp','element',0));
	md.calving.stress_threshold_floatingice(union(pos,shelf_pos)) = 3100; %DONE
	% STEENSBY_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==205), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE
	% STORM
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==208), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	%}}}
	% === NE Glaciers (DONE) === %{{{

	% AB_DRACHMANN_GLETSCHER_L_BISTRUP_BRAE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==1), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000;
	md.calving.stress_threshold_floatingice(pos) = 100;%200;
	% ADMIRALTY_TREFORK_KRUSBR_BORGJKEL_PONY
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==3), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% ADOLF_HOEL
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==4), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% BLSEBR_GAMMEL_HELLERUP_GLETSJER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==14), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000;
	% GERARD_DE_GEER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==53), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% HISINGER_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==72), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% JAETTEGLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==85), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, grounding line above sea level for now
	% NE_NONAME1
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==132), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NIOGHALVFJERDSFJORDEN
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==134), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %Shelf outside of domain
	shelf_pos = find(ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y,'exp/iceshelves/Nioghalvfjerdsfjorden_Shelf.exp','element',0));
	md.calving.stress_threshold_floatingice(union(pos,shelf_pos)) = 400; %800; %2500; %3100;
	% NORDENSKIOLD_NE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==140), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %2000; %DONE
	% PASSAGE_CHARPENTIER_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==151), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SORANERBRAEEN-EINAR_MIKKELSEN-HEINKEL-TVEGEGLETSCHER-PASTERZE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==201), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, grounding line above sea level for now
	% STORSTROMMEN
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==209), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000;
	shelf_pos = find(ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y,'exp/iceshelves/Drachmann_Storstrommen_Shelf.exp','element',0));
	md.calving.stress_threshold_floatingice(union(pos, shelf_pos)) = 50; %200; %2000;
	% WAHLENBERG_VIOLINGLETSJER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==248), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% WALTERSHAUSEN
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==249), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% WORDIE-VIBEKE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==250), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% ZACHARIAE_ISSTROM
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==253), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1500; %1000; % Shelf not within domain
	shelf_pos = find(ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y,'exp/iceshelves/Zachariae_Isstrom_Shelf.exp','element',0));
	md.calving.stress_threshold_floatingice(union(pos,shelf_pos)) = 100; %200; %400 %;3100;
	%}}}
	% === CE Glaciers (DONE) === %{{{

	% BORGGRAVEN
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==15), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% BREDEGLETSJER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==18), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, grounding line above sea level for now
	% CHARCOT
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==22), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% COURTAULD
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==24), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% DAUGAARD-JENSEN
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==28), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 400; 
	% DENDRITGLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==29), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 400; %800; %1500; %2000; %DONE
	% EIELSON_HARE_FJORD-ROLIGE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==34), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2500; %DONE
	% FREDERIKSBORG_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==40), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% F_GRAAE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==42), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% GEIKIE1
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==45), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% GEIKIE2
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==46), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% GEIKIE3
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==47), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% GEIKIE4
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==48), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% GEIKIE5
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==49), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% GEIKIE6
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==50), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %200; %600 %800; %1000; %DONE
	% GEIKIE7
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==51), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 200; %600 %1000; %DONE
	% GEIKIE_UNNAMED_VESTFORD_S
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==52), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 600; %1500; %2000; %3000; %DONE
	% KANGERLUSSUAQ
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==90), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100;
	md.calving.stress_threshold_floatingice(pos) = 425; %450; %500; %600; %400; %200;
	% KIV_STEENSTRUP_NODRE_BRAE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==96), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 600; %500; %300; %700; %1000; %3100;
	% KIV_STEENSTRUP_SONDRE_BRAE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==97), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 100; %300 %700; %1000; %3100
	% KOLVEGLETSJER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==104), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% KONG_CHRISTIAN
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==105), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2000; %3100 %DONE
	md.calving.stress_threshold_floatingice(pos) = 200; %3000
	% KRONBORG
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==107), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% KRUUSE_FJORD
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==108), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 500; %1000; %2000; %3100 %DONE
	% LAUBE_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==109), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% MAGGA_DAN_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==113), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NORDFJORD
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==141), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% POLARIC-DECEPTION_O_N
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==155), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 400; %200;
	% ROSENBORG
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==165), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SORGENFRI
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==202), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% SORTEBRAE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==203), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% STYRTE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==210), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SYDBR
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==215), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_DECEPTION_N
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==225), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_DECEPTION_O_CN_CS
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==226), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1400; %1500; %1000; %2000; %3100; %DONE
	% UNNAMED_KANGER_E
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==230), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_KANGER_W
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==231), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% UNNAMED_LAUBE_S
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==232), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_POLARIC_C
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==233), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_POLARIC_S
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==234), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_SORGENFRI_W
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==235), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_UUNARTIT_ISLANDS
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==237), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 900; %1000; %1200; %700; %2000; %3100;
	% VESTFJORD
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==247), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	%}}}
	% === SE Glaciers (DONE) === %{{{

	% ANORITUUP_KANGERLUA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==7), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 650; %700; %500; %1000; %2000 %3100;
	% APUSEERAJIK
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==8), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% APUSEERSERPIA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==9), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% AP_BERNSTOFF_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==10), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2000; %1000; %3000; %DONE
	md.calving.stress_threshold_floatingice(pos) = 150; %100; %200;
	% BRCKNER_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==17), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% BUSSEMAND
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==20), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% FENRISGLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==38), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 3100;
	% FIMBULGETLSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==39), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% GLACIERDEFRANCE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==54), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %500; %1000; %2000; %3100; %DONE
	% GRAULV
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==55), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% GYLDENLOVE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==56), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% GYLDENLOVE_S
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==57), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 800; %DONE, completely above sea level
	% GYLDENLOVE_SS
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==58), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 800; %DONE, completely above sea level
	% HEIMDAL_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==67), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% HEIM_GLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==68), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% HELHEIMGLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==69), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1800; %1500; %1150; %3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 800; %500; %1150; %3000;
	% HERLUF_TROLLE-KANGERLULUK-DANELL
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==71), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% IKERTIVAQ_M
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==75), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% IKERTIVAQ_N
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==76), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %2000; %3100; %DONE
	% IKERTIVAQ_NN
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==77), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %2000; %3100; %DONE
	% IKERTIVAQ_S
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==78), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 600; %DONE
	% KNUD-RASMUSSEN
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==99), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% KOGE_BUGT_C
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==100), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% KOGE_BUGT_N
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==101), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% KOGE_BUGT_S
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==102), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2000; %3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 3000;
	% KOGE_BUGT_SS
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==103), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% MAELKEVEJEN
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==112), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% MIDGARDGLETSCHER
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==118), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% MOGENS_HEINESEN_C
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==119), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% MOGENS_HEINESEN_N
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==120), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% MOGENS_HEINESEN_S
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==121), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% MOGENS_HEINESEN_SS_SSS
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==122), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NAPASORSUAQ_C_S
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==128), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 500; %1000; %2000; %3100; %DONE
	% NAPASORSUAQ_N
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==129), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	%md.calving.stress_threshold_floatingice(pos) = 3000;
	% NIGERTULUUP_KATTILERTARPIA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==133), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NONAME_IKERTIVAQ_N
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==135), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 300; %600; %800; %1000; %DONE
	% NONAME_IKERTIVAQ_S
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==136), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% PUISORTOQ_N
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==156), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; 
	md.calving.stress_threshold_floatingice(pos) = 1500; %1800; %2000; %3000; %200;
	% PUISORTOQ_S
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==157), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2400; %2300; %1800; %1300; %1000;
	md.calving.stress_threshold_floatingice(pos) = 175; %150; %200;
	% RIMFAXE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==162), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% SE_NONAME1
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==187), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SE_NONAME10
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==188), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SE_NONAME2
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==189), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% SE_NONAME4
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==190), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SE_NONAME5
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==191), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SE_NONAME6
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==192), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SE_NONAME7
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==193), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SE_NONAME8
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==194), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SE_NONAME9
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==195), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SKINFAXE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==200), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% SOUTHERN_TIP
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==204), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% TINGMIARMIUT_FJORD
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==216), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %2800; %2500; %3100;
	md.calving.stress_threshold_floatingice(pos) = 500; %1000; %2000; %150; %3000;
	% UMIIVIK_FJORD
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==221), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 600; %800; %1000; %3100;
	md.calving.stress_threshold_floatingice(pos) = 100; %150; %200;
	% UNNAMED_ANORITUUP_KANGERLUA_S
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==222), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_ANORITUUP_KANGERLUA_SS
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==223), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_DANELL_FJORD
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==224), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% UNNAMED_HERLUF_TROLLE_N
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==227), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %2800; %2500; %2000; %3100;
	% UNNAMED_HERLUF_TROLLE_S
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==228), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% UNNAMED_KANGERLULUK
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==229), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_SOUTH_DANELL_FJORD
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==236), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	%}}}
	% === SW Glaciers (DONE) === %{{{

	% AKULLERSUUP-QAMANAARSUUP
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==5), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %2000; %2500; %DONE
	% AVANNARLEQ-NIGERLIKASIK
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==11), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% AVANNARLEQ_N
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==12), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% EQALORUTSIT_KILLIIT_SERMIAT
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==35), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 150; %DONE, completely above sea level
	% FREDERIKSHABS-NAKKAASORSUAQ
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==41), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% ILORLIIT-SERMINNGUAQ
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==80), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% INUPPAAT_QUUAT
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==82), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% ISUNNGUATA-RUSSELL
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==84), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% KANGIATA_NUNAATA_SERMIA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==92), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2000; %DONE
	% KANGILINNGUATA_SERMIA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==94), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% KIATTUUT-QOOQQUP
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==95), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %2900; %2700; %2000; %3100; %DONE
	% NAAJAT_SERMIAT
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==126), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NARSAP_SERMIA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==130), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 420; %DONE, grounding line is above sea level for now
	% QAJUUTTAP_SERMIA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==158), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 2500; %3100; %DONE
	% QAJUUTTAP_SERMIA_N
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==159), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% QALERALLIT_SERMIAT
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==160), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SAQQAP-MAJORQAQ-SOUTHTERRUSSEL_SOUTHQUARUSSEL
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==167), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SERMEQ-KANGAASARSUUP
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==180), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SERMILIGAARSSUK_BRAE
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==185), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 800; %1000; %1500; %2000; %DONE
	% SERMILIK
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==186), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SIORALIK-ARSUK-QIPISAQQU
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==198), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SW_NONAME1
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==214), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UKAASORSUAQ
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==219), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% USULLUUP_SERMIA
	[tmp,pos] = ismember(intersect(find(md.miscellaneous.dummy.catchment_id==244), md.mesh.extractedelements), md.mesh.extractedelements); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, grounding line is above sea level for now
	%}}}
	%from kPa to Pa
	md.calving.stress_threshold_groundedice = md.calving.stress_threshold_groundedice*1000;
	md.calving.stress_threshold_floatingice = md.calving.stress_threshold_floatingice*1000;

	disp('	Adjusting ice front levelsets')
	md.levelset.spclevelset        = NaN(md.mesh.numberofvertices,1);
	disp('	--- Setting spclevelset for Petermann cliffs around shelf')
	cliff_pos = find(ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y,...
		'exp/spclevelset/Petermann_spcLevelset.exp','node',0));
	md.levelset.spclevelset(cliff_pos) = -1;
	disp('	--- Setting spclevelset for 79N cliffs around shelf')
	cliff_pos = find(ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y,...
		'exp/spclevelset/79N_spcLevelset.exp','node',0));
	md.levelset.spclevelset(cliff_pos) = -1;
   md.levelset.migration_max      = 10e5; %max possible ice front retreat rate (10 km/yr)

	savemodel(org,md);
end%}}}
if perform(org,['Greenland_ISMIP7Run_' climate_model '_' upper(scenario) '_' int2str(projection_start_time) '-' int2str(projection_end_time)]),% {{{

	md=loadmodel(org,['Greenland_ISMIP7Prep_' climate_model '_' upper(scenario)]);
	
	disp('   Updating start/end times');
	md.timestepping.time_step=0.01;
   md.timestepping.start_time=projection_start_time;
   md.timestepping.final_time=projection_end_time;
	md.settings.output_frequency=10;
	
	disp('   Removing unused SMB forcings');
	t = md.smb.mass_balance(end, :);                       % time stamps
	mask = (t >= projection_start_time) & (t < projection_end_time + 1);
	md.smb.mass_balance = md.smb.mass_balance(:, mask);

	md.verbose.solution = true;
	%md.verbose=verbose('all');

	md.transient.requested_outputs={'default','IceVolume','IceVolumeScaled','IceVolumeAboveFloatation','IceVolumeAboveFloatationScaled','MaskIceLevelset','MaskOceanLevelset','SmbMassBalance','GroundedArea','GroundedAreaScaled','FloatingArea','FloatingAreaScaled','TotalSmb','TotalSmbScaled', 'BasalforcingsGroundediceMeltingRate', 'TotalGroundedBmb', 'TotalGroundedBmbScaled', 'BasalforcingsFloatingiceMeltingRate','TotalFloatingBmb', 'TotalFloatingBmbScaled', 'TotalCalvingFluxLevelset', 'TotalCalvingMeltingFluxLevelset','Calvingratex','Calvingratey','CalvingMeltingrate','GroundinglineMassFlux'};
	
	md.cluster=cluster;

	md.toolkits.DefaultAnalysis=bcgslbjacobioptions();% biconjugate gradient with block Jacobi preconditioner
	md.settings.solver_residue_threshold = 1.e-4;
	md.masstransport.stabilization = 2;

	md.miscellaneous.name = org.steps(org.currentstep).string;
	md.cluster.interactive=interactive;
	md.settings.waitonlock=0;

	md=solve(md,'tr','runtimename',0,'loadonly',loadonly);

	if md.cluster.interactive == 1
		md=loadresultsfromcluster(md);
	end

	if (loadonly == 1) | (md.cluster.interactive == 1)
		savemodel(org,md);
	end

end%}}}

