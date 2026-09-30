steps=[5];
%Run Options{{{

% Cluster Options (for forward transients and inversions)
cluster = 1; % 0 for 'totten', 1 for 'andes', 2 for epica
clusternames = {'totten', 'andes', 'local'};
clustername = clusternames{cluster+1};

% Run parameters
max_run_time = 168; % hours
mem_alloc = 128; % in Gb
num_cpus = 64; % number of cores
interactive = 1;
loadonly = 0;

% Run start and end times
historical_start_time = 2007;
historical_end_time = 2015;
ocx_start_time = 2007;
ocx_end_time = 2026;
projection_start_time = 2015;
projection_end_time = 2301;

% Climate Models
climate_models = {'CESM2-WACCM', 'MRI-ESM2-0'};
climate_model = climate_models{2}
% Forecast Parameters
scenarios = {'ssp585', 'ssp370', 'ssp126', 'ctrl'};
%scenarios = {'ssp585', 'ssp126', 'ctrl'};
scenario = scenarios{4};
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

	md.verbose.solution = true;

	md=solve(md,'tr','runtimename',0,'loadonly',loadonly);

	if md.cluster.interactive == 1
		md=loadresultsfromcluster(md);
	end
	if (loadonly == 1) | (md.cluster.interactive == 1)
		savemodel(org,md);
	end
end%}}}
		
if perform(org,['Greenland_ISMIP7Prep_' climate_model '_' upper(scenario)]),% {{{

	md=loadmodel(org,['Greenland_TransientCalibration']);

	% Set the transient parameters
	md.transient.isthermal=0;
	md.transient.isstressbalance=1;
	md.transient.ismasstransport=1;
	md.transient.ismovingfront=1;
	md.transient.isgroundingline=1;
	md.groundingline.migration = 'SubelementMigration';
	md.groundingline.nomelt_under_lakes = 1;

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

	% === CW Glaciers (DONE) === %{{{
	% CW_NONAME1
	pos = find(md.miscellaneous.dummy.catchment_id==25); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% CW_NONAME2
	pos = find(md.miscellaneous.dummy.catchment_id==26); md.calving.stress_threshold_groundedice(pos) = 3100; %3000; %1500; %DONE
	% CW_NONAME3
	pos = find(md.miscellaneous.dummy.catchment_id==27); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% EQIP_SERMIA
	pos = find(md.miscellaneous.dummy.catchment_id==36); md.calving.stress_threshold_groundedice(pos) = 600; %700; %875; %DONE
	% JAKOBSHAVN_ISBRAE
	pos = find(md.miscellaneous.dummy.catchment_id==86); md.calving.stress_threshold_groundedice(pos) = 1450; % % % %1425;
	md.calving.stress_threshold_floatingice(pos) = 442; %445; %440; %430; %425; %400; %450; %460; %475; %500; %550; %450;  % % %
	side_pos = find(ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y,'exp/iceshelves/Jakobshavn_Side.exp','element',0));
	md.calving.stress_threshold_groundedice(side_pos) = 1650; %1600; % % %1500; %1400;
	md.calving.stress_threshold_floatingice(side_pos) = 372; %371; %372; %370; %375; %350; %300; %400; %200; %150; %100;
	% KANGERLUARSUUP_SERMIA
	pos = find(md.miscellaneous.dummy.catchment_id==89); md.calving.stress_threshold_groundedice(pos) = 3100; %2500;  %1500; %1000; %DONE
	% KANGERLUSSUUP_SERMERSUA
	pos = find(md.miscellaneous.dummy.catchment_id==91); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% KANGILERNGATA_SERMIA
	pos = find(md.miscellaneous.dummy.catchment_id==93); md.calving.stress_threshold_groundedice(pos) = 2000; %1200;  %800; %1500; %2000; %3000; %DONE
	md.calving.stress_threshold_floatingice(pos) = 300;
	% LILLE_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==111); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE
	% NORDENSKIOLD_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==139); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, grounding line above sea level for now
	% RINK_ISBRAE
	pos = find(md.miscellaneous.dummy.catchment_id==164); md.calving.stress_threshold_groundedice(pos) = 3000; %DONE
	md.calving.stress_threshold_floatingice(pos) = 3000; %DONE
	% SAQQARLIUP_ALANGORLIUP
	pos = find(md.miscellaneous.dummy.catchment_id==168); md.calving.stress_threshold_groundedice(pos) = 3100; %2500; %2000; %1200; %800;  %600; %800; %1000; %1500;%2000; %DONE
	md.calving.stress_threshold_floatingice(pos) = 500; %300;
	% SERMEQ_AVANNARLEQ
	pos = find(md.miscellaneous.dummy.catchment_id==181); md.calving.stress_threshold_groundedice(pos) = 400; %300;   %500; %600; %700; %780; %DONE
	% SERMEQ_AVANNARLEQ2
	pos = find(md.miscellaneous.dummy.catchment_id==182); md.calving.stress_threshold_groundedice(pos) = 3000; %DONE
	% SERMEQ_KUJALLEQ
	pos = find(md.miscellaneous.dummy.catchment_id==183); md.calving.stress_threshold_groundedice(pos) = 1500;     %2500; %DONE
	% SERMEQ_SILARLEQ
	pos = find(md.miscellaneous.dummy.catchment_id==184); md.calving.stress_threshold_groundedice(pos) = 2500; %3000;    %2000; %1000;%3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 500; %300;
	% STORE_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==207); md.calving.stress_threshold_groundedice(pos) = 3100; %2500; %1300; %DONE
	md.calving.stress_threshold_floatingice(pos) = 300;
	%}}}
	% === NW Glaciers (DONE) === %{{{

	% ALISON_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==6); md.calving.stress_threshold_groundedice(pos) = 795; %785; %775; %800; %825; %750; %900;     %650; %750; %DONE
	% BAMSE
	pos = find(md.miscellaneous.dummy.catchment_id==13); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% BOWDOIN
	pos = find(md.miscellaneous.dummy.catchment_id==16); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% CARLOS
	pos = find(md.miscellaneous.dummy.catchment_id==21); md.calving.stress_threshold_groundedice(pos) = 1000; %2000;  %3100; %DONE
	% CORNELL_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==23); md.calving.stress_threshold_groundedice(pos) = 600;    %300; %650 %700; %800; %DONE
	% DIEBITSCH
	pos = find(md.miscellaneous.dummy.catchment_id==30); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% DOCKER_SMITH_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==31); md.calving.stress_threshold_groundedice(pos) = 2500; %3100; %DONE
	% DOCKER_SMITH_GLETSCHER_W
	pos = find(md.miscellaneous.dummy.catchment_id==32); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% FARQUHAR_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==37); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% GABLE_MIRROR
	pos = find(md.miscellaneous.dummy.catchment_id==43); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% GADE-MORELL
	pos = find(md.miscellaneous.dummy.catchment_id==44); md.calving.stress_threshold_groundedice(pos) = 3000; %DONE
	% HARALD_MOLTKE_BRAE
	pos = find(md.miscellaneous.dummy.catchment_id==60); md.calving.stress_threshold_groundedice(pos) = 200; %500; %600; %700; %DONE
	md.calving.stress_threshold_floatingice(pos) = 100;
	% HART
	pos = find(md.miscellaneous.dummy.catchment_id==62); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% HART_W
	pos = find(md.miscellaneous.dummy.catchment_id==63); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% HAYES_GLETSCHER_M_SS
	pos = find(md.miscellaneous.dummy.catchment_id==64); md.calving.stress_threshold_groundedice(pos) = 2000; %3000; %DONE
	% HAYES_GLETSCHER_N_NN
	pos = find(md.miscellaneous.dummy.catchment_id==65); md.calving.stress_threshold_groundedice(pos) = 700; %800; %1200; %1500; %2000; %DONE
	% HEILPRIN_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==66); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% HELLAND
	pos = find(md.miscellaneous.dummy.catchment_id==70); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% HUBBARD
	pos = find(md.miscellaneous.dummy.catchment_id==73); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% ILLULLIP_SERMIA
	pos = find(md.miscellaneous.dummy.catchment_id==79); md.calving.stress_threshold_groundedice(pos) = 2000; %DONE
	% INNGIA_ISBRAE
	pos = find(md.miscellaneous.dummy.catchment_id==81); md.calving.stress_threshold_groundedice(pos) = 3100;   %2800; %2100; %1900; %1750; %DONE
	md.calving.stress_threshold_floatingice(pos) = 1000; %900; %700; %500; %300;
	% ISSUUARSUIT_SERMIA
	pos = find(md.miscellaneous.dummy.catchment_id==83); md.calving.stress_threshold_groundedice(pos) = 1100;   %800; %1000; %2000; %2500; %3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 400; %300;
	% KAKIVFAAT_SERMIAT
	pos = find(md.miscellaneous.dummy.catchment_id==88); md.calving.stress_threshold_groundedice(pos) = 1350; %1400; %1600; %1200; %DONE
	md.calving.stress_threshold_floatingice(pos) = 250; %300;
	% KJER_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==98); md.calving.stress_threshold_groundedice(pos) = 475; %450; %550; %400; %700;    %325; %350; %300; %450; %570; %630; %DONE
	% KONG_OSCAR_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==106); md.calving.stress_threshold_groundedice(pos) = 2000; %DONE
	% LEIDY-MARIE-SERMIARSUPALUK
	pos = find(md.miscellaneous.dummy.catchment_id==110); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% MEEHAN
	pos = find(md.miscellaneous.dummy.catchment_id==115); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% MEEHAN_W
	pos = find(md.miscellaneous.dummy.catchment_id==116); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% MELVILLE_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==117); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% MOHN_GLETSJER
	pos = find(md.miscellaneous.dummy.catchment_id==123); md.calving.stress_threshold_groundedice(pos) = 1000; %800;  %500; %1000; %1500; %2000; %DONE
	md.calving.stress_threshold_floatingice(pos) = 500; %300;
	% MORRIS_JESUP
	pos = find(md.miscellaneous.dummy.catchment_id==124); md.calving.stress_threshold_groundedice(pos) = 2500; %DONE
	% MORRIS_JESUP_W
	pos = find(md.miscellaneous.dummy.catchment_id==125); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NANSEN_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==127); md.calving.stress_threshold_groundedice(pos) = 850; %800;      %550; %600; %800; %1000; %1300; %DONE
	% NONAME_NORTH_OSCAR
	pos = find(md.miscellaneous.dummy.catchment_id==137); md.calving.stress_threshold_groundedice(pos) = 850; %900; %1000;  %800; %DONE
	% NORDENSKIOLD_GLESCHER_NW
	pos = find(md.miscellaneous.dummy.catchment_id==138); md.calving.stress_threshold_groundedice(pos) = 550; %600; %650; %700; %515; %500; %575; %DONE
	% NUNATAKASSAAP_SERMIA
	pos = find(md.miscellaneous.dummy.catchment_id==145); md.calving.stress_threshold_groundedice(pos) = 2000; %DONE
	% NW_NONAME1
	pos = find(md.miscellaneous.dummy.catchment_id==146); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NW_NONAME2
	pos = find(md.miscellaneous.dummy.catchment_id==147); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NW_NONAME3
	pos = find(md.miscellaneous.dummy.catchment_id==148); md.calving.stress_threshold_groundedice(pos) = 2500; %DONE
	% NW_NONAME4
	pos = find(md.miscellaneous.dummy.catchment_id==149); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% PITUGFIK
	pos = find(md.miscellaneous.dummy.catchment_id==154); md.calving.stress_threshold_groundedice(pos) = 3100; %2000; %1400; %1200; %800; %600; %400;  %200; %400; %700; %900; %1050;
	md.calving.stress_threshold_floatingice(pos) = 300;
	% QEQERTARSUUP_SERMIA
	pos = find(md.miscellaneous.dummy.catchment_id==161); md.calving.stress_threshold_groundedice(pos) = 800; %DONE
	% RINK_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==163); md.calving.stress_threshold_groundedice(pos) = 2000; %2500; %3100; %DONE
	% SAVISSUAQ
	pos = find(md.miscellaneous.dummy.catchment_id==169); md.calving.stress_threshold_groundedice(pos) = 600; %800;   %400; %DONE
	% SAVISSUAQ_UNNAMED1
	pos = find(md.miscellaneous.dummy.catchment_id==171); md.calving.stress_threshold_groundedice(pos) = 1000; %1100; %DONE
	% SAVISSUAQ_UNNAMED2
	pos = find(md.miscellaneous.dummy.catchment_id==172); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SAVISSUAQ_UNNAMED3
	pos = find(md.miscellaneous.dummy.catchment_id==173); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SAVISSUAQ_UNNAMED4
	pos = find(md.miscellaneous.dummy.catchment_id==174); md.calving.stress_threshold_groundedice(pos) = 3000; %DONE
	% SAVISSUAQ_W
	pos = find(md.miscellaneous.dummy.catchment_id==175); md.calving.stress_threshold_groundedice(pos) = 3100;     %2000; %DONE
	md.calving.stress_threshold_floatingice(pos) = 300;
	% SAVISSUAQ_WW
	pos = find(md.miscellaneous.dummy.catchment_id==176); md.calving.stress_threshold_groundedice(pos) = 1000; %2000; %3000; %DONE
	% SAVISSUAQ_WWW
	pos = find(md.miscellaneous.dummy.catchment_id==177); md.calving.stress_threshold_groundedice(pos) = 300; %350; %DONE
	% SAVISSUAQ_WWWW
	pos = find(md.miscellaneous.dummy.catchment_id==178); md.calving.stress_threshold_groundedice(pos) = 170; %160; %175; %115; %200; %400; %600;   %115; %DONE
	% SAVISSUAQ_WWWWW
	pos = find(md.miscellaneous.dummy.catchment_id==179); md.calving.stress_threshold_groundedice(pos) = 600; %900; %1100; %1150; %DONE
	% SHARP
	pos = find(md.miscellaneous.dummy.catchment_id==196); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SHARP_W
	pos = find(md.miscellaneous.dummy.catchment_id==197); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SIORARSUAQ
	pos = find(md.miscellaneous.dummy.catchment_id==199); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% STEENSTRUP-DIETRICHSON
	pos = find(md.miscellaneous.dummy.catchment_id==206); md.calving.stress_threshold_groundedice(pos) = 705; %800;   %650; %600; %800; %DONE
	% SUN
	pos = find(md.miscellaneous.dummy.catchment_id==211); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SUN_W
	pos = find(md.miscellaneous.dummy.catchment_id==212); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SVERDRUP_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==213); md.calving.stress_threshold_groundedice(pos) = 795; %785; %800; %725; %700; %725;   %525; %550; %575; %500; %700; %725; %DONE
	md.calving.stress_threshold_floatingice(pos) = 250; %300;
	% TRACY_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==217); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 350; %250;
	% TUGTO
	pos = find(md.miscellaneous.dummy.catchment_id==218); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% UMIAMMAKKU_ISBRAE
	pos = find(md.miscellaneous.dummy.catchment_id==220); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% UPERNAVIK_ISSTROM_C
	pos = find(md.miscellaneous.dummy.catchment_id==238); md.calving.stress_threshold_groundedice(pos) = 2700; %2500; %2200; %1750; %1500; %1250; %DONE
	% UPERNAVIK_ISSTROM_N
	pos = find(md.miscellaneous.dummy.catchment_id==239); md.calving.stress_threshold_groundedice(pos) = 1250; %1300; %1200;      %990; %DONE
	md.calving.stress_threshold_floatingice(pos) = 500; %300;
	% UPERNAVIK_ISSTROM_S
	pos = find(md.miscellaneous.dummy.catchment_id==240); md.calving.stress_threshold_groundedice(pos) = 3100;     %2000; %DONE
	md.calving.stress_threshold_floatingice(pos) = 400; %500; %300;
	% UPERNAVIK_ISSTROM_SS
	pos = find(md.miscellaneous.dummy.catchment_id==241); md.calving.stress_threshold_groundedice(pos) = 2500;  %1500; %DONE
	md.calving.stress_threshold_floatingice(pos) = 300;
	% USSING_BRAEER
	pos = find(md.miscellaneous.dummy.catchment_id==242); md.calving.stress_threshold_groundedice(pos) = 2500; %DONE
	% USSING_BRAEER_N
	pos = find(md.miscellaneous.dummy.catchment_id==243); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 300;
	% VERHOEFF
	pos = find(md.miscellaneous.dummy.catchment_id==245); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% VERHOEFF_W
	pos = find(md.miscellaneous.dummy.catchment_id==246); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% YNGVAR_NIELSEN_BRAE
	pos = find(md.miscellaneous.dummy.catchment_id==251); md.calving.stress_threshold_groundedice(pos) = 1100; %900; %700; %600;  %530; %575; %500; %400; %600; %DONE
	md.calving.stress_threshold_floatingice(pos) = 150; %200; %100; %300;
	% YNGVAR_NIELSEN_BRAE_W
	pos = find(md.miscellaneous.dummy.catchment_id==252); md.calving.stress_threshold_groundedice(pos) = 3100; %2500;    %1500; %DONE
	md.calving.stress_threshold_floatingice(pos) = 900; %700; %500; %300;
	%}}}
	% === NO Glaciers (DONE) === %{{{

	% ACADEMY
	pos = find(md.miscellaneous.dummy.catchment_id==2); md.calving.stress_threshold_groundedice(pos) = 3000; %DONE
	% BRIKKERNE_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==19); md.calving.stress_threshold_groundedice(pos) = 2500; %2000;  %3000; %DONE
	md.calving.stress_threshold_floatingice(pos) = 300;
	% DODGE
	pos = find(md.miscellaneous.dummy.catchment_id==33); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% HAGEN_BRAE
	pos = find(md.miscellaneous.dummy.catchment_id==59); md.calving.stress_threshold_groundedice(pos) = 1500;   %1001; %DONE
	% HARDER_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==61); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% HUMBOLDT_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==74); md.calving.stress_threshold_groundedice(pos) = 500;    %300; %350; %450; %600; %700; %800; %DONE
	md.calving.stress_threshold_floatingice(pos) = 100;   %50
	% JUNGERSEN_HENSON_NARAVANA
	pos = find(md.miscellaneous.dummy.catchment_id==87); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% MARIE_SOPHIE_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==114); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NEWMAN_BUGT
	pos = find(md.miscellaneous.dummy.catchment_id==131); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NO_NONAME1
	pos = find(md.miscellaneous.dummy.catchment_id==142); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NO_NONAME2
	pos = find(md.miscellaneous.dummy.catchment_id==143); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NO_NONAME3
	pos = find(md.miscellaneous.dummy.catchment_id==144); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% OSTENFELD_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==150); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE
	% PETERMANN_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==152); md.calving.stress_threshold_groundedice(pos) = 1000; %Domain does not cover shelf
	shelf_pos = find(ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y,'exp/iceshelves/Petermann_Shelf.exp','element',0));
	md.calving.stress_threshold_floatingice(union(pos,shelf_pos)) = 30; %50; %200; %DONE
	% PETERMANN_GLETSCHER_N
	pos = find(md.miscellaneous.dummy.catchment_id==153); md.calving.stress_threshold_groundedice(pos) = 1000; %3100; %1000;
	% RYDER_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==166); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	shelf_pos = find(ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y,'exp/iceshelves/Ryder_Shelf.exp','element',0));
	md.calving.stress_threshold_floatingice(union(pos,shelf_pos)) = 3100; %DONE
	% STEENSBY_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==205); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE
	% STORM
	pos = find(md.miscellaneous.dummy.catchment_id==208); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	%}}}
	% === NE Glaciers (DONE) === %{{{

	% AB_DRACHMANN_GLETSCHER_L_BISTRUP_BRAE
	pos = find(md.miscellaneous.dummy.catchment_id==1); md.calving.stress_threshold_groundedice(pos) = 900;  %1000; %DONE
	md.calving.stress_threshold_floatingice(pos) = 100;%200;
	% ADMIRALTY_TREFORK_KRUSBR_BORGJKEL_PONY
	pos = find(md.miscellaneous.dummy.catchment_id==3); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% ADOLF_HOEL
	pos = find(md.miscellaneous.dummy.catchment_id==4); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% BLSEBR_GAMMEL_HELLERUP_GLETSJER
	pos = find(md.miscellaneous.dummy.catchment_id==14); md.calving.stress_threshold_groundedice(pos) = 1000;
	% GERARD_DE_GEER
	pos = find(md.miscellaneous.dummy.catchment_id==53); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% HISINGER_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==72); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% JAETTEGLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==85); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, grounding line above sea level for now
	% NE_NONAME1
	pos = find(md.miscellaneous.dummy.catchment_id==132); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NIOGHALVFJERDSFJORDEN
	pos = find(md.miscellaneous.dummy.catchment_id==134); md.calving.stress_threshold_groundedice(pos) = 1000; %3100; %Shelf outside of domain (here I only apply the shelf value to the current shelf)
	shelf_pos = find(ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y,'exp/iceshelves/Nioghalvfjerdsfjorden_Shelf.exp','element',0));
	md.calving.stress_threshold_floatingice(shelf_pos) = 30; %50; %100;     %400; %800; %2500; %3100;
	% NORDENSKIOLD_NE
	pos = find(md.miscellaneous.dummy.catchment_id==140); md.calving.stress_threshold_groundedice(pos) = 3100; %2000; %DONE
	md.calving.stress_threshold_floatingice(union(pos,shelf_pos)) = 300;    %400; %800; %2500; %3100;
	% PASSAGE_CHARPENTIER_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==151); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SORANERBRAEEN-EINAR_MIKKELSEN-HEINKEL-TVEGEGLETSCHER-PASTERZE
	pos = find(md.miscellaneous.dummy.catchment_id==201); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, grounding line above sea level for now
	% STORSTROMMEN
	pos = find(md.miscellaneous.dummy.catchment_id==209); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, again only outlined shelf
	shelf_pos = find(ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y,'exp/iceshelves/Drachmann_Storstrommen_Shelf.exp','element',0));
	md.calving.stress_threshold_floatingice(shelf_pos) = 10; %20; %30;      %50; %200; %2000;
	% WAHLENBERG_VIOLINGLETSJER
	pos = find(md.miscellaneous.dummy.catchment_id==248); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% WALTERSHAUSEN
	pos = find(md.miscellaneous.dummy.catchment_id==249); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% WORDIE-VIBEKE
	pos = find(md.miscellaneous.dummy.catchment_id==250); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% ZACHARIAE_ISSTROM
	pos = find(md.miscellaneous.dummy.catchment_id==253); md.calving.stress_threshold_groundedice(pos) = 2000;     %1500; %1000; % Shelf not within domain
	shelf_pos = find(ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y,'exp/iceshelves/Zachariae_Isstrom_Shelf.exp','element',0));
	md.calving.stress_threshold_floatingice(union(pos,shelf_pos)) = 175; %200;     %100; %200; %400 %;3100; %DONE
	%}}}
	% === CE Glaciers (DONE) === %{{{

	% BORGGRAVEN
	pos = find(md.miscellaneous.dummy.catchment_id==15); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% BREDEGLETSJER
	pos = find(md.miscellaneous.dummy.catchment_id==18); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, grounding line above sea level for now
	% CHARCOT
	pos = find(md.miscellaneous.dummy.catchment_id==22); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% COURTAULD
	pos = find(md.miscellaneous.dummy.catchment_id==24); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% DAUGAARD-JENSEN
	pos = find(md.miscellaneous.dummy.catchment_id==28); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 400;
	% DENDRITGLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==29); md.calving.stress_threshold_groundedice(pos) = 1500; %1000;  %400; %800; %1500; %2000; %DONE
	% EIELSON_HARE_FJORD-ROLIGE
	pos = find(md.miscellaneous.dummy.catchment_id==34); md.calving.stress_threshold_groundedice(pos) = 3100;   %2500; %DONE
	% FREDERIKSBORG_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==40); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% F_GRAAE
	pos = find(md.miscellaneous.dummy.catchment_id==42); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% GEIKIE1
	pos = find(md.miscellaneous.dummy.catchment_id==45); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% GEIKIE2
	pos = find(md.miscellaneous.dummy.catchment_id==46); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% GEIKIE3
	pos = find(md.miscellaneous.dummy.catchment_id==47); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% GEIKIE4
	pos = find(md.miscellaneous.dummy.catchment_id==48); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% GEIKIE5
	pos = find(md.miscellaneous.dummy.catchment_id==49); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% GEIKIE6
	pos = find(md.miscellaneous.dummy.catchment_id==50); md.calving.stress_threshold_groundedice(pos) = 500; %600; %800;    %1000; %200; %600 %800; %1000; %DONE
	% GEIKIE7
	pos = find(md.miscellaneous.dummy.catchment_id==51); md.calving.stress_threshold_groundedice(pos) = 600;    %200; %600 %1000; %DONE
	% GEIKIE_UNNAMED_VESTFORD_S
	pos = find(md.miscellaneous.dummy.catchment_id==52); md.calving.stress_threshold_groundedice(pos) = 2000; %1500;  %600; %1500; %2000; %3000; %DONE
	% KANGERLUSSUAQ
	pos = find(md.miscellaneous.dummy.catchment_id==90); md.calving.stress_threshold_groundedice(pos) = 2300; %2000; %2500; %3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 350; %400; %300;   %%425; %450; %500; %600; %400; %200;
	% KIV_STEENSTRUP_NODRE_BRAE
	pos = find(md.miscellaneous.dummy.catchment_id==96); md.calving.stress_threshold_groundedice(pos) = 1000; %1500;     %600; %500; %300; %700; %1000; %3100; %DONE
	% KIV_STEENSTRUP_SONDRE_BRAE
	pos = find(md.miscellaneous.dummy.catchment_id==97); md.calving.stress_threshold_groundedice(pos) = 1800; %1500;  %100; %300 %700; %1000; %3100 %DONE
	md.calving.stress_threshold_floatingice(pos) = 250; %300;
	% KOLVEGLETSJER
	pos = find(md.miscellaneous.dummy.catchment_id==104); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% KONG_CHRISTIAN
	pos = find(md.miscellaneous.dummy.catchment_id==105); md.calving.stress_threshold_groundedice(pos) = 2000; %3100 %DONE
	md.calving.stress_threshold_floatingice(pos) = 200; %3000
	% KRONBORG
	pos = find(md.miscellaneous.dummy.catchment_id==107); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% KRUUSE_FJORD
	pos = find(md.miscellaneous.dummy.catchment_id==108); md.calving.stress_threshold_groundedice(pos) = 750; %650; %800;   % 500; %1000; %2000; %3100 %DONE
	% LAUBE_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==109); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% MAGGA_DAN_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==113); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NORDFJORD
	pos = find(md.miscellaneous.dummy.catchment_id==141); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% POLARIC-DECEPTION_O_N
	pos = find(md.miscellaneous.dummy.catchment_id==155); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 400; %200;
	% ROSENBORG
	pos = find(md.miscellaneous.dummy.catchment_id==165); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SORGENFRI
	pos = find(md.miscellaneous.dummy.catchment_id==202); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% SORTEBRAE
	pos = find(md.miscellaneous.dummy.catchment_id==203); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% STYRTE
	pos = find(md.miscellaneous.dummy.catchment_id==210); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SYDBR
	pos = find(md.miscellaneous.dummy.catchment_id==215); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_DECEPTION_N
	pos = find(md.miscellaneous.dummy.catchment_id==225); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_DECEPTION_O_CN_CS
	pos = find(md.miscellaneous.dummy.catchment_id==226); md.calving.stress_threshold_groundedice(pos) = 950; %1000; %800; %1200;    %1400; %1500; %1000; %2000; %3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 300;
	% UNNAMED_KANGER_E
	pos = find(md.miscellaneous.dummy.catchment_id==230); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_KANGER_W
	pos = find(md.miscellaneous.dummy.catchment_id==231); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 400;
	% UNNAMED_LAUBE_S
	pos = find(md.miscellaneous.dummy.catchment_id==232); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_POLARIC_C
	pos = find(md.miscellaneous.dummy.catchment_id==233); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_POLARIC_S
	pos = find(md.miscellaneous.dummy.catchment_id==234); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_SORGENFRI_W
	pos = find(md.miscellaneous.dummy.catchment_id==235); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_UUNARTIT_ISLANDS
	pos = find(md.miscellaneous.dummy.catchment_id==237); md.calving.stress_threshold_groundedice(pos) = 750; %700; %900; %800; %700;   %%900; %1000; %1200; %700; %2000; %3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 300;
	% VESTFJORD
	pos = find(md.miscellaneous.dummy.catchment_id==247); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	%}}}
	% === SE Glaciers (DONE) ===  %{{{

	% ANORITUUP_KANGERLUA
	pos = find(md.miscellaneous.dummy.catchment_id==7); md.calving.stress_threshold_groundedice(pos) = 1350; %1200; %1400; %1000; %1500;      %650; %700; %500; %1000; %2000 %3100;
	% APUSEERAJIK
	pos = find(md.miscellaneous.dummy.catchment_id==8); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% APUSEERSERPIA
	pos = find(md.miscellaneous.dummy.catchment_id==9); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% AP_BERNSTOFF_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==10); md.calving.stress_threshold_groundedice(pos) = 3000;   %2000; %1000; %3000; %DONE
	md.calving.stress_threshold_floatingice(pos) = 250; %150; %100; %200;
	% BRCKNER_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==17); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% BUSSEMAND
	pos = find(md.miscellaneous.dummy.catchment_id==20); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% FENRISGLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==38); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 3100;
	% FIMBULGETLSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==39); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% GLACIERDEFRANCE
	pos = find(md.miscellaneous.dummy.catchment_id==54); md.calving.stress_threshold_groundedice(pos) = 3100; %1000; %500; %1000; %2000; %3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 700; %500; %300;
	% GRAULV
	pos = find(md.miscellaneous.dummy.catchment_id==55); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% GYLDENLOVE
	pos = find(md.miscellaneous.dummy.catchment_id==56); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% GYLDENLOVE_S
	pos = find(md.miscellaneous.dummy.catchment_id==57); md.calving.stress_threshold_groundedice(pos) = 200; %400; %600;    %800; %DONE
	% GYLDENLOVE_SS
	pos = find(md.miscellaneous.dummy.catchment_id==58); md.calving.stress_threshold_groundedice(pos) = 800; %DONE, completely above sea level
	% HEIMDAL_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==67); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 300;
	% HEIM_GLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==68); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% HELHEIMGLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==69); md.calving.stress_threshold_groundedice(pos) = 2300; %2000;  %1800; %1500; %1150; %3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 700; %500;   %800; %500; %1150; %3000;
	% HERLUF_TROLLE-KANGERLULUK-DANELL
	pos = find(md.miscellaneous.dummy.catchment_id==71); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% IKERTIVAQ_M
	pos = find(md.miscellaneous.dummy.catchment_id==75); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% IKERTIVAQ_N
	pos = find(md.miscellaneous.dummy.catchment_id==76); md.calving.stress_threshold_groundedice(pos) = 1000; %2000; %3100; %DONE
	% IKERTIVAQ_NN
	pos = find(md.miscellaneous.dummy.catchment_id==77); md.calving.stress_threshold_groundedice(pos) = 1000; %2000; %3100; %DONE
	% IKERTIVAQ_S
	pos = find(md.miscellaneous.dummy.catchment_id==78); md.calving.stress_threshold_groundedice(pos) = 600; %DONE
	% KNUD-RASMUSSEN
	pos = find(md.miscellaneous.dummy.catchment_id==99); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% KOGE_BUGT_C
	pos = find(md.miscellaneous.dummy.catchment_id==100); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% KOGE_BUGT_N
	pos = find(md.miscellaneous.dummy.catchment_id==101); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% KOGE_BUGT_S
	pos = find(md.miscellaneous.dummy.catchment_id==102); md.calving.stress_threshold_groundedice(pos) = 2000; %3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 3000;
	% KOGE_BUGT_SS
	pos = find(md.miscellaneous.dummy.catchment_id==103); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% MAELKEVEJEN
	pos = find(md.miscellaneous.dummy.catchment_id==112); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% MIDGARDGLETSCHER
	pos = find(md.miscellaneous.dummy.catchment_id==118); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% MOGENS_HEINESEN_C
	pos = find(md.miscellaneous.dummy.catchment_id==119); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% MOGENS_HEINESEN_N
	pos = find(md.miscellaneous.dummy.catchment_id==120); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% MOGENS_HEINESEN_S
	pos = find(md.miscellaneous.dummy.catchment_id==121); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% MOGENS_HEINESEN_SS_SSS
	pos = find(md.miscellaneous.dummy.catchment_id==122); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NAPASORSUAQ_C_S
	pos = find(md.miscellaneous.dummy.catchment_id==128); md.calving.stress_threshold_groundedice(pos) = 500; %1000; %2000; %3100; %DONE
	% NAPASORSUAQ_N
	pos = find(md.miscellaneous.dummy.catchment_id==129); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	%md.calving.stress_threshold_floatingice(pos) = 3000;
	% NIGERTULUUP_KATTILERTARPIA
	pos = find(md.miscellaneous.dummy.catchment_id==133); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NONAME_IKERTIVAQ_N
	pos = find(md.miscellaneous.dummy.catchment_id==135); md.calving.stress_threshold_groundedice(pos) = 300; %600; %800; %1000; %DONE
	% NONAME_IKERTIVAQ_S
	pos = find(md.miscellaneous.dummy.catchment_id==136); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% PUISORTOQ_N
	pos = find(md.miscellaneous.dummy.catchment_id==156); md.calving.stress_threshold_groundedice(pos) = 2500; %2000;    %3100; %DONE
	% PUISORTOQ_S
	pos = find(md.miscellaneous.dummy.catchment_id==157); md.calving.stress_threshold_groundedice(pos) = 2000; %2400; %3100;      %2400; %2300; %1800; %1300; %1000;
	md.calving.stress_threshold_floatingice(pos) = 300;
	% RIMFAXE
	pos = find(md.miscellaneous.dummy.catchment_id==162); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% SE_NONAME1
	pos = find(md.miscellaneous.dummy.catchment_id==187); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SE_NONAME10
	pos = find(md.miscellaneous.dummy.catchment_id==188); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SE_NONAME2
	pos = find(md.miscellaneous.dummy.catchment_id==189); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% SE_NONAME4
	pos = find(md.miscellaneous.dummy.catchment_id==190); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SE_NONAME5
	pos = find(md.miscellaneous.dummy.catchment_id==191); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SE_NONAME6
	pos = find(md.miscellaneous.dummy.catchment_id==192); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SE_NONAME7
	pos = find(md.miscellaneous.dummy.catchment_id==193); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SE_NONAME8
	pos = find(md.miscellaneous.dummy.catchment_id==194); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SE_NONAME9
	pos = find(md.miscellaneous.dummy.catchment_id==195); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SKINFAXE
	pos = find(md.miscellaneous.dummy.catchment_id==200); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% SOUTHERN_TIP
	pos = find(md.miscellaneous.dummy.catchment_id==204); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% TINGMIARMIUT_FJORD
	pos = find(md.miscellaneous.dummy.catchment_id==216); md.calving.stress_threshold_groundedice(pos) = 3100; %2800; %2500; %3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 650; %700; %500;      %1000; %2000; %150; %3000;
	% UMIIVIK_FJORD
	pos = find(md.miscellaneous.dummy.catchment_id==221); md.calving.stress_threshold_groundedice(pos) = 650; %500; %700; %1000; %2000;    %%600; %800; %1000; %3100;
	% UNNAMED_ANORITUUP_KANGERLUA_S
	pos = find(md.miscellaneous.dummy.catchment_id==222); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_ANORITUUP_KANGERLUA_SS
	pos = find(md.miscellaneous.dummy.catchment_id==223); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_DANELL_FJORD
	pos = find(md.miscellaneous.dummy.catchment_id==224); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% UNNAMED_HERLUF_TROLLE_N
	pos = find(md.miscellaneous.dummy.catchment_id==227); md.calving.stress_threshold_groundedice(pos) = 3100; %2800; %2500; %2000; %3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 300;
	% UNNAMED_HERLUF_TROLLE_S
	pos = find(md.miscellaneous.dummy.catchment_id==228); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 300;
	% UNNAMED_KANGERLULUK
	pos = find(md.miscellaneous.dummy.catchment_id==229); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UNNAMED_SOUTH_DANELL_FJORD
	pos = find(md.miscellaneous.dummy.catchment_id==236); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 300;
	%}}}
	% === SW Glaciers (DONE) === %{{{

	% AKULLERSUUP-QAMANAARSUUP
	pos = find(md.miscellaneous.dummy.catchment_id==5); md.calving.stress_threshold_groundedice(pos) = 1000; %2000; %2500; %DONE
	% AVANNARLEQ-NIGERLIKASIK
	pos = find(md.miscellaneous.dummy.catchment_id==11); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% AVANNARLEQ_N
	pos = find(md.miscellaneous.dummy.catchment_id==12); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% EQALORUTSIT_KILLIIT_SERMIAT
	pos = find(md.miscellaneous.dummy.catchment_id==35); md.calving.stress_threshold_groundedice(pos) = 150; %DONE, completely above sea level
	% FREDERIKSHABS-NAKKAASORSUAQ
	pos = find(md.miscellaneous.dummy.catchment_id==41); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% ILORLIIT-SERMINNGUAQ
	pos = find(md.miscellaneous.dummy.catchment_id==80); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% INUPPAAT_QUUAT
	pos = find(md.miscellaneous.dummy.catchment_id==82); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% ISUNNGUATA-RUSSELL
	pos = find(md.miscellaneous.dummy.catchment_id==84); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% KANGIATA_NUNAATA_SERMIA
	pos = find(md.miscellaneous.dummy.catchment_id==92); md.calving.stress_threshold_groundedice(pos) = 2000; %DONE
	% KANGILINNGUATA_SERMIA
	pos = find(md.miscellaneous.dummy.catchment_id==94); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% KIATTUUT-QOOQQUP
	pos = find(md.miscellaneous.dummy.catchment_id==95); md.calving.stress_threshold_groundedice(pos) = 3100; %2900; %2700; %2000; %3100; %DONE
	md.calving.stress_threshold_floatingice(pos) = 700; %500;
	% NAAJAT_SERMIAT
	pos = find(md.miscellaneous.dummy.catchment_id==126); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% NARSAP_SERMIA
	pos = find(md.miscellaneous.dummy.catchment_id==130); md.calving.stress_threshold_groundedice(pos) = 1500; %1000; %420;
	% QAJUUTTAP_SERMIA
	pos = find(md.miscellaneous.dummy.catchment_id==158); md.calving.stress_threshold_groundedice(pos) = 3100;     % 2500; %3100;
	md.calving.stress_threshold_floatingice(pos) = 1000; %700; %500; %300;
	% QAJUUTTAP_SERMIA_N
	pos = find(md.miscellaneous.dummy.catchment_id==159); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% QALERALLIT_SERMIAT
	pos = find(md.miscellaneous.dummy.catchment_id==160); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SAQQAP-MAJORQAQ-SOUTHTERRUSSEL_SOUTHQUARUSSEL
	pos = find(md.miscellaneous.dummy.catchment_id==167); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SERMEQ-KANGAASARSUUP
	pos = find(md.miscellaneous.dummy.catchment_id==180); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SERMILIGAARSSUK_BRAE
	pos = find(md.miscellaneous.dummy.catchment_id==185); md.calving.stress_threshold_groundedice(pos) = 3100; %2500; %2000;      %800; %1000; %1500; %2000;
	md.calving.stress_threshold_floatingice(pos) = 1000; %700; %500; %300;
	% SERMILIK
	pos = find(md.miscellaneous.dummy.catchment_id==186); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SIORALIK-ARSUK-QIPISAQQU
	pos = find(md.miscellaneous.dummy.catchment_id==198); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% SW_NONAME1
	pos = find(md.miscellaneous.dummy.catchment_id==214); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, completely above sea level
	% UKAASORSUAQ
	pos = find(md.miscellaneous.dummy.catchment_id==219); md.calving.stress_threshold_groundedice(pos) = 3100; %DONE
	% USULLUUP_SERMIA
	pos = find(md.miscellaneous.dummy.catchment_id==244); md.calving.stress_threshold_groundedice(pos) = 1000; %DONE, grounding line is above sea level for now
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

	%% Add artificial ice on cliffs to prevent shelf from seperating from cliffs
	%disp('   Adjusting ice mask for Petermann cliffs around shelf')
	%cliff_pos_pet = find(ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y,...
	%   'exp/spclevelset/Petermann_spcLevelset.exp','node',0));
	%md.mask.ice_levelset(cliff_pos_pet) = -1;
	%disp('   Adjusting ice mask for 79N cliffs around shelf')
	%cliff_pos_79 = find(ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y,...
	%   'exp/spclevelset/79N_spcLevelset.exp','node',0));
	%md.mask.ice_levelset(cliff_pos_79) = -1;
	%%Tricky part here: we want to offset the mask by one element so that we don't end up with a cliff at the transition
	%md.mask.ice_levelset=reinitializelevelset(md,md.mask.ice_levelset);
	%
	%disp('   Setting spclevelset in cliff areas')
	%md.levelset.spclevelset        = NaN(md.mesh.numberofvertices,1);
	%mask  = md.mask.ice_levelset < 0;
	%spc_idx = union(cliff_pos_pet(mask(cliff_pos_pet)), cliff_pos_79(mask(cliff_pos_79)));
	%md.levelset.spclevelset(spc_idx) = -1;
	%md.levelset.migration_max      = 10e5; %max possible ice front retreat rate (10 km/yr)

	savemodel(org,md);
end%}}}
%for r = 1:length(climate_models)
%	climate_model = climate_models{r};
%	for s = 1:length(scenarios)
%		scenario = scenarios{s};
if perform(org,['Greenland_ISMIP7Run_' climate_model '_' upper(scenario) '_' int2str(projection_start_time) '-' int2str(projection_end_time) '_fix']),% {{{

	md=loadmodel(org,['Greenland_ISMIP7Prep_' climate_model '_' upper(scenario)]);
	md.groundingline.migration = 'SubelementMigration';

	disp('   Updating start/end times');
	md.timestepping.time_step=0.01;
	md.timestepping.start_time=projection_start_time;
	md.timestepping.final_time=projection_end_time;
	md.settings.output_frequency=10;

	disp('   Removing unused SMB and ocean forcings if relevant');
	t = md.smb.smbref(end, :);
	mask = (t >= projection_start_time) & (t < projection_end_time + 1);
	md.smb.smbref = md.smb.smbref(:, mask);
	t = md.smb.b_pos(end, :);
	mask = (t >= projection_start_time) & (t < projection_end_time + 1);
	md.smb.b_pos = md.smb.b_pos(:, mask);
	md.smb.b_neg = md.smb.b_neg(:, mask);

	t = md.frontalforcings.subglacial_discharge(end, :);
	mask = (t >= projection_start_time) & (t < projection_end_time + 1);
	md.frontalforcings.subglacial_discharge = md.frontalforcings.subglacial_discharge(:, mask);
	t = md.frontalforcings.thermalforcing(end, :);
	mask = (t >= projection_start_time) & (t < projection_end_time + 1);
	md.frontalforcings.thermalforcing = md.frontalforcings.thermalforcing(:, mask);

	md.verbose.solution = true;
	%		md.verbose=verbose('all');

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
%	end
%end
