## Submission Information

- **Contributor names, affiliations and emails:** Hugh Shields (nikolaus.h.shields.gr@dartmouth.edu), Mathieu Morlighem (mathieu.morlighem@dartmouth.edu)
- **Date of submission:** 09/30/2026
- **Ice Sheet Modeled / domain_id:** GrIS
- **Modeling group name / source_id:** Dartmouth
- **Ice Sheet Model Name / ism_id:** ISSM

---

## Initialization Methods

**1. Describe the initialization method used**, including assimilation of data, spin-up, or any other method.

>We use a transient calibration (see Badgeley et al. 2025) to initialize the friction coefficient used. The cost function used here minimizes an absolute velocity misfit, a logarithmic velocity misfit, and a regularization term to constrain the friction coefficient over a calibration period.

**2. What are the tuning targets, constraints and outputs of your initialization procedure?**

>The model is initialized to minimize the misfit between predicted and observed velocities on selected fast-flowing glaciers from InSAR (Joughin et al., 2021, NSIDC-0481) between 2008 and 2015 and monthly regional mosaics derived from SAR and optical imagery between 2015-2022 (Joughin et al., 2021, NSIDC-0731). Regularization in the cost function is determined using an L-curve.

**3. What procedures do you use to validate the initial conditions and/or the recent changes observed**, or tests that you apply to the final state of your initialization method?

>The transient calibration method inherently adjusts the friction coefficient to match observed velocities over the calibration period.

**4. What dataset is used for the SMB or ocean climatology?** Are calving fronts and ice margins held fixed or allowed to evolve during spin-up?

>During the calibration period, we use SMB from the Regional Atmospheric Climate Model v2.3p2 (RACMO, Noël et al. 2019). The calving fronts evolve following observed monthly ice front positions from Greene et al. (2024). Melting under shelves follows a linear melt parameterization with parameters from Holmes et al. (2025).

**5. Do you apply corrections, such as SMB corrections, during your initialization?** Please provide a description.

>No.

---

## Projections — GrIS Ice-Ocean Questions

**6. What is your typical spatial resolution near the calving fronts in your model?**

>500 m.

**7. What is your approach to representing calving?** What calving law or calving rate is used?

> We use the von Mises calving law (Morlighem et al. 2016).

**8. How is calving front migration treated?** Do you use any subgrid-scale schemes? How about ice shelf/grounding line migration?

> The calving front evolves via the level-set method (see Bondzio et al. 2016); likewise for the grounding line. Both use subgrid-scale parameterizations.

**9. What approach/parameterisation is used to estimate submarine melt of approximately vertical calving fronts?** What values do the parameters take? What about submarine melt of floating ice shelves? How are these submarine melt rates applied in your model (e.g. level-set or effective mass flux)? Do you use any subgrid-scale schemes in the implementation of these processes?

> We use the Rignot et al. (2016) submarine melt parameterization with inputs of the provided ocean thermal forcing, the provided subglacial discharge, and dynamically computed frontal areas. The melt rates are applied via the level-set, which is implemented using a subgrid-scale parameterization.

**10. Did you follow any calibration procedure for the representation of melt and calving?** If so, did you use observations, over what time period and over what spatial scale?

> We did not use a melt calibration procedure. The von Mises law was calibrated between 2007-2022 on a individual glacier catchment basis (catchments from  Mouginot et al. 2019) using the Greene et al. (2024) ice fronts.

**11. Did you use the provided forcing fields (Q/TF) or something you derived yourself?** For example, did you do any regridding of Q/TF to your model grid, and if so how? Do you retain multiple forcing values for a single glacier or somehow amalgamate them? For subglacial discharge, did you use the provided fields directly or calculate your own based on SMB/runoff?

> Q/TF were both regridded onto the native mesh and the melt was computed along the levelset using the local TF and the provided Q vlaue as well as the computed frontal area.

**12. Anything else relevant to understanding how you represented ocean forcing of the ice sheet in your model?**

> No.

---

## SMB Questions

**13. How did you implement SMB and what SMB corrections** (e.g., surface elevation feedback) do you apply in the experiments?

> SMB is implemented using the time-varying SMB anomaly from the ISMIP7 forcing. Surface elevation feedback are accounted for using the SMB gradients method from Helsen et al. (2025). See https://issmteam.github.io/ISSM-Documentation/using-issm/parameterization/smb.html for details.

**14. Anything else relevant to understanding how you represented atmospheric forcing of the ice sheet in your model?**

> No.

---

## GIA, Bedrock and Sea Level Questions

**15. Do you include bedrock adjustments?**

> No.

**16. What reference surface do you use for vertical heights** including bed elevation and surface elevation? This should be either a reference geoid or reference ellipsoid. (BedMachine reports elevations relative to EIGEN-6C4 geoid; Bedmap reports elevations relative to gl04c geoid.)

> EIGEN-6C4 geoid

**17. Explain the bedrock adjustment model and setup.** Include all chosen parameters. For ELRA, this includes relaxation time, flexural rigidity, and mantle density. For other models, include the assigned viscosity and density for each layer, and the lithospheric depth. Describe the spatial and temporal resolution and numerical methods employed. Note which of gravitational, rotational, and deformational effects are included, how Love numbers are calculated, and whether the model is incompressible or compressible.

> 

**18. Describe the spin-up or initialization procedure for the bedrock state**, detailing how it is prepared before the main experiment. Include the goal of the spin-up (e.g., reaching a steady state, representing adjustment to loading changes over a specified period, or another target configuration). If a relaxed bedrock state is defined in the setup, describe it.

>

**19. The protocol requests that far-field sea-level change is not included.** However, if it is, describe how it is applied.

>

---

## Other General Questions

**Describe how the ice-covered area was determined.** Were glaciers removed from the periphery? What is the target ice mask of the model, if any? How is the constraint enforced (e.g., during runtime per negative SMB, as a hard mask, or in postprocessing)?

>The domain includes peripheral glaciers. Ice covered area is determined from the ice levelset.

**Do you plan to participate in the Perturbed Parameter Ensemble (PPE) and/or Earth System Model experiments (ESM)?** Which experiments do you plan to participate in and how many experiments do you plan to perform?

>No.

**Summary paragraph** describing your model, for use in publications (see ISMIP6 group publications). Include references.

>The ice sheet is initialized to present conditions using transient calibration (Badgeley et al., 2025; Goldberg et al., 2026) to tune the basal friction coefficient. Calibration was performed at a regional scale and assimilated velocity data from InSAR between 2008 and 2015 (Joughin et al., 2021, NSIDC-0481) and monthly regional mosaics derived from SAR and optical imagery between 2015 and 2022 (Joughin et al., 2021, NSIDC-0731). The initial surface was interpolated from GIMP (Howat et al., 2014) and the bed topography from BedMachine v6 (Morlighem et al., 2025). The model uses the Shallow Shelf Approximation (SSA; MacAyeal, 1989) across the domain and a Budd-type friction law (Budd et al., 1979) with exponent one and effective pressure calculated assuming the bed is fully connected to the ocean. The depth-averaged viscosity is computed using Glen's flow law (Glen, 1955) with the standard exponent of three and ice hardness from the depth-averaged temperature field of the JPL-ISSM simulation for ISMIP6 (Goelzer et al., 2020). The mesh resolution varies between 1 and 10 km, with additional refinement to 500 m near the grounding line and in areas into which the grounding line retreats in projections. For the historical runs, calving front positions are prescribed from observations (Green et al., 2024). For the projection and observationally constrained experiments, calving is computed using the von Mises law (Morlighem et al., 2016), calibrated on a catchment-by-catchment basis to observed fronts between 2007 and 2022.

---

## Model Characteristic Table

Please complete the following table for your core experiments. This information will be used as a summary table in the publication. A separate spreadsheet will record information for non-core experiments (PPE and ESM). Please indicate your current plans, including the number of experiments where applicable.

| Characteristic | Main suite of experiments | Do you plan to change these as part of PPE? |
|---|---|---|
| Mesh discretization *(e.g. rectangular grid, Delaunay triangulation, ALE, Centroidal Voronoi tessellation)* | Delaunay triangulation (BAMG) | |
| Native Grid (horizontal and vertical) | adaptive (500-10000 m) | |
| Native Projection | EPSG 3031 | |
| Interpolation method to diagnostic grid | linear (P1) | |
| Time integration scheme; expected formal order of accuracy | semi-implicit, first-order | |
| Time Step |0.01 year | |
| Advection scheme, including numerical method and expected formal order of accuracy |streamline upwinding, first-order in space and time (Dias dos Santos et al., 2021)| |
| Ice Flow Mechanics, including numerical method |SSA (MacAyeal, 1989), FEM | |
| Ice Rheology | Glen with $n=3$ | |
| Basal Sliding | linear ($m=1$) Budd-type (Budd et al., 1979) | |
| Basal Hydrology | NA | |
| Advance and Retreat | level-set (Bondzio et al., 2016) | |
| Grounding Line: Determination, Parameterization | hydrostatic floatation criteria, sub-element parameterization, and levelset evolution (Seroussi and Morlighem, 2018) | |
| Calving | von Mises (Morlighem et al., 2016)| |
| Initial Surface Mass Balance | RACMO v2.3p2 (Noël et al., 2019)| |
| Do you include bedrock adjustment | No | |
| Year (or range of years) assigned to initial condition |2007 in historical runs | |
| Parameters for ice, ocean water and freshwater density ($ρ_i$, $ρ_o$, $ρ_w$), gravitational acceleration (g), etc. | $ρ_i=917$ kg/m^3, $ρ_o=1027$ kg/m^3, $ρ_w=1000$ kg/m^3, $g=9.81$ m/s^2 | |
| Variable in data request not included, and reason | `licalvf`, `lifmassbf`, `tendlicalvf`, and `tendlifmassbf` are omitted in the historical runs, as the calving front is prescribed from observation (calving and frontal melt are not separated)| |
| Number of days per year |365| |
| Other comments | | |

>See https://github.com/hughshields/ismip7greenland for run code.