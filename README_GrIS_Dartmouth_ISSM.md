## Submission Information

- **Contributor names, affiliations and emails:** Hugh Shields (nikolaus.h.shields.gr@dartmouth.edu), Mathieu Morlighem (mathieu.morlighem@dartmouth.edu)
- **Date of submission:** 07/14/2026
- **Ice Sheet Modeled / domain_id:** GrIS
- **Modeling group name / source_id:** Dartmouth
- **Ice Sheet Model Name / ism_id:** ISSM

---

## Initialization Methods

**1. Describe the initialization method used**, including assimilation of data, spin-up, or any other method.

>We use a transient calibration (see Badgeley et al. 2025) to initialize the friction coefficient used in the friction law (linear Budd). The cost function for this method includes an absolute velocity misfit, a logarithmic velocity misfit, and a regularization term.

**2. What are the tuning targets, constraints and outputs of your initialization procedure?**

>The model is initialized to minimize the misfit between predicted and observed velocities on selected fast-flowing glaciers from InSAR (Joughin et al., 2021, NSIDC-0481) between 2008 and 2015 and monthly regional mosaics derived from SAR and optical imagery between 2015-2022 (Joughin et al., 2021, NSIDC-0731).

**3. What procedures do you use to validate the initial conditions and/or the recent changes observed**, or tests that you apply to the final state of your initialization method?

>The transient calibration method inherently adjusts the friction coefficient to match recent observed velocities.

**4. What dataset is used for the SMB or ocean climatology?** Are calving fronts and ice margins held fixed or allowed to evolve during spin-up?

>During the calibration period, we use SMB from the Regional Atmospheric Climate Model v2.3p2 (RACMO, Noël et al. 2019). The calving fronts evolve following observed monthly ice front positions from Greene et al. (2024).

**5. Do you apply corrections, such as SMB corrections, during your initialization?** Please provide a description.

>No.

---

## Projections — GrIS Ice-Ocean Questions

**6. What is your typical spatial resolution near the calving fronts in your model?**

>500 m.

**7. What is your approach to representing calving?** What calving law or calving rate is used?

> We use the von Mises calving law (Morlighem et al. 2016).

**8. How is calving front migration treated?** Do you use any subgrid-scale schemes? How about ice shelf/grounding line migration?

> The calving front evolves via the level-set method (see Bondzio et al. 2016); likewise for the grounding line.

**9. What approach/parameterisation is used to estimate submarine melt of approximately vertical calving fronts?** What values do the parameters take? What about submarine melt of floating ice shelves? How are these submarine melt rates applied in your model (e.g. level-set or effective mass flux)? Do you use any subgrid-scale schemes in the implementation of these processes?

> Currently, we use the ISMIP6 miroc-esm-chem_rcp8.5 frontal forcing. This will be updated.

**10. Did you follow any calibration procedure for the representation of melt and calving?** If so, did you use observations, over what time period and over what spatial scale?

> We did not use a melt calibration procedure. The von Mises law was calibrated between 2007-2022 on a basin by basin basis (basins from  Mouginot et al. 2019) using the Greene et al. (2024) ice fronts.

**11. Did you use the provided forcing fields (Q/TF) or something you derived yourself?** For example, did you do any regridding of Q/TF to your model grid, and if so how? Do you retain multiple forcing values for a single glacier or somehow amalgamate them? For subglacial discharge, did you use the provided fields directly or calculate your own based on SMB/runoff?

> No.

**12. Anything else relevant to understanding how you represented ocean forcing of the ice sheet in your model?**

> No.

---

## SMB Questions

**13. How did you implement SMB and what SMB corrections** (e.g., surface elevation feedback) do you apply in the experiments?

> SMB is implemented as a fixed 1995–2014 RACMO2.3p2 reference climatology plus a time-varying anomaly from the ISMIP7 forcing. No surface elevation feedback or other elevation-dependent correction is applied to the SMB in this implementation.

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

>The domain did not include peripheral glaciers. Ice covered area is determined from the ice levelset.

**Do you plan to participate in the Perturbed Parameter Ensemble (PPE) and/or Earth System Model experiments (ESM)?** Which experiments do you plan to participate in and how many experiments do you plan to perform?

>No.

**Summary paragraph** describing your model, for use in publications (see ISMIP6 group publications). Include references.

>To be added in the future.

---

## Model Characteristic Table

Please complete the following table for your core experiments. This information will be used as a summary table in the publication. A separate spreadsheet will record information for non-core experiments (PPE and ESM). Please indicate your current plans, including the number of experiments where applicable.

| Characteristic | Main suite of experiments | Do you plan to change these as part of PPE? |
|---|---|---|
| Mesh discretization *(e.g. rectangular grid, Delaunay triangulation, ALE, Centroidal Voronoi tessellation)* | Delaunay triangulation (BAMG) | |
| Native Grid (horizontal and vertical) | adaptive (500-10000 m) | |
| Native Projection | EPSG 3031 | |
| Interpolation method to diagnostic grid | | |
| Time integration scheme; expected formal order of accuracy | | |
| Time Step |0.01 year | |
| Advection scheme, including numerical method and expected formal order of accuracy | | |
| Ice Flow Mechanics, including numerical method | | |
| Ice Rheology | glen n=3 | |
| Basal Sliding | linear budd | |
| Basal Hydrology | NA | |
| Advance and Retreat | | |
| Grounding Line: Determination, Parameterization | | |
| Calving | von Mises | |
| Initial Surface Mass Balance | | |
| Do you include bedrock adjustment | No | |
| Year (or range of years) assigned to initial condition | | |
| Parameters for ice, ocean water and freshwater density (ρ_i, ρ_o, ρ_w), gravitational acceleration (g), etc. | | |
| Variable in data request not included, and reason | | |
| Number of days per year | | |
| Other comments | | |

**References linked to table:**
>
