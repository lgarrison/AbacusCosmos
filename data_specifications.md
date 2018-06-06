---
layout: default
---

# Overview
The file formats for all the data products are documented below.
[Python code]({{ site.baseurl }}{% link code.md %}) is also provided to load the friends-of-friends and Rockstar halos,
but using it is completely optional.  If nothing else, it provides a concrete
example of how to read the files, in case you have any questions after reading the
documentation below.

Jump to:
- [Friends-of-friends Halos](#friends-of-friends-halos)
- [Rockstar Halos](#rockstar-halos)
- [Particle Subsamples](#particle-subsamples)
- [Power Spectra](#power-spectra)
- [Info directory](#info-directory)
- [Abacus Parameter and header files](#abacus-parameter-and-header-files)
- [Important parameters](#abacus-parameter-and-header-files)

# Friends-of-friends halos
## Directory structure
The friends of friends (FoF) halos are located in the `SimName_FoF_halos` subdirectory of each simulation.  An example directory structure follows:

<div markdown="1" class="tree">
- `AbacusCosmos_1100box_products`
  - `AbacusCosmos_1100box_00_products`
    - `AbacusCosmos_1100box_00_FoF_halos`
      - `info` (see [Info directory](#info-directory))
      - `z0.300`
        - `halos.tar.gz`: The main halo catalogs
        - `fof.cfg`: FoF configuration file
        - `header`: the Abacus parameter file describing this redshift slice
        - `particles.tar.gz`: a 10% subsample of particles inside halos
        - `particle_ids.tar.gz`: the corresponding particle IDs
        - `field_particles.tar.gz`: a 10% subsample of particle outside of halos
        - `field_ids.tar.gz`: the corresponding particle IDs
</div>
        
## filename: <code class="fn">halos_N</code>

These files form the main halo catalog (i.e. all halo information except for particle subsamples).  The
number of these files is arbitrary.
The file format is two integers, followed by the list of halo records:

| **Name** | **Type** | **Units** | **Description** |
|:---- |:---- |:----- |:----------- |
| `num_groups` | uint64_t | | Number of halo records in this file |
| `n_largest_subhalos` | uint64_t | | Number of subhalos records per halo |
| halos | `struct halo` | | The halos |

The halo records are unpacked C structs with the following fields:

| **Name** | **Type** | **Units** | **Description** |
|:---- |:---- |:----- |:----------- |
| `id` | int64_t | | Halo ID |
| `npstart` | uint64_t | | Particle subsample starting index in corresponding particles_N file |
| `npout` | uint32_t | | Number of subsampled particles |
| `N` | uint32_t | | Number of halo particles |
| `subhalo_N` | uint32_t[N_LARGEST_SUBHALOS] | | Number of particles in the few most massive subhalos |
| `x` | float[3] | Mpc/h | Center of mass of the halo |
| `v` | float[3] | km/s | Velocity of the center of mass |
| `sigmav` | float[3] | km/s | Velocity dispersion |
| `r25`, `r50`, `r75`, `r90` | floats | Mpc/h | Radial percentiles of the particle distribution |
| `vcirc_max` | float | km/s | The maximum circular velocity of the halo |
| `rvcirc_max` | float | Mpc/h | The radius at which `vcirc_max` occurs |
|:---- |:---- |:----- |:----------- |
|`subhalo_x` | float[3] | Mpc/h | Center of mass of the most massive subhalo |
| `subhalo_v` | float[3] | km/s | Velocity of the center of mass |
|:---- |:---- |:----- |:----------- |
| `subhalo_sigmav` | float[3] | km/s | Parent halo velocity dispersion, centered on the most massive subhalo |
| `subhalo_r25`, `subhalo_r50`, `subhalo_r75`, `subhalo_r90` | floats | Mpc/h | Parent halo radial percentiles, centered on the most massive subhalo |
| `subhalo_vcirc_max` | float | km/s | The maximum circular velocity of the parent halo, centered on the most massive subhalo |
| `subhalo_rvcirc_max` | float | Mpc/h | The radius at which `subhalo_vcirc_max` occurs |

## filename: <code class="fn">particles_N</code>
These files hold the particle subsamples for the halos in file `halos_N`.  The format is a list of particles stored as binary C structs.  Each struct has
the following format:

| **Name** | **Type** | **Units** | **Description** |
|:---- |:---- |:----- |:----------- |
| `pos` | float[3] | Mpc/h | The particle position |
| `vel` | float[3] | km/s | The particle velocity |

Particles can be associated with halos by the `npstart` and `npout` fields of the corresponding `halos_N` file.

The corresponding particle IDs can be found in `particle_ids_N`.  They appear in the same order as this `particles_N` file.

## filename: <code class="fn">particle_ids_N</code>
These files hold the particle IDs corresponding to the particle subsamples for the halos in file `halos_N`.  The format is a binary list of particle IDs:

| **Name** | **Type** | **Units** | **Description** |
|:---- |:---- |:----- |:----------- |
| `pid` | uint64_t |  | The particle ID |

Particle IDs appear in the same order as particle positions in the `particles_N` file.
Particles can be associated with halos by the `npstart` and `npout` fields of the corresponding `halos_N` file.

## filename: <code class="fn">field_particles_N</code>
These files hold the field particle subsamples.  Field particles are particles that do not fall in any FoF halo.  The format is a list of particles stored as binary C structs.  Each struct has
the following format:

| **Name** | **Type** | **Units** | **Description** |
|:---- |:---- |:----- |:----------- |
| `pos` | float[3] | Mpc/h | The particle position |
| `vel` | float[3] | km/s | The particle velocity |

The corresponding particle IDs can be found in `field_ids_N`.  They appear in the same order as this `field_particles_N` file.

## filename: <code class="fn">field_ids_N</code>
These files hold the particle IDs corresponding to the field particle subsamples.  Field particles are particles that do not fall in any FoF halo.  The format is a binary list of particle IDs:

| **Name** | **Type** | **Units** | **Description** |
|:---- |:---- |:----- |:----------- |
| `pid` | uint64_t |  | The particle ID |

Field particle IDs appear in the same order as field particle positions in the `field_particles_N` file.

## filename: <code class="fn">fof.cfg</code>
The parameters that FoF was invoked with.  The most important parameters are documented here:

| **Name** | **Type** | **Units** | **Description** |
|:---- |:---- |:----- |:----------- |
| `linklen` | float | interparticle spacing | Level 1 FoF linking length for finding halos |
| `linklen_L2` | float | interparticle spacing | Level 2 FoF linking length for finding subhalos |
| `n_largest_subhalos` | int | | How many subhalo masses to record in `subhalo_N` |
| `min_members` | int | Particle count | Minimum halo size to output |
| `n_block` | int |  | The number of blocks to divide the domain into.  Should equal the number of `halo_N` output files. |


## Units
All distances are comoving in the domain [-`BoxSize`/2, `BoxSize`/2), and all velocities are proper.  See the Rockstar unit notes for information about comparing FoF positions to Rockstar positions.


# Rockstar halos
The [Rockstar](https://bitbucket.org/gfcstanford/rockstar) halos are located in the `SimName_Rockstar_halos` subdirectory of each simulation.  An example directory structure follows:
<div markdown="1" class="tree">
- `AbacusCosmos_1100box_products`
  - `AbacusCosmos_1100box_00_products`
    - `AbacusCosmos_1100box_00_rockstar_halos`
      - `info` (see [Info directory](#info-directory))
      - `z0.300`
        - `halos.tar.gz`: the main halo catalogs
        - `rockstar.cfg`: the Rockstar configuration file
        - `header`: the Abacus parameter file describing this redshift slice
        - `particles.tar.gz`: a 10% subsample of particles inside halos
</div>

## filename: <code class="fn">halos_M.N.h5</code>
The Rockstar halos are stored in HDF5 files.  There is one dataset called `halos`.  The fields are identical to default Rockstar,
with the following modifications:
- We split `pos[6]` into `pos[3]` and `vel[3]`
- We add the `parent_id` field. For subhalos, this is the id of the parent halo. For parent halos, this is -1.
- We add the `subsamp_start` and `subsamp_len` fields, corresponding to the starting index and count of the halo particle subsample.
- We add `m_SO` and `alt_m_SO[4]` for the spherical overdensity masses.
- We add `N`, `alt_N[4]`, `N_SO`, and `alt_N_SO[4]`: the particle counts corresponding to the mass values.

The `header` file is also stored as an HDF5 dataset attribute, but this is purely for convenience; it contains the same information as the actual header file.

## filename: <code class="fn">particles_M.N.h5</code>
The particle subsamples corresponding to halos in `halos_M.N.h5`.  There is one dataset called `subsamples` (in `AbacusCosmos`) or `particles` (in `emulator_planck`).
The dataset contains the positions, velocities, and, in `AbacusCosmos`, the PIDs.

The `header` file is also stored as an HDF5 dataset attribute, but this is purely for convenience; it contains the same information as the actual header file.

## filename: <code class="fn">rockstar.cfg</code>
The Rockstar configuration file.  See the [Rockstar documentation](https://bitbucket.org/gfcstanford/rockstar) for definitions.

## Units
The mass definitions corresponding to the `m` and `alt_m` fields can be read from the `MASS_DEFINITION` fields in `rockstar.cfg`.  In almost all cases, we have left the mass definitions at the default, so `m` corresponds to the virial mass.

For convenience, we copy the Rockstar notes about units from the ASCII header here.  You may need to consult the [Rockstar documentation or source code](https://bitbucket.org/gfcstanford/rockstar)
for the exact definition of a quantity, however.
- Masses in Msun / h
- Positions in Mpc / h (comoving)
- Velocities in km / s (physical, peculiar)
- Halo Distances, Lengths, and Radii in kpc / h (comoving)
- Angular Momenta in (Msun/h) * (Mpc/h)*km/s (physical)
- Spins are dimensionless
- Total energy in (Msun/h) * (km/s)^2 (physical)
- idx, i_so, and i_ph are internal debugging quantities
- Np is an internal debugging quantity.

The Rockstar halo positions are in the domain [0, `BoxSize`).  Rockstar subsample particle positions may be negative if the halo is near the edge of the box.  Applying a periodic wrap (pos % BoxSize) will bring the particle positions back to the same domain as the halos.

A related issue is that Rockstar takes the Abacus particle positions modulo the box size, since the original particles span [-`BoxSize`/2, `BoxSize`/2) but Rockstar wants to work in [0, `BoxSize`).  FoF works in native Abacus units and does not do this conversion.  If you want to compare halos between Rockstar and FoF, you will need to apply the same unit conversion.  In the forward direction, one way to do this is as follows:

\\( \mathbf{x}\_\mathrm{Rockstar} = \mathbf{x}\_\mathrm{FoF} \bmod L \\)

In the backwards direction:

\\( \mathbf{x}\_\mathrm{FoF} = \mathbf{x}\_\mathrm{Rockstar} - L\, \mathrm{round}(\mathbf{x}\_\mathrm{Rockstar}/L)\\)

## SO halos
Rockstar can also produce spherical overdensity masses for the same set of halo centers that the main code uses. That is, no new halo finding is done, but new masses are computed that include all unbound particles.  We compute these masses and store them in the `m_SO` and `alt_m_SO` fields.

Caution should be used when interpreting SO halo masses of less than 100 particles.  About 2% of halos in this regime are so close to another halo that their spherical overdensities do not fall below the threshold before encompassing the larger halo.  To prevent small halos from artificially inflating their masses this way, Rockstar limits the SO search radius to 10% beyond the nominal halo radius (this is the `BCG2_R=1.1e-3` factor).  In these cases, the reported mass is not a "true" SO mass in the sense that the particle counting was truncated before the SO threshold was reached.  If nothing else, this is yet another reason not to trust small halos.


# Particle subsamples
Both the friends of friends and Rockstar catalogs contain 10% halo particle subsamples for use, e.g. in populating halos with HOD galaxies.
The subsample particles are selected based on their particle ID, so the set of "subsample-able" particles is the same across time slices, and
even between FoF and Rockstar.  The fact that the same particles are subsampled in every time slice (if they are in a halo) could allow one to
construct crude merger trees, for instance.

The FoF catalogs also include a 10% subsample of particles that fall outside of halos, called "field" particles.  Thus, the union of the halo particles
and the field particles is a uniform 10% subsample of the whole matter density field.

We do not include a field subsample from Rockstar, as it would increase the data volume significantly and require internal Rockstar changes.

# Power spectra
Power spectra from each redshift slice are located in the `SimName_power` subdirectory of each simulation.  An example directory structure follows:
<div markdown="1" class="tree">
- `AbacusCosmos_1100box_products`
  - `AbacusCosmos_1100box_00_products`
    - `AbacusCosmos_1100box_00_power`
      - `info` (see [Info directory](#info-directory))
      - `z0.300`
        - `power_nfft2048.csv`: The measured power spectrum
        - `header`: the Abacus parameter file describing this redshift slice
</div>
        
The power spectrum file `power_nfft2048.csv` is a comma separated values file with three columns:
wavenumber \\(k\\) [h/Mpc], power \\(P(k)\\) [(Mpc/h)^3], and number of modes \\(N_\mathrm{modes}\\).  The suffix `_nfft2048`
in the filename indicates that the power spectrum computation used a \\(2048^3\\) mesh.

We compute the matter power spectra by gridding the particles onto a mesh (\\(2048^3\\) or finer) with triangle-shaped cloud mass assignment.
We then Fourier transform the density field, convert the result to a power spectrum, de-convolve the TSC-aliased window function from
Jeong (2010), and bin in spherical annuli.  The number of mesh cells that fall into each annulus is recorded as N_modes in the third
column of the csv file.


# <code class="fn">info</code> directory
Every data product directory contains an `info` directory (alongside the redshift directories) with various input and output files.  An example `info` directory structure is shown here:
<div markdown="1" class="tree">
- `info`
    - `abacus.par`: the main Abacus parameter file
    - `abacus_params.par2`: level 2 parameter file that gets processed into the level 1 `abacus.par` file.  Useful for just seeing the parameters that vary among sims.
    - `camb_derived.out`: quantities that CAMB outputs as it runs
    - `camb_matterpower.dat`: the main CAMB matter power spectrum; this is the input power spectrum to the zeldovich-PLT code.
    - `camb_params.ini`: the CAMB input file
    - `camb_transfer_out.dat`: various transfer functions that CAMB produces
</div>

# Abacus parameter and header files
The parameter file `info/abacus.par` is the input configuration file to Abacus for the simulations.
It contains settings that do not vary with time, like the box size and `Omega_M`.

The header file `zZ.ZZZ/header` is an Abacus output file that describes the current redshift slice.
It contains the entire `abacus.par` file, followed by time-varying quantities, like the scale factor
and `OmegaNow_m`, the \\(\Omega_M\\) that an observer at this epoch would measure.

You will usually want to parse the `header` file for analysis tasks.  If you are working in Python,
[InputFile.py](https://github.com/lgarrison/AbacusCosmos/blob/master/AbacusCosmos/InputFile.py) can do this for you.

# Important parameters
We collect the most important simulation parameters here for convenience.  These are values can all be found in `info/abacus.par` or `zZ.ZZZ/header`.
If something isn't covered here, the [Abacus user guide](http://nbody.rc.fas.harvard.edu/public/AbacusCosmos/abacus_user_guide_public.pdf) contains detailed descriptions of all parameters.

| **Name** | **Type** | **Units** | **Description** |
|:---- |:---- |:----- |:----------- |
| `BoxSize` | float | Mpc/h | Comoving side length of the simulation cube |
| `H0` | float | km/s/Mpc | The Hubble constant at z=0 |
| `InitialRedshift` | float | | The starting redshift of the simulation |
| `NP` | int64 | | The number of particles in the simulation.  Equal to ppd^3. |
| `N_eff` | float | | The effective number of neutrino species.  Only used by CAMB, not Abacus. |
| `Omega_DE` | float | | Dark energy density parameter at z=0 |
| `Omega_M` | float | | Matter density parameter at z=0 |
| `SofteningLength` | float | Mpc/h in AbacusCosmos sims; Unit-box in emulator_planck sims | Plummer-equivalent softening length. Unit-box units can be converted to Mpc/h by multipling by BoxSize. |
| `ns` | float | | Scalar spectral index.  Only used by CAMB, not Abacus. |
| `ombh2` | float | | Physical baryon density parameter. Only used by CAMB, not Abacus. |
| `omch2` | float | | Physical CDM density parameter. Only used by CAMB, not Abacus. |
| `simga_8` | float | | Amplitude of density fluctuations at 8 Mpc/h.  Technically not read by CAMB or Abacus.  This is converted to `ZD_Pk_sigma` before the zeldovich-PLT IC code is invoked. |
| `w0` | float | | Dark energy equation of state parameter at z=0 |
| `wa` | float | | Dark energy equation of state evolution parameter: \\(w(a) = w_0 + (1-a)w_a\\) |
| `ppd` | float | | "Particles per dimension"; NP = ppd^3 |
| `ScaleFactor` | float | | The current scale factor a = 1 / (1 + z) |
| `BoxSizeMpc` and `BoxSizeHMpc` | floats | Mpc and Mpc/h | The comoving size of the box in Mpc and Mpc/h |
| `HubbleTimeGyr` and `HubbleTimeHGyr` | floats | Gyr and Gyr/h | The value of 1/`H0` in Gyr and Gyr/h; note that all times below are reported in units of `H0 = 1`. |
| `ParticleMassMsun` and `ParticleMassHMsun` | floats | M_sun and M_sun/h | The particle mass in M_sun and M_sun/h units |
| `VelZSpace_to_kms` | float | | Conversion factor of output particle velocities to km/s |
| `Redshift` | float | | The redshift z at this epoch; 1 + z = 1/a |
| `Time` | float | | Proper time in 1/`H0` units |
| `Growth` | float | | The linear growth function, normalized to \\(a\\) at early times (EdS) |
| `Growth_on_a` | float | | The linear growth function divided by \\(a\\), which is unity in Einstein-de Sitter. |
| `f_growth` | float | | The growth rate \\(d \ln D / d \ln a \\) |
| `w` | float | | The dark energy equation of state parameter at this epoch |
| `HubbleNow` | float | | The Hubble parameter at the current epoch, in H0 units, i.e., H(z)/H0 |
| `Htime` | float | | The product of the Hubble parameter and the current time, H(z)t(z). This 2/3 in Einstein-de Sitter. |
| `OmegaNow_m` | float | | The \\(\Omega_M\\) that an observer at this epoch would measure |
| `OmegaNow_DE` | float | | The \\(\Omega_{DE}\\) that an observer at this epoch would measure |
| `SofteningType` | string | | The softening technique employed.  See the Abacus Cosmos paper for softening details. |
