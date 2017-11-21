---
layout: default
---

# Overview
Simulations are organized into **sets**.  Sets have names like
`emulator_1100box_planck`, and simulations have numbers appended to
them, like `emulator_1100box_planck_00`.  A dash after the number (like 
`emulator_1100box_planck_00-1`) means this is a different
phase of the same cosmology; that is, an independent realization of the power spectrum.

# Simulation Sets

<!---
- ## AbacusCosmos_1100box
  - 40 boxes of size 1100 Mpc/h and particle mass \\(\sim 4\times 10^{10}\ M_\odot\\)
  - Different cosmologies centered on Planck 2015 with identical phases in the initial conditions
  - <http://nbody.rc.fas.harvard.edu/public/AbacusCosmos/AbacusCosmos_1100box_products/>

- ## AbacusCosmos_720box
  - Same as `AbacusCosmos_1100box` above, but at higher mass resolution (\\(\sim 1\times 10^{10}\ M_\odot\\)) and smaller box size (720 Mpc/h)
  - Note these are not "zoom-in" simulations of the larger boxes, but completely independent realizations of the power spectrum
  - <http://nbody.rc.fas.harvard.edu/public/AbacusCosmos/AbacusCosmos_720box_products/>

- ## emulator_planck_1100box
  - 21 total boxes divided into two sets:
    - 17 with identical fiducial cosmologies and different phases (useful for stacking volume)
    - 4 boxes with single-parameter deviations from the Planck cosmology (useful for measuring derivatives)
  - Box size 1100 Mpc/h and particle mass \\(\sim 4\times 10^{10}\ M_\odot\\)
  - <http://nbody.rc.fas.harvard.edu/public/AbacusCosmos/emulator_1100box_planck_products/>

- ## emulator_planck_720box
  - Same as `emulator_planck_1100box` above, but at higher mass resolution (\\(\sim 1\times 10^{10}\ M_\odot\\)) and smaller box size (720 Mpc/h)
  - Note these are not "zoom-in" simulations of the larger boxes, but completely independent realizations of the power spectrum
  - <http://nbody.rc.fas.harvard.edu/public/AbacusCosmos/emulator_720box_planck_products/>
-->
  
| **Simulation Set Name** | **# of sims** | **Box Size [\\(\mathrm{Mpc}/h\\)]** | **\\(N_\mathrm{part}\\)** | **Particle Mass [\\(M_\odot/h\\)]** | **Cosmologies** | **Initial Phases** | **Output Redshifts** | **Notes** | **Browse** |
|:-------------|--------------:|------------------------------------:|--------------------------:|:------------------------------------|:----------------|:-----------|:--------------|:----------|:-----------|
| **AbacusCosmos_1100box** | 40 | 1100 | \\(1440^3\\) | \\(\sim 4\times 10^{10}\\) | Latin Hypercube centered on Planck 2015 | Matched | 1.5, 1.0, 0.7, 0.5, 0.3 | | [Browse](http://nbody.rc.fas.harvard.edu/public/AbacusCosmos/AbacusCosmos_1100box_products/) |
| **AbacusCosmos_720box** | 40 | 720 | \\(1440^3\\) | \\(\sim 1\times 10^{10}\\) | Latin Hypercube centered on Planck 2015 | Matched | 1.5, 1.0, 0.7, 0.5, 0.3, 0.1 | Not zoom-in boxes of the above | [Browse](http://nbody.rc.fas.harvard.edu/public/AbacusCosmos/AbacusCosmos_720box_products/) |
| <span style="white-space: nowrap;"> **emulator_1100box_planck_00-{0..15}** </span> | 16 | 1100 | \\(1440^3\\) | \\(\sim 4\times 10^{10}\\) | Fiducial | Independent | 0.7, 0.57, 0.5, 0.3 | Volume-building boxes | [Browse](http://nbody.rc.fas.harvard.edu/public/AbacusCosmos/emulator_1100box_planck_products/) |
| **emulator_1100box_planck_{00..04}** | 5 | 1100 | \\(1440^3\\) | \\(\sim 4\times 10^{10}\\) | Fiducial + perturbations | Matched | 0.7, 0.57, 0.5, 0.3 | Derivative-measuring boxes | [Browse](http://nbody.rc.fas.harvard.edu/public/AbacusCosmos/emulator_1100box_planck_products/) |
| **emulator_720box_planck_00-{0..15}** | 16 | 720 | \\(1440^3\\) | \\(\sim 1\times 10^{10}\\) | Fiducial | Independent | 0.7, 0.57, 0.5, 0.3, 0.1 | Volume-building boxes | [Browse](http://nbody.rc.fas.harvard.edu/public/AbacusCosmos/emulator_720box_planck_products/) |
| **emulator_720box_planck_{00..04}** | 5 | 720 | \\(1440^3\\) | \\(\sim 1\times 10^{10}\\) | Fiducial + perturbations | Matched | 0.7, 0.57, 0.5, 0.3, 0.1 | Derivative-measuring boxes | [Browse](http://nbody.rc.fas.harvard.edu/public/AbacusCosmos/emulator_720box_planck_products/) |


# Cosmologies
## Fiducial Cosmology
The fiducial cosmology is taken from the [Planck 2015 cosmological parameters paper](https://arxiv.org/abs/1502.01589).
This is the cosmology of the `emulator_planck_00` sims.  The `AbacusCosmos` sims use 40 cosmologies scattered around this cosmology.

The cosmological parameters are listed here for convenience, but in analysis applications the cosmological parameters should always be read from
`abacus.par` or `camb_params.ini`.

| **Parameter** | **Value** |
|:---------|:--------|
| `ombh^2` | 0.02222 |
| `omcdmh^2` | 0.1199 |
| `omh^2` | 0.14212 |
| `w0` | -1.0 |
| `ns` | 0.9652 |
| `sigma_8` | 0.830 |
| `H0` | 67.26 |
| `N_eff` | 3.04 |
| `massless_neutrinos` | 3.04 |
| `omnuh2` | 0.0 |

## Fiducial + perturbations
The `emulator_planck_{00..04}` boxes provide a central cosmology and four perturbations: \\(\sigma_8 \pm 5\%\\) and \\(H_0 \pm 5\%\\).  Thus, these 5 boxes form a "plus sign" of cosmologies.

The perturbations are listed here.  Note that since the cosmologies are chosen in the space of physical densities, changes to \\(H_0\\) result in changes to \\(\Omega_m, \Omega_\Lambda\\).

| **Sim**      | \\(H_0\\)   | \\(\Omega_\mathrm{DE}\\)   | \\(\Omega_M\\)   | \\(\sigma_8\\)   |
|:-------------|------------:|---------------------------:|-----------------:|-----------------:|
| 00 (fiducial)    | 67.3    | 0.686      | 0.314     | 0.83      |
| 01    | 67.3    | 0.686          | 0.314         | **0.78**      |
| 02    | 67.3    | 0.686          | 0.314         | **0.88**      |
| 03    | **64.3** | **0.656**      | **0.344**     | 0.83      |
| 04    | **70.3** | **0.712**      | **0.288**     | 0.83         |

## AbacusCosmos Cosmologies
The `AbacusCosmos_1100box` and `AbacusCosmos_720box` sims use a set of 40 cosmologies chosen with a Latin hypercube algorithm centered on the Planck 2015 cosmology.
The cosmology of a given simulation can be read from the `info/abacus.par` file.
Below we show a corner plot representation of this 5-dimensional parameter space.  The blue square marks the fiducial cosmology.

<center>
<img src="{{ site.baseurl }}{% link cosmology_corner_plot.png %}" alt="cosmology corner plot" style="width: 50%;"/>
</center>