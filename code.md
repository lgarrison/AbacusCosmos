---
layout: default
---

# Overview
We provide example Python code to load the FoF and Rockstar halo catalogs and particle subsamples.
The main file is [Halos.py](https://github.com/lgarrison/AbacusCosmos/blob/master/AbacusCosmos/Halos.py);
we also provide [Halotools.py](https://github.com/lgarrison/AbacusCosmos/blob/master/AbacusCosmos/Halotools.py),
which is a thin wrapper around `Halos.py`
that puts the catalogs in [halotools](http://halotools.readthedocs.io/) format, where HODs can easily be applied.

All the code is available on the [AbacusCosmos GitHub](https://github.com/lgarrison/AbacusCosmos).
Please file an issue there if you find a problem with the code.

Using these interfaces is not necessary; all of the data formats are documented in
[Data Specifications]({{ site.baseurl }}{% link data_specifications.md %}).  But the example
code might help you get started using the catalogs faster or help with implementing your own
code to parse the catalogs.

The generally recommended usage is to clone the repository and add it to your `PYTHONPATH` environment variable.
This will allow you to import AbacusCosmos from anywhere.  Typically, this will involve
placing a command like the following in your `~/.bashrc` file:
```bash
export PYTHONPATH="/path/to/AbacusCosmos/:$PYTHONPATH"
```

[Yuan, Eisenstein, & Garrison](https://arxiv.org/abs/1802.10115) has also developed a decorated HOD package called GRAND-HOD that incorporates multiple HOD generalizations, such as assembly bias and velocity bias.  The [GRAND-HOD](https://github.com/SandyYuan/GRAND-HOD) repository contains details on its usage and code examples.

# Examples
The following examples are adapted from [this Jupyter notebook](https://github.com/lgarrison/AbacusCosmos/blob/master/AbacusCosmos_Python_Interface_Examples.ipynb).

<div markdown="1" class="jupyternb">
## [Halos.py](https://github.com/lgarrison/AbacusCosmos/blob/master/AbacusCosmos/Halos.py)


```python
from AbacusCosmos import Halos

import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np
```

### Load a single halo Rockstar catalog with particle subsamples:


```python
# The path to the catalog will be different on your system
cat = Halos.make_catalog_from_dir(dirname='/mnt/store2/bigsim_products/AbacusCosmos_1100box_products/AbacusCosmos_1100box_00_products/AbacusCosmos_1100box_00_rockstar_halos/z0.300',
                                  load_subsamples=True, load_pids=False)
halos = cat.halos
halo_1000 = halos[1000]
subsamples = cat.subsamples
```


```python
for field in sorted(halo_1000.dtype.fields):
    print field, ':', halo_1000[field]
```

    A : [ 53.04890823   5.39909792  -9.3987608 ]
    A2 : [ 39.85029221   6.90763903  -7.76677895]
    J : [ -1.96935156e+13   2.79465133e+13  -1.44266725e+13]
    N : 155
    N_SO : 983
    Voff : 29.6215
    Xoff : 22.5928
    alt_N : [135 121  91  29]
    alt_N_SO : [983 983 983  48]
    alt_m : [  5.03707704e+12   4.51471370e+12   3.39536301e+12   1.08203881e+12]
    alt_m_SO : [  3.66773849e+13   3.66773849e+13   3.66773849e+13   1.79096073e+12]
    b_to_a : 0.462127
    b_to_a2 : 0.381544
    bulkvel : [ -663.75305176    93.82868958 -1657.2800293 ]
    bullock_spin : 0.0713309
    c_to_a : 0.330447
    c_to_a2 : 0.274972
    child_r : 0.439964
    corevel : [ -682.87280273   108.45961761 -1700.54968262]
    desc : 0
    energy : -9.7564e+16
    flags : 25
    halfmass_radius : 154.012
    id : 298710
    kin_to_pot : 0.726075
    klypin_rs : 49.5365
    m : 5.78331e+12
    m_SO : 3.66774e+13
    m_pe_b : 5.29508e+12
    m_pe_d : 3.5073e+12
    mgrav : 4.7759e+12
    min_bulkvel_err : 984.373
    min_pos_err : 0.000250784
    min_vel_err : 1050.16
    n_core : 98
    num_child_particles : 206
    num_p : 206
    p_start : 73900406
    parent_id : 298727
    pos : [ 426.9362793    15.22821426  237.38449097]
    r : 375.432
    rs : 55.4336
    rvmax : 237.906
    spin : 0.0540608
    subsamp_len : 18
    subsamp_start : 40967
    vel : [ -682.87280273   108.45961761 -1700.54968262]
    vmax : 303.445
    vmax_r : 0.431498
    vrms : 326.531


The subsample particles should be near the halo location of `(426.9362793, 15.22821426, 237.38449097)`, and indeed they are:


```python
particles_for_halo_1000 = subsamples[halo_1000['subsamp_start']:halo_1000['subsamp_start']+halo_1000['subsamp_len']]
print particles_for_halo_1000['pos']
```

    [[ 426.91430664   15.24571419  237.90284729]
     [ 427.37713623   15.36857128  237.59571838]
     [ 427.23428345   15.09571362  237.58428955]
     [ 426.9442749    15.06714249  237.53713989]
     [ 426.88571167   15.22571373  237.47570801]
     [ 426.90570068   15.15857124  237.44285583]
     [ 426.75427246   15.26714325  237.44142151]
     [ 427.12286377   15.10857201  237.41999817]
     [ 427.16430664   15.28857136  237.35858154]
     [ 426.79858398   15.33000088  237.34571838]
     [ 426.95858765   15.26571369  237.34143066]
     [ 426.99429321   15.20714283  237.32713318]
     [ 426.67285156   15.55428505  237.28141785]
     [ 426.94857788   15.18714333  237.27999878]
     [ 427.09716797   15.38714218  237.17857361]
     [ 426.53857422   15.37285709  237.11143494]
     [ 426.64715576   15.35000038  236.78427124]
     [ 427.05999756   15.08571434  236.7285614 ]]


### Filter out subhalos:

Note that the above halo is actually a subhalo since its `parent_id` is not `-1`.  Let's see what fraction of Rockstar halos are actually subhalos:


```python
(halos['parent_id'] != -1).mean()
```




    0.094740065311804428



So about 10% are subhalos.  For some analyses you might want to only include top-level parent halos, since Rockstar halo masses always include substructure mass:


```python
print '# halos before subhalo filtering:', len(halos)
halos = halos[halos['parent_id'] == -1]
print '# halos after subhalo filtering:', len(halos)
```

    # halos before subhalo filtering: 9565499
    # halos after subhalo filtering: 8659263


### Load all redshifts of a FoF catalog:


```python
cats_by_z = Halos.make_catalogs(sim_name='emulator_1100box_planck_00',
                                products_dir='/mnt/alan1/lgarrison/bigsim_products/emulator_1100box_planck_products/',
                                redshifts='all', load_phases=False, load_subsamples=False, halo_type='FoF')
```

We can plot the evolution of the halo mass function with redshift:


```python
fig, ax = plt.subplots()
for z in sorted(cats_by_z.keys()):
    cat = cats_by_z[z]
    bin_edges, bin_centers, hist = Halos.halo_mass_hist(cat.halos[0]['N'])
    ax.loglog(bin_centers, hist, label='$z = {:.2f}$'.format(z))
ax.set_xlabel('Halo mass [# of particles]')
ax.set_ylabel('Number of halos')
ax.legend()
```




    <matplotlib.legend.Legend at 0x7fcc193639d0>




![png](/jupyter/AbacusCosmos_Python_Interface_Examples_files/AbacusCosmos_Python_Interface_Examples_17_1.png)


## [Halotools.py](https://github.com/lgarrison/AbacusCosmos/blob/master/AbacusCosmos/Halotools.py)


```python
from AbacusCosmos import Halotools
import halotools

import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np
```

### Load two different cosmologies at redshift 0.1:


```python
cats = Halotools.make_catalogs(sim_name='AbacusCosmos_720box', cosmologies=[0,1], redshifts=0.1,
                                products_dir='/mnt/store2/bigsim_products/AbacusCosmos_720box_products/',
                                phases=None, halo_type='Rockstar', load_halo_ptcl_catalog=False)
```

Generate mock galaxy catalogs for both:


```python
for cat in cats:
    # First apply an arbitrary mass cut to make the example run faster
    cat.halo_table = cat.halo_table[cat.halo_table['halo_N'] >= 100]
    # Make an approximate concentration column for mock purposes
    cat.halo_table['halo_conc'] = cat.halo_table['halo_rvir'] / cat.halo_table['halo_klypin_rs']
```


```python
from halotools.empirical_models import PrebuiltHodModelFactory
models = []
for cat in cats:
    model = PrebuiltHodModelFactory('zheng07', redshift=cats[0].redshift, concentration_key='halo_conc')
    model.populate_mock(cat)
    models += [model]
```

Compute the 2PCF on the mock galaxies:


```python
import Corrfunc
bins = np.logspace(-1,1.5,25)
bin_centers = (bins[:-1] + bins[1:])/2.
tpcfs = []
for model in models:
    gals = model.mock.galaxy_table
    results = Corrfunc.theory.xi(X=gals['x'], Y=gals['y'], Z=gals['z'],
                            boxsize=model.mock.BoxSize, nthreads=4,
                            binfile=bins)
    tpcfs += [results]
```


```python
plt.loglog(bin_centers, tpcfs[0]['xi'], label=cats[0].SimName)
plt.loglog(bin_centers, tpcfs[1]['xi'], label=cats[1].SimName)
plt.legend()
plt.xlabel(r'galaxy separation $s$ [Mpc/h]')
plt.ylabel(r'$\xi(s)$')
```




    <matplotlib.text.Text at 0x7fcc09fcaa50>




![png](/jupyter/AbacusCosmos_Python_Interface_Examples_files/AbacusCosmos_Python_Interface_Examples_27_1.png)


### Print the particle subsamples for the 1000th halo in the second cosmology:


```python
cats = Halotools.make_catalogs(sim_name='AbacusCosmos_720box', cosmologies=[0,1], redshifts=0.1,
                                products_dir='/mnt/store2/bigsim_products/AbacusCosmos_720box_products/',
                                phases=None, halo_type='Rockstar', load_halo_ptcl_catalog=True)
```


```python
halos = cats[1].halo_table
subsamples = cats[1].halo_ptcl_table
i = 999
```


```python
halos[i]
```




&lt;Row index=999&gt;
<table id="table140514219953936">
<thead><tr><th>halo_upid</th><th>halo_A2 [3]</th><th>halo_num_child_particles</th><th>halo_corevel [3]</th><th>halo_energy</th><th>halo_min_vel_err</th><th>halo_p_start</th><th>halo_desc</th><th>halo_rvmax</th><th>halo_rvir</th><th>halo_vx</th><th>halo_bullock_spin</th><th>halo_vmax_r</th><th>halo_min_pos_err</th><th>halo_m_pe_d</th><th>halo_b_to_a2</th><th>halo_m_pe_b</th><th>halo_child_r</th><th>halo_c_to_a2</th><th>halo_mgrav</th><th>halo_A [3]</th><th>halo_alt_m [4]</th><th>halo_mvir</th><th>halo_num_p</th><th>halo_subsamp_start</th><th>halo_N_SO</th><th>halo_spin</th><th>halo_subsamp_len</th><th>halo_b_to_a</th><th>halo_x</th><th>halo_z</th><th>halo_m_SO</th><th>halo_c_to_a</th><th>halo_n_core</th><th>halo_alt_N [4]</th><th>halo_min_bulkvel_err</th><th>halo_vrms</th><th>halo_Voff</th><th>halo_bulkvel [3]</th><th>halo_alt_N_SO [4]</th><th>halo_vmax</th><th>halo_alt_m_SO [4]</th><th>halo_y</th><th>halo_Xoff</th><th>halo_id</th><th>halo_halfmass_radius</th><th>halo_vy</th><th>halo_vz</th><th>halo_flags</th><th>halo_klypin_rs</th><th>halo_rs</th><th>halo_kin_to_pot</th><th>halo_N</th><th>halo_J [3]</th><th>halo_hostid</th></tr></thead>
<thead><tr><th>int64</th><th>float32</th><th>int64</th><th>float32</th><th>float32</th><th>float32</th><th>int64</th><th>int64</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>int64</th><th>int64</th><th>int32</th><th>float32</th><th>int64</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>int64</th><th>int32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>int32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>int64</th><th>float32</th><th>float32</th><th>float32</th><th>int64</th><th>float32</th><th>float32</th><th>float32</th><th>int32</th><th>float32</th><th>int64</th></tr></thead>
<tr><td>569974</td><td>14.8675 .. -12.6811</td><td>207</td><td>-152.748 .. -599.933</td><td>-7.82279e+15</td><td>751.386</td><td>32769850</td><td>0</td><td>0.0342442</td><td>0.202554</td><td>-152.748</td><td>0.0234053</td><td>0.267052</td><td>4.50487e-05</td><td>3.90289e+11</td><td>0.947171</td><td>1.33911e+12</td><td>0.000278079</td><td>0.597965</td><td>9.18326e+11</td><td>-1.91093 .. 25.4669</td><td>9.18326e+11 .. 3.90289e+11</td><td>1.44636e+12</td><td>207</td><td>50850</td><td>841</td><td>0.020859</td><td>17</td><td>0.655851</td><td>31.8827</td><td>145.677</td><td>9.6539e+12</td><td>0.514556</td><td>53</td><td>80 .. 34</td><td>751.386</td><td>203.216</td><td>18.2015</td><td>-152.748 .. -599.933</td><td>841 .. 46</td><td>195.099</td><td>9.6539e+12 .. 5.28038e+11</td><td>45.7946</td><td>13.8211</td><td>569976</td><td>0.059089</td><td>859.761</td><td>-599.933</td><td>25</td><td>0.0133996</td><td>0.0133996</td><td>0.708265</td><td>126</td><td>-5.92084e+11 .. 5.50002e+11</td><td>569974</td></tr>
</table>




```python
subsamples[halos['halo_subsamp_start'][i]:halos['halo_subsamp_start'][i] + halos['halo_subsamp_len'][i]]
```




&lt;Table length=17&gt;
<table id="table140514222685200" class="table-striped table-bordered table-condensed">
<thead><tr><th>y</th><th>x</th><th>vx</th><th>vy</th><th>vz</th><th>z</th></tr></thead>
<thead><tr><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th></tr></thead>
<tr><td>46.0763</td><td>31.5762</td><td>-469.43</td><td>693.177</td><td>-122.842</td><td>145.953</td></tr>
<tr><td>45.9182</td><td>31.9297</td><td>-362.675</td><td>732.662</td><td>-49.7216</td><td>145.876</td></tr>
<tr><td>45.7995</td><td>31.8661</td><td>-356.46</td><td>859.617</td><td>-504.528</td><td>145.8</td></tr>
<tr><td>46.2259</td><td>31.8782</td><td>-213.51</td><td>829.181</td><td>-650.768</td><td>145.767</td></tr>
<tr><td>45.6115</td><td>31.8745</td><td>-301.62</td><td>885.666</td><td>-222.102</td><td>145.736</td></tr>
<tr><td>45.4404</td><td>32.0737</td><td>-425.01</td><td>585.417</td><td>-160.407</td><td>145.709</td></tr>
<tr><td>46.2333</td><td>31.8689</td><td>-276.394</td><td>731.2</td><td>-451.882</td><td>145.694</td></tr>
<tr><td>45.8359</td><td>31.7978</td><td>-87.744</td><td>830.643</td><td>-365.6</td><td>145.69</td></tr>
<tr><td>45.8191</td><td>31.9025</td><td>-235.446</td><td>1061.7</td><td>-548.4</td><td>145.668</td></tr>
<tr><td>45.5405</td><td>31.9914</td><td>-490.818</td><td>774.615</td><td>-296.136</td><td>145.663</td></tr>
<tr><td>45.7911</td><td>31.8062</td><td>-242.667</td><td>549.771</td><td>-659.451</td><td>145.656</td></tr>
<tr><td>45.7752</td><td>31.9128</td><td>-315.33</td><td>913.086</td><td>-748.566</td><td>145.624</td></tr>
<tr><td>45.5246</td><td>31.7501</td><td>-677.274</td><td>952.845</td><td>-346.863</td><td>145.624</td></tr>
<tr><td>45.7658</td><td>31.8362</td><td>-153.552</td><td>686.871</td><td>-652.596</td><td>145.62</td></tr>
<tr><td>45.8537</td><td>31.9867</td><td>-225.21</td><td>884.752</td><td>-67.2704</td><td>145.575</td></tr>
<tr><td>45.8023</td><td>31.8829</td><td>50.727</td><td>793.809</td><td>-788.325</td><td>145.562</td></tr>
<tr><td>45.7901</td><td>31.9222</td><td>-227.586</td><td>1125.59</td><td>-667.677</td><td>145.341</td></tr>
</table>



</div>
