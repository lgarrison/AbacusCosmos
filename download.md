---
layout: default
---
# Overview
All simulation data products can be browsed at <http://nbody.rc.fas.harvard.edu/public/AbacusCosmos/>.

You'll probably want to use a utility like `wget` or `curl` to recursively download the catalogs.
The following tool will help you build a `wget` command to download the desired files.

We would like to thank Harvard Research Computing for their assistance in hosting the catalogs.

# Command-line Builder
<script src="{{ site.baseurl }}/assets/js/clipboard.min.js"></script>
<script>
new Clipboard('.btn');
</script>

<div id="regex-builder">
<div id="text-and-button">
<!-- Target -->
<textarea id="wget-command" value="" wrap="off" readonly></textarea>
<!-- Trigger -->
<br>
<button class="copy-button" type="button" data-clipboard-target="#wget-command">
    <img class="clippy" src="{{ site.baseurl }}/assets/images/clippy.svg" width="30" alt="Copy to clipboard">
</button>
</div>
<br>

<div class="checkbox-header">
<h2>Simulation Sets</h2>
<form name="sims" class="checkbox-group">
<!-- note: maybe could generate this from yml? -->
<label><input class="chk" type="checkbox" data-path="AbacusCosmos_1100box_products" checked /> AbacusCosmos_1100box </label>
<label><input class="chk" type="checkbox" data-path="AbacusCosmos_720box_products"/> AbacusCosmos_720box </label>
<label><input class="chk" type="checkbox" data-path="emulator_1100box_planck_products"/> emulator_1100box_planck </label>
<label><input class="chk" type="checkbox" data-path="emulator_720box_planck_products"/> emulator_720box_planck </label>
</form>
</div>

<div class="checkbox-header">
<h2>Data Products</h2>
<form name="products" class="checkbox-group">
<label><input class="chk" type="checkbox" data-product="power" checked /> Matter power spectrum </label>
<label><input id="fofchk" class="chk" type="checkbox" data-product="FoF_halos"/> Friends-of-friends halos</label>
<ul class="normal">
  <li><label><input name="fofchk_sub" class="subchk" type="checkbox" data-fn="halos.tar.gz"/>Halos</label></li>
  <li><label><input name="fofchk_sub" class="subchk" type="checkbox" data-fn="halo_subsamples.tar.gz"/>Halo particles</label></li>
  <li><label><input name="fofchk_sub" class="subchk" type="checkbox" data-fn="halo_subsample_ids.tar.gz"/>Halo particle IDs</label></li>
  <li><label><input name="fofchk_sub" class="subchk" type="checkbox" data-fn="field_subsamples.tar.gz"/>Field particles</label></li>
  <li><label><input name="fofchk_sub" class="subchk" type="checkbox" data-fn="field_subsample_ids.tar.gz"/>Field particle IDs</label></li>
</ul>
<label><input id="rockchk" class="chk" type="checkbox" data-product="rockstar_halos"/> Rockstar halos</label>
<ul class="normal">
  <li><label><input name="rockchk_sub" class="subchk" type="checkbox" data-fn="halos.tar.gz"/>Halos</label></li>
  <li><label><input name="rockchk_sub" class="subchk" type="checkbox" data-fn="halo_subsamples.tar.gz"/>Halo particles</label></li>
  <li><label><input name="rockchk_sub" class="subchk" type="checkbox" data-fn="halo_subsample_ids.tar.gz"/>Halo particle IDs</label></li>
</ul>
</form>
</div>

<div class="checkbox-header">
<h2>Redshifts</h2>
<form name="redshifts" class="checkbox-group" id="redshifts">
<label><input class="chk" type="checkbox" data-redshift="z0.100"/> z = 0.1</label>
<label><input class="chk" type="checkbox" data-redshift="z0.300"/> z = 0.3 </label>
<label><input class="chk" type="checkbox" data-redshift="z0.500" checked/> z = 0.5</label>
<label><input class="chk" type="checkbox" data-redshift="z0.570"/> z = 0.57</label>
<label><input class="chk" type="checkbox" data-redshift="z0.700"/> z = 0.7</label>
<label><input class="chk" type="checkbox" data-redshift="z1.000"/> z = 1.0</label>
<label><input class="chk" type="checkbox" data-redshift="z1.500"/> z = 1.5</label>
</form>
</div>

</div>  <!-- regex display -->


<script src="{{ site.baseurl }}/assets/js/regex_builder.js"></script>

# Extracting data products
The larger data products (halos and particles) are stored in compressed formats (`*.tar.gz` files) to save
disk space and reduce file transfer time. You will probably need to decompress them to perform analyses.  Use the command
```bash
tar -xzvf file_to_decompress.tar.gz
```
in each directory that contains a tarball.

To automate this process on a whole directory tree, the Unix `find` command may be useful:
```bash
find directory_tree_root -type f -name '*.tar.gz' -execdir tar -xzf {} \;
```
To save disk space, you can add `-delete` to that command to remove the tar file after extraction.
