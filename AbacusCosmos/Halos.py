"""
File:    Halos.py
Author:  Lehman Garrison
Website: http://lgarrison.github.io/AbacusCosmos

Overview
--------
This file provides an interface to FoF and Rockstar halo catalogs and particle subsamples
from the AbacusCosmos project.  Using this interface is completely optional; users who want to
interact directly with the catalogs can find the file formats documented here:
http://lgarrison.github.io/AbacusCosmos/data_specifications/

Note that Halotools.py provides a wrapper around the Halos.py interface to load catalogs into
Halotools format (http://halotools.readthedocs.io/), if you wish to use that package.


Primary Interfaces
------------------
The main function for loading catalogs is `make_catalogs()`, which can load a set of catalogs
from one or more simulations, possibly including multiple redshifts and IC phases.

`make_catalog_from_dir()` loads a single halo catalog; `make_catalogs()` is simply a wrapper
that calls this function many times to load a whole set of catalogs.

The rest of the functions are related to loading halo and particle catalogs from specific file
formats; for example, `read_halos_FoF()` loads FoF halos, and `read_halos_Rockstar()` loads
Rockstar halos.  Most of the time, there's no need to call these functions directly;
`make_catalog_from_dir()` will call the correct one based on the `halo_type` argument.


Examples
--------
See http://lgarrison.github.io/AbacusCosmos/code


Quick example
-------------
>>> from AbacusCosmos import Halos
>>> halo_catalogs = Halos.make_catalogs('emulator_1100box_planck_00',
                                        '/mnt/store2/bigsim_products/emulator_1100box_planck_products/',
                                        load_phases=True, load_subsamples=False,
                                        redshifts=[0.3, 0.5])
>>> cat_at_z_half = halo_catalogs[0.5]
>>> halos_3rd_phase = cat_at_z_half.halos[2]
>>> box = 1.*cat_at_z_half.header['BoxSize']

"""

import numpy as np
from glob import glob
import os
import os.path as path
from os.path import join as pjoin
import re
from itertools import izip, cycle
from collections import defaultdict, OrderedDict
import h5py
import numpy.lib.recfunctions as rfn
import seaborn
import astropy.cosmology as ac
import astropy.units as u
import matplotlib.pyplot as plt
import cycler

from InputFile import InputFile
from Reglob import reglob


def make_catalogs(sim_name, products_dir, label='', load_phases=False, load_subsamples=False, redshifts='all', suffix='', downsampled=1, halo_type='FoF', color=None, linestyle='-'):
    """
    Read halo catalogs for a simulation, possibly with multiple z and phases.
    This is the primary entry point for most use cases.
    
    The halo catalogs are guaranteed to have at least 'N', 'id', 'pos', 'vel'.
    The other fields are determined by which halo finder the catalogs came from.
    
    Parameters
    ----------
    sim_name: str
        The name of the simulation
    products_dir: str
        The products directory that contains the 'sim_name_products' directory
    label: str, optional
        Plotting/descriptive label
    load_phases: bool or int, optional
        Load other phases of this sim.  An int will load that many phases; `True` will load all.
    load_subsamples: bool, optional
        Load the halo particle subsamples
    halo_type: str, optional
        'FoF', 'Rockstar', or 'Rockstar_SO'
    redshifts: str, int, or list of int, optional
        'all', or a redshift, or a list of redshifts to load.  Default: 'all'
    suffix: str, optional
        Set this if the halo catalog has a suffix (e.g. _downsampled)
    downsampled: int, optional
        The downsample-per-dimension factor.
        Set this if the sim was downsampled before halo finding.
        This will increase the particle mass assumed when converting to
        particle counts.
    color: str, optional
        Color for plotting.  Default: next color in cycle
    linestyle: str, optional
        Linestyle for plotting. Default: '-'
        
    Returns
    -------
    cats_by_z: dict of Catalog
        Catalog for this sim indexed by redshift.
        cats_by_z[z].halos is a list N_phases long.
    """
    
    phase_regex = ''
    if load_phases:
        sim_name = re.sub(r'-\d+$', '', sim_name)
        sim_name += '{0}'
        phase_regex = '(-\d+)?'

    halo_type = halo_type.lower()
        
    tag = {'fof':r'_FoF_halos',
           'rockstar':r'_rockstar_halos',
           'rockstar_so':r'_rockstar_halos'}
    slice_pattern = pjoin(products_dir, sim_name + '_products', sim_name + tag[halo_type] + suffix, 'z*')
    
    # Read halo catalogs
    paths = reglob(slice_pattern, phase_regex)
    if type(load_phases) is int:
        paths = paths[:load_phases]
    if redshifts != 'all':
        if type(redshifts) != list:
            redshifts = [redshifts]
        paths = [p for p in paths if redshift_from_path(p) in redshifts]
    cats = [make_catalog_from_dir(path, downsampled=downsampled, load_subsamples=load_subsamples,
                                     halo_type=halo_type)
            for path in paths]
    assert len(cats) != 0, 'No catalogs found!'
    
    # Aggregate by redshift
    cats_by_z = defaultdict(list)
    for c in cats:
        cats_by_z[c.z].append(c)
    
    # Merge multiple phases into a single catalog
    # Halos is now a list of halo arrays, one per phase
    for z,cats in cats_by_z.iteritems():
        cats[0].halos = [c.halos for c in sorted(cats, key=lambda x:x.header.ZD_Seed)]
        if load_subsamples:
            cats[0].subsamples = [c.subsamples for c in sorted(cats, key=lambda x:x.header.ZD_Seed)]
        cats_by_z[z] = cats[0]
    cats_by_z = dict(cats_by_z)
            
    # Plot formatting
    if not color:
        color = color_cycle.next()
    for cat in cats_by_z.itervalues():
        cat.label = label if label else cat.header.SimName
        cat.linestyle = linestyle
        cat.color = color
        
    return cats_by_z


class Catalog:
    '''
    A container class that will be populated with a halo catalog by `make_catalog_from_dir()`
    '''
    pass


def make_catalog_from_dir(dirname, downsampled=1, load_subsamples=False, load_uniform_subsample=False, load_pids=False, halo_type='auto'):
    """
    Load a single halo catalog.  Use `vars(c)` to see the full list of loaded fields.  Some important ones are listed here.
    
    Some important fields:
    `c.halos`:             the halo catalog as a Numpy structured array
    `c.header`:            the InputFile containing information from the header (e.g. simulation and time slice parameters)
    `c.halo_finder_cfg`:   the InputFile containing information about how the halo finder was invoked
    `c.cosmology`:         an astropy.cosmology object representing the simulation cosmology
    `c.subsamples`:        a Numpy array containing the halo particle subsamples
    `c.uniform_subsample`: a Numpy array containing a uniform sample of particles
    
    Parameters
    ----------
    dirname: str
        A redshift slice directory containing a halo catalog
    downsampled: int, optional
        The downsample-per-dimension factor.
        Set this if the sim was downsampled before halo finding.
        This will increase the particle mass assumed when converting to
        particle counts.
    halo_type: str, optional
        'FoF', 'Rockstar', or 'Rockstar_SO'.  Default: 'auto' (try to guess from `dirname`)
    load_subsamples: bool, optional
        Load halo particle subsamples into `c.subsamples`
    load_uniform_subsample: bool, optional
        Load the uniform particle subsample into `c.uniform_subsample`
    load_pids: bool, optional
        For any subsamples, load the particle IDs into `c.subsample_pids`
        
    Returns
    -------
    c: Catalog
        Catalog containing the halos and some metadata
    """
    dirname = os.path.abspath(dirname)
    
    halo_type = halo_type.lower()
    hts = ['fof', 'rockstar', 'rockstar_so']
    if halo_type == 'auto':
        for ht in hts:
            if ht in dirname.lower():
                halo_type = ht
    if halo_type not in hts:
        raise ValueError(halo_type)
        
    # read the header to get the redshift
    header = InputFile(pjoin(dirname, 'header'))
    z = float(np.round(header.Redshift, decimals=3))

    # could extract from the dirname name to avoid roundoff error
    #z = redshift_from_path(dirname)

    c = Catalog()
    c.z = z
    c.header = header
    c.dirname = dirname
    c.halo_type = halo_type
    
    c.halo_finder_cfg = get_halo_finder_cfg(dirname, halo_type=halo_type)
    
    # A common pathology is that BoxSize is interpreted as an int by ParseHeader
    # We should endeavor to always write "50." instead of "50" in the .par files
    for field in ['BoxSize', 'InitialRedshift', 'ZD_PLT_target_z', 'wa']:
        setattr(c.header, field, float(getattr(c.header, field)))
    
    # Make common params accessible as cat.FIELD
    for k in vars(header):
        if k in important_header_fields:
            if k in header:
                setattr(c, k, header[k])
            
    # Presently, Flatw0waCDM is the most general cosmology that Abacus supports
    try:
        assert c.header.omnuh2 == 0.  # need to actually know the individual neutrino masses if non-zero
        Ob0 = c.header.ombh2 / (c.header.H0/100.)**2.
        Neff = c.header.N_eff
    except AttributeError:
        Ob0 = 0.0486  # Might not have saved the CAMB params.  Just use Planck params.
        Neff = 3.046
        pass
    c.cosmology = ac.Flatw0waCDM(H0=c.header.H0, Om0=c.header.Omega_M, w0=c.header.w0,
                                 wa=c.header.wa, Neff=Neff, m_nu=0.*u.eV, Ob0=Ob0)
    
    c.halos = halo_readers[halo_type](dirname, particle_mass=c.header.ParticleMassHMsun*downsampled**3, boxsize=c.BoxSize)
    
    if load_subsamples:
        c.subsamples = subsample_readers[halo_type](dirname, c.halos, load_pids=load_pids, boxsize=c.BoxSize)
    
    if load_uniform_subsample:
        c.uniform_subsample = uniform_subsample_readers[halo_type](dirname, load_pids=load_pids, boxsize=c.BoxSize)
    
    return c

def wrap_zero_centered(pos, box):
    """
    Wraps an array of positions in-place
    to the range [-box/2, box/2).
    
    Parameters
    ----------
    pos: ndarray
        The positions to wrap, already in the same units as box
    box: float
        The box edge length (box is zero-centered)
    """
    
    # numpy.round rounds to nearest even, which makes life harder
    pos -= box*np.round(pos/box)
        
    # Fix particles at box/2., which did not get rounded down because half rounds to zero
    pos[pos >= box/2.] -= box
    pos[pos < -box/2.] += box
    

def read_halos_FoF(dirname, **kwargs):
    """
    Loads friends-of-friends halos into a numpy array.
    
    Parameters
    ----------
    dirname: str
        The slice directory containing the halo catalog.
    halo_regex: str, optional
        The regex for the halo file number (i.e. {} for 'halos_{}').
    ret_fns: bool, optional
        Return the halo filenames that were read.
    boxsize: float, optional
        Wrap positions to zero-centered boxsize box.  Default: None
        
    Returns
    -------
    halos: ndarray
        An array of halos
    halo_fns: list of str, optional
        If `ret_fns`, returns the filenames of the halo files that were concatenated
        to build the catalog.
    """
    halo_regex = kwargs.get('halo_regex', r'\d+$')
    ret_fns = kwargs.get('ret_fns', False)
    boxsize = kwargs.get('boxsize', None)
    
    all_halos = []
    halo_fns = reglob(path.join(dirname, 'halos_{}'), halo_regex)
    for halo_fn in halo_fns:
        with open(halo_fn,"rb") as fp:
            # the file begins with two ints
            num_groups = np.fromfile(fp,dtype=np.uint64,count=1)
            n_largest_subhalos = np.fromfile(fp,dtype=np.uint64,count=1)
            assert (n_largest_subhalos,) == halo_dt_FoF()['subhalo_N'].shape, n_largest_subhalos

            halos = np.fromfile(fp, dtype=halo_dt_FoF())
            
            assert len(halos) == num_groups, (len(halos), num_groups, halo_fn)
            
            all_halos.append(halos)
    
    halos = np.concatenate(all_halos)
    
    if boxsize:
        wrap_zero_centered(halos['pos'], boxsize)
    
    if ret_fns:
        return halos, halo_fns
    return halos
    
    
def read_subsamples_FoF(dirname, halos, halo_regex=r'\d+', load_pids='auto', boxsize=None):
    """
    Load the particle subsamples for a set of halos.
    
    The input `halos` will have the `subsamp_*` index fields
    modified so that they correspond to indices in the returned
    subsample array.
    
    To get the subsample for an individual halo:
    subsamples[halos[i]['subsamp_start']:halos[i]['subsamp_start']+halos[i]['subsamp_len']]
    
    To split the returned array into individual subsample arrays:
    split_points = halos['subsamp_len'].cumsum()[:-1]
    subsamples = np.split(particles, split_points)
    
    Parameters
    ----------
    dirname: str
        The slice directory containing the halo catalog and subsamples
    halos: ndarray
        The FoF halo catalog from `read_halos_FoF()`.
        The halo catalog must be in its original order (i.e. the order
        on disk, because the subsamples on disk are in the same order)
    halo_regex: str, optional
        The regex for the halo file number (i.e. {} for 'halos_{}').
    load_pids: bool or str, optional
        Read the files containing the particle ids for the halo subsamples.
        Default of 'auto' loads them if present.
    boxsize: float, optional
        Wrap positions to zero-centered boxsize box.  Default: None
        
    Returns
    -------
    subsamples: ndarray
        An array holding the concatenated particle subsample positions
        and velocities, and possibly PIDs as well.
    """
    # Read all the subsample particles into a big array
    all_particles = []
    for ss_fn in reglob(pjoin(dirname, 'particles_{}'), halo_regex):
        all_particles.append(np.fromfile(ss_fn, dtype=[('pos',np.float32,3), ('vel',np.float32,3)]))
    particles = np.concatenate(all_particles)
    nfiles = len(all_particles)
    del all_particles
    
    pid_fns = reglob(pjoin(dirname, 'particle_ids_{}'), halo_regex)
    if load_pids == 'auto':
        load_pids = bool(pid_fns)
    
    if load_pids:
        all_pids = []
        for pid_fn in pid_fns:
            all_pids.append(np.fromfile(pid_fn, dtype=[('pid',np.uint64)]))
        pids = np.concatenate(all_pids)
        npidfiles = len(all_pids)
        assert npidfiles == nfiles
        del all_pids
        # add the 'pid' field to the existing particles array
        # somewhat inefficient compared to reading directly into a pre-allocated buffer
        particles = rfn.append_fields(particles, 'pid', pids['pid'], usemask=False)
        
    
    # Make sure we read the number of particles that the halos implied
    assert halos['subsamp_len'].sum() == len(particles), \
        "Read {} subsample particles, but the halo catalog expected {}".format(len(particles), halos['subsamp_len'].sum())
    
    njump = reindex_halo_subsamples(halos)
    assert njump == nfiles - 1, \
        "{} jumps in subsample indexing, but {} file splits.  Were the halos reordered?".format(njump, nfiles-1)
        
    if boxsize:
        wrap_zero_centered(particles['pos'], boxsize)
        
    return particles


def reindex_halo_subsamples(halos, subsamp_len_key='subsamp_len'):
    """
    If we concatenate halos and particles into big files/arrays, the "subsample start" indices
    in the halos table no longer correspond to the concatenated particle array.  But we can
    easily reconstruct the correct indices.
    
    Parameters
    ----------
    halos: ndarray
        The FoF halo catalog from `read_halos_FoF`
        The halo catalog must be in its original order (i.e. the order
        on disk, because the subsamples on disk are in the same order)
        
    Returns
    -------
    njump: int
        The number of discontinuities that were fixed in the subsample indexing.
        Should be equal to the number of file splits.
    """
    
    # The number of "discontinuities" in particle indexing should equal the number of files
    # This helps us make sure the halos were not reordered
    njump = (halos['subsamp_start'][:-1] + halos[subsamp_len_key][:-1] != halos['subsamp_start'][1:]).sum()
    
    # Now reindex the halo records
    halos['subsamp_start'][0] = 0
    halos['subsamp_start'][1:] = halos[subsamp_len_key].cumsum()[:-1]
    
    return njump
    

def read_uniform_subsample_FoF(dirname, halo_regex=r'\d+', load_pids='auto', boxsize=None):
    """
    Read the uniform particle subsample from the given catalog
    directory by concatenating field and halo particles.
    
    Parameters
    ----------
    dirname: str
        The catalog directory containing the `halo_*.field`
        and `halo_*.particles`.
    halo_regex: str, optional
        The regex for the halo file number (i.e. {} for 'halos_{}').
    load_pids: bool or str, optional
        Read the files containing the particle ids for the halo subsamples.
        Default of 'auto' loads them if present.
    boxsize: float, optional
        Wrap positions to zero-centered boxsize box.  Default: None
        
    Returns
    -------
    subsample: ndarray
        The uniform particle subsample
    """
    all_particles = []
    for ss_fn in reglob(pjoin(dirname, 'particles_{}'), halo_regex) + reglob(pjoin(dirname, 'field_particles_{}'), halo_regex):
        all_particles.append(np.fromfile(ss_fn, dtype=[('pos',np.float32,3), ('vel',np.float32,3)]))
    particles = np.concatenate(all_particles)
    nfiles = len(all_particles)
    del all_particles
    
    pid_fns = reglob(pjoin(dirname, 'particle_ids_{}'), halo_regex) + reglob(pjoin(dirname, 'field_ids_{}'), halo_regex)
    if load_pids == 'auto':
        load_pids = bool(pid_fns)
    
    if load_pids:
        all_pids = []
        for pid_fn in pid_fns:
            all_pids.append(np.fromfile(pid_fn, dtype=[('pid',np.uint64)]))
        pids = np.concatenate(all_pids)
        npidfiles = len(all_pids)
        assert npidfiles == nfiles
        del all_pids
        # add the 'pid' field to the existing particles array
        # somewhat inefficient compared to reading directly into a pre-allocated buffer
        particles = rfn.append_fields(particles, 'pid', pids['pid'], usemask=False)
    
    if boxsize:
        wrap_zero_centered(particles['pos'], boxsize)

    return particles
    
    
def read_halos_Rockstar(dirname, particle_mass=None, **kwargs):
    """
    Load Rockstar halos.  The format is the same as the Rockstar
    internal representation, except we add the 'N' and 'alt_N' particle counts,
    and equivalents for SO masses and particle counts.
    
    Our old format stored SO masses in separate .bin files;
    these will be detected and loaded automatically.
    
    Parameters
    ----------
    dirname: str
        The slice directory containing the Rockstar catalog
    particle_mass: float
        The particle mass to assume when converting Rockstar physical mass values
        into particle counts.
    boxsize: float, optional
        Wrap positions to zero-centered `boxsize` box.  Default is to read
        `boxsize` from the HDF5 headers.
        
    Returns
    -------
    halos: ndarray
        An array of halos
    """
    
    all_halos = []
    
    for h5halo_fn in reglob(pjoin(dirname, 'halos_{}.h5'),'.*'):
        # Load the hdf5 halos
        f = h5py.File(h5halo_fn)
        h5halos = f['halos'][:]
        
        if not particle_mass:
            particle_mass = f['halos'].attrs['ParticleMassHMsun']
        boxsize = kwargs.get('boxsize', f['halos'].attrs['BoxSize'])
        
        # Check that the mass is a whole number
        assert np.allclose(h5halos['N'] * particle_mass, h5halos['m'], rtol=1e-4)
        assert np.allclose(h5halos['N_SO'] * particle_mass, h5halos['m_SO'], rtol=1e-4)

        all_halos.append(h5halos)
        
    halos = np.concatenate(all_halos)
    
    #wrap_zero_centered(halos['pos'], boxsize)
    
    return halos
        
    
# Thin wrapper to put the SO masses in the main mass field
def read_halos_Rockstar_SO(*args, **kwargs):
    halos = read_halos_Rockstar(*args, **kwargs)
    halos = rfn.rename_fields(halos, {'m':'m_rockstar', 'alt_m':'alt_m_rockstar',
                                      'N':'N_rockstar', 'alt_N':'alt_N_rockstar'})
    halos = rfn.rename_fields(halos, {'m_SO':'m', 'alt_m_SO':'alt_m',
                                      'N_SO':'N', 'alt_N_SO':'alt_N'})
    return halos
    

def read_subsamples_Rockstar(dirname, halos, halo_regex=r'\d+\.\d+', load_pids=False, boxsize=None):
    """
    Load the particle subsamples for a set of halos.
    Wraps the particles to the zero-centered box.
    Note that the PIDs are in the same file as the positions,
    unlike FoF.
    
    Parameters
    ----------
    dirname: str
        The slice directory containing the halo catalog and subsamples
    halos: ndarray
        The Rockstar-based halo catalog
    load_pids: bool, optional
        Load the particle IDs for the subsampled particles.
        Default: False.
    boxsize: float, optional
        Wrap positions to zero-centered `boxsize` box.  Default is to read
        `boxsize` from the HDF5 headers.
        
    Returns
    -------
    subsamples: ndarray
        The particle subsamples
    """
    
    # Read all the subsample particles into a big array
    all_particles = []
    for ss_fn in reglob(pjoin(dirname, 'particles_{}.h5'), halo_regex):
        with h5py.File(ss_fn, 'r') as hfp:
            try:
                particles = hfp['particles'][:]
            except KeyError:
                particles = hfp['subsamples'][:]
            if not load_pids:
                particles = particles[['pos','vel']]
            if not boxsize:
                boxsize = particles.attrs['BoxSize']
        all_particles.append(particles)
    nfiles = len(all_particles)
    particles = np.concatenate(all_particles)
    del all_particles
    
    # Make sure we read the number of particles that the halos implied
    assert halos['subsamp_len'].sum() == len(particles), \
        "Read {} subsample particles, but the halo catalog expected {}".format(len(particles), halos['subsamp_len'].sum())
    
    njump = reindex_halo_subsamples(halos, subsamp_len_key='subsamp_len')
    assert njump == nfiles - 1, \
        "{} jumps in subsample indexing, but {} file splits.  Were the halos reordered?".format(njump, nfiles-1)
        
    #wrap_zero_centered(particles['pos'], boxsize)
        
    return particles
    
def get_halo_finder_cfg(dirname, halo_type='', cfg_fn=''):
    '''
    Parse the input parameter file for a halo finder.
    
    Parameters
    ----------
    dirname: str
        The catalog directory path
    halo_type: str, optional
        The halo finder.  One of ['fof', 'rockstar', 'rockstar_so', 'fof'].
        This is ignored if 'cfg_fn' is given.  Default: ''
    cfg_fn: str
        The name/globbing pattern of the parameter file to load.
        Overrides `halo_type`.  Default: ''
        
    Returns
    -------
    cfg: InputFile
        The parsed halo finder config.
    '''
    assert halo_type or cfg_fn
    
    default_fns = {'fof':'fof.cfg',
                   'rockstar':'rockstar.cfg',
                   'rockstar_so':'rockstar.cfg'}
    
    if not cfg_fn:
        halo_type = halo_type.lower()
        cfg_fn = default_fns[halo_type]
        
    param_fns = reglob(path.join(dirname, '{}'), cfg_fn)
    assert len(param_fns) == 1
    
    return InputFile(param_fns[0])


def redshift_from_path(path):
    z = re.findall(r'(?<=z).*\b', path)
    assert len(z) == 1
    z = float(z[0])
    return z


halo_readers = {'fof':read_halos_FoF,
                'rockstar':read_halos_Rockstar,
                'rockstar_so':read_halos_Rockstar_SO}
subsample_readers = {'fof':read_subsamples_FoF,
                     'rockstar':read_subsamples_Rockstar,
                     'rockstar_so':read_subsamples_Rockstar}
uniform_subsample_readers = {'fof':read_uniform_subsample_FoF,
                             'rockstar':None,  # not implemented
                             'rockstar_so':None}
    
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
color_cycle = cycler.cycle(colors)
del colors, prop_cycle
    
def halo_mass_hist(masses, bin_width_dex=.15, min_mass=None, max_mass=None, bin_edges=None):
    '''
    Take a single list of masses
    and compute the halo mass function
    
    Returns
    -------
    bin_edges, bin_centers, hist
    '''
    if not min_mass:
        min_mass = masses.min()
    if not max_mass:
        max_mass = masses.max()
        
    if bin_edges is None:
        nbins = int(round(np.log10(1.*max_mass/min_mass)/bin_width_dex))
        bin_edges = np.logspace(np.log10(min_mass),np.log10(max_mass), nbins)
        
    hist, be = np.histogram(masses, bins=bin_edges)
    bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.

    return bin_edges, bin_centers, hist
    
# A list of the header fields relevant for simulation analysis
# These will be available as cat.FIELD, while the rest will be in cat.header.FIELD
important_header_fields = ['SimName',
                           'BoxSize',
                           #'H0',  # included in cat.cosmology
                           'InitialRedshift',
                           'NP',
                           #'N_eff',
                           #'Omega_DE',
                           #'Omega_M',
                           #'Omega_K',
                           'SofteningLength',
                           'SofteningType',
                           'ZD_Seed',
                           'ns',
                           #'omnuh2',
                           #'ombh2',
                           #'omch2',
                           'sigma_8',
                           #'w0',
                           #'wa',
                           'w',
                           'ParticleMassHMsun',
                           'Growth',
                           'Growth_on_a',
                           'HubbleNow',
                           'OmegaNow_DE',
                           'OmegaNow_K',
                           'OmegaNow_m',
                           'Redshift',
                           'f_growth']

# dtypes are globally shared unless we create a new object every time
# so we require this be called as a function
halo_dt_FoF = lambda: np.dtype([("id",np.int64),("subsamp_start",np.uint64),("subsamp_len",np.uint32),
               ("N",np.uint32),("subhalo_N",np.uint32,4),("pos",np.float32,3),("vel",np.float32,3),("sigma_v",np.float32,3),
               ("r25",np.float32),("r50",np.float32),("r75",np.float32),("r90",np.float32),
               ("vcirc_max",np.float32),("rvcirc_max",np.float32),
               ("subhalo_pos",np.float32,3),("subhalo_vel",np.float32,3),("subhalo_sigma_v",np.float32,3),
               ("subhalo_r25",np.float32),("subhalo_r50",np.float32),("subhalo_r75",np.float32),("subhalo_r90",np.float32),
               ("subhalo_vcirc_max",np.float32),("subhalo_rvcirc_max",np.float32)], align=True)
