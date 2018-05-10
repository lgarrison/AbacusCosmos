"""
File:    Halotools.py
Author:  Lehman Garrison
Website: http://lgarrison.github.io/AbacusCosmos

Overview
--------
This file provides halotools (https://github.com/astropy/halotools) interfaces
for Abacus halo catalogs and particle samples. It is simply a wrapper around
Halos.py.


Primary Interfaces
------------------
The only class defined in this interface is `AbacusHaloCatalog`, which is a
subclass of the halotools `UserSuppliedHaloCatalog`.  Its purpose is to load
a single halo catalog (FoF or Rockstar) from a single sim at a single redshift.

The `make_catalogs()` function can be used to build `AbacusHaloCatalog`s for
multiple redshifts, cosmologies, and/or phases.

The `AbacusHaloCatalog` acts as a normal halotools halo catalog; that is,
halos are found in `cat.halo_table`, positions are wrapped to [0, Lbox), etc.
The uniform subsample, if requested, is stored in the standard `cat.ptcl_table`.
Since we also support halo particle subsamples, which halotools does not natively
support, we create the non-standard field `cat.halo_ptcl_table`.


Examples
--------
See http://lgarrison.github.io/AbacusCosmos/code


Quick Example
-------------
>>> from AbacusCosmos import Halotools
>>> halo_catalogs = Halotools.make_catalogs(sim_name='AbacusCosmos_1100box', cosmologies='all',
                                        redshifts=0.3, products_dir='/path/to/AbacusCosmos_1100box_products/',
                                        phases=None, halo_type='Rockstar')
>>> cat_for_cosmology_10 = halo_catalogs[10]
>>> halos = cat_for_cosmology_10.halo_table

"""

import halotools
from halotools.sim_manager import UserSuppliedHaloCatalog, UserSuppliedPtclCatalog
from halotools.empirical_models.phase_space_models.analytic_models.halo_boundary_functions import halo_mass_to_halo_radius

import Halos
from Reglob import reglob
import os
import numpy as np

def make_catalogs(sim_name, cosmologies, redshifts, products_dir, phases=None, halo_type='FoF', **cat_kwargs):
    """
    Load multiple catalogs for multiple redshifts, cosmologies, and/or phases into `AbacusHaloCatalog`s.
    This is just a convenience function; you can always load a catalog by calling the
    `AbacusHaloCatalog` constructor directly.
    
    A list of ints (or 'all') yields a list of catalogs, while a single number yields
    a bare catalog.  This applies to cosmologies, phases, and redshifts.
    
    If there are multiple redshifts, cosmologies, and/or phases, the lists are nested as:
    `catalogs[redshift][cosm][phase]`
    
    Parameters
    ----------
    sim_name: str
        The simulation name, without any phase or cosmology specifiers.
        e.g. 'emulator_1100box_Neff3', not 'emulator_1100box_Neff3_00-0'
    products_dir: str
        The directory containing the simulation product directories.
        e.g. 'emulator_1100box_Neff3_products', not 'emulator_1100box_Neff3_00_products'
    cosmologies: int, list of int, or 'all'
        One or more cosmologies to load.
    redshifts: float, list of float, or 'all'
        Redshifts to load
    phases: int, list of int, or 'all', optional
        One or more phases to load.  The default (None) loads the sim
        without any phase suffix like "-1".
    halo_type: str, optional
        'FoF', 'Rockstar', or 'Rockstar_SO'
    **cat_kwargs: dict, optional
        Keyword arguments to pass to `AbacusHaloCatalog`.
        
    Returns
    -------
    catalogs: `AbacusHaloCatalog`, or (nested) list(s) of `AbacusHaloCatalog`
        The catalogs(s).  If multiple cosmologies, phases, and/or redshifts
        were specified, then the lists are nested as
        `catalogs[redshift][cosmology][phase]`.
    """
    # Wrap all numeric inputs in lists
    return_flat_cosm = False
    all_cosm = cosmologies == 'all'
    if type(cosmologies) is int:
        cosmologies = [cosmologies]
        return_flat_cosm = True
    cstr = lambda c: '{:02d}'.format(c)
    
    return_flat_redshifts = False
    all_redshifts = redshifts == 'all'
    if type(redshifts) in (float, int):
        redshifts = [redshifts]
        return_flat_redshifts = True
    zstr = lambda z: '{:.3f}'.format(z)
    
    return_flat_phases = False
    all_phases = phases == 'all'
    if type(phases) is int or phases is None:
        phases = [phases]
        return_flat_phases = True
    pstr = lambda p: '' if p is None else '-{:d}'.format(p)
        
    tag = {'FoF':r'_FoF_halos',
           'Rockstar':r'_rockstar_halos',
           'Rockstar_SO':r'_rockstar_halos'}
    cat_dir_fmt = os.path.join(products_dir, sim_name + '_{cosm:s}{phase:s}_products', sim_name + '_{cosm:s}{phase:s}' + tag[halo_type], 'z{z:s}')
    
    def for_all_catalogs(action, redshifts=redshifts, cosmologies=cosmologies, phases=phases):
        if all_redshifts:
            redshifts = reglob(cat_dir_fmt, {'z':r'(\d+\.\d{3})', 'cosm':r'\d{2,}', 'phase':r'(?:-\d+)?'}, return_capture=1)
            redshifts = set(float(z) for z in redshifts)
            assert redshifts, "No catalogs found for any redshift. Search pattern was '{dir}'".format(dir=cat_dir_fmt)
        zcats = []
        for z in redshifts:
            if all_cosm:
                cat_dir = cat_dir_fmt.format(cosm='{cosm:s}', phase='{phase:s}', z=zstr(z))
                cosmologies = reglob(cat_dir, {'cosm':r'(\d{2,})', 'phase':r'(-\d+)?'}, return_capture=1)
                cosmologies = set(int(c) for c in cosmologies)
                assert cosmologies, "No cosmologies found for redshift {z}. Search pattern was '{dir}'".format(z=zstr(z), dir=cat_dir)
            ccats = []
            for c in cosmologies:
                if all_phases:
                    cat_dir = cat_dir_fmt.format(cosm=cstr(c), phase='{phase:s}', z=zstr(z))
                    phases = reglob(cat_dir, {'phase':r'(-\d+)?'}, return_capture=1)
                    phases = set((int(p[1:]) if p is not None else p) for p in phases)  # remove leading dash
                    assert phases, "No phases found for redshift {z}, cosmology {cosm}. Search pattern was '{dir}'".format(z=zstr(z), cosm=cstr(c), dir=cat_dir)
                pcats = []
                for p in phases:
                    cat_dir = cat_dir_fmt.format(cosm=cstr(c), phase=pstr(p), z=zstr(z))
                    result = action(cat_dir)
                    pcats.append(result)
                if return_flat_phases:
                    pcats = pcats[0]
                ccats.append(pcats)
            if return_flat_cosm:
                ccats = ccats[0]
            zcats.append(ccats)
        if return_flat_redshifts:
            zcats = zcats[0]
        return zcats
            
    # do quick validation that all directories exist
    def validate_dir(cat_dir):
        assert os.path.isdir(cat_dir), "Catalog directory '{}' does not exist".format(cat_dir)
    for_all_catalogs(validate_dir)
    
    def printcats(cat_dir):
        print cat_dir
    #for_all_catalogs(printcats)
    
    # now do the real catalog construction
    make_cat = lambda cat_dir: AbacusHaloCatalog(cat_dir, **cat_kwargs)
    zcats = for_all_catalogs(make_cat)

    return zcats
    

class AbacusHaloCatalog(UserSuppliedHaloCatalog):
    """
    Load an Abacus halo catalog in a format usable by Halotools.
    To load many catalogs, see `make_catalogs()`.
    """
    
    def __init__(self, path_to_cat, load_halo_ptcl_catalog=False, load_ptcl_catalog=False, load_pids='auto'):
        """
        Initialize an `AbacusHaloCatalog` by passing a halo catalog directory.
        To load many catalogs, see `make_catalogs()`.
        
        Note about FoF masses: Halotools depends on exact mass definitions
        (e.g. m500c) in several places.  Friends of friends does not use a
        strict overdensity definition, but we can convert our FoF linking length
        to an approximate average overdensity using the scaling relation of More+2011.
        Thus, the 'halo_m' field will appear as 'halo_mNNNm' in the halo table.
        
        Parameters
        ----------
        path_to_cat: str
            A directory containing a halo catalog
        load_halo_ptcl_catalog: bool, optional
            Load a list of `UserSuppliedPtclCatalog` (one for each halo)
            into `halocat.halo_ptcl_table`.
        load_ptcl_catalog: bool, optional
            Load a uniform particle subsample from the whole box
            into `halocat.ptcl_table`.  This is the usual meaning of
            a particle catalog in halotools (as opposed to halo subsamples).
        load_pids: bool or str, optional
            Load the particle IDs of the subsample catalogs.  Default
            of 'auto' loads them if present.
        """
        halocat = Halos.make_catalog_from_dir(path_to_cat, load_subsamples=load_halo_ptcl_catalog,
                                                           load_uniform_subsample=load_ptcl_catalog,
                                                           load_pids=load_pids)
        # Could convert units here, but would need to be careful to get every field, subsample, etc
        assert halocat.header.hMpc, "Positions are not in 1/h units! Halotools assumes this."
        
        # the important header fields will be stored as metadata
        metavars = vars(halocat)
        metadata = {k:metavars[k] for k in Halos.important_header_fields if k in metavars}
        metadata['redshift'] = halocat.header.Redshift
        metadata['Lbox'] = halocat.header.BoxSizeHMpc
        metadata['particle_mass'] = halocat.header.ParticleMassHMsun
        metadata['header'] = halocat.header
        metadata['halo_finder_cfg'] = halocat.halo_finder_cfg
        metadata['cosmology'] = halocat.cosmology
        
        halos = halocat.halos
        halos['pos'] %= halocat.header.BoxSize  # wrap to [0, Lbox)
        try:
            halos['subhalo_pos'] %= halocat.header.BoxSize
        except ValueError:
            pass  # no subhalo_pos field

        if halocat.halo_type == 'rockstar':
            # convert rockstar halo radii to Mpc/h from kpc/h
            kpc_fields = ['r', 'rs', 'klypin_rs', 'rvmax', 'child_r', 'halfmass_radius']
            for f in kpc_fields:
                try:
                    halocat.halos[f] /= 1000.
                except ValueError:
                    # halfmass_radius doesn't exist for the earlier catalogs because they came from an earlier version of Rockstar
                    pass
        
        # prepend 'halo_' to field names
        halos.dtype.names = [('halo_' + n) if not n.startswith('halo_') else n for n in halos.dtype.names]
        halo_catalog_columns = {n:halos[n] for n in halos.dtype.names}
        
        # Rockstar calls the primary mass field 'halo_m',
        # but halotools expects us to append a mass definition like "halo_mvir"
        # TODO: read the alternative mass definitions and modify the alt_m fields
        try:
            mdef = halocat.halo_finder_cfg['MASS_DEFINITION']
            mdef = mdef.lower().replace('b','m')
            halo_catalog_columns['halo_m{}'.format(mdef)] = halo_catalog_columns['halo_m']
            halo_catalog_columns['halo_r{}'.format(mdef)] = halo_catalog_columns['halo_r']
            del halo_catalog_columns['halo_m']
            del halo_catalog_columns['halo_r']
        except KeyError:
            pass  # not Rockstar
        
        # halotools requires that id is unique
        #assert len(np.unique(halo_catalog_columns['halo_id'])) == len(halo_catalog_columns['halo_id']), "halo_id not unique!"
        
        # break out x into x,y,z
        float3_columns = { 'halo_x':halo_catalog_columns['halo_pos'][:,0],  'halo_y':halo_catalog_columns['halo_pos'][:,1],  'halo_z':halo_catalog_columns['halo_pos'][:,2],
                          'halo_vx':halo_catalog_columns['halo_vel'][:,0], 'halo_vy':halo_catalog_columns['halo_vel'][:,1], 'halo_vz':halo_catalog_columns['halo_vel'][:,2]}
        halo_catalog_columns.update(float3_columns)
        del halo_catalog_columns['halo_pos']
        del halo_catalog_columns['halo_vel']
        
        # rename halo_parent_id to halo_upid, or create it if it does not exist
        try:
            halo_catalog_columns['halo_upid'] = halo_catalog_columns['halo_parent_id']
            del halo_catalog_columns['halo_parent_id']
        except KeyError:
            halo_catalog_columns['halo_upid'] = np.full_like(halo_catalog_columns['halo_id'], -1)  # no parent
            
        # if no mass-like variable exists, create it
        # this obviously isn't fool-proof, but it works for our FoF
        if not any(k.startswith('halo_m') for k in halo_catalog_columns):
            # Halotools depends on exact mass definitions in several places.
            # We can convert our FoF linking length to an approximate average overdensity
            # using the scaling relation of More+2011
            mu = lambda x:np.log(1.+x) - x/(1.+x)
            overdens = lambda b,c: 3*0.6529*b**-3*mu(c)*(1. + c)**2/c**2 - 1
            
            linklen = halocat.halo_finder_cfg['linklen']
            conc = 10.  # More+2011 take a fiducial concentration of 10, but this is arbitrary
            mdef = '{:d}m'.format(int(round(overdens(linklen, conc))))
            mass_key = 'halo_m{}'.format(mdef)
            radius_key = 'halo_r{}'.format(mdef)
            halo_catalog_columns[mass_key] = halo_catalog_columns['halo_N']*metadata['particle_mass']
            # Compute a spherical overdensity radius for this mass definition
            halo_catalog_columns[radius_key] = halo_mass_to_halo_radius(halo_catalog_columns[mass_key], halocat.cosmology, halocat.Redshift, mdef)
            
        
        # Prepare the uniform subsample
        if load_ptcl_catalog:
            u_ss = halocat.uniform_subsample
            u_ss['pos'] %= halocat.header.BoxSize  # wrap to [0, Lbox)
            ptcl_catalog_columns = { 'x':u_ss['pos'][:,0],  'y':u_ss['pos'][:,1],  'z':u_ss['pos'][:,2],
                                    'vx':u_ss['vel'][:,0], 'vy':u_ss['vel'][:,1], 'vz':u_ss['vel'][:,2]}
            if 'pid' in u_ss.dtype.names:
                ptcl_catalog_columns.update({'pid':u_ss['pid']})
            ptcl_catalog_columns.update(metadata)  # merge into one dict so we can pass as kwarg
            
            halo_catalog_columns['user_supplied_ptclcat'] = UserSuppliedPtclCatalog(**ptcl_catalog_columns)
        
        halo_catalog_columns.update(metadata)  # merge into one dict so we can pass as kwarg
        super(AbacusHaloCatalog, self).__init__(**halo_catalog_columns)
        halotools.utils.add_halo_hostid(self.halo_table)
        
        # Manually bind the halo subsamples
        if load_halo_ptcl_catalog:
            subsamples = halocat.subsamples
            subsamples['pos'] %= halocat.header.BoxSize  # wrap to [0, Lbox)
            ptcl_catalog_columns = { 'x':subsamples['pos'][:,0],  'y':subsamples['pos'][:,1],  'z':subsamples['pos'][:,2],
                                    'vx':subsamples['vel'][:,0], 'vy':subsamples['vel'][:,1], 'vz':subsamples['vel'][:,2]}
            if 'id' in subsamples.dtype.names:
                ptcl_catalog_columns.update({'pid':subsamples['id']})
            if 'pid' in subsamples.dtype.names:
                ptcl_catalog_columns.update({'pid':subsamples['pid']})
            ptcl_catalog_columns.update(metadata)  # merge into one dict so we can pass as kwarg
            
            self.halo_ptcl_table = UserSuppliedPtclCatalog(**ptcl_catalog_columns).ptcl_table
