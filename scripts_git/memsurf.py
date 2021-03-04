import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import resource

import MDAnalysis as mda
import MDAnalysis.transformations as trans
from MDAnalysis.analysis.leaflet import LeafletFinder

import memsurfer
from memsurfer import *
memsurfer.utils.create_logger(3,1,0,'','')

sigma = 20             # standard deviation?
rho_type = 3           # 1 = geodesic, 2 = 2D, 3 = 3D
flipflop_offset = 11   # Distance to membrane to determine flipflopping
knbrs = 36             # Amount of neighbours for surface pnormals
boundary_layer = 0.5   # gives boundary layer for membranes

outputname = sys.argv[1]

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.float32):
            return str(obj)
        return json.JSONEncoder.default(self, obj)

class mem_analysis:
    """
    Class for analysis of membranes from .data and .xtc files
    """

    def __init__(self, prefix, pre_out=None):
        """
        Reads .data and .xtc into universe
        Defines static atomgroups based on type
        Reads parameter values (nsteps, nlipids)
        Initiates bbox and leaflets based on first frame
        Initiates dictionary out for storage of data
        Writes lipid type to out['lip']
        """
        if not pre_out:
            self.pre_out = prefix
        else:
            self.pre_out = pre_out
        
        self.u = mda.Universe(prefix + '.data', prefix + '.xtc')
        self.ref = mda.Universe(prefix + '.data')
        self.nsteps = len(self.u.trajectory)
        self.ts = 0

        b = self.u.trajectory[0].dimensions[0:3]
        self.bbox = np.array([0, 0, 0, b[0], b[1], b[2]]).reshape(2,3)
        
        # Define atomgroups
        self.atomgroups = {}
        self.atomgroups['nc3'] = self.u.select_atoms('type 1', updating=True)
        self.atomgroups['dopc'] = self.u.select_atoms('type 2', updating=True)
        self.atomgroups['dogl'] = self.u.select_atoms('type 9', updating=True)
        self.atomgroups['tails'] = self.u.select_atoms('type 7', updating=True)
        self.atomgroups['heads'] = self.u.select_atoms('type 2 9', updating=True)
        self.dnlipids = len(self.atomgroups['heads'])

        self.membranes = {}

        # Set leaflets
        self.leaflets = LeafletFinder(self.u, self.atomgroups['dopc'], cutoff=15, pbc=True)
        sizes = np.array(list(self.leaflets.sizes().values()))
        l = [self.leaflets.groups(i) for i in np.argsort(-sizes)[:2]]
        if l[0].center_of_geometry()[2] > l[1].center_of_geometry()[2]:
            self.atomgroups['ul_dopc'], self.atomgroups['ll_dopc'] = l[0], l[1]
        else:
            self.atomgroups['ul_dopc'], self.atomgroups['ll_dopc'] = l[1], l[0]

        # Make out dict
        self.out = {'lip':[],'r':[], 's': [], 'j': [], 'k': [], 'd': [], 'rho': [], 'l': [], 'apl': [], 'd_nc3':[]}
        self.out_box = []

        for key in self.out.keys():
            self.out[key] = [[] for i in range(len(self.atomgroups['heads']))]
            
        # Add lipid types to out['lip']
        for head in self.atomgroups['heads']:
            if head in self.atomgroups['dopc']:
                self.out['lip'][head.resid - 1] = 'dopc'
            else:
                self.out['lip'][head.resid - 1] = 'dogl'

    def update_membranes(self):
        """
        Computes all membranes 
        Updates dynamic atomgroups
        """

        self.membranes['ul_dopc'] = Membrane(self.atomgroups['ul_dopc'].atoms.positions, labels=self.atomgroups['ul_dopc'].atoms.resids, periodic=True, bbox=self.bbox, boundary_layer=0.5)
        self.membranes['ul_dopc'].fit_points_to_box_xy()
        self.membranes['ul_dopc'].compute_pnormals(knbrs=36)
        self.membranes['ul_dopc'].compute_approx_surface(exactness_level=8)
        self.membranes['ul_dopc'].compute_membrane_surface()
        self.membranes['ul_dopc'].compute_properties()

        self.membranes['ll_dopc'] = Membrane(self.atomgroups['ll_dopc'].atoms.positions, labels=self.atomgroups['ll_dopc'].atoms.resids, periodic=True, bbox=self.bbox, boundary_layer=0.5)
        self.membranes['ll_dopc'].fit_points_to_box_xy()
        self.membranes['ll_dopc'].compute_pnormals(knbrs=36)
        self.membranes['ll_dopc'].compute_approx_surface(exactness_level=8)
        self.membranes['ll_dopc'].compute_membrane_surface()
        self.membranes['ll_dopc'].compute_properties()
        
        Membrane.compute_thickness(self.membranes['ul_dopc'], self.membranes['ll_dopc'])
        
        # Generate NC3 membranes
        self.atomgroups['ul_nc3'] = self.u.select_atoms('resid ' + ' '.join(map(str, self.atomgroups['ul_dopc'].atoms.resids))).intersection(self.atomgroups['nc3'])
        self.atomgroups['ll_nc3'] = self.u.select_atoms('resid ' + ' '.join(map(str, self.atomgroups['ll_dopc'].atoms.resids))).intersection(self.atomgroups['nc3'])

        self.membranes['ul_nc3'] = Membrane(self.atomgroups['ul_nc3'].atoms.positions, labels=self.atomgroups['ul_nc3'].atoms.resids, periodic=True, bbox=self.bbox, boundary_layer=0.5)
        self.membranes['ul_nc3'].fit_points_to_box_xy()
        self.membranes['ul_nc3'].compute_pnormals(knbrs=36)
        self.membranes['ul_nc3'].compute_approx_surface(exactness_level=8)
        self.membranes['ul_nc3'].compute_membrane_surface()
        self.membranes['ul_nc3'].compute_properties()

        self.membranes['ll_nc3'] = Membrane(self.atomgroups['ll_nc3'].atoms.positions, labels=self.atomgroups['ll_nc3'].atoms.resids, periodic=True, bbox=self.bbox, boundary_layer=0.5)
        self.membranes['ll_nc3'].fit_points_to_box_xy()
        self.membranes['ll_nc3'].compute_pnormals(knbrs=36)
        self.membranes['ll_nc3'].compute_approx_surface(exactness_level=8)
        self.membranes['ll_nc3'].compute_membrane_surface()
        self.membranes['ll_nc3'].compute_properties()

        # Compute all distances from dogl to closest membrane and determine upper / lower / flipflop
        self.membranes['dogl'] = Membrane.compute(self.atomgroups['dogl'].atoms.positions, self.atomgroups['dogl'].atoms.resids, periodic=True, bbox=self.bbox)
        
        Membrane.compute_thickness(self.membranes['dogl'], self.membranes['ul_nc3'], mtype='exact')
        self.ul_dogl_dist = np.column_stack((self.membranes['dogl'].labels, self.membranes['dogl'].properties['thickness']))
        ul_dogl_labels = np.array([int(lip[0]) for lip in self.ul_dogl_dist if lip[1] < flipflop_offset])
        if len(ul_dogl_labels) > 0:
            self.atomgroups['ul_dogl'] = self.u.select_atoms('resid ' + ' '.join(map(str, ul_dogl_labels))).intersection(self.atomgroups['dogl'])
        else:
            self.atomgroups['ul_dogl'] = self.u.atoms[[]]
        
        Membrane.compute_thickness(self.membranes['dogl'], self.membranes['ll_nc3'], mtype='exact')
        self.ll_dogl_dist = np.column_stack((self.membranes['dogl'].labels, self.membranes['dogl'].properties['thickness']))
        ll_dogl_labels = np.array([int(lip[0]) for lip in self.ll_dogl_dist if lip[1] < flipflop_offset])
        if len(ll_dogl_labels) > 0:
            self.atomgroups['ll_dogl'] = self.u.select_atoms('resid ' + ' '.join(map(str, ll_dogl_labels))).intersection(self.atomgroups['dogl'])
        else:
            self.atomgroups['ll_dogl'] = self.u.atoms[[]]

        # Generate ff_dogl membrane from lipids not in ul and ll
        ff_dogl_labels = np.array([lip for lip in self.membranes['dogl'].labels if lip not in np.union1d(ul_dogl_labels, ll_dogl_labels)])
        if len(ff_dogl_labels) > 0:
            self.atomgroups['ff_dogl'] = self.u.select_atoms('resid ' + ' '.join(map(str, ff_dogl_labels))).intersection(self.atomgroups['dogl'])
        else:
            self.atomgroups['ff_dogl'] = self.u.atoms[[]]

        # Generate combines membranes         
        self.atomgroups['ul'] = self.atomgroups['ul_dopc'] + self.atomgroups['ul_dogl']
        self.atomgroups['ll'] = self.atomgroups['ll_dopc'] + self.atomgroups['ll_dogl']

        self.membranes['ul'] = Membrane(self.atomgroups['ul'].atoms.positions, labels=self.atomgroups['ul'].atoms.resids, periodic=True, bbox=self.bbox, boundary_layer=0.5)
        self.membranes['ul'].fit_points_to_box_xy()
        self.membranes['ul'].compute_pnormals(knbrs=36)
        self.membranes['ul'].compute_approx_surface(exactness_level=8)
        self.membranes['ul'].compute_membrane_surface()
        self.membranes['ul'].compute_properties()
        
        self.membranes['ll'] = Membrane(self.atomgroups['ll'].atoms.positions, labels=self.atomgroups['ll'].atoms.resids, periodic=True, bbox=self.bbox, boundary_layer=0.5)
        self.membranes['ll'].fit_points_to_box_xy()
        self.membranes['ll'].compute_pnormals(knbrs=36)
        self.membranes['ll'].compute_approx_surface(exactness_level=8)
        self.membranes['ll'].compute_membrane_surface()
        self.membranes['ll'].compute_properties()
        
        # Compute thickness and densities
        self.membranes['ul'].compute_density(rho_type, sigma, 'density', True, self.membranes['ul'].labels)
        self.membranes['ll'].compute_density(rho_type, sigma, 'density', True, self.membranes['ll'].labels)
        Membrane.compute_thickness(self.membranes['ul_nc3'], self.membranes['ll_nc3'])

        for key in self.membranes:
            self.membranes[key].write_all(self.pre_out + '/' + key + '_' + str(self.ts))

    def write_to_out(self):
    		"""
    		Save all lipids to self.out. Switch case for lipid type and u/l/f
    		"""
        for head in self.atomgroups['heads']:
            if head in self.atomgroups['ul_dopc']:
                i_ul_dopc = np.where(self.membranes['ul_dopc'].labels==head.resid)[0][0]
                i_ul_nc3 = np.where(self.membranes['ul_nc3'].labels==head.resid)[0][0]
                i_ul = np.where(self.membranes['ul'].labels==head.resid)[0][0]
                
                self.out['l'][head.resid - 1].append('u')
                self.out['r'][head.resid - 1].append(head.position)

                self.out['d'][head.resid - 1].append('na')
                self.out['d_nc3'][head.resid - 1].append(self.membranes['ul_nc3'].properties['thickness'][i_ul_nc3])

                self.out['rho'][head.resid - 1].append(self.membranes['ul'].properties['density'][i_ul])
                self.out['j'][head.resid - 1].append(self.membranes['ul'].memb_smooth.mean_curv[i_ul])
                self.out['k'][head.resid - 1].append(self.membranes['ul'].memb_smooth.gaus_curv[i_ul])
                self.out['apl'][head.resid - 1].append(self.membranes['ul'].memb_smooth.pareas[i_ul])

            elif head in self.atomgroups['ll_dopc']:
                i_ll_dopc = np.where(self.membranes['ll_dopc'].labels==head.resid)[0][0]
                i_ll_nc3 = np.where(self.membranes['ll_nc3'].labels==head.resid)[0][0]
                i_ll = np.where(self.membranes['ll'].labels==head.resid)[0][0]
                
                self.out['l'][head.resid - 1].append('l')
                self.out['r'][head.resid - 1].append(head.position)

                self.out['d'][head.resid - 1].append('na')
                self.out['d_nc3'][head.resid - 1].append(self.membranes['ll_nc3'].properties['thickness'][i_ll_nc3])

                self.out['rho'][head.resid - 1].append(self.membranes['ll'].properties['density'][i_ll])
                self.out['j'][head.resid - 1].append(self.membranes['ll'].memb_smooth.mean_curv[i_ll])
                self.out['k'][head.resid - 1].append(self.membranes['ll'].memb_smooth.gaus_curv[i_ll])
                self.out['apl'][head.resid - 1].append(self.membranes['ll'].memb_smooth.pareas[i_ll])

            elif head in self.atomgroups['ul_dogl']:
                i_ul = np.where(self.membranes['ul'].labels==head.resid)[0][0]
                
                self.out['l'][head.resid - 1].append('u')
                self.out['r'][head.resid - 1].append(head.position)

                d_ul = self.ul_dogl_dist[self.ul_dogl_dist[:,0]==head.resid, 1]
                d_ll = self.ll_dogl_dist[self.ll_dogl_dist[:,0]==head.resid, 1]
                self.out['d'][head.resid - 1].append([d_ul, d_ll])
                
                self.out['d_nc3'][head.resid - 1].append('na')

                self.out['rho'][head.resid - 1].append(self.membranes['ul'].properties['density'][i_ul])
                self.out['j'][head.resid - 1].append(self.membranes['ul'].memb_smooth.mean_curv[i_ul])
                self.out['k'][head.resid - 1].append(self.membranes['ul'].memb_smooth.gaus_curv[i_ul])
                self.out['apl'][head.resid - 1].append(self.membranes['ul'].memb_smooth.pareas[i_ul])

            elif head in self.atomgroups['ll_dogl']:
                i_ll = np.where(self.membranes['ll'].labels==head.resid)[0][0]
                
                self.out['l'][head.resid - 1].append('l')
                self.out['r'][head.resid - 1].append(head.position)

                d_ul = self.ul_dogl_dist[self.ul_dogl_dist[:,0]==head.resid, 1]
                d_ll = self.ll_dogl_dist[self.ll_dogl_dist[:,0]==head.resid, 1]
                self.out['d'][head.resid - 1].append([d_ul, d_ll])
                self.out['d_nc3'][head.resid - 1].append('na')

                self.out['rho'][head.resid - 1].append(self.membranes['ll'].properties['density'][i_ll])
                self.out['j'][head.resid - 1].append(self.membranes['ll'].memb_smooth.mean_curv[i_ll])
                self.out['k'][head.resid - 1].append(self.membranes['ll'].memb_smooth.gaus_curv[i_ll])
                self.out['apl'][head.resid - 1].append(self.membranes['ll'].memb_smooth.pareas[i_ll])

            elif head in self.atomgroups['ff_dogl']:
                self.out['l'][head.resid - 1].append('f')
                self.out['r'][head.resid - 1].append(head.position)

                d_ul = self.ul_dogl_dist[self.ul_dogl_dist[:,0]==head.resid, 1]
                d_ll = self.ll_dogl_dist[self.ll_dogl_dist[:,0]==head.resid, 1]
                self.out['d'][head.resid - 1].append([d_ul, d_ll])

                self.out['d_nc3'][head.resid - 1].append('na')

                self.out['rho'][head.resid - 1].append('na')
                self.out['j'][head.resid - 1].append('na')
                self.out['k'][head.resid - 1].append('na')
                self.out['apl'][head.resid - 1].append('na')

    def write_to_json(self):
        with open(self.pre_out + '.json', 'w') as f:
            json.dump(self.out, f, cls=NumpyEncoder)

    def run(self, trange=None):
        """
        Loops over trajectory
        Update ul_dopc and ll_dopc every n steps (dopc's generally don't flipflop)
        Call update_membranes to recalculate membranes
        Call write_out to add properties to self.out
        """
        if not trange:
            trange = (0, self.nsteps)
        elif len(trange) == 1:
            trange = (trange[0], self.nsteps)

        print(trange)
        print(len(self.u.trajectory))
        
        for ts in self.u.trajectory[trange[0] : trange[1] : 1]:
            print(ts.frame)
            self.ts = ts.frame
            b = ts.dimensions[0:3]
            self.bbox = np.array([0,0,0, b[0], b[1], b[2]]).reshape(2,3)
            self.out_box.append(self.bbox)

            self.update_membranes()
            self.write_to_out()

        self.write_to_json()

sim0 = mem_analysis(sys.argv[1], pre_out=outputname)
sim0.run()
