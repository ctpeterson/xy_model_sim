""" External modules """
# Numpy for number crunching operations
import numpy as np

""" Cluster class """
class XYCluster:
    """Cluster update class.
    
    Class implements cluster update according to algorithm
    in Phys. Rev. Lett. 62, 361 ("Collective Monte Carlo 
    Updating for Spin Systems", Ulli Wolf)

    Attributes:
    """
    def __init__(self):
        """ Give user information """
        # Print outout to user
        print('Cluster update class initialized.')
        
        # Return nothing
        return None

    """ Private methods """
    def __grow_cluster(self, site):
        """Grows cluster

        Recursively grows cluster. Recursion implemented as "depth-first"
        according to pseudocode found in lecture series by Kari Rummukainen: 
        -  https://www.mv.helsinki.fi/home/rummukai/simu/

        Attributes:
           lat_vec_site (array): Lattice vector at this site
           lat_vec_nhbr (array): Lattice vector at neighboring site
           dt_prd_st_rfl (float): Dot product of current site with rfl. dir.
           dt_prd_nghbr_rfl (flat): Dot product of neighbor site with rfl. dir.
           spin_prod (float): Product of dot products
           min_of_change (float): Minimum of change of energy and zero
           prob_acc (float): Probability of reflecting neighboring spine
           refl_spin (array): Reflected spin of neighbor           
        """

        """ Define a few convient quantities """
        # Temporarily store lattice vector of site
        lat_vec_site = self.lattice[site].lat_vec

        # Temporarily store dot product
        dt_prd_st = np.dot(lat_vec_site, self.refl_dir)

        """ Flip value at site """
        # Define reflected spin
        refl_spin = lat_vec_site - 2. * dt_prd_st * self.refl_dir
                
        # Save reflected spin
        self.lattice[site].set_new_lat_vec(refl_spin)

        # Add site to cluster
        self.cluster_sites.append(site)

        """ Add neighbors to cluster """
        # Get neighbors
        nghbrs = self.lattice[site].lat_neighbors

        # Cycle through neighbors
        for nghbr_site in [nghbr for nghbr in nghbrs if nghbr not in self.cluster_sites]:
            """ Calculate a few things """
            # Temporarily store lattice vector of neighbor
            lat_vec_nhbr = self.lattice[nghbr_site].lat_vec

            # Temporarily store dot product
            dt_prd_nghbr = np.dot(lat_vec_nhbr, self.refl_dir)
            
            # Get product of spin dot products
            spin_prod = dt_prd_st * dt_prd_nghbr
            
            # Get minimum of 0 and change in energy
            min_of_change = min(0., -2. * self.J * spin_prod)

            """ Do acc./rej. step """
            # Generate random number between 0 and 1
            rand_num = self._rand_zero_to_one()
            
            # Define probability of acceptance
            prob_acc = 1. - np.exp(min_of_change)

            # Check if change is to be accepted
            if (prob_acc >= rand_num):
                # Grow cluster more
                self.__grow_cluster(nghbr_site)
                
        # Return nothing
        return None

    """ Public methods """
    def cluster_update(self, site = None, save = False,
                       prnt = False, log = False):
        """Implement cluster update
        
        Implement cluster update

        Attributes:
           site (int): Random site to start cluster
           angle (float): Random angle to define reflection
        """
        
        """ Initial tasks """
        # Print separator if requested
        print(25 * '-.' + '\n') if prnt is True else None
        
        # Sample a random site
        site = self._rand_site() if site is None else site
        
        # Sample a random angle
        angle = self._rand_angle()

        """ Get direction of reflection and grow cluster """
        # Initialize list of cluster sites
        self.cluster_sites = []

        # Create lattice vector representing reflection
        self.refl_dir = np.array([np.cos(angle), np.sin(angle)])
        
        # Grow cluster
        self.__grow_cluster(site)

        """ Update lattice and save if necessary """
        # Update configuration number
        self.conf_num += 1
        
        # Check if save is true
        if save is True:
            # Save angles to location
            self._save_conf('cluster', prnt)

        """ Save or print extra info if necessary """
        # Check if prnt is True
        if prnt is True:
            # Print information about configuration
            print('Info about conf. num. ' + str(self.conf_num) + ':')

            # Print information about ensemble name
            print('Ensemble name:', self.ens_name.strip('_'))
            
            # Print information about site and random angle
            print('site, angle, alg = ' + str(site) + ',', 
                  str(angle) + ', cluster')

            # Print cluster size
            print('Cluster size:', len(self.cluster_sites))

            # Create final separator
            print('\n' + 25 * '-.')

        """ Option to log info to be added """
        
        # Return nothing
        return None
