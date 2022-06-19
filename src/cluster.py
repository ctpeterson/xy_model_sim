""" External modules """
# Numpy for number crunching operations
import numpy as np

""" Cluster class """
class XYCluster(object):
    """Do cluster upadate
    
    Do cluster update

    Attributes:
    """
    def __init__(self):
        # Return nothing
        return None

    """ Private methods """
    def __grow_cluster(self, site):
        """Grows cluster

        Recursively grows cluster

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
        dt_prd_st_rfl = np.dot(lat_vec_site, self.refl_dir)

        """ Flip value at site """
        # Define reflected spin
        refl_spin = lat_vec_site - 2. * dt_prd_st_rfl * self.refl_dir
                
        # Save reflected spin
        self.lattice[site].set_new_lat_vec(refl_spin)
        
        """ Add neighbors to cluster """
        # Cycle through neighbors
        for nghbr_site in self.lattice[site].lat_neighbors:
            """ Calculate a few things """
            # Temporarily store lattice vector of neighbor
            lat_vec_nhbr = self.lattice[nghbr_site].lat_vec

            # Temporarily store dot product
            dt_prd_nghbr_rfl = np.dot(lat_vec_nhbr, self.refl_dir)
            
            # Get product of spin dot products
            spin_prod = dt_prd_st_rfl * dt_prd_nghbr_rfl
            
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
    def cluster_update(self, save = False, prnt = False, log = False):
        """Implement cluster update
        
        Implement cluseter update

        Attributes:
           site (int): Random site to start cluster
           angle (float): Random angle to define reflection
        """
        
        """ Initial tasks """
        # Print separator if requested
        print(50 * '-.' + '\n') if prnt is True else None
        
        # Sample a random site
        site = self._rand_site()
        
        # Sample a random angle
        angle = self._rand_angle()

        """ Get direction of reflection and grow cluster """
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

            # Create final separator
            print('\n' + 50 * '-.')

        """ Option to log info to be added """
        
        # Return nothing
        return None
