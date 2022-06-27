""" Import modules """
from sim_imports import *

""" Lattice Class """
class XYLattice(object):
    """Class defining XY lattice

    Class defining XY lattice sites and their attributes

    Attributes:
        site (int): Location of lattice site
        value (float): Value of spin at lattice site
        lat_vec (array): Array of vector components
        lat_dirs (array): Array of neighbors (counter-clockwise)
    """
    def __init__(self, angle, site, N):
        # Save value at site
        self.angle = angle

        # Set lattice vector
        self.lat_vec = np.array([np.cos(self.angle), np.sin(self.angle)])
        
        # Save site location
        self.site = site
        
        # Initialize neighbors
        self.__neighbors(N)

        # Return nothing
        return None

    """ Private Methods """
    def __set_x(self, N):
        """Set neighbor in x

        Set neighbor in x-direction

        Attributes:
           N (int): Linear dimension of lattice
        """

        """ Figure out neighbors in positive x-direction """
        # Perform check to impose periodic BC's
        if ( (self.site + 1) % N ) == 0:
            # Set index of next neighbor on period
            self.next_X = self.site - (N - 1)
        else: # Set to next lattice site
            # Set index of next neighbor
            self.next_X = self.site + 1

        """ Figure out neighbors in negative x-direction """
        # Perform check to impose periodic BC's
        if (self.site % N == 0):
            # Set index of neighbor behind on period
            self.last_X = self.site + (N - 1)
        else: # Set to next lattice site
            # Set index of neighbor behind
            self.last_X = self.site - 1
            
        # Return nothing
        return None
    
    def __set_y(self, N):
        """Set neighbor in y

        Set neighbor in y-direction

        Attributes:
           N (int): Linear dimension of lattice
        """

        """ Figure out neighbors in positive y-direction """
        # Perform check to impose periodic BC's
        if (0 <= self.site <= N - 1):
            # Set index of neighbor on period
            self.next_Y = self.site + N*(N - 1)
        else: # Set to next lattice site
            # Set index of neighbor
            self.next_Y = self.site - N

        """ Figure out neighbros in negative y-direction """
        if (N*(N - 1) <= self.site <= N*N - 1):
            # Set index of last neighbor in period
            self.last_Y = self.site - N*(N - 1)
        else: # Set to next lattice site
            # Set index of last neighbor
            self.last_Y = self.site + N
        
        # Return nothing
        return None

    def __neighbors(self, N):
        """Saves index of neighbors

        Saves index of neighbors

        Attributes
            ndim (int): Number of dimensions
            dr_mthds (array): Array of dir. sttng. mthds.
        """

        # Set number of dimensions
        ndim = 2

        # Define direction methods
        dr_mthds = [self.__set_x, self.__set_y]

        # Set directions
        [set_dir(N) for set_dir in dr_mthds]

        # Set array of counter-clockwise directions
        self.lat_neighbors = [self.next_X, self.next_Y,
                              self.last_X, self.last_Y]

        # Return nothing
        return None    
    
    """ Public methods """
    def set_new_lat_vec(self, lat_vec):
        """Sets new lattice vector
        
        Sets new lattice vector, figures 
        out angle

        Attributes:
        """

        # Set lattice vector
        self.lat_vec = lat_vec

        # Set new angle
        self.angle = np.mod(np.arctan(lat_vec[-1]/lat_vec[0]), 2. * np.pi)

        # Return nothing
        return None
