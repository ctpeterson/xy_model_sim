""" Import modules """
from sim_imports import *

""" XY measurement class """
class XYMeas:
    """Measure observables for XY model

    Measure observables for XY model

    Attributes:
    """
    def __init__(self):
        """ Give user information """
        # Print outout to user
        print('Measurement class initialized.')
        
        # Return nothing
        return None

    def meas_loc_energ(self, site):
        """Measure dot product with neighbors
        
        Measure dot product with neighbors

        Attributes:
        """

        # Get next x and y ind
        next_X = self.lattice[site].next_X; next_Y = self.lattice[site].next_Y;
        
        # Get dot with neighbor in positive x
        dot_x = np.dot(self.lattice[site].lat_vec, self.lattice[next_X].lat_vec)
        
        # Get dot with neighbor in positive y
        dot_y = np.dot(self.lattice[site].lat_vec, self.lattice[next_Y].lat_vec)
        
        # Return local energy
        return 2. - dot_x - dot_y
    
    def energ(self):
        """Measure energy

        Measure energy for this configuration

        Attributes:
        """

        # Calculate energy
        self.energy = self.J * sum(self.meas_loc_energ(site) for site in range(self.size))
        
        # Return Nothing
        return None

    def magn(self):
        """Measure mean squared magnetization

        Measure mean suared magnetization

        Attributes:
        """

        # Calculate mean squared magnetization
        self.mag = sum(obj.lat_vec for obj in self.lattice) / self.size
        
        # Return nothing
        return None
