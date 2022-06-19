""" Import modules """
from sim_imports import *

""" XY Model Class """
class XYSimulation(XYCluster):
    """Performs XY simulation

    Class implements simulation of Ising model

    Attributes:
        N (int): Linear dimension of lattice
        J (float): Value of spin-spin coupling
    """
    def __init__(self, J, N, config = None, seeds = None,
                 start = None, alg = 'cluster',
                 load_loc = './', save_loc = './'):
        """ Initialize class """
        # Initialize cluster update
        super().__init__()
        
        """ Set initial variables """
         # Set coupling
        self.J = J
        
        # Set linear dimension of lattice
        self.N = N

        # Set size of lattice
        self.size = self.N * self.N
        
        # Set default algorithm
        self.alg = alg
        
        # Define name for ensemble
        self.ens_name = 'xyl' + str(self.N) + 't' + str(self.N)
        self.ens_name += 'J' + str(self.J) + '_'

        # Define place to load configrations/rng states
        self.load_loc = load_loc

        # Define place to save configurations/rng states
        self.save_loc = save_loc

        """ Print out some information """
        # Create first separator
        print('\n' + 50 * '--' + '\n')

        # Tell user ensemble name
        print('Ensemble name:', self.ens_name.strip('_'))

        # Tell user default updating algorithm
        print('Default updates: ' + self.alg)

        # Tell user where configs will be loaded from
        print('Loading dir:', self.load_loc) if config is not None else None

        # Tell user where configurations will be saved to
        print('Saving dir:', self.save_loc)
        
        """ Set lattice and rng up """
        # Check if seeds is None
        if seeds is None:
            # Set default seeds
            self.seeds = {'sites' : 'def_sts', 'angles' : 'def_angs',
                          'probabilities' : 'def_prbs', 'start' : 'def_start'}
        else: # Otherwise, set to user specification
            # Set to specification
            self.seeds = seeds
            
        # Check if start is None
        if start is None:
            # Set to default start
            start = {'start' : 'hot'}
            
        # Set lattice up
        self.__setup_lattice(config, start)

        # Create last separator
        print('\n' + 50 * '--' + '\n')
        
        # Return nothing
        return None

    """ Protected methods """
    def _rand_site(self):
        """Sample lattice sites

        Sample lattice sites

        Attributes:
        """
        # Return lattice site
        return self._site_rng.randint(0, self.size - 1)

    def _rand_angle(self):
        """Sample angle
        
        Sample angle from uniform distribution

        Attributes:
        """
        # Return angle
        return self._angl_rng.uniform(0., 2. * np.pi)

    def _rand_zero_to_one(self):
        """Sample randum number for test
        
        Sample random number between 0 and 1

        Attributes:
        """
        # Return random number between zero and one
        return self._prob_rng.uniform(0., 1.)

    def _init_rngs(self):
        """Initializes RNG random states

        Initializes RNG random states

        Attributes:
        """

        """ Initialize RNG's """
        # Initialize RNG for sites
        self._site_rng = np.random.RandomState(seed = self.seeds['sites'])
        
        # Initialize RNG for angles
        self._angl_rng = np.random.RandomState(seed = self.seeds['angles'])

        # Initialize RNG for MC
        self._prob_rng = np.random.RandomState(seed = self.seeds['probabilities'])

        """ Extra information """
        # Cycle through seed keys
        for key in self.seeds.keys():
            # Print out name of seed
            print('Initializd ' + key + ' state w/',
                  self.seeds[key] + ' string')
        
        # Return nothing
        return None
    
    def _save_conf(self, alg, prnt):
        """Save configuration

        Save configuration

        Attributes:
           alg (str): MC algorithm that generated configuration
        """

        """ General information """
        # Create configuration name
        full_conf_name = self.ens_name + alg + '.' + str(self.conf_num)
        
        """ First, save information about the lattice """
        # Define lattice name
        lat_name = self.save_loc + full_conf_name + '.lat'
        
        # Open file
        with open(lat_name, 'wb+') as out_file:
            # Get angles of configuration
            angles = [site.angle for site in self.lattice]

            # Save angles to file
            pickle.dump(angles, out_file)

            # Tell user what you did
            print('Saved lattice to ' + lat_name) if prnt is True else None

        """ Now, save information about RNG state """
        # Define rng name
        rng_name = self.save_loc + full_conf_name + '.rng'
        
        # Get state of site rng
        site_state = self._site_rng.get_state()

        # Get state of angle rng
        angl_state = self._angl_rng.get_state()

        # Get state of prob rng
        prob_state = self._prob_rng.get_state()

        # Open file to set rng states
        with open(rng_name, 'wb+') as out_file:
            # Put states in numpy array
            state_arr = [site_state, angl_state, prob_state]

            # Save array of states to rng file
            pickle.dump(state_arr, out_file)

            # Tell user what you did
            print('Saved rng states to ' + rng_name) if prnt is True else None
            
        # Return nothing
        return None

    def _get_conf(self):
        """Grabs configuration

        Grabs configuration

        Attributes:
        """

        """ General information """
        # Create configuration name
        full_conf_name = self.ens_name + alg + '.' + str(self.conf_num)
        
        """ Grab lattice configuration """
        # Define lattice name
        lat_name = self.load_loc + full_conf_name + '.lat'
        
        # Open file containing information about lattice
        with open(lat_name, 'rb') as in_file:
            # Get angles
            angles = pickle.load(in_file)
            
            # Reconstruct lattice
            self.lattice = [XYLattice(angle, site, self.N)
                            for site, angle in enumerate(angles)]

            # Tell user what you did
            print('Grabbed lattice from ' + lat_name)

        """ Grab information on RNG state """
        # Define rng name
        rng_name = self.load_loc + full_conf_name + '.rng'

        # Initialize RNG states (a bit redundant)
        self._init_rngs()
        
        # Open file containing information about lattice
        with open(rng_name, 'rb') as in_file:
            # Define state array
            state_arr = pickle.load(in_file)

            # Set state for site rng
            self._site_rng.set_state(state_arr[0])

            # Set state for angle rng
            self._angl_rng.set_state(state_arr[1])

            # Set state for probability rng
            self._prob_rng.set_state(state_arr[-1])

            # Tell user what you did
            print('Grabbed rng states from ' + rng_name)
            
        # Return nothing
        return None

    def _setup_init_lattice(self, start):
        """Setup lattice

        Sets lattice up with particular start

        Attributes:
           start (dict): Dictionary containing information about start
        """

        """ Initialize lattice """
        # Check start
        if start['start'] == 'hot':
            # Create rng state for start
            start_rng = np.random.RandomState(seed = self.seeds['start'])
            
            # Set lattice with hot start
            self.lattice = [XYLattice(start_rng.uniform(0., 2. * np.pi),
                                      site, self.N)
                            for site in range(self.size)]

            # Tell user what you did
            print('Initialized lattice with hot start')
        elif start['start'] == 'cold':
            # Set lattice with cold start
            self.lattice = [XYLattice(start['angle'], site, self.N)
                            for site in range(self.size)]

            # Tell user what you did
            print('Initialized lattice with hot start')
            
        """ Last housekeeping """
        # Initialize RNG states
        self._init_rngs()

        # Set configuration number
        self.conf_num = 0
        
        # Return nothing
        return None

    """ Private methods """
    def __setup_lattice(self, config, start):
        """Set lattice up
        
        Set lattice up. General information:
           - If config is an integer, this method will look for the
             corresponding lattice and rng state and load it up. 
             If config is None, program will default to creating 
             a new lattice with information specified by "start"
             dictionary.
           - The start dictionary contains information about the
             starting configuration. If config['start'] is 'hot', 
             then the lattice will contain random angles
             drawn from 0 to 2*pi. If config['start'] is 'cold',
             then lattice will contain same value at each site,
             and the value give is specified by start['angle'].
 
        Attributes:
           config (int or None): Configuration to start off on
           start (dict): Dictionary containing information about start
        """

        # Check if configuration is none
        if config is None:
            # Setup lattice
            self._setup_init_lattice(start)
        else: # Otherwise, pick up last location
            # Load configuration
            self._get_conf()
        
        # Return nothing
        return None
