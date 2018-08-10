from blues.moves import Move
import numpy as np

class RandomLigandTranslationMove(Move):
    """Move that provides methods for calculating properties on the
    object 'model' (i.e ligand) being perturbed in the NCMC simulation.
    Current methods calculate the object's atomic masses and center of masss.
    Calculating the object's center of mass will get the positions
    and total mass.
    Ex.
        from blues.move import RandomLigandTranslationMove
        ligand = RandomLigandTranslationMove(structure, 'LIG')
        ligand.resname : string specifying the residue name of the ligand
        ligand.calculateProperties()
    Attributes
    ----------
    ligand.topology : openmm.topology of ligand
    ligand.atom_indices : list of atom indicies of the ligand.
    #Dynamic attributes that must be updated with each iteration
    ligand.positions : np.array of ligands positions
    """

    def __init__(self, structure, resname='LIG', step=0.1, random_state=None):
        """Initialize the model.
        Parameters
        ----------
        resname : str
             String specifying the resiue name of the ligand.
        step: float
             Distance to translate in a given direction for each move in angstroms.
        structure: parmed.Structure
             ParmEd Structure object of the relevant system to be moved.
        random_state : integer or numpy.RandomState, optional
             The generator used for random numbers. If an integer is given, it fixes the seed. Defaults to the global numpy random number generator.
        """
        self.structure = structure
        self.resname = resname
        self.step = step
        self.random_state = random_state
        self.atom_indices = self.getAtomIndices(structure, self.resname)
        self.topology = structure[self.atom_indices].topology
        self.positions = structure[self.atom_indices].positions

    def getAtomIndices(self, structure, resname):
        """
        Get atom indices of a ligand from ParmEd Structure.
        Arguments
        ---------
        resname : str
             String specifying the resiue name of the ligand.
        structure: parmed.Structure
             ParmEd Structure object of the atoms to be moved.
        Returns
        -------
        atom_indices : list of ints
             list of atoms in the coordinate file matching lig_resname
        """
#        TODO: Add option for resnum to better select residue names
        atom_indices = []
        topology = structure.topology
        for atom in topology.atoms():
             if str(resname) in atom.residue.name:
                  atom_indices.append(atom.index)
        return atom_indices

    def _random_three_vector(self):
        """
        Generates a random 3D unit vector (direction) with a uniform spherical distribution
        Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution

        Adapted from Gist here: https://gist.github.com/andrewbolster/10274979,
        Credit Andrew Bolster
        """
        phi = np.random.uniform(0,np.pi*2)
        costheta = np.random.uniform(-1,1)

        theta = np.arccos(costheta)
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        return np.array([x,y,z])

    def move(self, context):
        """Function that performs a random translation to the ligand.

        Parameters
        ----------
        context: simtk.openmm.Context object
             Context containing the positions to be moved.
        Returns
        -------
        context: simtk.openmm.Context object
             The same input context, but whose positions were changed by this function.
        """
        positions = context.getState(getPositions=True).getPositions(asNumpy=True)
        self.positions = positions[self.atom_indices]

        vec = self._random_three_vector()

        """
        self.center_of_mass = self.getCenterOfMass(self.positions, self.masses)
        reduced_pos = self.positions - self.center_of_mass

        # Define random rotational move on the ligand
        rand_quat = mdtraj.utils.uniform_quaternion(size=None, random_state=self.random_state)
        rand_rotation_matrix = mdtraj.utils.rotation_matrix_from_quaternion(rand_quat)
        #multiply lig coordinates by rot matrix and add back COM translation from origin
        rot_move = np.dot(reduced_pos, rand_rotation_matrix) * positions.unit + self.center_of_mass
        """

        positions = positions * positions.unit + vec*self.step

        context.setPositions(positions)
        positions = context.getState(getPositions=True).getPositions(asNumpy=True)
        self.positions = positions[self.atom_indices]
        return context