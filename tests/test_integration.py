# coding=utf-8
"""Integration tests"""

import unittest
import os
import filecmp

import horton
import icecube.grid
import icecube.cube
import iceube.surfaces

class IntegrationTests(unittest.TestCase):
    def test_cube(self):
        
        """
        test = 'PNA_in_CCl4'
        tests_dir = '{0}'.format(os.path.dirname(__file__))
        project = pyframe.Project(work_dir='{0}'.format(tests_dir))
        system = pyframe.MolecularSystem(input_file='{0}/{1}/{1}.pdb'.format(tests_dir, test), bond_threshold=1.15)
        core = system.get_fragments_by_name(names=['PNA'])
        system.set_core_region(core)
        solvent = system.get_fragments_by_name(names=['TET'])
        system.add_region(name='solvent', fragments=solvent, use_standard_potentials=True, standard_potential_model='SEP')
        project.create_embedding_potential(system)
        project.write_core(system)
        self.assertTrue(os.path.isfile('{0}/{1}/{1}.mol'.format(tests_dir, test)))
        self.assertTrue(filecmp.cmp('{0}/{1}/{1}.mol'.format(tests_dir, test), '{0}/{1}/{1}.mol.ref'.format(tests_dir, test)))
        os.remove('{0}/{1}/{1}.mol'.format(tests_dir, test))
        project.write_potential(system)
        self.assertTrue(os.path.isfile('{0}/{1}/{1}.pot'.format(tests_dir, test)))
        self.assertTrue(filecmp.cmp('{0}/{1}/{1}.pot'.format(tests_dir, test), '{0}/{1}/{1}.pot.ref'.format(tests_dir, test)))
        os.remove('{0}/{1}/{1}.pot'.format(tests_dir, test))
        """
