import sys
import os
import shutil
import warnings



class Pele_env_Builder(object):
		"""
			Base class wher the needed pele environment
			is build by creating folders and files
		"""

		def __init__(self, input, folders, files, forcefield, template, rotamers, pele_dir):
			self.input = input
                        self.folders = folders
                        self.files = files
			self.forcefield = forcefield
			self.template = template
			self.rotamers = rotamers
			self.pele_dir = pele_dir
			self.templates = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__)), "PeleTemplates"))

		def folder_levels(self):
			"""
				Create pele folders
			"""
			
			for folder in self.folders:
				self.create_dir(self.pele_dir, folder)
			

		def file_dist(self):
			"""
				Copy control rotamer
				and template files
			"""
			#paths
			self.template_dir = os.path.join(
				self.pele_dir, "DataLocal/Templates/{}/HeteroAtoms/".format(self.forcefield))
			self.rotamers_dir = os.path.join(
				self.pele_dir, "DataLocal/LigandRotamerLibs")

			#actions
                        for file in self.files:
                            self.copy(file, self.pele_dir)
			shutil.move(self.template, self.template_dir)
			shutil.move(self.rotamers, self.rotamers_dir)





		def create_dir(self, base_dir, extension=None):
			"""
				Class Method to manage
				directory creation only if that
				ones doesn't exist

				Location:
					base_dir+extension
					or base_dir if extension is None
			"""
			if extension:				
				path = os.path.join(base_dir, extension)
 				if os.path.isdir(path):
 					warnings.warn("Directory {} already exists.".format(path),RuntimeWarning)
 				else:
					os.makedirs(path)
			else:
				if os.path.isdir(base_dir):
					warnings.warn("Directory {} already exists.".format(base_dir), RuntimeWarning)
				else:
					os.makedirs(base_dir)

		def copy(self, standard, destination, user=None):
			if user:
				shutil.copy(user, os.path.join(self.pele_dir, standard))
			else:
				shutil.copy(standard, self.pele_dir)
                        return os.path.join(self.pele_dir,standard)

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def set_pele_env(system,  folders, files, forcefield, template, rotamers_file, pele_dir):
	pele_env = Pele_env_Builder(system, folders, files,  forcefield, template, rotamers_file, pele_dir)
	pele_env.folder_levels()
	pele_env.file_dist()