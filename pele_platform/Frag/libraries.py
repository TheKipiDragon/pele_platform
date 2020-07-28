from rdkit import Chem
import yaml, glob, argparse, os
from string import Template

#InBond="C6-H7"
InBond="$BOND"

def growing_sites(fragment):
	"""
	Function that obtains the possible starting points to grow a ligand. 
	Entrance: File to obtain the possible bonds of a protein to commence the growth
	Exit: List of strings which represent possible starting bonds
        """
	#breakpoint()
	bonds = []
	if (fragment[-3:] == "sdf"):
		moleculas = Chem.SDMolSupplier(fragment, removeHs=False)
		for mol in moleculas:
			for a in mol.GetAtoms():
				if (a is None or a.GetSymbol() == 'H'): continue
				#print(a.GetSymbol())
				#[print(neigh.GetAtomicNum()) for neigh in a.GetNeighbors()]
				#breakpoint()
				for neigh in a.GetNeighbors():
					if (neigh.GetAtomicNum() == 1): 
						ind1 = int(a.GetIdx()) + 1
						ind2 = int(neigh.GetIdx()) + 1 
						#print(a.GetSymbol() + str(ind1) + '-' + neigh.GetSymbol() + str(ind2))
						bonds.append(a.GetSymbol() + str(ind1) + '-' + neigh.GetSymbol() + str(ind2))
	else:
		mol = Chem.MolFromPDBFile(fragment, removeHs=False)
		for a in mol.GetAtoms():
				if (a is None or a.GetSymbol() == 'H'): continue
				#print(a.GetSymbol())
				#[print(neigh.GetAtomicNum()) for neigh in a.GetNeighbors()]
				#breakpoint()
				for neigh in a.GetNeighbors():
					if (neigh.GetAtomicNum() == 1): 
						ind1 = int(a.GetIdx()) + 1
						ind2 = int(neigh.GetIdx()) + 1 
						#print(a.GetSymbol() + str(ind1) + '-' + neigh.GetSymbol() + str(ind2))
						bonds.append(a.GetSymbol() + str(ind1) + '-' + neigh.GetSymbol() + str(ind2))
	return bonds
	

def input_file(fragment_pdb, core_atom, growing_site, output="input.conf", createFile=True):
	"""
	Function that creates the pele.conf file to for the frag_pele module
	Entrance: 
		-Fragment_pdb: The file of the protein to be grown.
		-Core_atom: The starting bond of the fragment to grow.
		-Growing_site: list of bonds in which the growth will begin
		-Output: file to save the configuration.
		-CreateFile: Parameter used to create a new pele.conf file, or append into an existing one. 
        """
	if createFile:
		with open(output, "w+") as file: output = write_File(file, False, fragment_pdb, core_atom, growing_site, output)
	else:
		with open(output, "a+") as file: output = write_File(file, True, fragment_pdb, core_atom, growing_site, output)
	return output

def write_File(file, isAppending, fragment_pdb, core_atom, growing_site, output):
	StrForFile = ""
	for bond in growing_site:
		if (isAppending):
			#StrForFile += ("\n" + fragment_pdb + " " + core_atom + " " + bond + "\n")
			StrForFile += ("\n" + fragment_pdb + " $BOND " + bond + "\n")
			isAppending = False
		else:
			#StrForFile += (fragment_pdb + " " + core_atom + " " + bond + "\n")
			StrForFile += (fragment_pdb + " $BOND " + bond + "\n")
	s = Template(StrForFile)
	s = s.safe_substitute(BOND=core_atom)
	file.write(s[:-1])
	return file

def parse_Args():
	parser = argparse.ArgumentParser()
	parser.add_argument("--user_bond", help="Bond to start the growth",  type=str, required=True)
	#parser.add_argument("--library", help="Indicates if a folder with pdb files will be used, if null, then the .yaml config will be used",  type=str)
	#parser.add_argument("--UseLibrary", help="Indicates if a folder with pdb files will be used",  type=bool, required=True)
	args = parser.parse_args()
	return args.user_bond, args.UseLibrary

#UseLibrary
def main(bondUser, frag_library):
	
	path = os.path.join(frag_library, "Core,ligand*.pdb")
	#bondUser, library = parse_Args()
	#dc = yaml.safe_load(open("input.yaml","r"))
		#print (path)
	firstTime=True
	for file in glob.glob(path):
		listaBonds = growing_sites(file)
		if (firstTime): 
			output = input_file(file, bondUser, listaBonds)
			firstTime=False 
		else: output = input_file(file, bondUser, listaBonds, createFile=False)
	"""	
	else:	
		pdbFile = dc["frag_core"]
		sdfFile = dc["frag_ligands"]
		listaBonds = growing_sites(sdfFile)
		input_file(pdbFile, InBond, listaBonds)
	"""
	return output.name	#Returns the input.conf file
	
	
if __name__ == "__main__":
	dc = yaml.safe_load(open("input.yaml","r"))
	frag_library = dc["frag_library"]
	bondUser = dc["frag_library_core"]
	#bondUser, UseLibrary = parse_Args()
	main(bondUser, frag_library)

