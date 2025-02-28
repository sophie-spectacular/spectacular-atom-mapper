from rdkit import Chem
from preprocess import consistent_mapping


def unmap_smiles(smiles: str) -> str:
    """
    Remove atom mapping from a SMILES string (reaction or molecule).
    
    :param smiles: SMILES string of a molecule or reaction.
    :return: Unmapped SMILES string.
    """
    parts = smiles.split('>>')
    
    unmapped_parts = []
    for part in parts:
        molecules = part.split('.')
        unmapped_molecules = []
        for mol in molecules:
            mol_obj = Chem.MolFromSmiles(mol)
            if mol_obj:
                for atom in mol_obj.GetAtoms():
                    atom.SetAtomMapNum(0)
                unmapped_molecules.append(Chem.MolToSmiles(mol_obj))
            else:
                unmapped_molecules.append(mol)
        unmapped_parts.append('.'.join(unmapped_molecules))
    
    return '>'.join(unmapped_parts)

from rdkit import Chem
from collections import defaultdict

def unmap_reactants(reactants):
    mols = [Chem.MolFromSmiles(smi) for smi in reactants.split('.')]
    for mol in mols:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
    return '.'.join(Chem.MolToSmiles(mol) for mol in mols)

def get_multi(file_path):
    reaction_groups = defaultdict(list)
    rank_map = {}
    
    with open(file_path, "r") as file:
        for line in file:
            reaction, metadata = line.rsplit(" ", 1)
            rank = float(metadata.split(" ")[-1])
            reactants, products = reaction.split(">>")
            unmarked_reactants = unmap_reactants(reactants)
            reaction_groups[unmarked_reactants].append(line.strip())
            rank_map[unmarked_reactants] = rank
    
    sorted_groups = sorted(reaction_groups.items(), key=lambda x: rank_map[x[0]], reverse=True)
    
    return [group for _, group in sorted_groups]


def group_reactions(input_file, output_file):
    sets = get_multi(input_file)
    print(sets)
    print([len(set) for set in sets])
    with open(output_file, 'w') as file1:
        for set in sets:
            amalgamated = consistent_mapping(set)[2]
            file1.write(amalgamated+"\n")
    