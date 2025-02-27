from rdkit import Chem
import re
import pickle
import pandas as pd
from tqdm import tqdm
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') 

def remove_atom_map(mol, isotope=False):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
        if isotope: atom.SetIsotope(0)
    return mol

def numbering(mol, mapping_number):
    numbering_dict = {}
    for index, atom in enumerate(mol.GetAtoms()):
        numbering_dict[atom.GetIntProp('molAtomMapNumber')] = index + 1 + mapping_number
    return numbering_dict, mapping_number+index+1

def get_changed_bonds(rxn_smi):
    reactants = Chem.MolFromSmiles(rxn_smi.split('>')[0])
    products  = Chem.MolFromSmiles(rxn_smi.split('>')[2])

    conserved_maps = [a.GetAtomMapNum() for a in products.GetAtoms() if a.HasProp('molAtomMapNumber')]
    bond_changes = set() # keep track of bond changes

    # Look at changed bonds
    bonds_prev = {}
    for bond in reactants.GetBonds():
        nums = sorted(
            [bond.GetBeginAtom().GetAtomMapNum(),
             bond.GetEndAtom().GetAtomMapNum()])
        if (nums[0] not in conserved_maps) and (nums[1] not in conserved_maps): continue
        bonds_prev['{}~{}'.format(nums[0], nums[1])] = bond.GetBondTypeAsDouble()
    bonds_new = {}
    for bond in products.GetBonds():
        nums = sorted(
            [bond.GetBeginAtom().GetAtomMapNum(),
             bond.GetEndAtom().GetAtomMapNum()])
        bonds_new['{}~{}'.format(nums[0], nums[1])] = bond.GetBondTypeAsDouble()

    for bond in bonds_prev:
        if bond not in bonds_new:
            bond_changes.add((bond.split('~')[0], bond.split('~')[1], 0.0)) # lost bond
        else:
            if bonds_prev[bond] != bonds_new[bond]:
                bond_changes.add((bond.split('~')[0], bond.split('~')[1], bonds_new[bond])) # changed bond
    for bond in bonds_new:
        if bond not in bonds_prev:
            bond_changes.add((bond.split('~')[0], bond.split('~')[1], bonds_new[bond]))  # new bond

    return bond_changes


def check_non_continuous(bond_changes):
    atom_set = set()
    for bond in bond_changes:
        atom_set.add(bond[0])
        atom_set.add(bond[1])
    if len(atom_set) <= len(bond_changes) + 1:  # Allowing 1 disconnection between edits
        return False
    else:
        return True


def check_invalid_bonds(bond_changes):
    for bond in bond_changes:
        if bond[0] == bond[1]:
            return True
    return False


# def numbering(mol, mapping_number):
#     numbering_dict = {}
#     for index, atom in enumerate(mol.GetAtoms()):
#         numbering_dict[atom.GetIntProp('molAtomMapNumber')] = index + 1 + mapping_number
#     return numbering_dict, mapping_number+index+1


def consistent_mapping(reaction_list):
    rxn_dict={}
    mapping_number=0
    product_smi=[]
    reactant_smi=[]
    for idx, rxn_smiles in enumerate(reaction_list):
        split = '.'.join([rxn_smiles.split('>')[0], rxn_smiles.split('>')[1]]) if rxn_smiles.split('>')[1] else rxn_smiles.split('>')[0]
        mols_reactants = Chem.MolFromSmiles(split)
        mols_products=Chem.MolFromSmiles(rxn_smiles.split('>')[-1])
        
        print(numbering(mols_products, mapping_number))

        if idx==0:
            numbering_dict, mapping_number=numbering(mols_reactants, mapping_number)
#             print(numbering_dict)

        else:
            addition_rxn_number=dict()            
            # Avoid the same atom map number. I assume each rxn has different set of mapping.
            for atom in mols_reactants.GetAtoms():
                addition_rxn_number[atom.GetIntProp('molAtomMapNumber')] = atom.GetIntProp('molAtomMapNumber')+mapping_number
                atom.SetIntProp('molAtomMapNumber', atom.GetIntProp('molAtomMapNumber')+mapping_number)
            for atom in mols_products.GetAtoms():
                atom.SetIntProp('molAtomMapNumber', addition_rxn_number[atom.GetIntProp('molAtomMapNumber')])
                
        # Make dictionary for checking shared reactants        
        used_rmol=[]
        for mol_id, mol in enumerate(Chem.GetMolFrags(mols_reactants, asMols=True)):
#             print(idx, mol_id, Chem.MolToSmiles(mol))

            mol_dict = {'mapped_SMILES': Chem.MolToSmiles(mol, canonical=True)  ,
                        'unmapped_SMILES': Chem.MolToSmiles(remove_atom_map(mol)),
                        'used': False,
                        'conversion': numbering_dict}
    
            # Check reactants whether they are used or not
            notused_rmol=[k for k, v in rxn_dict.items() if v['unmapped_SMILES']==mol_dict['unmapped_SMILES'] and not v['used']]
#             used_rmol=[k for k, v in rxn_dict.items() if v['unmapped_SMILES']==mol_dict['unmapped_SMILES'] and v['used']]
            if idx==0:                
                rxn_dict[(idx, mol_id)]=mol_dict
            else:
                if notused_rmol!=[]:
                    rxn_dict[notused_rmol[0]]['used']=True  # It should be turned off after this for loop
                    used_rmol.append(notused_rmol[0])
                    mol_dict['used']=True
#                     print(notused_rmol)
                    present_mol=Chem.MolFromSmiles(rxn_dict[notused_rmol[0]]['mapped_SMILES'])
                    mol=Chem.MolFromSmiles(mol_dict['mapped_SMILES'])

                    present_order_list=list(present_mol.GetSubstructMatches(mol)[0])

                    to_be_change_dict={}
                    for atom in mol.GetAtoms():
                        to_be_change_dict[atom.GetIdx()]=atom.GetAtomMapNum()
                    ref_dict={}
                    for atom in present_mol.GetAtoms():
                        ref_dict[atom.GetIdx()]=atom.GetAtomMapNum()
                    changing_dict={}

                    for k, v in to_be_change_dict.items():
                        changing_dict[v] = ref_dict[present_order_list[k]]
                    mol_dict['mapped_SMILES']=rxn_dict[notused_rmol[0]]['mapped_SMILES']
                    mol_dict['conversion'] = changing_dict
                    rxn_dict[(idx, mol_id)]=mol_dict
                else:
                    mol=Chem.MolFromSmiles(mol_dict['mapped_SMILES'])
                    new_numbering_dict, mapping_number=numbering(mol, mapping_number)
                    numbering_dict.update(new_numbering_dict)

                    for atom in mol.GetAtoms():
                        atom.SetIntProp('molAtomMapNumber', numbering_dict[atom.GetIntProp('molAtomMapNumber')])

                    mol_dict['mapped_SMILES']=Chem.MolToSmiles(mol, canonical=True)
                    mol_dict['conversion']=numbering_dict
                    rxn_dict[(idx, mol_id)]=mol_dict 
        if idx!=0:            
            changing_dict = {}
            for value in rxn_dict.values():
                if value['conversion']:
                    changing_dict.update(value['conversion'])
            for atom in mols_products.GetAtoms():
                atom.SetIntProp('molAtomMapNumber', changing_dict[atom.GetIntProp('molAtomMapNumber')])
        else: pass
        
        product_smi.append(Chem.MolToSmiles(mols_products, isomericSmiles=False))
        
        for used_mol_id in used_rmol:
            rxn_dict[used_mol_id]['used']=False 
#     print(product_smi)
    for value in rxn_dict.values():
        reactant_smi.append(value['mapped_SMILES'])
        
    reactant_smi=list(set(reactant_smi))
    reactant_side='.'.join(reactant_smi)
      
    
    ### Revised 
    rmol_2=Chem.MolFromSmiles(reactant_side)
    numbering_dict2, mapping_number2=numbering(rmol_2, 0)
    reactant_smi2=[]
    product_smi2=[]
    
    rmol2=Chem.MolFromSmiles(reactant_side)
    for atom in rmol2.GetAtoms():
        atom.SetIntProp('molAtomMapNumber', numbering_dict2[atom.GetIntProp('molAtomMapNumber')])
    reactant_side2=Chem.MolToSmiles(rmol2)
    
    for psmi2 in product_smi:
        pmol2=Chem.MolFromSmiles(psmi2)
        for atom in pmol2.GetAtoms():
            atom.SetIntProp('molAtomMapNumber', numbering_dict2[atom.GetIntProp('molAtomMapNumber')])
        product_smi2.append(Chem.MolToSmiles(pmol2))
    
    ## Revised end
#     product_side='.'.join(product_smi)  
#     reactant_side2='.'.join(reactant_smi)
    product_side2='.'.join(product_smi2)  
    
    reaction_smi='>>'.join([reactant_side2, product_side2])
    
    return reactant_side2, product_smi2, reaction_smi



def check_non_continuous(bond_changes):
    atom_set = set()
    for bond in bond_changes:
        atom_set.add(bond[0])
        atom_set.add(bond[1])
    if len(atom_set) <= len(bond_changes) + 1:  # Allowing 1 disconnection between edits
        return False
    else:
        return True


def check_invalid_bonds(bond_changes):
    for bond in bond_changes:
        if bond[0] == bond[1]:
            return True
    return False


def decrease_indices_correctly(chem_list, removed_index):
    updated_chem_list = []
    for chem in chem_list:
        updated_chem = ""
        i = 0
        while i < len(chem):
            if chem[i] == ':' and i + 1 < len(chem) and chem[i+1].isdigit():
                # ':' 발견 후 숫자 시작 지점 찾기
                start = i + 1
                end = start
                while end < len(chem) and chem[end].isdigit():
                    end += 1
                
                # 해당 인덱스 추출 및 조건에 따른 처리
                index = int(chem[start:end])
                if index > removed_index:
                    index -= 1  # 인덱스 감소
                updated_chem += chem[i:start] + str(index)
                i = end - 1
            else:
                updated_chem += chem[i]
            i += 1
        updated_chem_list.append(updated_chem)
    return updated_chem_list

def remove_explicit_hydrogen(chemical_structure):
    split_chemicals = chemical_structure.split('.')
    removed_index = None

    # '[H:'로 시작하는 원소의 인덱스 찾기
    for chem in split_chemicals:
        if chem.startswith('[H'):
            removed_index = int(''.join(filter(str.isdigit, chem)))
            split_chemicals.remove(chem)

    if removed_index is not None:
        return '.'.join(split_chemicals)
    else:
        # '[H:'로 시작하는 원소가 없으면 원본 반환
        return chemical_structure
    
def check_validity(line):

    # 원자들과 그 숫자를 추출하기 위한 정규 표현식 사용
#     rxn = line.split()
    r, p = line.split('>>')
    # 원자와 숫자를 추출
    atoms_with_numbers = re.findall(r'\[([^\]:]+):(\d+)\]', r)
    # 원자 숫자만 추출해서 정수로 변환
    numbers = [int(num) for _, num in atoms_with_numbers]
    # 수소 원자가 있는지 확인
    
    if checkConsecutive(numbers):
        return True
    else: return False
    
def checkConsecutive(l):
    return sorted(l) == list(range(min(l), max(l)+1))


def remapping(rxn):
    r, p = rxn.split('>>')
    mols_reactants=Chem.MolFromSmiles(r)
    numbering_dict, _ =numbering(mols_reactants, 0)
    
    
    rmol2=Chem.MolFromSmiles(r)
    for atom in rmol2.GetAtoms():
        atom.SetIntProp('molAtomMapNumber', numbering_dict[atom.GetIntProp('molAtomMapNumber')])
    reactant_side2=Chem.MolToSmiles(rmol2)
    
    pmol2=Chem.MolFromSmiles(p)
    for atom in pmol2.GetAtoms():
        atom.SetIntProp('molAtomMapNumber', numbering_dict[atom.GetIntProp('molAtomMapNumber')])
    product_smi2=Chem.MolToSmiles(pmol2)
    
    reaction_smi='>>'.join([reactant_side2, product_smi2])
    
    return reaction_smi

# def process_file(fpath, keep_invalid):
#     non_continuous_reactions, five_plus_active_bond_reactions, no_active_bonds, explicit_H  = 0, 0, 0, 0
#     six_active_bond_reactions, seven_active_bond_reactions = 0, 0
#     invalid_bonds, too_many_atoms, invalid_reactions, invalid_numbering = 0, 0, 0, 0
#     valid_reactions = 0
    
#     with open(fpath, 'rb') as f:
#         fid_in=pickle.load(f)
        
#     with open(fpath + '.proc', 'w') as fid_out:
#         for i in tqdm(range(len(fid_in))):
#             retain = True            
#             total_bond_changes=[]
#             reactant = fid_in['reactant'][i]
#             for product in fid_in['product_list'][i]:
#                 new_reactant = remove_explicit_hydrogen(reactant)
#                 new_rxn_smi = '>>'.join([new_reactant, product])
#                 rxn_smi=remapping(new_rxn_smi)
#                 bond_changes = get_changed_bonds(rxn_smi)
#                 if bond_changes == set():
#                     no_active_bonds += 1
#                     retain = False
#                 elif check_invalid_bonds(bond_changes):
#                     invalid_bonds += 1
#                     retain = False
#                 elif not check_validity(rxn_smi):
#                     invalid_numbering += 1
#                     retain = False
#                 elif check_non_continuous(bond_changes):
#                     non_continuous_reactions += 1
#                     retain = False
#                 elif '[H' in rxn_smi.split('>')[0]:
#                     explicit_H += 1
#                     retain = False
#                 elif ':200]' in rxn_smi.split('>')[0]:
#                     too_many_atoms += 1
#                     retain = False
#                 elif len(bond_changes) > 5:
#                     five_plus_active_bond_reactions += 1
#                     if len(bond_changes) == 6:
#                         six_active_bond_reactions += 1
#                     elif len(bond_changes) == 7:
#                         seven_active_bond_reactions += 1
#                     retain = False

#                 if not retain:
#                     invalid_reactions += 1
#                     if not keep_invalid:
#                         continue
#                 else: 
#                     valid_reactions+=1
#                 fid_out.write('{} {}\n'.format(rxn_smi, ';'.join(['{}-{}-{}'.format(x[0], x[1], x[2]) for x in bond_changes])))

#     print(f'Finished processing {fpath}. \n')
#     print(f'In total, there were {invalid_reactions} invalid reactions: \n')
#     print(f'In total, there were {valid_reactions} valid reactions: \n')
#     print(f'* {no_active_bonds} reactions involved no active bond \n')
#     print(f'* {invalid_bonds} reactions involved one or more invalid bonds \n')
#     print(f'* {invalid_numbering} reactions involved a non-continuous atom mapping \n')
#     print(f'* {non_continuous_reactions} reactions involved a non-continuous series of active bonds \n')
#     print(f'* {explicit_H} reactions involved explicit hydrogens \n')
#     print(f'* {too_many_atoms} reactions involved reaction smiles which are too long \n')
#     print(f'* {five_plus_active_bond_reactions} reactions involved more than five active bonds '
#           f'({six_active_bond_reactions} involved 6 and {seven_active_bond_reactions} involved 7 active bonds) \n')
#     if not keep_invalid:
#         print('\nThe invalid reactions have been filtered out. \n')

def process_file(fpath, keep_invalid):
    non_continuous_reactions, five_plus_active_bond_reactions, no_active_bonds, explicit_H  = 0, 0, 0, 0
    six_active_bond_reactions, seven_active_bond_reactions = 0, 0
    invalid_bonds, too_many_atoms, invalid_reactions, invalid_numbering = 0, 0, 0, 0
    valid_reactions = 0
    output = []
    
    turn_off = False
    
    with open(fpath, 'rb') as f:
        fid_in=pickle.load(f)
        
    for i in tqdm(range(len(fid_in))):
        retain = True            
        total_bond_changes=[]
        rxn_label = []
        for rxn_smi in fid_in['consistent_rxnsmi'][i]:
#             try:
#                 new_reactant = remove_explicit_hydrogen(reactant)
#                 new_rxn_smi = '>>'.join([new_reactant, product])
#                 rxn_smi=remapping(new_rxn_smi)
#             except Exception as e:
#                 invalid_numbering += 1
#                 retain = False
#                 continue
            bond_changes = get_changed_bonds(rxn_smi)
            if bond_changes == set():
                no_active_bonds += 1
                retain = False
            elif check_invalid_bonds(bond_changes):
                invalid_bonds += 1
                retain = False
            elif not check_validity(rxn_smi):
                invalid_numbering += 1
                retain = False
            elif check_non_continuous(bond_changes):
                non_continuous_reactions += 1
                retain = False
            elif '[H' in rxn_smi.split('>')[0]:
                explicit_H += 1
                retain = False
            elif ':200]' in rxn_smi.split('>')[0]:
                too_many_atoms += 1
                retain = False
            elif len(bond_changes) > 5:
                five_plus_active_bond_reactions += 1
                if len(bond_changes) == 6:
                    six_active_bond_reactions += 1
                elif len(bond_changes) == 7:
                    seven_active_bond_reactions += 1
                retain = False
            if not retain:
                invalid_reactions += 1
                if not keep_invalid:
                    continue
            else: 
                valid_reactions+=1
            rxn_label.append('{} {}\n'.format(rxn_smi, ';'.join(['{}-{}-{}'.format(x[0], x[1], x[2]) for x in bond_changes])))
        if rxn_label:
            output_format= {'consistent_rxnsmi': fid_in['consistent_rxnsmi'][i],
                'processed_rxnsmi': fid_in['processed_rxnsmi'][i],
                'rxn_class': fid_in['rxn_class'][i],
                'rank': fid_in['rank'][i],
                'solvent': fid_in['solvent'][i],
                'rxnid': fid_in['rxnid'][i],
                'orgn_smiles': fid_in['orgn_smiles'][i],
                'WLN input': rxn_label
                }
            output.append(output_format)
        if turn_off:
            break
    output_df = pd.DataFrame(output)      
    with open(f'{fpath}.proc', 'wb') as handle:
        pickle.dump(output_df, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print(f'Finished processing {fpath}. \n')
    print(f'In total, there were {invalid_reactions} invalid reactions: \n')
    print(f'In total, there were {valid_reactions} valid reactions: \n')
    print(f'* {no_active_bonds} reactions involved no active bond \n')
    print(f'* {invalid_bonds} reactions involved one or more invalid bonds \n')
    print(f'* {invalid_numbering} reactions involved a non-continuous atom mapping \n')
    print(f'* {non_continuous_reactions} reactions involved a non-continuous series of active bonds \n')
    print(f'* {explicit_H} reactions involved explicit hydrogens \n')
    print(f'* {too_many_atoms} reactions involved reaction smiles which are too long \n')
    print(f'* {five_plus_active_bond_reactions} reactions involved more than five active bonds '
          f'({six_active_bond_reactions} involved 6 and {seven_active_bond_reactions} involved 7 active bonds) \n')
    if not keep_invalid:
        print('\nThe invalid reactions have been filtered out. \n')
        
# def process_file_modifiedWLN(fpath, keep_invalid):
#     non_continuous_reactions, five_plus_active_bond_reactions, no_active_bonds, explicit_H  = 0, 0, 0, 0
#     six_active_bond_reactions, seven_active_bond_reactions = 0, 0
#     invalid_bonds, too_many_atoms, invalid_reactions, invalid_numbering = 0, 0, 0, 0
#     valid_reactions = 0
#     output = []
    
#     turn_off = False
    
#     with open(fpath, 'rb') as f:
#         fid_in=pickle.load(f)
        
#     for i in tqdm(range(len(fid_in))):
#         retain = True            
#         total_bond_changes=[]
#         reactant = fid_in['reactant'][i]
#         rxn_label = []
#         product_list = [product for product in fid_in['product_list'][i]]
#         product = '.'.join(product_list)
 
#             try:
#                 new_reactant = remove_explicit_hydrogen(reactant)
#                 new_rxn_smi = '>>'.join([new_reactant, product])
#                 rxn_smi=remapping(new_rxn_smi)
#             except Exception as e:
#                 invalid_numbering += 1
#                 retain = False
#             bond_changes = get_changed_bonds(rxn_smi)
#             if bond_changes == set():
#                 no_active_bonds += 1
#                 retain = False
#             elif check_invalid_bonds(bond_changes):
#                 invalid_bonds += 1
#                 retain = False
#             elif not check_validity(rxn_smi):
#                 invalid_numbering += 1
#                 retain = False
#             elif check_non_continuous(bond_changes):
#                 non_continuous_reactions += 1
#                 retain = False
#             elif '[H' in rxn_smi.split('>')[0]:
#                 explicit_H += 1
#                 retain = False
#             elif ':200]' in rxn_smi.split('>')[0]:
#                 too_many_atoms += 1
#                 retain = False
#             elif len(bond_changes) > 5:
#                 five_plus_active_bond_reactions += 1
#                 if len(bond_changes) == 6:
#                     six_active_bond_reactions += 1
#                 elif len(bond_changes) == 7:
#                     seven_active_bond_reactions += 1
#                 retain = False
#             if not retain:
#                 invalid_reactions += 1
#                 if not keep_invalid:
#                     continue
#             else: 
#                 valid_reactions+=1
#             rxn_label.append('{} {}\n'.format(rxn_smi, ';'.join(['{}-{}-{}'.format(x[0], x[1], x[2]) for x in bond_changes])))
#         if rxn_label:
#             output_format= {'reactant': fid_in['reactant'][i],
#                 'product_list': fid_in['product_list'][i],
#                 'processed_rxnsmi': fid_in['processed_rxnsmi'][i],
#                 'rxn_class': fid_in['rxn_class'][i],
#                 'rank': fid_in['rank'][i],
#                 'solvent': fid_in['solvent'][i],
#                 'rxnid': fid_in['rxnid'][i],
#                 'orgn_smiles': fid_in['orgn_smiles'][i],
#                 'WLN input': rxn_label
#                 }
#             output.append(output_format)
#         if turn_off:
#             break
#     output_df = pd.DataFrame(output)      
#     with open(f'{fpath}.proc', 'wb') as handle:
#         pickle.dump(output_df, handle, protocol=pickle.HIGHEST_PROTOCOL)

#     print(f'Finished processing {fpath}. \n')
#     print(f'In total, there were {invalid_reactions} invalid reactions: \n')
#     print(f'In total, there were {valid_reactions} valid reactions: \n')
#     print(f'* {no_active_bonds} reactions involved no active bond \n')
#     print(f'* {invalid_bonds} reactions involved one or more invalid bonds \n')
#     print(f'* {invalid_numbering} reactions involved a non-continuous atom mapping \n')
#     print(f'* {non_continuous_reactions} reactions involved a non-continuous series of active bonds \n')
#     print(f'* {explicit_H} reactions involved explicit hydrogens \n')
#     print(f'* {too_many_atoms} reactions involved reaction smiles which are too long \n')
#     print(f'* {five_plus_active_bond_reactions} reactions involved more than five active bonds '
#           f'({six_active_bond_reactions} involved 6 and {seven_active_bond_reactions} involved 7 active bonds) \n')
#     if not keep_invalid:
#         print('\nThe invalid reactions have been filtered out. \n')
        
     