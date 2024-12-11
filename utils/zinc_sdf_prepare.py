from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
import pandas as pd
import argparse
import os

# Function to optimize 3D structure
def optimize_mol(mol: Chem.Mol):
    """
    Optimize the 3D structure of a molecule using RDKit's UFF (Universal Force Field).
    """
    # Add hydrogens to the molecule
    mol = Chem.AddHs(mol)
    # Generate 3D coordinates for the molecule
    AllChem.EmbedMolecule(mol)
    # Optimize the structure using UFF (Universal Force Field)
    AllChem.UFFOptimizeMolecule(mol)
    return mol

# Function to convert SMILES to SDF
def convert_smiles_to_sdf(file_path: str, sdf_file:str, add_properties: bool=True):
    """
    Load a file containing SMILES strings and convert them to an SDF file with 3D optimization.

    :param csv_file: input csv file with a column of "SMILES" strings.
    :param output_file: Path to the output SDF file.
    """

    # load file
    smiles_df = pd.read_csv(file_path, sep=" ", header=None, names=["SMILES", "ZINC_ID"])
    #smiles_df = pd.read_csv(file_path)

    # Create an SDF writer to write the molecules to an output file
    writer = rdmolfiles.SDWriter(sdf_file)

    # Iterate over each csv file

    for _, row in smiles_df.iterrows():
        smile = row['SMILES']  # Assuming the column name with SMILES is 'SMILES'
        
        # Convert SMILES to molecule
        mol = Chem.MolFromSmiles(smile)
        
        if mol is None:
            print(f"Error: Invalid SMILES string '{smile}'")
            continue
        
        # Optimize the molecule
        optimized_mol = optimize_mol(mol)
        
        if add_properties:
            # Optionally add any properties (like SMILES) to the molecule in the SDF
            optimized_mol.SetProp('SMILES', smile)

            # Add ZINC ID as a property
            optimized_mol.SetProp('ZINC_ID', row['ZINC_ID'])

        # Write the optimized molecule to the SDF file
        writer.write(optimized_mol)
    
    # Close the writer
    writer.close()
    print(f"Conversion complete. SDF saved to: {sdf_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", "-i", type=str, required=True, help="Path to the input CSV file with SMILES strings.")
    parser.add_argument("--output_sdf_file", "-o", type=str, required=True, help="Path to the output SDF file.")
    parser.add_argument("--add_properties", default=False, action="store_true", help="Add properties (like SMILES) to the SDF.")
    args = parser.parse_args()

    convert_smiles_to_sdf(args.input_file, args.output_sdf_file, args.add_properties)   
