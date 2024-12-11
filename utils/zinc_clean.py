"""
Author: Yonglan Liu
Date: 2024-12-10
This script is to process the ZINC database, which includes
1. remove duplicates
2. check validity of SMILES strings
3. check druglikeness of SMILES strings: 
    a. Lipinski’s rule of 5 
    b. Quantitative Estimate of Druglikeness (QED score) with a cutoff of 0.5
The final output is a CSV file with the cleaned and filtered ZINC molecules."""

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

import argparse

# Function to read ZINC database
def load_zinc_file(file_path: str) -> pd.DataFrame:
    """
    Download ZINC database from https://files.docking.org/zinc20-ML/
    code: wget https://files.docking.org/zinc20-ML/ZINC20-ML_smiles.tar.gz
    extract: tar -xvzf ZINC20-ML_smiles.tar.gz
    """
    df = pd.read_csv(file_path, sep=" ", header=None, names=["SMILES", "ZINC_ID"]) 
    # test, load 100 molecules
    # df = pd.read_csv(file_path, sep=" ", header=None, names=["SMILES", "ZINC_ID"], nrows=100)
    return df

# Function to delte replicates
def remove_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop_duplicates(subset='SMILES', keep='first')

# Funtion to check validity of SMILES
def check_validity(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    # Check if molecule is valid
    if mol is None:
        return False  # Invalid molecule
    # Check if molecule is connected properly (single molecule, no disconnected fragments)
    if mol.GetNumAtoms() == 0 or not mol.GetNumBonds():
        return False  # Invalid molecule
    return True

# Function to filter out undesirable molecules
def determine_druglikeness(smiles: str) -> bool:
    """Calculates the specified property for the input smiles 
    and returns a boolean 
    1. Lipinski’s rule of 5 compliance
        Molecule weight <= 500 Dalton 
        LogP <= 5
        Number of hydrogen donors <= 5
        Number of hydrogen acceptors <= 10
    2. Quantitative Estimate of Druglikeness (QED score) with a cutoff of 0.5
    """
    mol = Chem.MolFromSmiles(smiles)
    if ((Descriptors.ExactMolWt(mol) <= 500) and 
        (Chem.Crippen.MolLogP(mol) <= 5) and 
        (Descriptors.NumHDonors(mol) <= 5) and 
        (Descriptors.NumHAcceptors(mol) <= 10) and 
        (Chem.QED.qed(mol) >= 0.5)):
        return True
    else:
        return False

def prepare_zinc(file_path: str) -> pd.DataFrame:
    df = load_zinc_file(file_path)
    df = remove_duplicates(df)
    df["validity"] = df["SMILES"].apply(check_validity)
    df_valid = df[df["validity"] == True]
    df_valid["druglikeness"] = df_valid["SMILES"].apply(determine_druglikeness)
    df_final = df_valid[df_valid["druglikeness"] == True]
    return df_final.drop(columns=["validity", "druglikeness"])


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_zinc_file", "-i", type=str, required=True, help="Path to ZINC file")
    parser.add_argument("--output_zinc_file", '-o', type=str, required=True, help="Path to output ZINC file")
    args = parser.parse_args()
    df_final = prepare_zinc(args.input_zinc_file)
    df_final.to_csv(args.output_zinc_file, index=False, header=False, sep=" ")
