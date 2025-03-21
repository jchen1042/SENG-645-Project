import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

# Initialize ChEMBL client
molecule_client = new_client.molecule

# Prepare list of compound IDs
chembl_ids = [f"CHEMBL{i}" for i in range(1, 100)]

# Store all data here
data = []

for chembl_id in chembl_ids:
    print(f"Fetching data for {chembl_id}...")

    try:
        # Get compound data
        mol_data = molecule_client.get(chembl_id)

        # Get canonical SMILES
        smiles = mol_data.get("molecule_structures", {}).get("canonical_smiles", None)
        if not smiles:
            print(f"  ❌ SMILES not found for {chembl_id}")
            continue

        # Convert to RDKit Mol
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            print(f"  ❌ Invalid SMILES for {chembl_id}")
            continue

        # Calculate properties
        mol_weight = Descriptors.MolWt(mol)
        alogp = Crippen.MolLogP(mol)
        polar_surface_area = rdMolDescriptors.CalcTPSA(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        cx_logp = mol_data.get("molecule_properties", {}).get("cx_logp", "Not available")

        # Append to data list
        data.append({
            "ChEMBL_ID": chembl_id,
            "SMILES": smiles,
            "Molecular_Weight": mol_weight,
            "ALogP": alogp,
            "Polar_Surface_Area": polar_surface_area,
            "Rotatable_Bonds": rotatable_bonds,
            "CX_LogP": cx_logp
        })

    except Exception as e:
        print(f"  ⚠️ Error with {chembl_id}: {e}")
        continue

# Save to CSV
df = pd.DataFrame(data)
df.to_csv("chembl_compound_properties.csv", index=False)
print("\n✅ Saved as chembl_compound_properties.csv")






