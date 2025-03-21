# Step 1: Import libraries
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import joblib

# Step 2: Load dataset
df = pd.read_csv("chembl_compound_properties.csv")  # Update to your filename

# Step 3: Clean & prepare data
df = df[df["CX_LogP"] != "Not available"]
df["CX_LogP"] = pd.to_numeric(df["CX_LogP"], errors='coerce')
df.dropna(inplace=True)

# Step 4: Convert SMILES to molecular descriptors
def smiles_to_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return [
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            Descriptors.TPSA(mol),
            Descriptors.NumRotatableBonds(mol),
            Descriptors.NumHAcceptors(mol),
            Descriptors.NumHDonors(mol)
        ]
    return None

df["features"] = df["SMILES"].apply(smiles_to_descriptors)
df = df[df["features"].notnull()]  # Remove invalid SMILES

# Step 5: Create feature and target matrices
X = np.array(df["features"].tolist())
target_columns = ["Molecular_Weight", "ALogP", "Polar_Surface_Area", "Rotatable_Bonds", "CX_LogP"]
y = df[target_columns]

# Step 6: Train & save model for each target property
models = {}
for col in target_columns:
    print(f"\n📊 Training for: {col}")
    y_col = y[col]

    X_train, X_test, y_train, y_test = train_test_split(X, y_col, test_size=0.2, random_state=42)

    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)

    preds = model.predict(X_test)
    print(f"  RMSE: {np.sqrt(mean_squared_error(y_test, preds)):.3f}")
    print(f"  R² Score: {r2_score(y_test, preds):.3f}")

    joblib.dump(model, f"{col.lower().replace(' ', '_')}_rf_model.pkl")
    models[col] = model
    print(f"✅ Model saved: {col.lower().replace(' ', '_')}_rf_model.pkl")

# Step 7: Prediction function
def predict_properties_from_smiles(smiles):
    desc = smiles_to_descriptors(smiles)
    if desc is None:
        return "Invalid SMILES"

    predictions = {}
    for col in target_columns:
        model = joblib.load(f"{col.lower().replace(' ', '_')}_rf_model.pkl")
        predictions[col] = round(model.predict([desc])[0], 3)
    return predictions

# Step 8: Example usage
test_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
print("\n🔮 Prediction for:", test_smiles)
print(predict_properties_from_smiles(test_smiles))


