from flask import Flask, request, jsonify, render_template
import joblib
from rdkit import Chem
from rdkit.Chem import Descriptors

app = Flask(__name__)

# Descriptor function
def smiles_to_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.TPSA(mol),
        Descriptors.NumRotatableBonds(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.NumHDonors(mol)
    ]

# Properties you're predicting
properties = ["Molecular_Weight", "ALogP", "Polar_Surface_Area", "Rotatable_Bonds", "CX_LogP"]

@app.route('/')
def index():
    return render_template("index.html")

@app.route('/predict', methods=['POST'])
def predict():
    data = request.get_json()
    smiles = data.get("smiles", "")
    descriptors = smiles_to_descriptors(smiles)

    if descriptors is None:
        return jsonify({"error": "Invalid SMILES"}), 400

    predictions = {}
    for prop in properties:
        model = joblib.load(f"{prop.lower().replace(' ', '_')}_rf_model.pkl")
        predictions[prop] = round(model.predict([descriptors])[0], 3)

    return jsonify(predictions)

if __name__ == '__main__':
    app.run(debug=True)
