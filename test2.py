from flask import Flask, render_template, request, redirect, session, url_for, flash
import mysql.connector
from werkzeug.security import generate_password_hash, check_password_hash
from rdkit import Chem
from rdkit.Chem import Descriptors
import joblib
import numpy as np

app = Flask(__name__)
app.secret_key = 'your_secret_key_here'

# Database connection - update with your AWS credentials
db = mysql.connector.connect(
    host="database-1.cvakiyk6o3wb.us-east-2.rds.amazonaws.com",     # e.g., "your-db.us-east-1.rds.amazonaws.com"
    user="admin",
    password="ChemicalCompound4!",
    database="userInfo"
)
cursor = db.cursor()

# Load ML models using joblib
alogp_model = joblib.load('model/alogp_rf_model.pkl')
cxlogp_model = joblib.load('model/cx_logp_rf_model.pkl')
mw_model = joblib.load('model/molecular_weight_rf_model.pkl')
psa_model = joblib.load('model/polar_surface_area_rf_model.pkl')
rot_model = joblib.load('model/rotatable_bonds_rf_model.pkl')

# Feature extraction from SMILES
def extract_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.TPSA(mol),
        Descriptors.NumRotatableBonds(mol),
        Descriptors.HeavyAtomCount(mol),
        0
    ]

@app.route('/')
def index():
    return redirect('/login')

@app.route('/register', methods=['GET', 'POST'])
def register():
    if request.method == 'POST':
        email = request.form['email']
        password = generate_password_hash(request.form['password'])

        cursor.execute("SELECT * FROM users WHERE email = %s", (email,))
        if cursor.fetchone():
            flash('Email already registered')
            return redirect('/register')

        cursor.execute("INSERT INTO users (email, password_hash) VALUES (%s, %s)", (email, password))
        db.commit()
        flash('Registered successfully, please log in')
        return redirect('/login')
    return render_template('register.html')

@app.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        email = request.form['email']
        password = request.form['password']

        cursor.execute("SELECT * FROM users WHERE email = %s", (email,))
        user = cursor.fetchone()

        if user and check_password_hash(user[2], password):  # user[2] is password_hash
            session['user'] = email
            return redirect('/home')
        else:
            flash('Invalid credentials')
    return render_template('login.html')

@app.route('/home', methods=['GET', 'POST'])
def home():
    if 'user' not in session:
        return redirect('/login')

    predictions = None
    error = None

    if request.method == 'POST':
        smiles = request.form['smiles']
        features = extract_features(smiles)

        if features:
            features = np.array(features).reshape(1, -1)
            predictions = {
                'ALogP': alogp_model.predict(features)[0],
                'CxLogP': cxlogp_model.predict(features)[0],
                'Molecular Weight': mw_model.predict(features)[0],
                'Polar Surface Area': psa_model.predict(features)[0],
                'Rotatable Bonds': rot_model.predict(features)[0],
            }
        else:
            error = "Invalid SMILES notation."

    return render_template('home.html', user=session['user'], predictions=predictions, error=error)

@app.route('/logout')
def logout():
    session.pop('user', None)
    return redirect('/login')

if __name__ == '__main__':
    app.run(debug=True)