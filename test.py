from flask import Flask, render_template, request, redirect, url_for
from flask_sqlalchemy import SQLAlchemy
from flask_login import LoginManager, UserMixin, login_user, login_required, logout_user, current_user
from werkzeug.security import generate_password_hash, check_password_hash
from pymongo import MongoClient
from bson.objectid import ObjectId
import pickle
import joblib
import numpy as np

# Flask setup
app = Flask(__name__)
app.secret_key = 'secretkey123'

# ✅ MySQL (AWS RDS)
app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+pymysql://admin:ChemicalCompound4%21@database-1.cvakiyk6o3wb.us-east-2.rds.amazonaws.com:3306/myapp'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

# ✅ MongoDB Atlas for search history
mongo_uri = "mongodb+srv://myuser:ChemicalCompound4!@chemicalcompoundhistory.6hdzk.mongodb.net/?retryWrites=true&w=majority&appName=ChemicalCompoundHistory"
mongo_client = MongoClient(mongo_uri)
mongo_db = mongo_client["ChemicalCompoundHistory"]
history_collection = mongo_db["search_history"]

# Load ML models using joblib
alogp_model = joblib.load('model/alogp_rf_model.pkl')
cxlogp_model = joblib.load('model/cx_logp_rf_model.pkl')
mw_model = joblib.load('model/molecular_weight_rf_model.pkl')
psa_model = joblib.load('model/polar_surface_area_rf_model.pkl')
rot_model = joblib.load('model/rotatable_bonds_rf_model.pkl')

# ✅ Login setup
login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = "login"

# ✅ User table model
class User(UserMixin, db.Model):
    id = db.Column(db.Integer, primary_key=True)
    email = db.Column(db.String(100), unique=True, nullable=False)
    password = db.Column(db.String(255), nullable=False)

# Load user for session
@login_manager.user_loader
def load_user(user_id):
    return User.query.get(int(user_id))

# Root route redirects to login
@app.route("/")
def index():
    return redirect(url_for("login"))

# ✅ Register page
@app.route("/register", methods=["GET", "POST"])
def register():
    if request.method == "POST":
        email = request.form["email"]
        password = request.form["password"]

        existing_user = User.query.filter_by(email=email).first()
        if existing_user:
            return "User already exists!"

        hashed_password = generate_password_hash(password)
        new_user = User(email=email, password=hashed_password)
        db.session.add(new_user)
        db.session.commit()
        return redirect(url_for("login"))

    return render_template("register.html")

# ✅ Login page
@app.route("/login", methods=["GET", "POST"])
def login():
    if request.method == "POST":
        email = request.form["email"]
        password = request.form["password"]

        user = User.query.filter_by(email=email).first()
        if user and check_password_hash(user.password, password):
            login_user(user)
            return redirect(url_for("home"))
        else:
            return "Invalid credentials!"

    return render_template("login.html")

# ✅ Logout
@app.route("/logout")
@login_required
def logout():
    logout_user()
    return redirect(url_for("login"))

# ✅ Home (prediction) page
@app.route("/home")
@login_required
def home():
    return render_template("index.html")

# ✅ Predict route
@app.route("/predict", methods=["POST"])
@login_required
def predict():
    try:
        smiles = request.form["smiles"]
        features = compute_descriptors(smiles)
        prediction = model.predict([features])[0]

        # Save prediction to MongoDB
        history_collection.insert_one({
            "user_email": current_user.email,
            "smiles": smiles,
            "prediction": prediction
        })

        return render_template("index.html", prediction=prediction, smiles=smiles)

    except Exception as e:
        print("Prediction error:", e)
        return render_template("index.html", prediction="Prediction error")

# ✅ History page
@app.route("/history")
@login_required
def history():
    records = list(history_collection.find({"user_email": current_user.email}))
    return render_template("history.html", records=records)

# ✅ SMILES descriptor computation
def compute_descriptors(smiles):
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return [
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            Descriptors.TPSA(mol),
            Descriptors.NumRotatableBonds(mol),
            Descriptors.MolMR(mol)
        ]
    else:
        raise ValueError("Invalid SMILES")

# ✅ Run app
if __name__ == "__main__":
    db.create_all()
    app.run(debug=True)