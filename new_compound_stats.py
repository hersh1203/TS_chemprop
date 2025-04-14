import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, QED
from SAscorer import sascorer  # Make sure this points to Ertl's script

# === Load CSV ===
input_csv = "/Users/harshmallow/Desktop/chemprop_ensemble_eval_results.csv"
df = pd.read_csv(input_csv)

# === Function to compute SA, Lipinski, and QED ===
def compute_all(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None, None, None
    
    # SA Score
    sa_score = sascorer.calculateScore(mol)
    
    # Lipinski
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    violations = sum([
        mw >= 500,
        logp >= 5,
        hbd > 5,
        hba > 10
    ])
    lipinski_pass = violations == 0

    # QED
    qed_score = QED.qed(mol)
    
    return sa_score, violations, lipinski_pass, qed_score

# === Apply to all SMILES ===
results = df["smiles"].apply(compute_all)

# === Assign columns ===
df[["SA_Score", "Lipinski_Violations", "Lipinski_Pass", "QED"]] = pd.DataFrame(results.tolist(), index=df.index)

# === Save to CSV ===
output_file = "chemprop_with_SA_Lipinski_QED.csv"
df.to_csv(output_file, index=False)

print(f"✅ SA + Lipinski + QED analysis complete. Output saved to:\n{output_file}")
