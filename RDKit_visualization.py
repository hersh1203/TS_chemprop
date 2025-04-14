from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd

# Load your data
df = pd.read_csv("/Users/harshmallow/Desktop/chemprop_ensemble_eval_results.csv")

# Sort by score descending
df_sorted = df.sort_values("score", ascending=False)

# Pick top N molecules
top_n = 12
top_df = df_sorted.head(top_n)

# Convert to RDKit molecules
mols = [Chem.MolFromSmiles(smile) for smile in top_df["smiles"]]

# Create legends with Rank, Name, short SMILES, and Score
legends = [
    f"#{i+1}\n{row['Name']}\n{row['smiles'][:25]}...\nScore: {row['score']:.3f}"
    for i, row in top_df.iterrows()
]

# Draw grid of molecules
img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(300, 300), legends=legends)
img.show()
