from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

candidates = {
    "1. Aminoguanidinium 5-aminotetrazolate": "NNC(N)=[NH2+].Nc1nnn[n-]1",
    "2. Guanidinium dinitramide": "NC(N)=[NH2+].O=[N+]([O-])[N-][N+](=O)[O-]",
    "3. Diaminoguanidinium dicyanamide": "N#C[N-]C#N.NNC(=[NH2+])NN",
    "4. BMIM 5-aminotetrazolate": "CCCC[n+]1ccn(C)c1.Nc1nnn[n-]1",
    "5. N-Methylmorpholinium dicyanamide": "C[NH+]1CCOCC1.N#C[N-]C#N",
    "6. Hydroxylammonium dicyanamide": "N#C[N-]C#N.[NH3+]O",
    "7. 1,5-Diaminotetrazolium nitrocyanamide": "N#C[N-][N+](=O)[O-].Nc1[nH+]nnn1N",
    "8. Choline dicyanamide": "C[N+](C)(C)CCO.N#C[N-]C#N",
    "9. 1-Ethylnicotinium nitrocyanamide": "CC[n+]1cccc(C(=O)O)c1.N#C[N-][N+](=O)[O-]",
    "10. 4-Amino-triazolium 5-aminotetrazolate": "Nc1nnn[n-]1.Nn1cn[nH+]c1",
}

# Structural-alert catalogs available in RDKit (fast red-flag screening)
params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
catalog = FilterCatalog(params)

print(f"{'Candidate':45s} {'MW':>7s} {'cLogP':>7s} {'TPSA':>6s} {'HBD':>4s} {'HBA':>4s}  Structural alerts (whole ion pair)")
for name, smi in candidates.items():
    mol = Chem.MolFromSmiles(smi)
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    matches = catalog.GetMatches(mol)
    alert_names = [m.GetDescription() for m in matches] if matches else ["none"]
    print(f"{name:45s} {mw:7.1f} {logp:7.2f} {tpsa:6.1f} {hbd:4d} {hba:4d}  {'; '.join(alert_names)}")