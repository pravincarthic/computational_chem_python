from openbabel import pybel
import numpy as np
import os

def min_dist_check(mol, label="Molecule"):
    """Calculates and prints the minimum pairwise interatomic distance."""
    coords = np.array([a.coords for a in mol.atoms])
    
    # Vectorized pairwise distance calculation (avoids slow Python loops)
    diff = coords[:, None, :] - coords[None, :, :]
    dists = np.linalg.norm(diff, axis=-1)
    np.fill_diagonal(dists, np.inf)
    min_d = dists.min()
    
    if min_d > 1.4:
        status = "✓ OK"
    elif min_d > 1.0:
        status = "⚠️ SHORT"
    else:
        status = "❌ CLASH"
        
    print(f"  {label} min dist: {min_d:.3f} Å  {status}")
    return min_d

def process_smiles(smi_string, output_file="compound.xyz"):
    """Processes a SMILES string, generates 3D geometries, and saves to XYZ."""
    
    # Dynamically split fragments by '.' (standard SMILES separator)
    smiles_frags = [s.strip() for s in smi_string.split(".") if s.strip()]
    all_atoms_data = []
    current_x_offset = 0.0
    total_atom_count = 0

    print(f"Detected {len(smiles_frags)} fragment(s).")

    for i, frag_smi in enumerate(smiles_frags):
        print(f"\nProcessing fragment {i+1}...")
        mol = pybel.readstring("smi", frag_smi)
        
        # MMFF94 is typically much better than UFF for strained silicon/nitrogen cages
        mol.make3D(forcefield="mmff94", steps=5000)
        min_dist_check(mol, f"Fragment {i+1}")
        
        # VERSION-SAFE TRICK: Use OpenBabel's internal writer to get symbols/coords
        # This completely avoids common AttributeError bugs with OBElementTable across different conda versions
        xyz_output = mol.write("xyz").splitlines()
        
        # XYZ format: line 0 = count, line 1 = comment, line 2+ = atoms
        atom_lines = xyz_output[2:]
        
        frag_coords = []
        parsed_atoms = []
        for line in atom_lines:
            parts = line.split()
            if len(parts) < 4: 
                continue
            symbol = parts[0]
            x, y, z = map(float, parts[1:4])
            frag_coords.append([x, y, z])
            parsed_atoms.append((symbol, x, y, z))
            
        if not frag_coords: 
            continue
        
        # Calculate bounding box to prevent fragments from overlapping
        coords_np = np.array(frag_coords)
        min_x = coords_np[:, 0].min()
        max_x = coords_np[:, 0].max()
        
        # Shift this fragment so it starts right after the previous one
        shift = current_x_offset - min_x
        
        for symbol, x, y, z in parsed_atoms:
            all_atoms_data.append(f"{symbol:<2} {x + shift:>12.6f} {y:>12.6f} {z:>12.6f}\n")
            total_atom_count += 1
            
        # Update offset: width of the current fragment + a safe 10.0 Å buffer
        current_x_offset += (max_x - min_x) + 10.0

    # Write the final combined XYZ file using Python's native file handling
    with open(output_file, "w") as f:
        f.write(f"{total_atom_count}\nGenerated combined system\n")
        f.writelines(all_atoms_data)
        
    print(f"\nSuccess: Saved {total_atom_count} atoms to {output_file}")
    print("\nNext steps:")
    print(f"  xtb {output_file} --opt --gfn 2 --verbose")

# --- Execution ---
if __name__ == "__main__":
    # You can safely swap this string out for single molecules or multi-fragment salts (e.g., "C.C")
    user_smiles = "N1=NN2N3N4N5N=NN=NN6N(N7N(N=N1)N1N7N3N21)N1N6N5N41 |c:0,6,8,14|"
    process_smiles(user_smiles)

# Optional: os.system("xtb compound.xyz --opt --gfn 2 --verbose") followed by os.system("xtb xtbopt.xyz --opt --gfn 2 --verbose --cycles 500") to run the optimization immediately after generation.