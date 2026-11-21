# Predicting Monopropellant Decomposition Products
### Note: RocketCEA already does most of the work, however this is just hypothesising and explaining why.

## Core Principle

**Decomposition products form based on thermodynamic stability, not kinetics.**

Atoms preferentially form the **lowest-energy compounds** they can make:
- N prefers N₂ (not N atoms, NO, or N₂O)
- F prefers HF (not F₂ or other F compounds)
- C prefers CO₂ (not CO or C)
- H prefers H₂O (not H₂ or H)

## Hierarchy of Stability (What Bonds Form First)

```
Most stable → Least stable

HF  >>>  H₂O  >>>  CO₂  >>>  N₂  >>>  CO  >>>  O₂  >>>  H₂

Why?
HF:   ΔH_f = -272 kJ/mol  (massive negative = super stable)
H₂O:  ΔH_f = -242 kJ/mol  (very stable)
CO₂:  ΔH_f = -393 kJ/mol  (very stable)
N₂:   ΔH_f = 0 (by definition, but diatomic is stable)
CO:   ΔH_f = -110 kJ/mol  (less stable than CO₂)
O₂:   ΔH_f = 0 (by definition)
H₂:   ΔH_f = 0 (by definition)
```

**Translation:** When decomposing, atoms grab partners in this order to maximize energy release.

---

## Decision Tree for Predicting Products

### **Step 1: Account for Fluorine First**

**Rule: ALL fluorine → HF**

Why? HF is the most exothermic compound possible (-272 kJ/mol per HF). Fluorine will steal hydrogen from anywhere to form HF.

```
Example: C₆H₁₂F₁₂N₃₂O₂₄
├─ 12 F atoms → 12 HF (uses 12 H)
├─ Remaining: 6C, 0H, 32N, 24O
└─ Continue to step 2...
```

**Edge case:** If F > H, then some HF forms, rest stays as F₂ (unlikely).

---

### **Step 2: Account for Nitrogen**

**Rule: ALL nitrogen → N₂**

N₂ is extremely stable. Even in O-rich conditions, N₂ forms preferentially over NOₓ.

```
Example (continued):
├─ 32 N atoms → 16 N₂
├─ Remaining: 6C, 0H, 24O
└─ Continue to step 3...
```

**Exception:** Only in extreme O-excess (rare for my compounds) might NO form, but assume N₂.

---

### **Step 3: Account for Carbon**

**Check oxygen availability:**

**If O ≥ 2C:**
- All C → CO₂
- Remaining O → O₂

Example: 6C + 24O available
- 6C needs 12O → 6CO₂ ✓
- 24 - 12 = 12O left → 6O₂

**If C ≤ O < 2C:**
- Some C → CO₂, some C → CO
- Example: 6C + 15O available
  - Make 7.5 CO₂ uses 15O (but I only have 6C)
  - So: 6C + 12O → 6CO₂ (uses all 6C, 12O)
  - 3O left → 1.5 O₂

**If O < C:**
- Most C → CO (limited by O)
- Leftover C → elemental C (solid, not ideal for rockets)
- Example: 6C + 5O available
  - 5O allows 5CO
  - 1C left (doesn't form anything, stays as C)

For my compounds: **They're all O-balanced or O-rich**, so assume CO₂.

---

### **Step 4: Account for Hydrogen (After F is removed)**

**After F has taken H for HF:**

**If H remaining & O available:**
- H → H₂O (with O)

Example: Originally 12H, 12 used for HF → 0H remaining
- If I had 12H and only 6F:
  - 6H → 6HF
  - 6H remaining → 3H₂O

**If H remaining but no O left:**
- H → H₂

This is rare for my compounds (O-rich), so usually → H₂O.

---

### **Step 5: Leftover Oxygen**

**After all C, H, N accounted for:**
- Remaining O → O₂

---

## Systematic Algorithm

```
INPUT: Molecular formula (e.g., C₆H₁₂F₁₂N₃₂O₂₄)

1. Extract element counts:
   n_C = 6,  n_H = 12,  n_F = 12,  n_N = 32,  n_O = 24

2. Fluorine → HF:
   HF_formed = n_F = 12
   n_H_remaining = n_H - n_F = 12 - 12 = 0
   
3. Nitrogen → N₂:
   N2_formed = n_N / 2 = 32 / 2 = 16
   
4. Carbon → CO₂ (or CO if O-poor):
   O_needed_for_all_CO2 = 2 * n_C = 2 * 6 = 12
   O_available = n_O = 24
   Since O_available > O_needed:
   CO2_formed = n_C = 6
   O_remaining = 24 - 12 = 12
   
5. Hydrogen → H₂O (if H remaining):
   H2O_formed = n_H_remaining / 2 = 0 / 2 = 0
   H_leftover = 0
   
6. Oxygen → O₂:
   O2_formed = O_remaining / 2 = 12 / 2 = 6
   
RESULT: 16 N₂ + 12 HF + 6 CO₂ + 6 O₂

Verify atom balance:
  N: 16×2 = 32 ✓
  H: 12×1 = 12 ✓
  F: 12×1 = 12 ✓
  C: 6×1 = 6 ✓
  O: 6×2 + 6×2 = 12 + 12 = 24 ✓
```

---


## Edge Cases: H-Limited Compounds

For **C₄HF₃N₁₉O₁₂** and **C₇F₄N₂₈O₁₈** (not enough H for all F):

**Option 1: Conservative (ignore F)**
- Don't form HF if H is insufficient
- Compute scenario without HF
- Use for comparison

**Option 2: Alternative (some HF, rest F₂)**
- Form HF with available H
- Remaining F → F₂ (unlikely but thermodynamically possible)
- Compute scenario with F₂

**Recommended:** Compute BOTH scenarios for these compounds. Report the one with lower (more negative/favorable) ΔH as "most likely."

---

## Validation Checklist

For each product set, verify:

```
For C₆H₁₂F₁₂N₃₂O₂₄ → 16 N₂ + 12 HF + 6 CO₂ + 6 O₂:

Nitrogen count:
  Reactant: 32 N
  Products: 16 N₂ × 2 = 32 N ✓

Hydrogen count:
  Reactant: 12 H
  Products: 12 HF × 1 = 12 H ✓

Fluorine count:
  Reactant: 12 F
  Products: 12 HF × 1 = 12 F ✓

Carbon count:
  Reactant: 6 C
  Products: 6 CO₂ × 1 = 6 C ✓

Oxygen count:
  Reactant: 24 O
  Products: 6 CO₂ × 2 + 6 O₂ × 2 = 12 + 12 = 24 O ✓

All balanced! ✓
```

---

## For QE Calculations

Once we have products, we need to compute:

1. **E(monopropellant)** - Already have from my .xyz files ✓
2. **E(N₂)** - Single N₂ molecule
3. **E(HF)** - Single HF molecule
4. **E(CO₂)** - Single CO₂ molecule
5. **E(H₂O)** - Single H₂O molecule
6. **E(O₂)** - Single O₂ molecule
7. **E(CO)** - Single CO molecule (if needed)
8. **E(H₂)** - Single H₂ molecule (if needed)
9. **E(F₂)** - Single F₂ molecule (rare)

Then: ΔH_decomp = Σ(E_products) - E(monopropellant)

Most negative ΔH = most exothermic = thermodynamically favored.

# Predicted Products and Ratios will be there in the **CEA output**
