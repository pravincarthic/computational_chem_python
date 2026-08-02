# Predicting Monopropellant Decomposition Products: Chemistry Guide

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

**Exception:** Only in extreme O-excess (rare for your compounds) might NO form, but assume N₂.

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
  - Make 7.5 CO₂ uses 15O (but you only have 6C)
  - So: 6C + 12O → 6CO₂ (uses all 6C, 12O)
  - 3O left → 1.5 O₂

**If O < C:**
- Most C → CO (limited by O)
- Leftover C → elemental C (solid, not ideal for rockets)
- Example: 6C + 5O available
  - 5O allows 5CO
  - 1C left (doesn't form anything, stays as C)

For your compounds: **They're all O-balanced or O-rich**, so assume CO₂.

---

### **Step 4: Account for Hydrogen (After F is removed)**

**After F has taken H for HF:**

**If H remaining & O available:**
- H → H₂O (with O)

Example: Originally 12H, 12 used for HF → 0H remaining
- If you had 12H and only 6F:
  - 6H → 6HF
  - 6H remaining → 3H₂O

**If H remaining but no O left:**
- H → H₂

This is rare for your compounds (O-rich), so usually → H₂O.

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
   
RESULT: 16 N₂ + 12 HF + 6 CO₂ + 6 H₂O + 6 O₂
Wait, we had 0 H remaining, so 0 H₂O. Let me recalculate...

Actually:
RESULT: 16 N₂ + 12 HF + 6 CO₂ + 6 O₂

Verify atom balance:
  N: 16×2 = 32 ✓
  H: 12×1 = 12 ✓
  F: 12×1 = 12 ✓
  C: 6×1 = 6 ✓
  O: 6×2 + 6×2 = 12 + 12 = 24 ✓
```

---

## Your 10 Compounds: Predicted Products

Apply the algorithm to each:

### **C₃H₃N₂₁O₁₀**

```
1. Extract: C=3, H=3, F=0, N=21, O=10

2. F→HF: 0 HF, H_remaining = 3

3. N→N₂: 21/2 = 10.5 N₂

4. C→CO₂: 3C needs 6O, have 10O → 3 CO₂, O_left = 4

5. H→H₂O: 3H → 1.5 H₂O, O_left = 4 - 1.5 = 2.5

6. O→O₂: 2.5O → 1.25 O₂

PRODUCTS: 10.5 N₂ + 3 CO₂ + 1.5 H₂O + 1.25 O₂
```

### **C₂N₂₂O₁₆**

```
1. Extract: C=2, H=0, F=0, N=22, O=16

2. F→HF: 0 HF

3. N→N₂: 22/2 = 11 N₂

4. C→CO₂: 2C needs 4O, have 16O → 2 CO₂, O_left = 12

5. H→H₂O: 0H → 0 H₂O

6. O→O₂: 12O → 6 O₂

PRODUCTS: 11 N₂ + 2 CO₂ + 6 O₂
```

### **C₂N₂₈O₂₀**

```
1. Extract: C=2, H=0, F=0, N=28, O=20

2. F→HF: 0 HF

3. N→N₂: 28/2 = 14 N₂

4. C→CO₂: 2C needs 4O, have 20O → 2 CO₂, O_left = 16

5. H→H₂O: 0H → 0 H₂O

6. O→O₂: 16O → 8 O₂

PRODUCTS: 14 N₂ + 2 CO₂ + 8 O₂
```

### **C₄H₈F₄N₁₆O₁₆**

```
1. Extract: C=4, H=8, F=4, N=16, O=16

2. F→HF: 4 HF, H_remaining = 8-4 = 4

3. N→N₂: 16/2 = 8 N₂

4. C→CO₂: 4C needs 8O, have 16O → 4 CO₂, O_left = 8

5. H→H₂O: 4H → 2 H₂O, O_left = 8-2 = 6

6. O→O₂: 6O → 3 O₂

PRODUCTS: 8 N₂ + 4 HF + 4 CO₂ + 2 H₂O + 3 O₂
```

### **C₄H₈F₅N₂₃O₁₆**

```
1. Extract: C=4, H=8, F=5, N=23, O=16

2. F→HF: 5 HF, H_remaining = 8-5 = 3

3. N→N₂: 23/2 = 11.5 N₂

4. C→CO₂: 4C needs 8O, have 16O → 4 CO₂, O_left = 8

5. H→H₂O: 3H → 1.5 H₂O, O_left = 8-1.5 = 6.5

6. O→O₂: 6.5O → 3.25 O₂

PRODUCTS: 11.5 N₂ + 5 HF + 4 CO₂ + 1.5 H₂O + 3.25 O₂
```

### **C₄HF₃N₁₉O₁₂**

```
1. Extract: C=4, H=1, F=3, N=19, O=12

2. F→HF: 3 HF, H_remaining = 1-3 = -2 (PROBLEM!)
   → Not enough H to make 3 HF. Only 1 H available.
   → Only 1 HF forms, 2 F left over
   
   For this case: Some F can't find H, forms F₂ or stays unaccounted
   → Assume 1 HF only (most likely scenario)
   → H_remaining = 0, F_unaccounted = 2
   
   OR assume some flexibility: this is edge case, compute both scenarios

   Scenario A (F→HF maximally):
   1 HF, H_remaining = 0, F₂ = 1
   
3. N→N₂: 19/2 = 9.5 N₂

4. C→CO₂: 4C needs 8O, have 12O → 4 CO₂, O_left = 4

5. H→H₂O: 0H → 0 H₂O

6. O→O₂: 4O → 2 O₂

PRODUCTS (Scenario A): 9.5 N₂ + 1 HF + 4 CO₂ + 2 O₂ + 1 F₂
(F₂ is unlikely, so maybe this compound has issues, or compute differently)

Better approach: This might decompose to:
9.5 N₂ + 4 CO₂ + 2 O₂ (ignore F, treat as leftover)
OR compute scenarios with/without HF formation
```

### **C₅H₈F₃N₂₁O₁₆**

```
1. Extract: C=5, H=8, F=3, N=21, O=16

2. F→HF: 3 HF, H_remaining = 8-3 = 5

3. N→N₂: 21/2 = 10.5 N₂

4. C→CO₂: 5C needs 10O, have 16O → 5 CO₂, O_left = 6

5. H→H₂O: 5H → 2.5 H₂O, O_left = 6-2.5 = 3.5

6. O→O₂: 3.5O → 1.75 O₂

PRODUCTS: 10.5 N₂ + 3 HF + 5 CO₂ + 2.5 H₂O + 1.75 O₂
```

### **C₆H₁₂F₁₂N₃₂O₂₄**

```
1. Extract: C=6, H=12, F=12, N=32, O=24

2. F→HF: 12 HF, H_remaining = 12-12 = 0

3. N→N₂: 32/2 = 16 N₂

4. C→CO₂: 6C needs 12O, have 24O → 6 CO₂, O_left = 12

5. H→H₂O: 0H → 0 H₂O

6. O→O₂: 12O → 6 O₂

PRODUCTS: 16 N₂ + 12 HF + 6 CO₂ + 6 O₂
```

### **C₇F₄N₂₈O₁₈**

```
1. Extract: C=7, H=0, F=4, N=28, O=18

2. F→HF: No H available! 
   → 0 HF, but 4 F leftover (forms F₂? Stays unaccounted?)
   → Likely scenario: F stays bonded to C or just ignore
   → Conservative: Compute without F products

3. N→N₂: 28/2 = 14 N₂

4. C→CO₂: 7C needs 14O, have 18O → 7 CO₂, O_left = 4

5. H→H₂O: 0H → 0 H₂O

6. O→O₂: 4O → 2 O₂

PRODUCTS (ignoring F): 14 N₂ + 7 CO₂ + 2 O₂

Alternative if F₂ forms:
14 N₂ + 7 CO₂ + 2 O₂ + 2 F₂

Recommend: Compute both scenarios
```

---

## Summary Table: Predicted Products

| Compound | Scenario | Products |
|----------|----------|----------|
| C₃H₃N₂₁O₁₀ | Standard | 10.5 N₂ + 3 CO₂ + 1.5 H₂O + 1.25 O₂ |
| C₂N₂₂O₁₆ | O-rich | 11 N₂ + 2 CO₂ + 6 O₂ |
| C₂N₂₈O₂₀ | O-rich | 14 N₂ + 2 CO₂ + 8 O₂ |
| C₄H₈F₄N₁₆O₁₆ | With HF | 8 N₂ + 4 HF + 4 CO₂ + 2 H₂O + 3 O₂ |
| C₄H₈F₅N₂₃O₁₆ | With HF | 11.5 N₂ + 5 HF + 4 CO₂ + 1.5 H₂O + 3.25 O₂ |
| C₄HF₃N₁₉O₁₂ | H-limited | See note (edge case) |
| C₅H₈F₃N₂₁O₁₆ | With HF | 10.5 N₂ + 3 HF + 5 CO₂ + 2.5 H₂O + 1.75 O₂ |
| C₆H₁₂F₁₂N₃₂O₂₄ | With HF | 16 N₂ + 12 HF + 6 CO₂ + 6 O₂ |
| C₇F₄N₂₈O₁₈ | H-limited | See note (edge case) |

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

**Recommendation:** Compute BOTH scenarios for these compounds. Report the one with lower (more negative/favorable) ΔH as "most likely."

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

## For Your QE Calculations

Once you have products, you need to compute:

1. **E(monopropellant)** - Already have from your .xyz files ✓
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

EOF
