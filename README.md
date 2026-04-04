# Computational Chemistry with Python ⚛️🐍 #
This repository is where I "dabble" with the invisible forces that hold the universe together—specifically, the chemical bonds in rocket propellants. I'm currently self-learning how to use Python to simulate molecular energy so I can understand why some fuels (like Hydrazine) are so much more "exciting" (and energetic) than others.

**🧪 What’s happening here?**

I’m moving away from just reading about chemistry in textbooks and trying to actually compute it. I self-read about bond energy in chemistry and I am still learning the theory behind it.
Most people use Kaggle for Machine Learning, but I’ve been "repurposing" Kaggle’s free GPUs as a cloud-based supercomputer to run physics simulations. It is fun to get free GPU resources and simulation environment in Kaggle.

Current Focus:
- Potential Energy Surfaces (PES): Using the Hessian Matrix (from my Multivariable Calculus self-studies) to find local minima and stable molecular geometries.
- Quantum ESPRESSO: Someone said this is good to simulate. I am still learning what is Density Functional Theory (DFT) calculation to "see" electron density.

Because the computation of full DFT for screening large molecules is very high, I found a pipeline using xtb. By automating the conversion from Moldraw-generated SMILES to GFN2-xTB optimized geometries via Python, I found a quicker workflow to analyze bond vibration modes / energies

If we want to go to Mars or build better aircraft, we need to understand the chemistry of propulsion at an atomic level.
