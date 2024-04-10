# QuantumPSP
Solving Project Scheduling Problem instances using D-Wave's Quantum Annealing

WORK IN PROGRESS

# Preparation
Download datasets from https://solutionsupdate.ugent.be/rcpsp and put the folders inside the data/ folder or use the pre-downloaded instances.

# Installation
```bash
pip install dwave-ocean-sdk
pip install tqdm
```

# How to run 
```bash
python main.py \
    --instance data/DC2/DC2/npv75/rcpspdc364.rcp # for example \
    --time_limit 60 # optional, tells the D-wave hybrid annealer how much time to use on the problem.
                    # Increase if the results are unsatisfactory. Watch out, you will run out of time quicker!
```
