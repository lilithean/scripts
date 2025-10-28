# a2F-utilities

This directory contains post-processing scripts for analyzing **electron-phonon coupling (EPC)** data and constructing the **Eliashberg spectral function α2F(ω)** from *ab initio* outputs.
The tools automate the evaluation of mode-resolved coupling constants λ, spectral densities α2F(ω), and derived superconducting parameters such as ωₗₙ and Tc.

---

### Scripts

#### `clgo.py`

Calculates **λ**, **α2F(ω)**, and related quantities from `matdyn.x` and `elph.gamma.*` files produced by **Quantum ESPRESSO**.

**Inputs**

* `freq.gp` — phonon frequency file
* `elph.gamma.N` — electron-phonon linewidths for each q-point
* `lambda.dat` — Fermi-level DOS and smearing parameters

**Outputs**

* `a2F_band.*` — α2F(ω) per band
* `gamma_band.*` — phonon linewidths
* `elph.lambda.*`, `elph.lo.*`, `elph.A.*` — mode-resolved EPC strength and auxiliary data

Use this script to generate α2F(ω) distributions and visualize coupling “bubbles” along the phonon branches.

---

#### `int_a2F.py`

Integrates **α2F(ω)** from a Quantum ESPRESSO `a2f.dos` file and evaluates cumulative functions:

[
\lambda(\omega),\quad \omega_{\ln}(\omega),\quad \langle\omega^2\rangle,\quad T_c
]

**Outputs**

* `int_lambda.dat`, `int_wln.dat`, `int_a2f.dat`, `int_gamma.dat` — integrated results
* `a2F.dat`, `ALPHA2F.OUT` — processed α2F(ω) data for Eliashberg calculations

Run as:

```bash
python int_a2F.py a2f.dos
```

Optionally provide μ* as a second argument (default = 0.10).

---

### Citation

If you use or adapt these scripts, please cite:

> **X. Wang et al.**, *Superconductivity in Dilute Hydrides of Ammonia under Pressure*,
> *J. Phys. Chem. Lett.* **2024**, DOI: [10.1021/acs.jpclett.4c01223](https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.4c01223)

---

### Notes

* Both scripts are written in pure Python 3 and require only standard libraries (`math`, `sys`, `numpy`, `scipy` for optional plotting).
* The results are fully compatible with Quantum ESPRESSO’s `matdyn.x` and `elph.x` outputs.
