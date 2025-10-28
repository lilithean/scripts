# trivex2.4.py — Ternary Convex Hull Plotter

`trivex2.4.py` is a Python script for constructing and visualizing **ternary convex-hull diagrams** from formation enthalpy data.
It supports both 2D and 3D plotting modes, contour visualization, and distance-to-hull analysis for identifying stable and metastable compounds in multicomponent systems.

---

### Features

* Generates convex-hull diagrams from user-provided compositions and enthalpies.
* Supports **2D contour** and **3D surface** visualization.
* Calculates **vertical distance to the convex hull**, identifying stability regions.
* Accepts **custom color themes** and unit specifications.
* Outputs publication-quality figures (`png`, `pdf`, `eps`, `svg`, etc.) via Matplotlib.

---

### Usage

```bash
python trivex2.4.py input.txt
```

Example input file:

```
pressure 0
color spring
file hull.pdf
contour 15
3d
component
Ca1 S1 H1 meV/atom
3   3   10   -72.631
4   4   2    -88.923
7   7   8    -15.376
```

---

### Outputs

* 2D or 3D convex-hull diagram (`.pdf`, `.png`, etc.)
* Printed table of compositions, enthalpies, and vertical distances from the hull
* Optional animation (planned feature)

---

### Citation

If you use or adapt this script, please cite:

> **X. Wang et al.**, *Structure, stability, and superconductivity of N-doped lutetium hydrides at kbar pressures*,
> *Phys. Rev. B* **108**, 014511 (2023).
> DOI: [10.1103/PhysRevB.108.014511](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.108.014511)

---

### Notes

* Requires Python ≥3.7, with `numpy`, `matplotlib`, and `scipy`.
* Licensed under *Xiaoyu’s You Happy Jiu OK To Public License*:
  “**You just HAPPILY DO WHATEVER YOU WANT TO.**”
