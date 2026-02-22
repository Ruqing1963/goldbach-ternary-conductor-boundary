# The Ternary Conductor Boundary

**Why Conductor Rigidity Is Specific to the Binary Goldbach Problem**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

This repository accompanies the paper:

> R. Chen, *The Ternary Conductor Boundary: Why Conductor Rigidity Is Specific to the Binary Goldbach Problem*, February 2026.

We investigate whether the conductor rigidity framework established for the binary Goldbach conjecture (genus 2, GSp(4)) extends to the ternary problem $N = p_1 + p_2 + p_3$ (genus 3, GSp(6)). By computing the **true** discriminant of the genus-3 Frey curve, we show that it does not: the algebraic identity $p^2 - q^2 = (p-q)N$ has no ternary analogue, and the Band Shifting Law ceases to hold ($R^2 = 0.0002$).

### Key Results

| Property | Binary ($g = 2$) | Ternary ($g = 3$) |
|----------|:---:|:---:|
| $N$ in discriminant | Independent factor | Via $(N - p_k)$ only |
| Static conduit | $\operatorname{rad}_{\mathrm{odd}}(M)$ | Does not exist |
| BSL $R^2$ vs. $\xi$ | 0.997 | **0.0002** |
| PPP–CCC mean gap | > 0.5, stable | 0.10 ± 0.07 |

### Figures

| | |
|---|---|
| ![Geometric](figures/fig_geometric.png) | ![Gap](figures/fig_gap.png) |
| **Figure 1.** Conductor decomposition | **Figure 2.** PPP–CCC gap stability |

## Repository Structure

```
├── README.md
├── LICENSE
├── .gitignore
├── paper/
│   ├── Ternary_Conductor_Boundary.tex
│   └── Ternary_Conductor_Boundary.pdf
├── figures/
│   ├── fig_geometric.pdf / .png
│   └── fig_gap.pdf / .png
└── scripts/
    ├── ternary_geometric.py    # Main analysis (true discriminant)
    └── ternary_corrected.py    # Buggy-vs-correct comparison
```

## Quick Start

```bash
python3 scripts/ternary_geometric.py
```

**Dependencies:** Python ≥ 3.10, NumPy ≥ 1.24, Matplotlib ≥ 3.7.

## Series Context

This is **Paper #11** in the Titan Project conductor rigidity series:

| # | Paper | Repository |
|---|-------|-----------|
| 7 | The Goldbach Mirror | [goldbach-mirror-conductor-rigidity](https://github.com/Ruqing1963/goldbach-mirror-conductor-rigidity) |
| 8 | The Goldbach Mirror II | [goldbach-mirror-II-geometric-foundations](https://github.com/Ruqing1963/goldbach-mirror-II-geometric-foundations) |
| 9 | The Algebraic Vacuum | [goldbach-algebraic-vacuum-zero-ramification](https://github.com/Ruqing1963/goldbach-algebraic-vacuum-zero-ramification) |
| 10 | Dynamic Stability | [goldbach-dynamic-stability-band-shifting](https://github.com/Ruqing1963/goldbach-dynamic-stability-band-shifting) |
| **11** | **Ternary Conductor Boundary** | **this repository** |

## Citation

```bibtex
@misc{chen2026ternaryboundary,
  author       = {Ruqing Chen},
  title        = {The Ternary Conductor Boundary: Why Conductor Rigidity
                  Is Specific to the Binary Goldbach Problem},
  year         = {2026},
  howpublished = {\url{https://github.com/Ruqing1963/goldbach-ternary-conductor-boundary}}
}
```

## License

[MIT](LICENSE)
