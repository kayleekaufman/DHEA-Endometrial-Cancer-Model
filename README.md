# DHEA-driven endometrial tumour growth model

MATLAB code for: Basson, H. & Kaufman, K. (2026). *Non-Linear Growth Dynamics of Endometrial Adenocarcinoma Under Exogenous Influx of Dehydroepiandrosterone (DHEA)*. Biomathematics 374, Department of Mathematical Sciences, Stellenbosch University.

## Model

A three-equation ODE model linking an oral DHEA dose *S* to local oestradiol and tumour cell growth:

```
dD_S/dt = q*S - deltaD*D_S                     (supplement-derived tissue DHEA)
dE/dt   = Phi*(Db + D_S) - deltaE*E            (local oestradiol)
dN/dt   = r*H(E)*N*(1 - N/L) - dT*N            (tumour cells)
H(E)    = E^n / (K^n + E^n)                    (Hill receptor activation)
```

## Requirements

MATLAB R2019b or later. No additional toolboxes are needed (uses `ode45` and `fzero`).

## How to run

1. Download this repository (green **Code** button → **Download ZIP**) and unzip it.
2. Open MATLAB and set the current folder to the unzipped folder.
3. Open `BassonKaufmanDHEAEndometrialCancer.m` and press **Run**.

The script prints E_crit (≈ 0.520 nM) and S_crit (≈ 189.7 mg/day) to the Command Window and saves seven figures as `Fig1.png` – `Fig7.png`.

## What each section produces

| Code section | Output | Report |
|---|---|---|
| 1. Parameters | All parameter values with units and sources | Table 1 |
| 2. Model and thresholds | E_crit and S_crit | Sections 3.3–3.4 |
| Figure 1 | Hormone dynamics after a dose | Figure 2 |
| Figure 2 | Tumour trajectories at 0.8–1.4 × S_crit | Figure 3 |
| Figure 3 | Graphical check of E_crit and S_crit | Figure 4 |
| Figure 4 | Bifurcation diagram | Figure 5 |
| Figures 5–6 | Sensitivity of S_crit to Φ, K and d | Figures 6–7 |
| 8. Section 3 checks | Equilibria, eigenvalues and S_crit verified numerically (Figure 7) | Section 3 |
