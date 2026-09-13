# RACE

## Radial–Angular Correlation Descriptor for Atomistic Machine Learning

**RACE** is a modular and physically interpretable atomistic descriptor based on explicit
radial, angular, and joint radial–angular correlations.

RACE decomposes each local atomic environment into three complementary components:

- **RDF** — radial distribution functions;
- **ADF** — angular distribution functions;
- **F(r, θ)** — an explicit joint radial–angular correlation map.

Chemical identity is introduced through configurable elemental weighting functions
rather than mandatory element-pair channels.

In its single-channel scalar-weighted form, the dimensionality of RACE is independent
of the number of chemical species.

---

# Quick start

Clone the repository:

```bash
git clone https://github.com/tromer-unb/RACE.git
cd RACE
```

Install the required dependencies:

```bash
pip install -r requirements.txt
```

Then run RACE from the **root directory of the repository**:

```bash
python src/RACE_descriptor.py
```

By default, the program reads all `.cif` files located in:

```text
structures/
```

For example:

```text
RACE/
├── structures/
│   ├── 1.cif
│   ├── 2.cif
│   ├── 3.cif
│   ├── 4.cif
│   └── ...
│
├── src/
│   └── RACE_descriptor.py
│
├── results/
│   └── descriptors.csv
│
├── requirements.txt
├── LICENSE
└── README.md
```

The input directory is defined near the beginning of the Python script:

```python
STRUCTURE_FOLDER = "structures"
```

Therefore,

```bash
python src/RACE_descriptor.py
```

reads the structures from:

```text
structures/*.cif
```

and computes one structure-level RACE representation for each input structure.

The recommended output configuration for this repository is:

```python
STRUCTURE_FOLDER = "structures"
OUTPUT_FILE = "results/descriptors.csv"
```

With this configuration, the workflow is:

```text
structures/*.cif
        │
        ▼
src/RACE_descriptor.py
        │
        ▼
local atomic environments
        │
        ▼
RDF + ADF + F(r,θ)
        │
        ▼
mean + standard deviation
        │
        ▼
results/descriptors.csv
```

Each row of the CSV corresponds to one structure, and the first column retains the
original CIF filename.

> **Important:** run the script from the repository root when using relative paths such
> as `STRUCTURE_FOLDER = "structures"`.

---

## Using another structure directory

To calculate RACE descriptors for another dataset, replace:

```python
STRUCTURE_FOLDER = "structures"
```

with another relative directory, for example:

```python
STRUCTURE_FOLDER = "my_structures"
```

giving:

```text
RACE/
├── my_structures/
│   ├── 1.cif
│   ├── 2.cif
│   └── ...
```

An absolute path can also be used:

```python
STRUCTURE_FOLDER = "/path/to/my/structures"
```

The output filename can similarly be modified:

```python
OUTPUT_FILE = "results/my_descriptors.csv"
```

---

## Input structure ordering

The current implementation identifies the leading integer in each filename:

```python
def numeric_sort_key(filename):
    m = re.match(r"(\d+)", filename)
    return int(m.group(1)) if m else float("inf")
```

This ensures that numerically named structures such as:

```text
1.cif
2.cif
3.cif
10.cif
11.cif
100.cif
```

are processed numerically rather than using ordinary lexicographic ordering:

```text
1.cif
10.cif
100.cif
11.cif
2.cif
...
```

For maximum reproducibility, filenames beginning with a numerical index are
recommended.

For example:

```text
1.cif
2.cif
3.cif
...
150.cif
```

Files that do **not** begin with a number receive the same fallback sorting key in the
current implementation. Therefore, numerical prefixes are recommended when the exact
processing order is important.

---

> **Associated publication**
>
> Raphael M. Tromer, Jhionathan de Lima, Cristiano Francisco Woellner,
> Paola Ferreira Barbosa, Roberto Ventura Santos, and Luiz Antonio Ribeiro Junior  
>
> **RACE: A Modular and Interpretable Radial–Angular Correlation Descriptor for
> Atomistic Machine Learning**
>
> *Journal of Chemical Theory and Computation*  
> Accepted for publication.
>
> DOI and final bibliographic information will be added when available.

---

# Overview

RACE represents a local atomic environment using three explicit numerical objects:

1. a radial distribution function (**RDF**);
2. an angular distribution function (**ADF**);
3. a joint radial–angular correlation map **F(r, θ)**.

The descriptor is designed so that these components can be inspected, visualized,
removed, or concatenated independently.

This explicit organization makes it possible to examine which type of structural
information is contributing to a machine-learning representation without reformulating
a complete density expansion.

---

# Local atomic environment

For a central atom \(i\), RACE considers neighboring atoms \(j\) within a cutoff radius

\[
r_{ij} < R_c.
\]

The current implementation uses:

```python
R_CUT = 4.5
```

Neighbor environments are generated using the ASE `NeighborList` implementation.

The code uses:

```python
cutoffs = [R_CUT / 2] * n

nl = NeighborList(
    cutoffs,
    self_interaction=False,
    bothways=True
)
```

Because ASE defines neighbors through overlap of the atomic cutoff spheres, assigning
a radius of \(R_c/2\) to each atom produces a pair cutoff corresponding to \(R_c\).

The descriptor subsequently evaluates the actual distance and applies the explicit
smooth cutoff function, so contributions at or beyond \(R_c\) are discarded.

Periodic-image offsets returned by ASE are included when reconstructing neighboring
atomic positions, allowing periodic CIF structures to be treated correctly.

---

# Smooth cutoff

Neighbor contributions are smoothly attenuated as they approach the cutoff radius.

For

\[
r < R_c,
\]

RACE uses

\[
w_c(r)
=
\frac{1}{2}
\left[
\cos\left(
\frac{\pi r}{R_c}
\right)+1
\right].
\]

For

\[
r \ge R_c,
\]

\[
w_c(r)=0.
\]

In the implementation:

```python
def smooth_cutoff(r):

    if r >= R_CUT:
        return 0.0

    return 0.5 * (np.cos(np.pi * r / R_CUT) + 1.0)
```

The smooth cosine cutoff prevents an abrupt discontinuity when a neighbor approaches
the boundary of the local environment.

---

# Chemical weighting

Chemical information is introduced in the current implementation using atomic numbers.

For atoms \(i\) and \(j\),

\[
w_{ij}
=
w_c(r_{ij})
\sqrt{Z_i Z_j},
\]

where \(Z_i\) and \(Z_j\) are atomic numbers.

The corresponding code is:

```python
def chem_weight(Zi, Zj):

    return np.sqrt(Zi * Zj)
```

followed by:

```python
w = wc * chem_weight(Zi, Zj)
```

Atomic number is intended as a simple and universally available default chemical
encoding.

It should **not** be interpreted as a unique physical metric of chemical similarity.

The RACE formulation can also accommodate other elemental maps, including quantities
such as:

- electronegativity;
- covalent radius;
- valence-electron count;
- multiple chemical channels.

In the associated study, different chemical maps were evaluated explicitly and no
single scalar map was found to be optimal for all datasets.

---

# Radial distribution function

For each central atom, the radial component accumulates Gaussian contributions from
neighboring atoms:

\[
\mathrm{RDF}^{(i)}_p
=
\sum_j
w_{ij}
\exp
\left[
-\frac{(r_{ij}-\mu_p)^2}
{2\sigma_r^2}
\right].
\]

The radial grid is defined by:

```python
R_MIN = 1.0
R_CUT = 4.5

N_R = 16
SIGMA_R = 0.08
```

with centers constructed using:

```python
r_centers = np.linspace(R_MIN, R_CUT, N_R)
```

Thus, 16 radial basis centers are uniformly distributed between 1.0 and 4.5 Å.

---

# Angular distribution function

For each central atom, pairs of neighbors define an angle

\[
\theta_{jil}.
\]

The angular channel is represented using Gaussian functions centered on a fixed angular
grid:

\[
\mathrm{ADF}^{(i)}_q
=
\sum_{j<l}
w_{ij}w_{il}
\exp
\left[
-\frac{(\theta_{jil}-\nu_q)^2}
{2\sigma_\theta^2}
\right].
\]

The implementation uses:

```python
N_A = 16
SIGMA_A_DEG = 7.0
```

and:

```python
a_centers = np.deg2rad(
    np.linspace(0, 180, N_A)
)
```

Therefore, the angular basis contains 16 centers distributed between:

```text
0° and 180°
```

with a Gaussian width of:

```text
7°
```

---

# Joint radial–angular correlation map

The third component of RACE explicitly couples radial and angular information.

For every neighbor pair, the code evaluates both:

- the radial positions of the two neighbors;
- the angle formed around the central atom.

The resulting two-dimensional map is

\[
F^{(i)}(r,\theta).
\]

In discretized form:

\[
F^{(i)}_{p,q}
=
\sum_{j<l}
w_{ij}w_{il}
\left[
K_r(r_{ij};\mu_p)
+
K_r(r_{il};\mu_p)
\right]
K_\theta(\theta_{jil};\nu_q).
\]

This component records **where radial and angular motifs occur together**.

For example, two structures can exhibit similar radial distributions and similar
angle distributions while differing in which particular distances are associated with
particular angles.

The joint map can therefore contain structural information that is not explicit in
separate RDF and ADF histograms.

However, its usefulness is dataset dependent and it should not be interpreted as a
universally beneficial addition.

---

# Default RACE parameters

The default parameter set used in the associated study is:

| Parameter | Symbol | Code variable | Default |
|---|---:|---|---:|
| Cutoff radius | \(R_c\) | `R_CUT` | 4.5 Å |
| Minimum radial center | \(r_{\min}\) | `R_MIN` | 1.0 Å |
| Radial centers | \(N_R\) | `N_R` | 16 |
| Radial Gaussian width | \(\sigma_r\) | `SIGMA_R` | 0.08 Å |
| Angular centers | \(N_A\) | `N_A` | 16 |
| Angular Gaussian width | \(\sigma_\theta\) | `SIGMA_A_DEG` | 7° |
| Local normalization | — | `NORMALIZE_LOCAL` | `True` |
| Structure-level mean normalization | — | `NORMALIZE_GLOBAL` | `True` |
| Numerical tolerance | — | `EPS` | \(10^{-12}\) |
| Structure aggregation | — | — | mean + standard deviation |

The geometric parameters were held fixed across the benchmark datasets rather than
individually optimized for each system.

---

# Discussion of the parameters

## Cutoff radius — `R_CUT`

```python
R_CUT = 4.5
```

The cutoff determines the spatial extent of each local atomic environment.

The value of **4.5 Å** was chosen as a fixed compromise between:

- including the nearest-neighbor environment;
- including additional local coordination shells;
- retaining meaningful local structural information;
- limiting the number of neighbor pairs entering the angular and joint calculations.

This final point is particularly relevant because angular calculations involve pairs
of neighbors around every central atom.

The number of such combinations can increase rapidly when the cutoff is enlarged.

Therefore, increasing `R_CUT` can provide information from more distant coordination
shells but also increases computational cost.

The value 4.5 Å should **not** be interpreted as universally optimal for every material
or target property.

---

## Minimum radial center — `R_MIN`

```python
R_MIN = 1.0
```

`R_MIN` defines the location of the first radial Gaussian center.

It is important to distinguish `R_MIN` from a hard lower-distance cutoff.

The code does **not** discard a neighbor simply because:

\[
r < R_{\min}.
\]

Instead, its distance is evaluated against the Gaussian basis whose first center is at
1.0 Å.

Therefore:

```python
R_MIN = 1.0
```

defines the radial basis grid rather than defining an excluded region of the atomic
environment.

---

## Number of radial centers — `N_R`

```python
N_R = 16
```

`N_R` determines the radial resolution of the descriptor.

A larger value provides a finer sampling of radial space.

However, the joint component has dimensions:

\[
N_R \times N_A.
\]

Therefore increasing `N_R` affects not only the RDF block but also the much larger
joint radial–angular block.

For this reason, increasing radial resolution has a direct dimensional and computational
cost.

---

## Radial Gaussian width — `SIGMA_R`

```python
SIGMA_R = 0.08
```

The radial Gaussian width controls how strongly a neighbor contributes to radial
centers around its actual interatomic distance.

A smaller value produces sharper radial features.

A larger value produces broader and smoother radial features.

Thus there is a trade-off between:

```text
small sigma_R
    → higher localization

large sigma_R
    → stronger smoothing
```

The value used in the study is:

```text
0.08 Å
```

and was kept fixed across the datasets.

---

## Number of angular centers — `N_A`

```python
N_A = 16
```

`N_A` determines the angular resolution.

The centers span:

\[
0^\circ \leq \theta \leq 180^\circ.
\]

Increasing `N_A` gives a finer angular grid but also enlarges the joint block because:

\[
D_{\mathrm{joint}}
\propto N_R N_A.
\]

---

## Angular Gaussian width — `SIGMA_A_DEG`

```python
SIGMA_A_DEG = 7.0
```

This parameter controls the angular smoothing.

A small angular width produces sharper discrimination between neighboring angular
motifs.

A larger width produces smoother distributions and greater overlap between nearby
angular centers.

The default value is:

```text
7°
```

---

# Why the parameters were kept fixed

The default geometric parameters were intentionally kept fixed across the benchmark
systems rather than optimized separately for each dataset.

This avoids introducing dataset-specific descriptor tuning into comparisons between
representations.

The parameter set should therefore be interpreted as a **common reference
configuration**, not as a claim that these values are optimal for all possible
applications.

For a new material class or target property, parameter sensitivity can be investigated
by varying:

```text
R_CUT
N_R
SIGMA_R
N_A
SIGMA_A_DEG
```

while keeping the train/test evaluation protocol fixed.

---

# Local normalization

The current implementation optionally normalizes each local block:

```python
NORMALIZE_LOCAL = True
```

using L1 normalization:

```python
def safe_l1_normalize(x):

    s = np.sum(x)

    if s > EPS:
        return x / s

    return x
```

Therefore, for every atomic environment, the following components are normalized
independently:

```text
RDF
ADF
F(r,θ)
```

provided their total contribution is greater than the numerical tolerance.

This reduces sensitivity to the absolute magnitude of each block and emphasizes its
distribution over radial or angular space.

---

# Structure-level aggregation

After calculating the local descriptor for every atom, RACE converts the set of local
representations into one fixed-size structure-level vector.

For each descriptor component, the code evaluates:

\[
\text{mean}
\]

and:

\[
\text{standard deviation}.
\]

For example:

```python
rdf_mean = rdf_stack.mean(axis=0)
rdf_std  = rdf_stack.std(axis=0)
```

and analogously for the ADF and joint components.

The mean describes the average local environment in the structure.

The standard deviation provides information about the heterogeneity of local
environments.

This is particularly useful when structures contain inequivalent atomic sites or
structural disorder.

---

# Global normalization

The current implementation also contains:

```python
NORMALIZE_GLOBAL = True
```

When enabled, the **mean blocks** are L1-normalized after structure-level aggregation:

```python
rdf_mean
adf_mean
F_mean
```

The corresponding standard-deviation blocks:

```python
rdf_std
adf_std
F_std
```

are not subjected to this second normalization.

---

# Descriptor dimensionality

For one chemical channel, the core RACE descriptor has dimensionality:

\[
D_{\mathrm{RACE}}
=
2N_R + 2N_A + 2N_RN_A.
\]

The factor of two appears because each local component is represented at the structure
level by:

```text
mean + standard deviation
```

Using:

```text
N_R = 16
N_A = 16
```

gives:

### RDF

\[
2N_R = 32
\]

features.

### ADF

\[
2N_A = 32
\]

features.

### Joint map

\[
2N_RN_A
=
2\times16\times16
=
512
\]

features.

Therefore:

\[
D_{\mathrm{RACE}}
=
32+32+512
=
576.
\]

| Component | Dimension |
|---|---:|
| RDF mean + std | 32 |
| ADF mean + std | 32 |
| Joint mean + std | 512 |
| **Core RACE** | **576** |

---

# Output columns

The output CSV begins with:

```text
filename
```

followed by the descriptor components.

## RDF

```text
rdf_mean_00
rdf_std_00
rdf_mean_01
rdf_std_01
...
rdf_mean_15
rdf_std_15
```

## ADF

```text
adf_mean_00
adf_std_00
...
adf_mean_15
adf_std_15
```

## Joint radial–angular map

The \(16\times16\) matrix is flattened before storage.

Therefore:

```text
joint_mean_0000
joint_std_0000

joint_mean_0001
joint_std_0001

...

joint_mean_0255
joint_std_0255
```

represent the 256 joint grid positions for the mean and standard deviation.

---

# Auxiliary structure information

The current Python implementation additionally writes:

```text
n_atoms
z_mean
z_std
z_min
z_max
```

These quantities are useful global metadata, but they are **not part of the
576-dimensional core RACE descriptor** defined above.

Consequently, the current script produces:

```text
576 core RACE features
+ 5 auxiliary structure quantities
= 581 numerical columns
```

plus:

```text
filename
```

for structure identification.

When comparing descriptor dimensionalities with the associated publication, the
five auxiliary quantities should therefore be distinguished from the core RACE
representation.

---

# Vacancy-defective graphene example

The structures currently provided in this repository correspond to
vacancy-defective graphene.

The associated benchmark was constructed from graphene supercells initially containing:

```text
288 atoms
```

Vacancies were introduced at nominal concentrations of:

```text
5%
10%
15%
```

with:

```text
50 configurations per vacancy concentration
```

giving a total of:

```text
150 structures
```

The default RACE parameters were used:

```python
R_CUT = 4.5

R_MIN = 1.0
N_R = 16
SIGMA_R = 0.08

N_A = 16
SIGMA_A_DEG = 7.0
```

The corresponding descriptor table is stored in:

```text
results/descriptors.csv
```

---

# Interpretation for the graphene vacancy dataset

Because the vacancy-defective graphene structures contain only carbon atoms,

\[
Z_\mathrm{C}=6.
\]

Therefore, every C–C chemical pair has:

\[
\sqrt{Z_iZ_j}
=
\sqrt{6\times6}
=
6.
\]

The chemical weight is consequently the same for every carbon-carbon pair.

With block normalization enabled, differences among these descriptors are therefore
primarily associated with structural changes caused by vacancies, including:

- changes in local coordination;
- removal of neighboring atoms;
- modification of radial coordination shells;
- changes in angular distributions;
- changes in radial–angular correlations;
- variation in the distribution of inequivalent local environments.

This makes vacancy-defective graphene a useful example for illustrating the geometric
content of RACE without additional chemical-species complexity.

---

# Physical interpretation

## RDF

The RDF describes where neighboring atoms are located radially around each central
atom.

Sharp features indicate characteristic interatomic-distance environments.

For ordered crystalline systems these can correspond to distinct coordination shells.

---

## ADF

The ADF describes local angular geometry.

For example, characteristic angular environments include approximately:

```text
sp2 carbon       → 120°
sp3 carbon       → 109.5°
```

The angular component can therefore distinguish local coordination motifs that may have
similar radial distributions.

---

## Joint map

The map:

\[
F(r,\theta)
\]

connects radial and angular information explicitly.

Instead of asking only:

```text
Which distances occur?
```

or:

```text
Which angles occur?
```

the joint map also asks:

```text
Which distances occur together with which angles?
```

This provides a directly inspectable representation of radial–angular correlations.

---

# Dimensionality and chemical complexity

A central characteristic of the single-channel scalar-weighted form of RACE is that its
descriptor dimension does not depend on the number of chemical species.

For fixed:

```text
N_R
N_A
```

the dimension remains:

\[
2N_R+2N_A+2N_RN_A.
\]

Introducing more chemical species therefore does not automatically introduce separate
element-pair blocks.

However, this should not be confused with complete dimension independence.

RACE dimensionality still depends on **geometric resolution**, particularly because:

\[
D_{\mathrm{joint}}
=
2N_RN_A.
\]

Increasing either `N_R` or `N_A` therefore increases the descriptor size.

---

# Species-independent dimensionality

For a single scalar chemical map:

```text
1 species
5 species
10 species
20 species
```

can all be represented using the same geometric feature dimension, provided `N_R`
and `N_A` remain unchanged.

This differs organizationally from conventional explicit species-channel
representations.

However, compressed and alchemical formulations of other descriptors can likewise
control chemical-channel growth.

RACE should therefore be described specifically as having **species-independent
descriptor dimensionality in its scalar-weighted form**, rather than as having a
unique or universal chemical-scaling advantage.

---

# Computational considerations

RACE explicitly evaluates neighbor pairs for the ADF and joint map.

Therefore, its computational requirements are influenced by:

- number of atoms;
- local coordination;
- cutoff radius;
- number of radial centers;
- number of angular centers.

A larger cutoff increases the number of neighbors and consequently the number of
neighbor pairs.

Likewise, increasing:

```text
N_R
N_A
```

increases the amount of work needed to construct the radial–angular representation.

The associated benchmarks showed that the current Python implementation has low
peak memory during descriptor generation, but it is not the fastest descriptor
generator among the tested implementations.

The current implementation therefore emphasizes transparency and direct manipulation
of the descriptor components rather than maximum runtime optimization.

---

# Scope and limitations

RACE is a finite-resolution descriptor.

It is not claimed to be a complete or systematically improvable many-body basis.

Important considerations include:

- the optimal descriptor parameters can depend on the application;
- the joint block has substantially more features than the RDF and ADF blocks;
- the incremental value of the joint block is dataset dependent;
- the joint block can contain redundant information in finite-data regimes;
- scalar chemical weighting does not fully represent global composition;
- atomic-number weighting is a convenient default rather than a universal measure
  of chemical similarity;
- structure-level mean/std aggregation can dilute strongly localized environments;
- the present Python implementation prioritizes interpretability over maximum speed.

For heterogeneous materials databases, explicit global composition information may be
useful in addition to local RACE structural features.

---

# Relationship to other atomistic descriptors

RACE is intended as an interpretable and configurable structural representation.

Its contribution is not that radial correlations, angular correlations, or chemical
compression are themselves new concepts.

Instead, RACE organizes these quantities into explicitly accessible components:

```text
RDF
ADF
F(r,θ)
chemical map
```

that can be modified and inspected independently.

RACE should therefore be regarded as a complement to established representations such
as:

- SOAP;
- compressed SOAP;
- ACSF;
- ACE;
- learned graph and equivariant representations.

It is not presented as a universal replacement for these methods.

---

# Modifying the descriptor

The main parameters are grouped near the beginning of:

```text
src/RACE_descriptor.py
```

For example:

```python
R_CUT = 4.5

N_R = 16
SIGMA_R = 0.08
R_MIN = 1.0

N_A = 16
SIGMA_A_DEG = 7.0

NORMALIZE_LOCAL = True
NORMALIZE_GLOBAL = True
```

This makes it straightforward to perform parameter-sensitivity studies.

For example, a researcher interested in a longer structural range could modify:

```python
R_CUT = 6.0
```

while keeping the other parameters fixed.

Alternatively, the radial or angular resolution could be changed using:

```python
N_R = ...
N_A = ...
```

When changing these values, remember that the core descriptor dimensionality becomes:

\[
D_{\mathrm{RACE}}
=
2N_R+2N_A+2N_RN_A.
\]

---

# Requirements

The current implementation uses:

```text
Python
NumPy
pandas
ASE
```

Install the dependencies using:

```bash
pip install -r requirements.txt
```

---

# Repository structure

```text
RACE/
├── src/
│   └── RACE_descriptor.py
│
├── structures/
│   ├── 1.cif
│   ├── 2.cif
│   ├── 3.cif
│   └── ...
│
├── results/
│   └── descriptors.csv
│
├── requirements.txt
├── LICENSE
└── README.md
```

### `src/`

Python implementation of the RACE descriptor.

### `structures/`

Atomic structures used as input.

The current repository provides vacancy-defective graphene structures in CIF format.

### `results/`

Precomputed descriptor results corresponding to the structures distributed with the
repository.

---

# Reproducing the included descriptor table

From the repository root:

```bash
pip install -r requirements.txt
```

make sure the script contains:

```python
STRUCTURE_FOLDER = "structures"
OUTPUT_FILE = "results/descriptors.csv"
```

and run:

```bash
python src/RACE_descriptor.py
```

The program will process the structures and generate:

```text
results/descriptors.csv
```

Each row corresponds to one input CIF file.

---

# Citation

If you use RACE in published work, please cite:

```text
Raphael M. Tromer, Jhionathan de Lima, Cristiano Francisco Woellner,
Paola Ferreira Barbosa, Roberto Ventura Santos, and
Luiz Antonio Ribeiro Junior,

"RACE: A Modular and Interpretable Radial–Angular Correlation
Descriptor for Atomistic Machine Learning",

Journal of Chemical Theory and Computation.
Accepted for publication.
```

The DOI and final bibliographic information will be added when available.

---

# Data availability

The Python implementation of RACE, representative structures, and precomputed
descriptor output are available at:

```text
https://github.com/tromer-unb/RACE
```

---

# License

See the `LICENSE` file distributed with this repository.
