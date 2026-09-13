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

Run RACE from the **root directory of the repository**:

```bash
python src/RACE_descriptor.py
```

By default, the script reads all `.cif` files located in:

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

The structure directory is defined near the beginning of the Python script:

```python
STRUCTURE_FOLDER = "structures"
```

Therefore, when the program is executed from the repository root:

```bash
python src/RACE_descriptor.py
```

it searches for:

```text
structures/*.cif
```

Each CIF file is read using ASE, and one structure-level RACE representation is calculated
for each structure.

The input directory can be changed directly in the source code.

For example:

```python
STRUCTURE_FOLDER = "my_structures"
```

or an absolute path can be supplied:

```python
STRUCTURE_FOLDER = "/path/to/my/structures"
```

---

## Output file

The output filename is controlled by:

```python
OUTPUT_FILE = "descriptor_RACE.csv"
```

To write the calculated descriptors directly to the `results/` directory used in this
repository, use:

```python
OUTPUT_FILE = os.path.join("results", "descriptors.csv")
```

and create the output directory automatically before saving:

```python
os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
```

The complete workflow is then:

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

Each row of the resulting CSV corresponds to one input structure.

The first column stores the original CIF filename, allowing each descriptor vector to
be mapped directly to the corresponding atomic structure.

> **Important:** when relative paths such as `STRUCTURE_FOLDER = "structures"` are used,
> run the script from the root directory of the repository.

---

# Input structure ordering

The current implementation determines the processing order using the leading integer
in each filename:

```python
def numeric_sort_key(filename):
    m = re.match(r"(\d+)", filename)
    return int(m.group(1)) if m else float("inf")
```

Therefore:

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

For maximum reproducibility, numerically indexed filenames are recommended:

```text
1.cif
2.cif
3.cif
...
150.cif
```

Files that do not begin with a number receive the same fallback sorting key in the
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
>
> Accepted for publication.
>
> DOI and final bibliographic information will be added when available.

---

# Overview

RACE represents an atomic environment using three explicit numerical components:

1. a radial distribution function (**RDF**);
2. an angular distribution function (**ADF**);
3. a joint radial–angular correlation map **F(r, θ)**.

These components can be inspected, visualized, removed, or concatenated independently.

The descriptor therefore provides an explicit decomposition of structural information
into radial, angular, and coupled radial–angular contributions.

---

# Local atomic environment

For a central atom `i`, RACE considers neighboring atoms `j` whose interatomic distance
satisfies:

```math
r_{ij} < R_c.
```

The default implementation uses:

```python
R_CUT = 4.5
```

where the cutoff is expressed in ångström.

Neighbor environments are generated using the ASE `NeighborList` implementation:

```python
cutoffs = [R_CUT / 2] * n

nl = NeighborList(
    cutoffs,
    self_interaction=False,
    bothways=True
)
```

The actual interatomic distance is subsequently calculated explicitly, and the smooth
cutoff ensures that contributions at or beyond `R_c` vanish.

Periodic-image offsets returned by ASE are included when reconstructing neighboring
atomic positions, allowing periodic structures stored in CIF format to be treated
correctly.

---

# Smooth cutoff

Neighbor contributions are smoothly attenuated as they approach the cutoff radius.

For `r < R_c`:

```math
w_c(r)
=
\frac{1}{2}
\left[
\cos\left(
\frac{\pi r}{R_c}
\right)
+
1
\right].
```

For `r >= R_c`:

```math
w_c(r)=0.
```

The corresponding implementation is:

```python
def smooth_cutoff(r):

    if r >= R_CUT:
        return 0.0

    return 0.5 * (np.cos(np.pi * r / R_CUT) + 1.0)
```

The cosine cutoff avoids an abrupt discontinuity in the descriptor when an atom
approaches the boundary of the local environment.

---

# Chemical weighting

Chemical identity is introduced through a scalar atomic-number map.

For a central atom `i` and a neighbor `j`, the current implementation uses:

```math
w_{ij}
=
w_c(r_{ij})
\sqrt{Z_i Z_j}.
```

Here, `Z_i` and `Z_j` are the corresponding atomic numbers.

The chemical part is implemented as:

```python
def chem_weight(Zi, Zj):

    return np.sqrt(Zi * Zj)
```

and combined with the smooth cutoff through:

```python
w = wc * chem_weight(Zi, Zj)
```

Atomic number should be regarded as a simple and universally available default
chemical map.

It should **not** be interpreted as a unique physical measure of chemical similarity.

The RACE framework can accommodate alternative elemental maps such as:

- electronegativity;
- covalent radius;
- valence-electron count;
- multiple concatenated chemical channels.

The associated study shows that no single scalar chemical map is optimal for every
dataset.

---

# Radial distribution function

For each central atom `i`, the radial component is constructed by accumulating Gaussian
contributions from neighboring atoms:

```math
\mathrm{RDF}^{(i)}_p
=
\sum_j
w_{ij}
\exp
\left[
-\frac{
(r_{ij}-\mu_p)^2
}{
2\sigma_r^2
}
\right].
```

Here:

- `mu_p` is the `p`-th radial center;
- `sigma_r` is the radial Gaussian width;
- `w_ij` contains the smooth cutoff and chemical weighting.

The radial grid is defined by:

```python
R_MIN = 1.0
R_CUT = 4.5

N_R = 16
SIGMA_R = 0.08
```

with centers generated using:

```python
r_centers = np.linspace(R_MIN, R_CUT, N_R)
```

Therefore, the default descriptor uses 16 radial Gaussian centers uniformly distributed
between 1.0 and 4.5 Å.

---

# Angular distribution function

For a central atom `i`, every pair of neighbors `j` and `l` defines an angle:

```math
\theta_{jil}.
```

The angular component is:

```math
\mathrm{ADF}^{(i)}_q
=
\sum_{j<l}
w_{ij}w_{il}
\exp
\left[
-\frac{
(\theta_{jil}-\nu_q)^2
}{
2\sigma_\theta^2
}
\right].
```

Here:

- `nu_q` is the `q`-th angular center;
- `sigma_theta` is the angular Gaussian width.

The implementation uses:

```python
N_A = 16
SIGMA_A_DEG = 7.0
```

with:

```python
a_centers = np.deg2rad(
    np.linspace(0, 180, N_A)
)
```

Thus, 16 angular centers span the interval from 0° to 180°.

---

# Joint radial–angular correlation map

The third RACE component explicitly couples radial and angular information.

The joint representation is:

```math
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
```

The radial kernel is:

```math
K_r(r;\mu_p)
=
\exp
\left[
-\frac{
(r-\mu_p)^2
}{
2\sigma_r^2
}
\right].
```

The angular kernel is:

```math
K_\theta(\theta;\nu_q)
=
\exp
\left[
-\frac{
(\theta-\nu_q)^2
}{
2\sigma_\theta^2
}
\right].
```

The joint map therefore retains information about **which radial and angular motifs occur
together**.

Two structures may exhibit similar RDF and ADF distributions while differing in the
association between particular interatomic distances and particular bond angles.

The joint component explicitly exposes these correlations.

Its usefulness, however, is dataset dependent and it should not be regarded as a
universally beneficial addition.

---

# Default RACE parameters

The default parameter set used in the associated study is:

| Parameter | Symbol | Code variable | Default |
|---|---|---|---:|
| Cutoff radius | `R_c` | `R_CUT` | 4.5 Å |
| Minimum radial center | `r_min` | `R_MIN` | 1.0 Å |
| Number of radial centers | `N_R` | `N_R` | 16 |
| Radial Gaussian width | `sigma_r` | `SIGMA_R` | 0.08 Å |
| Number of angular centers | `N_A` | `N_A` | 16 |
| Angular Gaussian width | `sigma_theta` | `SIGMA_A_DEG` | 7° |
| Local normalization | — | `NORMALIZE_LOCAL` | `True` |
| Structure-level mean normalization | — | `NORMALIZE_GLOBAL` | `True` |
| Numerical tolerance | — | `EPS` | `1e-12` |
| Structure aggregation | — | — | mean + standard deviation |

The geometric parameters were held fixed across the benchmark datasets rather than
optimized independently for every system.

---

# Discussion of the parameters

## Cutoff radius — `R_CUT`

```python
R_CUT = 4.5
```

The cutoff radius determines the spatial extent of each local atomic environment.

The value of **4.5 Å** was selected as a fixed compromise between:

- including the nearest-neighbor environment;
- retaining additional local coordination information;
- capturing structural information beyond only the first bond shell;
- limiting the number of neighbor pairs entering the angular and joint calculations.

This last point is particularly important because the angular and joint components
operate on pairs of neighbors around each central atom.

If a central atom has `m` neighbors, the number of unique neighbor pairs is:

```math
N_{\mathrm{pairs}}
=
\frac{m(m-1)}{2}.
```

Consequently, increasing the cutoff radius can significantly increase the number of
angular evaluations.

A larger cutoff can provide information from more distant coordination shells, but at
additional computational cost.

The value 4.5 Å should therefore **not** be interpreted as universally optimal for every
material or target property.

It was used as a common reference value in the associated study to avoid
dataset-specific descriptor tuning.

---

## Minimum radial center — `R_MIN`

```python
R_MIN = 1.0
```

`R_MIN` specifies the location of the first radial Gaussian center.

It is important to distinguish `R_MIN` from a hard lower-distance cutoff.

The implementation does **not** discard a neighbor simply because its distance is less
than `R_MIN`.

Instead, all valid neighbors inside `R_CUT` contribute to the Gaussian basis.

Thus, `R_MIN` determines the starting point of the radial basis grid, rather than a
forbidden interatomic-distance region.

The radial centers span:

```math
R_{\min}
\leq
\mu_p
\leq
R_c.
```

For the default configuration:

```math
1.0\ \mathrm{\AA}
\leq
\mu_p
\leq
4.5\ \mathrm{\AA}.
```

---

## Number of radial centers — `N_R`

```python
N_R = 16
```

`N_R` controls the radial resolution of the representation.

Increasing `N_R` produces a finer sampling of radial space.

However, its effect on dimensionality extends beyond the RDF because the joint map
contains:

```math
N_R N_A
```

grid positions.

For fixed `N_A`, the joint contribution has dimension:

```math
D_{\mathrm{joint}}
=
2N_RN_A.
```

Therefore, increasing radial resolution also increases the size and computational cost
of the joint radial–angular representation.

---

## Radial Gaussian width — `SIGMA_R`

```python
SIGMA_R = 0.08
```

The radial Gaussian width controls how strongly a neighbor contributes to basis centers
around its actual interatomic distance.

A smaller value produces narrower and more localized radial features.

A larger value produces broader and smoother distributions.

Conceptually:

```text
smaller SIGMA_R
        ↓
sharper radial localization

larger SIGMA_R
        ↓
stronger radial smoothing
```

The default value is:

```math
\sigma_r
=
0.08\ \mathrm{\AA}.
```

This value was kept fixed across the benchmark datasets rather than optimized separately
for every system.

---

## Number of angular centers — `N_A`

```python
N_A = 16
```

`N_A` controls the angular resolution of the descriptor.

The angular centers span:

```math
0^\circ
\leq
\theta
\leq
180^\circ.
```

A larger `N_A` allows finer sampling of angular space.

However, increasing `N_A` also enlarges the joint radial–angular block because:

```math
D_{\mathrm{joint}}
=
2N_RN_A.
```

---

## Angular Gaussian width — `SIGMA_A_DEG`

```python
SIGMA_A_DEG = 7.0
```

The angular Gaussian width determines the smoothing applied to bond-angle information.

The default value is:

```math
\sigma_\theta
=
7^\circ.
```

A smaller angular width provides sharper discrimination between nearby angular motifs.

A larger width increases overlap between neighboring angular basis functions and
produces smoother angular distributions.

---

# Why the parameters were kept fixed

The geometric RACE parameters were intentionally held fixed across the benchmark
datasets in the associated study.

They were not individually optimized for every material system.

This reduces the possibility that descriptor comparisons are dominated by
dataset-specific hyperparameter tuning.

The default configuration should therefore be interpreted as a **common reference
parameterization**, rather than as a claim that these values are universally optimal.

For a new application, parameter sensitivity can be investigated by varying:

```text
R_CUT
R_MIN
N_R
SIGMA_R
N_A
SIGMA_A_DEG
```

while maintaining a controlled training and validation protocol.

---

# Local normalization

The implementation optionally performs local L1 normalization:

```python
NORMALIZE_LOCAL = True
```

using:

```python
def safe_l1_normalize(x):

    s = np.sum(x)

    if s > EPS:
        return x / s

    return x
```

When enabled, each local block is normalized independently:

```text
RDF
ADF
F(r,θ)
```

provided that its sum is greater than the numerical tolerance.

The normalization is:

```math
\mathbf{x}_{\mathrm{norm}}
=
\frac{\mathbf{x}}
{\sum_k x_k}.
```

It is applied when:

```math
\sum_k x_k
>
\varepsilon.
```

The numerical tolerance is:

```math
\varepsilon
=
10^{-12}.
```

Local normalization reduces dependence on the absolute magnitude of each block and
emphasizes the distribution of intensity over radial and angular space.

---

# Structure-level aggregation

After calculating a local descriptor for every atom, the collection of local
representations is converted into one fixed-length structure-level vector.

For each descriptor component, the implementation calculates both the atomic mean and
standard deviation.

For example:

```python
rdf_mean = rdf_stack.mean(axis=0)
rdf_std = rdf_stack.std(axis=0)
```

and analogously for the ADF and joint blocks.

The mean:

```math
\left\langle
\mathbf{x}^{(i)}
\right\rangle_i
```

describes the average local environment.

The standard deviation:

```math
\operatorname{std}_i
\left[
\mathbf{x}^{(i)}
\right]
```

provides information about variation among the local atomic environments.

This is particularly useful for structures containing inequivalent sites, disorder,
defects, surfaces, or other forms of local heterogeneity.

---

# Global normalization

The implementation also contains:

```python
NORMALIZE_GLOBAL = True
```

When enabled, the structure-level **mean blocks** are L1-normalized:

```text
rdf_mean
adf_mean
F_mean
```

The standard-deviation blocks:

```text
rdf_std
adf_std
F_std
```

are retained without this second normalization.

---

# Descriptor dimensionality

For one chemical channel, the dimensionality of the **core RACE representation** is:

```math
D_{\mathrm{RACE}}
=
2N_R
+
2N_A
+
2N_RN_A.
```

The factor of two appears because each component contributes both:

```text
mean
+
standard deviation
```

at the structure level.

Using the default values:

```text
N_R = 16
N_A = 16
```

the contributions are:

## RDF contribution

```math
D_{\mathrm{RDF}}
=
2N_R
=
2(16)
=
32.
```

## ADF contribution

```math
D_{\mathrm{ADF}}
=
2N_A
=
2(16)
=
32.
```

## Joint radial–angular contribution

```math
D_{\mathrm{joint}}
=
2N_RN_A.
```

Therefore:

```math
D_{\mathrm{joint}}
=
2(16)(16)
=
512.
```

Finally:

```math
D_{\mathrm{RACE}}
=
32
+
32
+
512
=
576.
```

| Component | Dimension |
|---|---:|
| RDF mean + standard deviation | 32 |
| ADF mean + standard deviation | 32 |
| Joint mean + standard deviation | 512 |
| **Core RACE** | **576** |

Thus, for `N_R = N_A = 16`:

```math
D_{\mathrm{RACE}} = 576.
```

---

# Output columns

The CSV begins with:

```text
filename
```

followed by the descriptor features.

## RDF features

```text
rdf_mean_00
rdf_std_00
rdf_mean_01
rdf_std_01
...
rdf_mean_15
rdf_std_15
```

The RDF block contains:

```math
2N_R
=
32
```

features.

---

## ADF features

```text
adf_mean_00
adf_std_00
...
adf_mean_15
adf_std_15
```

The ADF block contains:

```math
2N_A
=
32
```

features.

---

## Joint radial–angular features

The `N_R × N_A` joint matrix is flattened before being written to the CSV.

For the default configuration:

```math
N_RN_A
=
16(16)
=
256.
```

Therefore, the mean joint block contains 256 values and the standard-deviation block
contains another 256:

```text
joint_mean_0000
joint_std_0000

joint_mean_0001
joint_std_0001

...

joint_mean_0255
joint_std_0255
```

The complete joint component therefore contains:

```math
2(256)
=
512
```

features.

---

# Auxiliary structure information

The current Python implementation additionally appends five global quantities:

```text
n_atoms
z_mean
z_std
z_min
z_max
```

These quantities can be useful as structural and chemical metadata, but they are
**not part of the 576-dimensional core RACE representation** defined above.

Therefore, the current implementation produces:

```math
576 + 5 = 581
```

numerical quantities:

```text
576 core RACE features
+ 5 auxiliary structure quantities
= 581 numerical quantities
```

The CSV additionally contains:

```text
filename
```

for structure identification.

When comparing descriptor dimensionality with the associated publication, these five
auxiliary quantities should therefore be distinguished from the core RACE descriptor.

---

# Vacancy-defective graphene example

The structures currently provided in this repository correspond to vacancy-defective
graphene.

The benchmark associated with the study was generated from graphene supercells
initially containing:

```text
288 atoms
```

Nominal vacancy concentrations of:

```text
5%
10%
15%
```

were considered.

For each concentration:

```text
50 independent configurations
```

were retained.

The total number of structures was therefore:

```math
3(50)
=
150.
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

No vacancy-specific optimization of the geometric parameters was performed.

The corresponding precomputed descriptor table is available in:

```text
results/descriptors.csv
```

---

# Interpretation for vacancy-defective graphene

The vacancy-defective graphene structures contain only carbon atoms.

For carbon:

```math
Z_{\mathrm{C}}
=
6.
```

Therefore, every C–C pair has the same chemical factor:

```math
\sqrt{Z_i Z_j}
=
\sqrt{6(6)}
=
6.
```

Consequently, the scalar chemical weight is identical for all C–C pairs.

After block normalization, differences among the descriptor vectors are therefore
primarily associated with structural changes generated by the vacancies.

These changes include:

- modification of local coordination;
- removal of neighboring atoms;
- changes in first- and higher-neighbor distance distributions;
- changes in bond-angle distributions;
- changes in radial–angular correlations;
- variation in the distribution of inequivalent local environments.

This makes vacancy-defective graphene a useful example for illustrating the geometric
content of RACE without the additional complexity introduced by multiple chemical
species.

---

# Physical interpretation

## RDF

The RDF describes how neighboring atoms are distributed as a function of distance from
a central atom.

Sharp radial features can reflect characteristic interatomic distances and coordination
shells.

Ordered crystalline environments generally produce more localized radial features than
strongly disordered environments.

---

## ADF

The ADF describes local angular geometry.

Typical carbon environments include approximately:

```text
sp2 carbon → 120°
sp3 carbon → 109.5°
```

Angular information can therefore distinguish environments having similar characteristic
bond lengths but different local coordination geometries.

---

## Joint radial–angular map

The joint map:

```math
F(r,\theta)
```

connects radial and angular information explicitly.

Instead of asking only:

```text
Which distances occur?
```

or:

```text
Which angles occur?
```

the joint representation additionally asks:

```text
Which distances occur together with which angles?
```

The joint map therefore provides a directly inspectable representation of
distance-angle correlations.

---

# Dimensionality and chemical complexity

A key characteristic of single-channel scalar-weighted RACE is that its feature
dimension does not depend on the number of chemical species.

For fixed `N_R` and `N_A`:

```math
D_{\mathrm{RACE}}
=
2N_R
+
2N_A
+
2N_RN_A
```

remains unchanged when the number of distinct chemical elements increases.

For example:

```text
1 species
5 species
10 species
20 species
```

can all be represented with the same single-channel geometric feature dimension when
the RACE grid is unchanged.

This is possible because chemical identity enters through a scalar elemental map rather
than mandatory separate species-pair blocks.

However, species-independent dimensionality should not be confused with dimensionality
that is independent of geometric resolution.

The joint block scales as:

```math
D_{\mathrm{joint}}
=
2N_RN_A.
```

Increasing either `N_R` or `N_A` therefore increases the total number of features.

---

# Species-independent dimensionality

For the default geometric resolution:

```text
N_R = 16
N_A = 16
```

the single-channel core RACE dimension remains:

```text
576 features
```

independently of the number of chemical species represented by the scalar chemical map.

This is an organizational characteristic of RACE.

It should not be interpreted as a claim that chemical compression is unique to RACE.

Compressed SOAP, alchemical mappings, and other low-dimensional chemical
representations can also control species-channel growth.

---

# Computational considerations

RACE explicitly evaluates pairs of neighbors when constructing the ADF and joint
radial–angular components.

The computational cost is therefore influenced by:

- number of atoms;
- number of neighbors per atom;
- cutoff radius;
- radial resolution;
- angular resolution.

If a central atom has `m` neighbors, the number of distinct neighbor pairs is:

```math
N_{\mathrm{pairs}}
=
\frac{m(m-1)}{2}.
```

Increasing the cutoff radius can therefore substantially increase the number of
neighbor-pair operations.

The joint grid additionally scales with:

```math
N_RN_A.
```

The computational benchmarks in the associated study show that the present Python
implementation has low peak memory during descriptor generation.

However, it is not the fastest descriptor generator among the evaluated
implementations.

The current implementation should therefore be understood as emphasizing:

- transparency;
- physical interpretability;
- direct access to individual descriptor components;
- ease of modification;

rather than maximum runtime optimization.

---

# Modifying the descriptor

The principal parameters are grouped near the beginning of:

```text
src/RACE_descriptor.py
```

as:

```python
R_CUT = 4.5

N_R = 16
SIGMA_R = 0.08
R_MIN = 1.0

N_A = 16
SIGMA_A_DEG = 7.0

NORMALIZE_LOCAL = True
NORMALIZE_GLOBAL = True

EPS = 1e-12
```

This makes parameter-sensitivity studies straightforward.

For example, a longer-range local environment could be tested with:

```python
R_CUT = 6.0
```

A finer radial grid could be tested with:

```python
N_R = 24
```

and a finer angular grid with:

```python
N_A = 24
```

When changing the resolution, the core descriptor dimension should be recalculated using:

```math
D_{\mathrm{RACE}}
=
2N_R
+
2N_A
+
2N_RN_A.
```

For example, with:

```text
N_R = 24
N_A = 24
```

the dimension becomes:

```math
D_{\mathrm{RACE}}
=
2(24)
+
2(24)
+
2(24)(24)
=
1248.
```

Thus, increasing geometric resolution can substantially increase the descriptor
dimension because the joint block grows as the product of the radial and angular
resolutions.

---

# Scope and limitations

RACE is a finite-resolution atomistic descriptor.

It is not claimed to be a complete or systematically improvable many-body basis.

Important considerations include:

- descriptor parameters may require assessment for new applications;
- the joint radial–angular block is much larger than the separate RDF and ADF blocks;
- the incremental predictive value of the joint block is dataset dependent;
- the joint block can introduce redundant information in finite-data regimes;
- scalar chemical weighting does not fully represent global composition;
- atomic-number weighting is a convenient default rather than a universal metric of
  chemical similarity;
- structure-level mean/std pooling can dilute strongly localized environments;
- the present Python implementation prioritizes transparency over maximum generation
  speed.

For chemically broad materials databases, explicit global composition information may
be useful in addition to the local RACE structural descriptor.

---

# Relationship to other atomistic descriptors

RACE is intended as a transparent and configurable atomistic representation.

Its contribution is not that radial correlations, angular correlations, or chemical
compression are individually new concepts.

Instead, RACE organizes structural information into explicitly accessible components:

```text
RDF
ADF
F(r,θ)
chemical map
```

that can be inspected, modified, removed, or combined independently.

RACE should therefore be regarded as a complement to established representations such
as:

- SOAP;
- compressed SOAP;
- ACSF;
- ACE;
- learned graph representations;
- learned equivariant representations.

It is not presented as a universal replacement for these methods.

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

Contains the Python implementation of the RACE descriptor.

### `structures/`

Contains input atomic structures in CIF format.

The current repository includes vacancy-defective graphene structures.

### `results/`

Contains precomputed descriptors corresponding to the structures distributed with the
repository.

---

# Reproducing the included descriptor table

Install the dependencies:

```bash
pip install -r requirements.txt
```

Use:

```python
STRUCTURE_FOLDER = "structures"
OUTPUT_FILE = os.path.join("results", "descriptors.csv")
```

and ensure that the output directory is created:

```python
os.makedirs("results", exist_ok=True)
```

Then execute from the repository root:

```bash
python src/RACE_descriptor.py
```

The workflow is:

```text
structures/
     │
     ▼
RACE_descriptor.py
     │
     ├── RDF
     ├── ADF
     └── F(r,θ)
             │
             ▼
     mean + standard deviation
             │
             ▼
results/descriptors.csv
```

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

https://github.com/tromer-unb/RACE

---

# License

See the `LICENSE` file distributed with this repository.
