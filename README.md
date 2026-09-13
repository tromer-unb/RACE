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

Therefore, when the program is executed from the repository root,

```bash
python src/RACE_descriptor.py
```

it searches for:

```text
structures/*.cif
```

Each CIF file is read using ASE and one structure-level RACE representation is
calculated for each structure.

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

In the current implementation, the output filename is controlled by:

```python
OUTPUT_FILE = "descriptor_RACE.csv"
```

With this setting, the CSV is written to the repository root.

To reproduce the organization used in this repository, the recommended setting is:

```python
OUTPUT_FILE = os.path.join("results", "descriptors.csv")
```

and, before writing the CSV, the output directory can be created automatically using:

```python
os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
```

The resulting workflow is:

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

Each row in the CSV corresponds to one structure.

The first column stores the original input filename, allowing every descriptor vector
to be mapped back to its corresponding atomic structure.

> **Important:** when relative paths such as `STRUCTURE_FOLDER = "structures"` are
> used, run the script from the root directory of the repository.

---

# Input structure ordering

The current implementation determines the processing order from the leading integer
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

are processed numerically rather than using standard lexicographic ordering:

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
current implementation. Therefore, their relative ordering is not explicitly defined
by `numeric_sort_key`.

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

For a central atom $i$, RACE considers neighboring atoms $j$ whose interatomic
distance satisfies:

$$
r_{ij} < R_c.
$$

The default implementation uses:

```python
R_CUT = 4.5
```

where the cutoff is expressed in ångström.

Neighbor environments are generated using the ASE `NeighborList` class.

The code uses:

```python
cutoffs = [R_CUT / 2] * n

nl = NeighborList(
    cutoffs,
    self_interaction=False,
    bothways=True
)
```

The actual interatomic distance is subsequently evaluated explicitly, and the smooth
cutoff function ensures that contributions with $r \geq R_c$ vanish.

Periodic-image offsets returned by ASE are included when reconstructing neighbor
positions, which allows periodic structures stored in CIF format to be treated.

---

# Smooth cutoff

Neighbor contributions are smoothly attenuated as they approach the cutoff radius.

For $r < R_c$:

$$
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
$$

For $r \geq R_c$:

$$
w_c(r)=0.
$$

The corresponding implementation is:

```python
def smooth_cutoff(r):

    if r >= R_CUT:
        return 0.0

    return 0.5 * (np.cos(np.pi * r / R_CUT) + 1.0)
```

The cosine cutoff avoids an abrupt discontinuity in the descriptor as an atom approaches
the boundary of the local environment.

---

# Chemical weighting

Chemical identity is introduced through a scalar atomic-number map.

For a central atom $i$ and neighbor $j$, the current implementation uses:

$$
w_{ij}
=
w_c(r_{ij})
\sqrt{Z_i Z_j},
$$

where $Z_i$ and $Z_j$ are the corresponding atomic numbers.

The chemical part is implemented as:

```python
def chem_weight(Zi, Zj):

    return np.sqrt(Zi * Zj)
```

and combined with the distance-dependent cutoff through:

```python
w = wc * chem_weight(Zi, Zj)
```

Atomic number should be regarded as a simple and universally available default
chemical map.

It is **not** intended to represent a unique physical metric of chemical similarity.

The RACE framework can also accommodate alternative elemental maps such as:

- electronegativity;
- covalent radius;
- valence-electron count;
- multiple concatenated chemical channels.

The associated study shows that no single scalar chemical map is optimal for every
dataset.

---

# Radial distribution function

For each central atom $i$, the radial component is constructed by accumulating Gaussian
contributions from neighboring atoms:

$$
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
$$

Here:

- $\mu_p$ is the $p$-th radial center;
- $\sigma_r$ is the Gaussian radial width;
- $w_{ij}$ contains the smooth cutoff and chemical weighting.

The radial grid is defined by:

```python
R_MIN = 1.0
R_CUT = 4.5

N_R = 16
SIGMA_R = 0.08
```

and the centers are generated as:

```python
r_centers = np.linspace(R_MIN, R_CUT, N_R)
```

Therefore, the default descriptor uses 16 radial Gaussian centers uniformly distributed
between 1.0 and 4.5 Å.

---

# Angular distribution function

For a central atom $i$, each pair of neighbors $j$ and $l$ defines an angle:

$$
\theta_{jil}.
$$

The angular component is:

$$
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
$$

Here:

- $\nu_q$ is the $q$-th angular center;
- $\sigma_\theta$ is the angular Gaussian width.

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

$$
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
$$

The radial kernel is:

$$
K_r(r;\mu_p)
=
\exp
\left[
-\frac{
(r-\mu_p)^2
}{
2\sigma_r^2
}
\right],
$$

and the angular kernel is:

$$
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
$$

The joint map therefore retains information about **which radial and angular motifs
occur together**.

For example, two structures can have similar independent RDF and ADF distributions
while differing in the association between particular bond lengths and particular
bond angles.

The joint component provides an explicit numerical representation of these correlations.

Its usefulness, however, is dataset dependent and it should not be regarded as a
universally beneficial addition.

---

# Default RACE parameters

The default parameter set used in the associated study is:

| Parameter | Symbol | Code variable | Default |
|---|---|---|---:|
| Cutoff radius | $R_c$ | `R_CUT` | 4.5 Å |
| Minimum radial center | $r_{\min}$ | `R_MIN` | 1.0 Å |
| Number of radial centers | $N_R$ | `N_R` | 16 |
| Radial Gaussian width | $\sigma_r$ | `SIGMA_R` | 0.08 Å |
| Number of angular centers | $N_A$ | `N_A` | 16 |
| Angular Gaussian width | $\sigma_\theta$ | `SIGMA_A_DEG` | 7° |
| Local normalization | — | `NORMALIZE_LOCAL` | `True` |
| Structure-level mean normalization | — | `NORMALIZE_GLOBAL` | `True` |
| Numerical tolerance | — | `EPS` | $10^{-12}$ |
| Structure aggregation | — | — | mean + standard deviation |

The geometric parameters were kept fixed across the benchmark datasets rather than
being independently optimized for every system.

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

The final point is particularly relevant for RACE because the angular and joint
components operate on pairs of neighbors around each central atom.

If the number of neighbors of an atom is $m$, the number of unique neighbor pairs is:

$$
\frac{m(m-1)}{2}.
$$

Consequently, increasing the cutoff can significantly increase the number of angular
evaluations.

A larger cutoff may include useful information from more distant coordination shells,
but this comes at additional computational cost.

The value 4.5 Å should therefore **not** be interpreted as universally optimal for every
material or target property.

It was used as a common reference value in the associated study in order to avoid
dataset-specific descriptor tuning.

---

## Minimum radial center — `R_MIN`

```python
R_MIN = 1.0
```

`R_MIN` defines the position of the first radial Gaussian center.

It is important to distinguish `R_MIN` from a hard lower cutoff.

The implementation does **not** reject a neighbor simply because:

$$
r < R_{\min}.
$$

Instead, all valid neighbors inside $R_c$ contribute to the radial Gaussian functions.

`R_MIN` therefore controls the starting position of the radial basis grid rather than
defining a forbidden interatomic-distance region.

The radial centers are:

$$
\mu_p
\in
[R_{\min},R_c].
$$

For the default configuration:

$$
1.0\ \mathrm{\AA}
\leq
\mu_p
\leq
4.5\ \mathrm{\AA}.
$$

---

## Number of radial centers — `N_R`

```python
N_R = 16
```

`N_R` controls the radial resolution of the descriptor.

Increasing `N_R` produces a finer sampling of radial space.

However, the effect on descriptor dimensionality is greater than simply adding more RDF
features because the joint component contains:

$$
N_RN_A
$$

grid positions.

Therefore, increasing radial resolution also increases the size of the joint
radial–angular representation.

For fixed $N_A$:

$$
D_{\mathrm{joint}}
=
2N_RN_A.
$$

Thus, radial resolution is directly connected to both descriptor size and computational
cost.

---

## Radial Gaussian width — `SIGMA_R`

```python
SIGMA_R = 0.08
```

The radial Gaussian width controls the degree of smoothing applied to interatomic
distances.

A smaller $\sigma_r$ gives narrower Gaussian functions and therefore sharper radial
localization.

A larger $\sigma_r$ gives broader features and stronger smoothing.

Conceptually:

```text
smaller sigma_r
        ↓
sharper radial resolution

larger sigma_r
        ↓
broader and smoother radial distributions
```

The default value used in the study is:

$$
\sigma_r = 0.08\ \mathrm{\AA}.
$$

This value was kept fixed across the benchmark systems rather than optimized separately
for each dataset.

---

## Number of angular centers — `N_A`

```python
N_A = 16
```

`N_A` controls the angular resolution.

The angular centers span:

$$
0^\circ
\leq
\theta
\leq
180^\circ.
$$

A larger `N_A` allows finer sampling of the angular distribution.

However, because the joint block contains $N_RN_A$ positions, increasing `N_A` also
increases the dimensionality of the joint representation.

For fixed $N_R$:

$$
D_{\mathrm{joint}}
=
2N_RN_A.
$$

---

## Angular Gaussian width — `SIGMA_A_DEG`

```python
SIGMA_A_DEG = 7.0
```

The angular Gaussian width determines the smoothing applied to bond-angle information.

The default value is:

$$
\sigma_\theta = 7^\circ.
$$

A smaller angular width gives sharper discrimination between nearby angular motifs.

A larger width increases overlap between neighboring angular basis functions and
produces smoother angular distributions.

---

# Why the parameters were kept fixed

The geometric RACE parameters were intentionally held fixed across the benchmark
datasets in the associated study.

They were not individually optimized for each material system.

This avoids introducing dataset-specific descriptor tuning when comparing RACE with
alternative atomistic representations.

The default configuration should therefore be interpreted as a **common reference
parameterization**, rather than as a claim that these values are universally optimal.

For a new application, parameter sensitivity may be investigated by varying:

```text
R_CUT
R_MIN
N_R
SIGMA_R
N_A
SIGMA_A_DEG
```

while maintaining a controlled model-training and validation protocol.

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

When enabled, the three blocks of each local atomic environment are normalized
independently:

```text
RDF
ADF
F(r,θ)
```

provided that their total contribution is larger than the numerical tolerance.

The normalization transforms a vector $\mathbf{x}$ according to:

$$
\mathbf{x}_{\mathrm{norm}}
=
\frac{\mathbf{x}}
{\sum_k x_k},
$$

when:

$$
\sum_k x_k > \varepsilon.
$$

The numerical tolerance used by the code is:

$$
\varepsilon = 10^{-12}.
$$

Local normalization reduces dependence on the absolute magnitude of each block and
emphasizes how its intensity is distributed over radial and angular space.

---

# Structure-level aggregation

After computing a local descriptor for every atom, the local representations are
converted into a fixed-length structure-level vector.

For each descriptor component, the implementation calculates both the atomic mean and
standard deviation.

For the RDF:

```python
rdf_mean = rdf_stack.mean(axis=0)
rdf_std = rdf_stack.std(axis=0)
```

and similarly for the ADF and joint blocks.

Conceptually:

$$
\left\langle
\mathbf{x}^{(i)}
\right\rangle_i
$$

describes the average local environment, while:

$$
\operatorname{std}_i
\left[
\mathbf{x}^{(i)}
\right]
$$

describes the variation among atomic environments in the structure.

The standard-deviation component is particularly useful for distinguishing structures
containing heterogeneous or inequivalent atomic environments.

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

$$
D_{\mathrm{RACE}}
=
2N_R
+
2N_A
+
2N_RN_A.
$$

The factor of two appears because each local descriptor component contributes both:

```text
mean
+
standard deviation
```

at the structure level.

Using the default values:

$$
N_R = 16
$$

and:

$$
N_A = 16,
$$

the individual contributions are as follows.

## RDF contribution

$$
D_{\mathrm{RDF}}
=
2N_R
=
2(16)
=
32.
$$

## ADF contribution

$$
D_{\mathrm{ADF}}
=
2N_A
=
2(16)
=
32.
$$

## Joint radial–angular contribution

$$
D_{\mathrm{joint}}
=
2N_RN_A.
$$

Therefore:

$$
D_{\mathrm{joint}}
=
2(16)(16)
=
512.
$$

Finally:

$$
D_{\mathrm{RACE}}
=
32
+
32
+
512
=
576.
$$

| Component | Dimension |
|---|---:|
| RDF mean + standard deviation | 32 |
| ADF mean + standard deviation | 32 |
| Joint mean + standard deviation | 512 |
| **Core RACE** | **576** |

Thus, for $N_R=N_A=16$:

$$
\boxed{
D_{\mathrm{RACE}}=576
}
$$

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

There are:

$$
2N_R = 32
$$

RDF columns.

---

## ADF features

```text
adf_mean_00
adf_std_00
...
adf_mean_15
adf_std_15
```

There are:

$$
2N_A = 32
$$

ADF columns.

---

## Joint radial–angular features

The $N_R \times N_A$ joint matrix is flattened before being written to the CSV.

For the default configuration:

$$
N_RN_A
=
16\times16
=
256.
$$

Therefore the mean block contains 256 values and the standard-deviation block contains
another 256 values:

```text
joint_mean_0000
joint_std_0000

joint_mean_0001
joint_std_0001

...

joint_mean_0255
joint_std_0255
```

The joint component therefore contains:

$$
2(256)=512
$$

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

These variables are useful structural and chemical metadata, but they are **not part of
the 576-dimensional core RACE representation** defined above.

Therefore, the current implementation produces:

$$
576 + 5 = 581
$$

numerical columns, consisting of:

```text
576 core RACE features
+ 5 auxiliary structure quantities
= 581 numerical quantities
```

The CSV additionally contains the non-numerical identification column:

```text
filename
```

When comparing descriptor dimensionality with the associated publication, the five
auxiliary quantities should therefore be distinguished from the core RACE descriptor.

---

# Vacancy-defective graphene example

The structures currently provided in this repository correspond to vacancy-defective
graphene.

The benchmark associated with the study was generated from graphene supercells
initially containing:

$$
288\ \text{atoms}.
$$

Nominal vacancy concentrations of:

```text
5%
10%
15%
```

were considered.

For each vacancy concentration:

```text
50 independent configurations
```

were retained, giving:

$$
3\times50 = 150
$$

structures in total.

The default RACE parameters were retained for this dataset:

```python
R_CUT = 4.5

R_MIN = 1.0
N_R = 16
SIGMA_R = 0.08

N_A = 16
SIGMA_A_DEG = 7.0
```

The corresponding precomputed descriptor table is stored in:

```text
results/descriptors.csv
```

---

# Interpretation for vacancy-defective graphene

The vacancy-defective graphene structures contain only carbon atoms.

For carbon:

$$
Z_{\mathrm{C}}=6.
$$

Therefore, for every C–C pair:

$$
\sqrt{Z_iZ_j}
=
\sqrt{6\times6}
=
6.
$$

The scalar chemical factor is consequently identical for all C–C pairs.

After block normalization, differences among the RACE vectors are therefore dominated
by changes in the structural environment rather than by changes in chemical identity.

Vacancies can modify:

- local coordination numbers;
- first- and higher-neighbor distance distributions;
- bond-angle distributions;
- radial–angular correlations;
- the distribution of inequivalent local atomic environments.

This makes vacancy-defective graphene a useful example for illustrating the geometric
information encoded by RACE without the additional complexity of multiple chemical
species.

---

# Physical interpretation

## RDF

The RDF describes how neighboring atoms are distributed as a function of distance from
a central atom.

Sharp radial features can reflect characteristic bond lengths and coordination shells.

For an ordered crystalline environment, the radial distribution generally exhibits more
localized features than in a strongly disordered structure.

---

## ADF

The ADF describes the local angular geometry.

For example, characteristic carbon environments include approximately:

```text
sp2 carbon → 120°
sp3 carbon → 109.5°
```

Thus, angular information can distinguish environments that may contain similar bond
lengths but different coordination geometry.

---

## Joint radial–angular map

The joint map:

$$
F(r,\theta)
$$

connects radial and angular information explicitly.

Rather than asking only:

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

The joint map therefore makes specific distance-angle relationships directly
inspectable.

---

# Dimensionality and chemical complexity

A key characteristic of single-channel scalar-weighted RACE is that its feature
dimension is independent of the number of chemical species.

For fixed $N_R$ and $N_A$:

$$
D_{\mathrm{RACE}}
=
2N_R
+
2N_A
+
2N_RN_A
$$

remains unchanged as the number of distinct elements increases.

For example, with the same geometric parameters:

```text
1 species
5 species
10 species
20 species
```

all use the same single-channel RACE dimensionality.

This is possible because chemical identity enters through a scalar elemental map rather
than through separate mandatory species-pair feature blocks.

However, species-independent dimensionality should not be confused with dimensionality
that is independent of descriptor resolution.

The joint block scales as:

$$
D_{\mathrm{joint}}
=
2N_RN_A.
$$

Therefore, increasing either the radial or angular resolution increases the total number
of features.

---

# Species-independent dimensionality

For the default geometric resolution:

$$
N_R=N_A=16,
$$

the scalar-weighted core RACE dimension remains:

$$
D_{\mathrm{RACE}}=576
$$

regardless of whether the dataset contains one or several chemical species.

This is an organizational characteristic of the scalar-weighted RACE representation.

It should not be interpreted as a claim that chemical compression is unique to RACE.

Compressed SOAP and other alchemical or low-dimensional chemical representations can
also reduce species-channel growth.

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

If a central atom has $m$ neighbors, the number of distinct neighbor pairs is:

$$
N_{\mathrm{pairs}}
=
\frac{m(m-1)}{2}.
$$

Thus, increasing the cutoff radius can substantially increase the amount of work required
for the angular and joint terms.

The joint grid additionally scales with:

$$
N_RN_A.
$$

The associated computational benchmarks showed that the present Python implementation
has low peak memory during descriptor generation.

However, it was not the fastest descriptor generator among the evaluated
implementations.

The current implementation should therefore be understood as emphasizing:

- transparency;
- direct access to the descriptor components;
- interpretability;
- ease of modification;

rather than maximum runtime optimization.

---

# Modifying the descriptor

The main parameters are grouped near the beginning of:

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

This organization makes it straightforward to perform parameter-sensitivity studies.

For example, a larger local environment could be investigated using:

```python
R_CUT = 6.0
```

A finer radial representation could be generated using:

```python
N_R = 24
```

and a finer angular grid using:

```python
N_A = 24
```

When modifying the grid resolution, the descriptor dimensionality must be recalculated
using:

$$
D_{\mathrm{RACE}}
=
2N_R
+
2N_A
+
2N_RN_A.
$$

For example, if:

$$
N_R=N_A=24,
$$

then:

$$
D_{\mathrm{RACE}}
=
2(24)
+
2(24)
+
2(24)(24)
=
1248.
$$

Thus, increasing geometric resolution can substantially increase the descriptor
dimension because of the joint block.

---

# Scope and limitations

RACE is a finite-resolution atomistic descriptor.

It is not claimed to be a complete or systematically improvable many-body basis.

Important considerations include:

- descriptor parameters may require assessment for new applications;
- the joint radial–angular block is much larger than the independent RDF and ADF
  blocks;
- the incremental predictive value of the joint block is dataset dependent;
- the joint block can introduce redundant information in finite-data regimes;
- scalar chemical weighting does not fully encode global composition;
- atomic-number weighting is a convenient default rather than a universal chemical
  similarity metric;
- structure-level mean/std pooling can dilute strongly localized environments;
- the current Python implementation prioritizes transparency over maximum generation
  speed.

For chemically broad materials databases, explicit global composition information may
be useful in addition to the local RACE descriptor.

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

It is not presented as a universal replacement for these approaches.

---

# Requirements

The current implementation uses:

```text
Python
NumPy
pandas
ASE
```

Install the dependencies with:

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

Contains precomputed descriptor output corresponding to the structures distributed
with the repository.

---

# Reproducing the included descriptor table

Install the dependencies:

```bash
pip install -r requirements.txt
```

For the repository organization shown above, set:

```python
STRUCTURE_FOLDER = "structures"
OUTPUT_FILE = os.path.join("results", "descriptors.csv")
```

and make sure the output directory exists:

```python
os.makedirs("results", exist_ok=True)
```

Then execute from the repository root:

```bash
python src/RACE_descriptor.py
```

The calculation follows:

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
