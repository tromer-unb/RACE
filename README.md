# RACE

## Radial–Angular Correlation Descriptor for Atomistic Machine Learning

**RACE** is a modular and physically interpretable atomistic descriptor based on explicit
radial, angular, and joint radial–angular correlations.

The descriptor decomposes each local atomic environment into three components:

- **RDF** — radial distribution functions;
- **ADF** — angular distribution functions;
- **F(r, θ)** — an explicit joint radial–angular correlation map.

Chemical identity is introduced through configurable elemental weighting functions rather
than mandatory element-pair channels.

In its single-channel scalar-weighted form, the dimensionality of RACE is independent
of the number of chemical species.

---

## Quick start

After cloning the repository and installing the dependencies, RACE can be executed directly
from the root directory:

```bash
python src/RACE_descriptor.py
```

By default, the script reads all CIF structures located in:

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

The input directory is defined near the beginning of `RACE_descriptor.py`:

```python
STRUCTURE_FOLDER = "structures"
```

Therefore, running

```bash
python src/RACE_descriptor.py
```

from the repository root causes the program to search for:

```text
structures/*.cif
```

Each CIF structure is read using ASE, its RACE descriptor is calculated, and one row is
generated for each structure.

Numerically named files are processed according to their leading integer, so that:

```text
1.cif
2.cif
3.cif
10.cif
100.cif
```

are processed in numerical order.

### Using a different structure directory

If the CIF files are stored somewhere else, change:

```python
STRUCTURE_FOLDER = "structures"
```

to the desired directory.

For example:

```python
STRUCTURE_FOLDER = "my_structures"
```

for:

```text
RACE/
├── my_structures/
│   ├── sample_1.cif
│   ├── sample_2.cif
│   └── ...
```

or use an absolute path:

```python
STRUCTURE_FOLDER = "/path/to/my/structures"
```

### Output file

The current script defines the output file through:

```python
OUTPUT_FILE = "descriptor_RACE.csv"
```

If you want the calculated descriptors to be written directly to the `results/`
directory distributed with this repository, use:

```python
OUTPUT_FILE = "results/descriptors.csv"
```

The recommended configuration is therefore:

```python
STRUCTURE_FOLDER = "structures"
OUTPUT_FILE = "results/descriptors.csv"
```

With this configuration, running:

```bash
python src/RACE_descriptor.py
```

performs the complete workflow:

```text
structures/*.cif
        │
        ▼
src/RACE_descriptor.py
        │
        ▼
RACE descriptor calculation
        │
        ▼
results/descriptors.csv
```

The first column of the resulting CSV contains the original CIF filename, allowing every
descriptor vector to be associated with its corresponding atomic structure.

> **Important:** run the command from the root directory of the repository. Since
> `STRUCTURE_FOLDER = "structures"` is a relative path, running the script from inside
> `src/` would make Python search for `src/structures/` instead.

---

> **Associated paper**
>
> R. M. Tromer, J. de Lima, C. F. Woellner, P. F. Barbosa,
> R. V. Santos, and L. A. Ribeiro Junior  
> **RACE: A Modular and Interpretable Radial–Angular Correlation Descriptor for Atomistic Machine Learning**  
> *Journal of Chemical Theory and Computation*  
> Accepted for publication.

---
