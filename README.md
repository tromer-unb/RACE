# RACE Descriptor

RACE (Radial–Angular Correlation Expansion) is an atomistic descriptor with species-independent dimensionality, designed for scalable and memory-efficient representations of chemically diverse datasets with many distinct elements.

---

## Overview

RACE encodes atomic environments using:
- Radial distribution functions (RDF)
- Angular distribution functions (ADF)
- Joint radial–angular correlations

Unlike descriptors such as SOAP, the dimensionality of RACE does **not depend on the number of chemical species**, making it particularly suitable for:
- high-entropy materials  
- multi-component systems  
- large materials databases  
- defective systems such as vacancy-containing structures  

---

## Installation

Install dependencies:

```bash
pip install -r requirements.txt
