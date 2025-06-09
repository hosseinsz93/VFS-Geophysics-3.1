# VFS-Geophysics 3.1

[![OSTI Record](https://img.shields.io/badge/OSTI-10.15473%2F1997004-blue)](https://www.osti.gov/biblio/1997004)

## Overview

**VFS-Geophysics** is a three-dimensional (3D) incompressible Navier-Stokes solver based on the Curvilinear Immersed Boundary (CURVIB) method. The CURVIB method is a sharp interface type of immersed boundary (IB) method, enabling the simulation of fluid flows in the presence of geometrically complex moving bodies. VFS-Geophysics is designed for high-fidelity simulations of geophysical and engineering flows, including wind and marine hydrokinetic (MHK) turbine simulations and other energy-related applications.

This repository contains the source code and documentation for **VFS-Geophysics 3.1**.

---

## Features

- 3D incompressible Navier-Stokes solver
- Curvilinear Immersed Boundary (CURVIB) method
- Sharp interface immersed boundary for complex/moving geometries
- Suitable for wind turbine, MHK turbine, and energy application simulations
- Developed through extensive research at the Computational Hydrodynamics and Biofluids Laboratory (Prof. Fotis Sotiropoulos)

---

## Documentation & Reference

For detailed documentation and the user manual, please refer to the official OSTI record:

- [VFS-Geophysics on OSTI.gov](https://www.osti.gov/biblio/1997004)

---

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/hosseinsz93/VFS-Geophysics-3.1.git
   cd VFS-Geophysics-3.1
   ```

2. **Dependencies:**
   - [List dependencies here, e.g., MPI, Fortran/C compilers, libraries, etc.]
   - (Please refer to the [manual](https://www.osti.gov/biblio/1997004) for detailed requirements.)

3. **Build:**
   ```bash
   # Example build instructions
   make
   ```

---

## Usage

Usage instructions depend on your simulation goals (e.g., turbine, channel flow, etc.). Please consult the [manual](https://www.osti.gov/biblio/1997004) for full details, input file descriptions, and example setups.

Typical workflow:
1. Prepare input files and geometry definitions.
2. Configure simulation parameters in the provided configuration files.
3. Run the solver:
   ```bash
   ./vfs-geophysics input_file.in
   ```
4. Post-process results using your preferred tools.

---

## Applications

- Wind turbine simulation
- Marine and hydrokinetic (MHK) turbine simulation
- Environmental and geophysical flow modeling
- Energy systems and renewable energy research

---

## Citation

If you use this code in your research, please cite:

```bibtex
@misc{osti_1997004,
  author       = {Khosronejad, Ali and Zhang, Zexia and Yang, Xiaolei and Santoni, Christian and Borazjani, Iman and Calderer, Antoni and Kang, Seokkoo and Gilmanov, Anvar and Le, Trung},
  title        = {Virtual Flow Solver - Geophysics: A 3D Incompressible Navier-Stokes Solver},
  annote       = {Virtual Flow Solver - Geophysics (VFS-Geophysics) is a three-dimensional (3D) incompressible Navier-Stokes solver based on the Curvilinear Immersed Boundary (CURVIB) method. The CURVIB is a sharp interface type of immersed boundary (IB) method that enables the simulation of fluid flows in the presence of geometrically complex moving bodies. The CURVIB method can be applied to wind/MHK turbine simulations and energy applications. VFS-Geophysics is the result of many years of research work by several graduate students, post-docs, and research associates that have been involved in the Computational Hydrodynamics and Biofluids Laboratory directed by Professor Fotis Sotiropoulos. The preparation of the present manual has been supported by the U.S. Department of Energy (DE-EE 0005482).},
  doi          = {10.15473/1997004},
  url          = {https://www.osti.gov/biblio/1997004},
  place        = {United States},
  year         = {2023},
  month        = {07}
}
```

---

## License



---

## Acknowledgments

VFS-Geophysics is the result of many years of research by graduate students, post-docs, and research associates at the Computational Hydrodynamics and Biofluids Laboratory, directed by Professor Fotis Sotiropoulos. The development has been supported by the U.S. Department of Energy (DE-EE 0005482).

---

## Contact

For questions, bug reports, or contributions, please open an issue or contact the repository maintainer.
