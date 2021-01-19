# DataDrivenTurbulenceModels
OpenFOAM's Data-Driven Turbulence Models.

These models were developed in OpenFOAM-4.x and OpenFOAM-7

Models are used to correct RANS turbulent flows by using Machine Learning predicted quantities or direct injection of DNS
or other high-fidelity simulations'data (e.g. LES).

Models' source terms:
- **RStress** - Reynolds Stress Tensor R
  - Directly injects the deviatoric part of ***R*** into the momentum balance
- **tForce** - Reynolds Force Vector t
  - Based on the work by Cruz et al. (2019) available on https://doi.org/10.1016/j.compfluid.2019.104258
  - Directly injects the vector ***t*** into the momentum balance
- **nutGammaR** - Symmetric tensor Gamma
  - Injects the tensor ***Gamma*** into a Reynolds Stress Transport Equation, which produces an ***R*** whose deviatoric part is injected into the momentum balance. A turbulent viscosity nut from the baseline RANS simulation is required.
- **nutGammat** - gamma vector
  - Injects the vector ***gamma*** into a Reynolds Force Vector Transport Equation, which produces a ***t*** that is injected into the momentum balance. A turbulent viscosity nut from the baseline RANS simulation is required.

Models were constructed using OF's *ShihQuadraticKE* turbulence model.

To include the library in your OF installation use the command:
1) Pull the repository, preferably into your $WM_PROJECT_USER_DIR
2) Go to the directory where you copied the repository's content
3) Use the command `wmake libso`
4) To use the models it's necessary to include the line below into your simulation's controlDict:
  `libs ("libmyDataDrivenRASModels.so");`
5) Change the turbulence model in `constant/turbulenceProperties` into one of the 4 models of this library.
  
# References
- Cruz, Matheus A., et al. "The use of the Reynolds force vector in a physics informed machine learning approach for predictive turbulence modeling." Computers & Fluids 192 (2019): 104258. https://doi.org/10.1016/j.compfluid.2019.104258
