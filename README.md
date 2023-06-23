# MagmaChamberModel

Thermo mechanical open system magma chamber model includes mixed CO2 and H2O magmatic volatile phase (MVP)

This magma chamber evolution box model is an expansion of the model presented in Degruyter and Huber (2014)and is written for Matlab.
## Model Overview
  - Governing equations conserve total mass, H2O mass, CO2 mass, and enthalpy
  - Tracks how MVP volume fraction, MVP composition (CO2 mole fraction), pressure, and temperature change over time
  - Options for simulating basaltic or rhyolitic magma chambers (including melting curve and volatile solubility relationships for basaltic or rhyolitic compositions)
  - Chamber undergoes conductive heat loss to crust, magma recharge, and eruptions (if crital overpressure is reached)
  - Chamber hosted in viscoelastic crust

## Inputs
Example runcode files included in usr directory. Specify  magma composition, total magma co2 and h2o content, inital chamber volume, and recharge rate

## Outputs 
- Run results saved as .mat file
- Overview timeseries for run
  
