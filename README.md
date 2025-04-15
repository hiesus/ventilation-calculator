# Bathroom Ventilation System Analysis

This repository contains tools for analyzing and designing ventilation systems for multi-floor facility bathrooms. The code calculates required airflow rates, duct velocities, pressure losses, and recommends appropriate centrifugal fans based on the system requirements.

## Project Overview

The ventilation system analysis covers:

- Calculation of required air changes per hour (ACH) for each bathroom
- Analysis of airflow velocities through grilles and ducts
- Calculation of pressure losses throughout the duct system
- Recommendations for appropriate centrifugal fans
- Visualization of the ventilation system layout
- Interactive adjustment of parameters to optimize the system

## Installation

### Local Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/ventilation-analysis.git
cd ventilation-analysis
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

### Running on Binder

You can run this project directly in your browser without installing anything using Binder:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/yourusername/ventilation-analysis/main?filepath=ventilation_analysis.ipynb)

## Usage

### Jupyter Notebook

The main analysis is contained in the `ventilation_analysis.ipynb` notebook. Open it with Jupyter:

```bash
jupyter notebook ventilation_analysis.ipynb
```

The notebook provides an interactive interface to:
- Run the ventilation system analysis
- View improvement suggestions
- See fan recommendations
- Visualize the system layout
- Interactively adjust parameters

### Python Script

You can also run the analysis directly using the Python script:

```bash
python ventilation_calculator.py
```

This will:
1. Analyze the ventilation system under the worst-case climate conditions
2. Generate suggestions for improvements
3. Recommend appropriate fans
4. Save results to CSV files
5. Generate a system diagram

## System Description

The ventilation system analyzed in this project consists of:

### System 1: Ground Floor & Cafeteria
- **PB_Bano_C**: 25m³, 2 grilles (4"x11")
- **PB_Bano_D**: 18m³, 3 grilles (6"x11")
- **Cafetin_Bano_C**: 6.5m³, 1 grille (7.5"x10")
- **Cafetin_Bano_D**: 9m³, 2 grilles (7.5"x11")

### System 2: Floor 1 & Floor 2
- **Piso1_Bano_C**: 10m³, 3 grilles (6"x14")
- **Piso1_Bano_D**: 10m³, 3 grilles (6"x14")
- **Piso2_Bano_C**: 10m³, 3 grilles (4"x11")
- **Piso2_Bano_D**: 10m³, 3 grilles (4"x11")
- **Ducto_Principal_Piso2**: Main duct (11"x11")

## Climate Conditions

The analysis considers the following climate conditions:
- **Average**: 25°C, 60% RH
- **Extreme Hot & Humid**: 28°C, 93% RH
- **Extreme Cold & Dry**: 16°C, 30% RH

The facility is located at an altitude of 850m above sea level.

## Output Files

The analysis generates the following output files:
- `ventilation_results.csv`: Detailed results for each bathroom
- `fan_recommendations.csv`: Fan recommendations for each system
- `ventilation_system_diagram.png`: Visual representation of the ventilation system

## Standards and Guidelines

The analysis follows these general guidelines:
- Minimum grille velocity: 2.5 m/s
- Maximum grille velocity: 4.0 m/s
- Minimum duct velocity: 3.0 m/s
- Maximum duct velocity: 10.0 m/s
- Minimum air changes per hour for bathrooms: 10 ACH

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

Startup information:

This is the ventilation of a facility bathrooms:

Planta Baja Bano C; 25m3, 2 rejillas 4"x11" conectadas a ducto 6"x11" de 4m que termina en expansión brusca a ducto11"x11" de 2m con codo 90 grados  11"x11" con cambio de sección de 11"x11" a sección circular de 12.5" y ducto circular vertical de 6m de 12.5" con codo 90 grados para conectar con ventilador en techo.
Planta Baja Baño D; 18m3, 3 rejillas 6x11 conectadas por codos 90 grados a ducto 6"x11" de 4m conectado a ducto Planta Baja Bano C.

Cafetín Bano C: 6.5m3, 1 rejilla 7.5"x10", ducto 10" 1,44m con codo 90 grados 10" y conexión con ducto 12.5" ascendente de PB.
Cafetín Bano D: 9m3, 2 rejillas 7.5"x11", ducto 10" 2m conecta con ducto 10" Cafetín Baño C.

Piso 1 Bano C: 3 10m3, rejillas 6x14, ducto 3.20m 11x6, Codo 90 11x6 con tubería ascendente 11x6 de 3.5m y conexión en T con tubería principal piso 2.
Piso 1 Bano D: 3 10m3, rejillas 6x14, ducto 3.20m 11x6, Codo 90 11x6 con tubería ascendente 11x6 de 3.5m  y conexión en T con tubería principal piso 2.

Piso 2 Bano C; 10m3, 3 rejillas 4x11, ducto 3.20m 11x6, conexión en T con Ducto Principal piso 2.
Piso 2 Bano D; 10m3, 3 rejillas 4x11, ducto 3.20m 11x6, conexión en T con Ducto Principal piso 2
Ducto Principal piso 2: 11x11 12m largo con terminal cuadrado-circular 12” diámetro con conexión directa al ventilador en techo.

Planta baja y cafetín comparten el mismo ventilador en el techo.
Piso1 y Piso 2 comparte el mismo ventilador en el techo.

las condiciones climáticas promedio son 25C y 60% HR, las extremas son 28C, 93% HR, 16C 30%. La altitud es de 850msnm.

evaluate the worst climatic case and adjust changes per hour to allow that smalest velocity be in range of standards. suggest grid changes if necessary.

Calcula ventiladores centrífugos comerciales adecuados para cada instalación sin dejar que las velocidades mínimas estén fuera de norma.  Crea un código para guardar en git y utilizar en mybinder.org

