# Bathroom Ventilation System Calculator

This repository contains a Python-based tool for calculating ventilation requirements for bathroom facilities. It analyzes duct systems, calculates pressure drops, and recommends appropriate centrifugal fans.

## Features

- Calculates required airflow based on room volume and air changes per hour
- Analyzes pressure drops in complex duct systems
- Checks duct velocities against recommended standards
- Recommends suitable centrifugal fans
- Visualizes results with charts and graphs
- Exports data to CSV files for further analysis

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/ventilation-calculator.git
cd ventilation-calculator
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Using Jupyter Notebook

The easiest way to use this tool is through the provided Jupyter notebook:

```bash
jupyter notebook ventilation_analysis.ipynb
```

This will open the notebook in your browser where you can run the analysis step by step.

### Using Python Script

Alternatively, you can run the main script directly:

```bash
python -c "import ventilation_calculator as vc; report = vc.calculate_ventilation_system(); vc.print_report(report); vc.visualize_system(report); vc.export_to_csv(report)"
```

## Online Usage with MyBinder

You can run this tool online without installing anything by using MyBinder:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/yourusername/ventilation-calculator/main?filepath=ventilation_analysis.ipynb)

Click the Binder badge above to launch the notebook in your browser.

## Project Structure

- `ventilation_calculator.py`: Main module containing all calculation functions
- `ventilation_analysis.ipynb`: Jupyter notebook for interactive analysis
- `requirements.txt`: List of Python dependencies
- `README.md`: This file

## Customization

You can modify the following parameters in the code to match your specific requirements:

- `CHANGES_PER_HOUR`: Number of air changes per hour (default: 12)
- `ALTITUDE`: Altitude in meters above sea level (default: 850)
- `TEMP_AVG`: Average temperature in °C (default: 25)
- `RH_AVG`: Average relative humidity in % (default: 60)

## License

This project is licensed under the MIT License - see the LICENSE file for details.


This is the ventilation startup specifications:

Planta Baja Bano C; 25m3, 2 rejillas 4"x11" conectadas a ducto 6"x11" de 4m que termina en expansión brusca a ducto11"x11" de 2m con codo 90 grados  11"x11" con cambio de sección de 11"x11" a sección circular de 12.5" y ducto circular vertical de 6m de 12.5" con codo 90 grados para conectar con ventilador en techo.
Planta Baja Baño D; 18m3, 3 rejillas 6x11 conectadas por codos 90 grados a ducto 6"x11" de 4m conectado a ducto Planta Baja Bano C.

Cafetín Bano C: 6.5m3 1 rejilla 7.5"x10", ducto 10" 1,44m con codo 90 grados 10" y conexión con ducto 12.5" ascendente de PB.
Cafetín Bano D: 9m3 2 rejillas 7.5"x11", ducto 10" 2m conecta con ducto 10" Cafetín Baño C.

Piso 1 Bano C: 3 rejillas 6x14, ducto 3.20m 11x6, Codo 90 11x6 con tubería ascendente 11x6 de 3.5m y conexión en T con tubería principal piso 2.
Piso 1 Bano D: 3 rejillas 6x14, ducto 3.20m 11x6, Codo 90 11x6 con tubería ascendente 11x6 de 3.5m  y conexión en T con tubería principal piso 2.

Piso 2 Bano C; 3 rejillas 4x11, ducto 3.20m 11x6, conexión en T con Ducto Principal piso 2.
Piso 2 Bano D; 3 rejillas 4x11, ducto 3.20m 11x6, conexión en T con Ducto Principal piso 2
Ducto Principal piso 2: 11x11 12m largo con terminal cuadrado-circular 12” diámetro con conexión directa al ventilador en techo.

Planta baja y cafetín comparten el mismo ventilador en el techo.
Piso1 y Piso 2 comparte el mismo ventilador en el techo.

las condiciones climáticas promedio son 25C y 60% HR, las extremas son 28C, 93% HR, 16C 30%. La altitud es de 850msnm.

Use 12 changes per hour. Do only 1 phyton file

Calcula ventiladores centrífugos comerciales adecuados para cada instalación sin dejar que las velocidades mínimas estén fuera de norma


