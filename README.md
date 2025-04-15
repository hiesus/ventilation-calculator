# Ventilation System Analysis Tool

This tool analyzes ventilation systems for bathroom facilities, calculates required airflow, pressure losses, and recommends appropriate centrifugal fans.

## Features

- Calculates required airflow based on room volumes and recommended air changes
- Analyzes duct systems and calculates pressure losses
- Ensures minimum velocities are within standards
- Recommends appropriate centrifugal fans
- Generates comprehensive reports and visualizations

## Usage

### Online with MyBinder

1. Visit [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/YOUR_USERNAME/ventilation-analysis/main)
2. Open and run the `ventilation_analysis


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

