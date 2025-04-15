import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Constants and conversion factors
INCH_TO_M = 0.0254
FT_TO_M = 0.3048
CFM_TO_M3S = 0.000471947
PA_TO_INWG = 0.00402
M3S_TO_CFM = 1/CFM_TO_M3S

# Climate conditions
climate_conditions = {
    'average': {'temperature': 25, 'humidity': 60},
    'extreme_hot_humid': {'temperature': 28, 'humidity': 93},
    'extreme_cold_dry': {'temperature': 16, 'humidity': 30}
}

# Altitude
altitude = 850  # meters above sea level

# Air density calculation based on temperature, humidity and altitude
def air_density(temp_c, humidity, altitude_m):
    # Barometric pressure at altitude (Pa)
    p0 = 101325  # sea level standard pressure, Pa
    g = 9.80665  # gravitational acceleration, m/s²
    M = 0.0289644  # molar mass of dry air, kg/mol
    R = 8.31447  # universal gas constant, J/(mol·K)
    T0 = 288.15  # sea level standard temperature, K

    # Temperature lapse rate, K/m
    L = 0.0065

    # Barometric pressure at altitude
    p = p0 * (1 - L * altitude_m / T0) ** (g * M / (R * L))

    # Saturation vapor pressure (Pa)
    temp_k = temp_c + 273.15
    es = 611.2 * np.exp(17.67 * temp_c / (temp_c + 243.5))

    # Actual vapor pressure (Pa)
    e = humidity / 100 * es

    # Density of moist air (kg/m³)
    Rd = 287.058  # specific gas constant for dry air, J/(kg·K)
    Rv = 461.495  # specific gas constant for water vapor, J/(kg·K)

    density = (p - e) / (Rd * temp_k) + e / (Rv * temp_k)

    return density

# Room data
rooms = {
    'PB_Bano_C': {'volume': 25, 'grids': [{'size': (4, 11), 'count': 2}]},
    'PB_Bano_D': {'volume': 18, 'grids': [{'size': (6, 11), 'count': 3}]},
    'Cafetin_Bano_C': {'volume': 6.5, 'grids': [{'size': (7.5, 10), 'count': 1}]},
    'Cafetin_Bano_D': {'volume': 9, 'grids': [{'size': (7.5, 11), 'count': 2}]},
    'Piso1_Bano_C': {'volume': 10, 'grids': [{'size': (6, 14), 'count': 3}]},
    'Piso1_Bano_D': {'volume': 10, 'grids': [{'size': (6, 14), 'count': 3}]},
    'Piso2_Bano_C': {'volume': 10, 'grids': [{'size': (4, 11), 'count': 3}]},
    'Piso2_Bano_D': {'volume': 10, 'grids': [{'size': (4, 11), 'count': 3}]}
}

# Duct system data
duct_systems = {
    'System1': {  # PB and Cafetin
        'rooms': ['PB_Bano_C', 'PB_Bano_D', 'Cafetin_Bano_C', 'Cafetin_Bano_D'],
        'sections': [
            {'name': 'PB_C_main', 'shape': 'rectangular', 'size': (6, 11), 'length': 4, 'connected_to': 'PB_C_expansion'},
            {'name': 'PB_C_expansion', 'shape': 'rectangular', 'size': (11, 11), 'length': 2, 'connected_to': 'PB_C_elbow'},
            {'name': 'PB_C_elbow', 'shape': 'elbow', 'size': (11, 11), 'angle': 90, 'connected_to': 'PB_C_transition'},
            {'name': 'PB_C_transition', 'shape': 'transition', 'size_in': (11, 11), 'size_out': (12.5, 12.5), 'connected_to': 'PB_C_vertical'},
            {'name': 'PB_C_vertical', 'shape': 'circular', 'diameter': 12.5, 'length': 6, 'connected_to': 'PB_C_final_elbow'},
            {'name': 'PB_C_final_elbow', 'shape': 'elbow', 'diameter': 12.5, 'angle': 90, 'connected_to': 'fan1'},

            {'name': 'PB_D_main', 'shape': 'rectangular', 'size': (6, 11), 'length': 4, 'connected_to': 'PB_C_main'},

            {'name': 'Cafetin_C_main', 'shape': 'circular', 'diameter': 10, 'length': 1.44, 'connected_to': 'Cafetin_C_elbow'},
            {'name': 'Cafetin_C_elbow', 'shape': 'elbow', 'diameter': 10, 'angle': 90, 'connected_to': 'PB_C_vertical'},

            {'name': 'Cafetin_D_main', 'shape': 'circular', 'diameter': 10, 'length': 2, 'connected_to': 'Cafetin_C_main'}
        ]
    },
    'System2': {  # Piso 1 and Piso 2
        'rooms': ['Piso1_Bano_C', 'Piso1_Bano_D', 'Piso2_Bano_C', 'Piso2_Bano_D'],
        'sections': [
            {'name': 'Piso1_C_main', 'shape': 'rectangular', 'size': (11, 6), 'length': 3.2, 'connected_to': 'Piso1_C_elbow'},
            {'name': 'Piso1_C_elbow', 'shape': 'elbow', 'size': (11, 6), 'angle': 90, 'connected_to': 'Piso1_C_vertical'},
            {'name': 'Piso1_C_vertical', 'shape': 'rectangular', 'size': (11, 6), 'length': 3.5, 'connected_to': 'Piso2_main'},

            {'name': 'Piso1_D_main', 'shape': 'rectangular', 'size': (11, 6), 'length': 3.2, 'connected_to': 'Piso1_D_elbow'},
            {'name': 'Piso1_D_elbow', 'shape': 'elbow', 'size': (11, 6), 'angle': 90, 'connected_to': 'Piso1_D_vertical'},
            {'name': 'Piso1_D_vertical', 'shape': 'rectangular', 'size': (11, 6), 'length': 3.5, 'connected_to': 'Piso2_main'},

            {'name': 'Piso2_C_main', 'shape': 'rectangular', 'size': (11, 6), 'length': 3.2, 'connected_to': 'Piso2_main'},
            {'name': 'Piso2_D_main', 'shape': 'rectangular', 'size': (11, 6), 'length': 3.2, 'connected_to': 'Piso2_main'},

            {'name': 'Piso2_main', 'shape': 'rectangular', 'size': (11, 11), 'length': 12, 'connected_to': 'Piso2_transition'},
            {'name': 'Piso2_transition', 'shape': 'transition', 'size_in': (11, 11), 'size_out': (12, 12), 'connected_to': 'fan2'}
        ]
    }
}

# Ventilation standards
min_air_changes = {
    'bathroom': 8,  # Air changes per hour for bathrooms
    'recommended': 12  # Recommended air changes per hour
}

# Velocity standards (m/s)
velocity_standards = {
    'min_duct': 2.5,  # Minimum velocity in ducts
    'max_duct': 10.0,  # Maximum velocity in ducts
    'min_grid': 1.0,   # Minimum velocity at grids
    'max_grid': 2.5    # Maximum velocity at grids
}

def calculate_required_airflow(room_volume, air_changes):
    """Calculate required airflow in m³/s based on room volume and air changes per hour"""
    return room_volume * air_changes / 3600  # Convert from m³/h to m³/s

def calculate_grid_area(grid_size):
    """Calculate grid area in m² from dimensions in inches"""
    width, height = grid_size
    return (width * INCH_TO_M) * (height * INCH_TO_M)

def calculate_duct_area(duct_size_or_diameter, shape='rectangular'):
    """Calculate duct area in m² from dimensions in inches"""
    if shape == 'rectangular':
        width, height = duct_size_or_diameter
        return (width * INCH_TO_M) * (height * INCH_TO_M)
    elif shape == 'circular':
        diameter = duct_size_or_diameter
        return np.pi * ((diameter * INCH_TO_M) / 2) ** 2
    else:
        raise ValueError(f"Unknown shape: {shape}")

def calculate_hydraulic_diameter(duct_size, shape='rectangular'):
    """Calculate hydraulic diameter in m from dimensions in inches"""
    if shape == 'rectangular':
        width, height = duct_size
        width_m, height_m = width * INCH_TO_M, height * INCH_TO_M
        return 4 * (width_m * height_m) / (2 * (width_m + height_m))
    elif shape == 'circular':
        return duct_size * INCH_TO_M
    else:
        raise ValueError(f"Unknown shape: {shape}")

def calculate_velocity(airflow, area):
    """Calculate velocity in m/s from airflow (m³/s) and area (m²)"""
    return airflow / area

def calculate_pressure_loss(airflow, section, density):
    """Calculate pressure loss in Pa for a duct section"""
    if section['shape'] == 'rectangular':
        area = calculate_duct_area(section['size'], 'rectangular')
        dh = calculate_hydraulic_diameter(section['size'], 'rectangular')
        velocity = calculate_velocity(airflow, area)

        # Friction factor (Colebrook equation approximation)
        roughness = 0.00015  # m, typical for galvanized steel
        reynolds = velocity * dh * density / 1.81e-5  # Reynolds number

        # Use Swamee-Jain equation for friction factor
        if reynolds > 2300:
            f = 0.25 / (np.log10(roughness/(3.7*dh) + 5.74/reynolds**0.9))**2
        else:
            f = 64 / reynolds

        # Pressure loss due to friction
        pressure_loss = f * (section['length'] / dh) * (density * velocity**2 / 2)

    elif section['shape'] == 'circular':
        area = calculate_duct_area(section['diameter'], 'circular')
        dh = calculate_hydraulic_diameter(section['diameter'], 'circular')
        velocity = calculate_velocity(airflow, area)

        # Friction factor
        roughness = 0.00015  # m
        reynolds = velocity * dh * density / 1.81e-5

        if reynolds > 2300:
            f = 0.25 / (np.log10(roughness/(3.7*dh) + 5.74/reynolds**0.9))**2
        else:
            f = 64 / reynolds

        pressure_loss = f * (section['length'] / dh) * (density * velocity**2 / 2)

    elif section['shape'] == 'elbow':
        if 'diameter' in section:
            area = calculate_duct_area(section['diameter'], 'circular')
            velocity = calculate_velocity(airflow, area)
        else:
            area = calculate_duct_area(section['size'], 'rectangular')
            velocity = calculate_velocity(airflow, area)

        # Loss coefficient for 90° elbow
        k = 0.3 if section['angle'] == 90 else 0.2

        pressure_loss = k * (density * velocity**2 / 2)

    elif section['shape'] == 'transition':
        area_in = calculate_duct_area(section['size_in'], 'rectangular')
        area_out = calculate_duct_area(section['size_out'], 'rectangular' if len(section['size_out']) == 2 else 'circular')

        velocity_in = calculate_velocity(airflow, area_in)
        velocity_out = calculate_velocity(airflow, area_out)

        # Loss coefficient for expansion or contraction
        if area_out > area_in:  # Expansion
            k = (1 - area_in/area_out)**2
        else:  # Contraction
            k = 0.5 * (1 - area_out/area_in)

        pressure_loss = k * (density * velocity_out**2 / 2)

    else:
        raise ValueError(f"Unknown section shape: {section['shape']}")

    return pressure_loss

def analyze_ventilation_system():
    """Analyze the ventilation system and recommend changes"""
    results = {}

    # Use the worst climatic case for calculations
    worst_case = 'extreme_hot_humid'  # Typically the hot and humid condition is worst for ventilation
    temp = climate_conditions[worst_case]['temperature']
    humidity = climate_conditions[worst_case]['humidity']
    # Calculate air density for the worst case
    density = air_density(temp, humidity, altitude)

    # Analyze each system
    for system_name, system in duct_systems.items():
        system_results = {
            'rooms': {},
            'total_airflow': 0,
            'pressure_loss': 0,
            'sections': {}
        }

        # Calculate required airflow for each room
        for room_name in system['rooms']:
            room = rooms[room_name]

            # Start with minimum air changes and adjust if needed
            air_changes = min_air_changes['bathroom']

            # Calculate required airflow
            required_airflow = calculate_required_airflow(room['volume'], air_changes)

            # Calculate total grid area
            total_grid_area = 0
            for grid in room['grids']:
                grid_area = calculate_grid_area(grid['size']) * grid['count']
                total_grid_area += grid_area

            # Calculate velocity at grids
            grid_velocity = calculate_velocity(required_airflow, total_grid_area)

            # If grid velocity is below minimum, increase air changes
            if grid_velocity < velocity_standards['min_grid']:
                # Calculate minimum air changes needed to meet velocity standard
                min_required_airflow = velocity_standards['min_grid'] * total_grid_area
                min_air_changes = min_required_airflow * 3600 / room['volume']

                # Update air changes and recalculate
                air_changes = max(air_changes, min_air_changes)
                required_airflow = calculate_required_airflow(room['volume'], air_changes)
                grid_velocity = calculate_velocity(required_airflow, total_grid_area)

            # If grid velocity is above maximum, suggest increasing grid size
            grid_size_increase = None
            if grid_velocity > velocity_standards['max_grid']:
                # Calculate required grid area
                required_grid_area = required_airflow / velocity_standards['max_grid']

                # Calculate increase factor
                increase_factor = required_grid_area / total_grid_area

                # Suggest increasing grid size
                grid_size_increase = {
                    'current_area': total_grid_area,
                    'required_area': required_grid_area,
                    'increase_factor': increase_factor
                }

            # Store room results
            system_results['rooms'][room_name] = {
                'volume': room['volume'],
                'air_changes': air_changes,
                'required_airflow': required_airflow,
                'grid_area': total_grid_area,
                'grid_velocity': grid_velocity,
                'grid_size_increase': grid_size_increase
            }

            # Add to total system airflow
            system_results['total_airflow'] += required_airflow

        # Analyze duct sections
        for section in system['sections']:
            # Calculate area
            if section['shape'] == 'rectangular':
                area = calculate_duct_area(section['size'], 'rectangular')
            elif section['shape'] == 'circular':
                area = calculate_duct_area(section['diameter'], 'circular')
            elif section['shape'] == 'elbow':
                if 'diameter' in section:
                    area = calculate_duct_area(section['diameter'], 'circular')
                else:
                    area = calculate_duct_area(section['size'], 'rectangular')
            elif section['shape'] == 'transition':
                area = calculate_duct_area(section['size_in'], 'rectangular')

            # Calculate velocity
            velocity = calculate_velocity(system_results['total_airflow'], area)

            # Calculate pressure loss
            pressure_loss = calculate_pressure_loss(system_results['total_airflow'], section, density)

            # Store
