import math
import numpy as np
import matplotlib.pyplot as plt

# Constants
AIR_CHANGES_PER_HOUR = 12
ALTITUDE = 850  # meters above sea level
STANDARD_TEMP = 25  # °C
STANDARD_RH = 0.60  # 60% relative humidity

# Air density correction for altitude
def calculate_air_density(altitude, temp, rh):
    # Standard air density at sea level (kg/m³)
    rho_0 = 1.225
    
    # Altitude correction (simplified)
    pressure_ratio = (1 - 2.25577e-5 * altitude) ** 5.25588
    temp_ratio = (273.15 + temp) / 288.15
    
    # Basic density calculation
    rho = rho_0 * pressure_ratio / temp_ratio
    
    # Humidity effect (simplified)
    humidity_factor = 1 - 0.378 * rh * 0.0315  # approximation of humidity effect
    
    return rho * humidity_factor

# Function to convert dimensions from inches to meters
def inch_to_meter(inches):
    return inches * 0.0254

# Function to calculate area in square meters
def calculate_area(width_inch, height_inch):
    width_m = inch_to_meter(width_inch)
    height_m = inch_to_meter(height_inch)
    return width_m * height_m

# Function to calculate equivalent diameter for rectangular ducts
def equivalent_diameter(width_inch, height_inch):
    width_m = inch_to_meter(width_inch)
    height_m = inch_to_meter(height_inch)
    return 1.3 * ((width_m * height_m) ** 0.625) / ((width_m + height_m) ** 0.25)

# Function to calculate required airflow based on room volume and air changes
def calculate_airflow(volume_m3, air_changes_per_hour):
    # Convert to m³/h
    airflow_m3h = volume_m3 * air_changes_per_hour
    # Convert to CFM (Cubic Feet per Minute)
    airflow_cfm = airflow_m3h * 0.589
    return airflow_m3h, airflow_cfm

# Function to calculate velocity in a duct
def calculate_velocity(airflow_m3h, area_m2):
    # Convert m³/h to m³/s
    airflow_m3s = airflow_m3h / 3600
    return airflow_m3s / area_m2

# Function to calculate friction loss in straight ducts
def calculate_friction_loss(velocity, equiv_diameter, length, friction_factor=0.02):
    # Darcy-Weisbach equation for pressure loss
    return friction_factor * (length / equiv_diameter) * (velocity ** 2) / 2

# Function to calculate pressure loss due to fittings
def calculate_fitting_loss(velocity, k_factor):
    # Dynamic pressure loss equation
    return k_factor * (velocity ** 2) / 2

# Function to calculate total pressure loss
def calculate_total_pressure_loss(straight_losses, fitting_losses):
    return sum(straight_losses) + sum(fitting_losses)

# K-factors for common fittings
def get_k_factor(fitting_type):
    k_factors = {
        'elbow_90': 0.3,       # 90° elbow
        'tee': 1.0,            # T-connection
        'entry': 0.5,          # Entry from room to duct
        'exit': 1.0,           # Exit from duct
        'sudden_expansion': 0.8,  # Sudden expansion
        'transition': 0.2      # Section change
    }
    return k_factors.get(fitting_type, 0.5)  # Default if not found

# Define bathroom spaces
bathrooms = {
    'PB_Bano_C': {'volume': 25, 'grilles': 2, 'grille_size': (4, 11)},
    'PB_Bano_D': {'volume': 18, 'grilles': 3, 'grille_size': (6, 11)},
    'Cafeteria_Bano_C': {'volume': 6.5, 'grilles': 1, 'grille_size': (7.5, 10)},
    'Cafeteria_Bano_D': {'volume': 9, 'grilles': 2, 'grille_size': (7.5, 11)},
    'Piso1_Bano_C': {'volume': 25, 'grilles': 3, 'grille_size': (6, 14)},
    'Piso1_Bano_D': {'volume': 25, 'grilles': 3, 'grille_size': (6, 14)},
    'Piso2_Bano_C': {'volume': 20, 'grilles': 3, 'grille_size': (4, 11)},
    'Piso2_Bano_D': {'volume': 20, 'grilles': 3, 'grille_size': (4, 11)}
}

# Calculate airflow requirements for each bathroom
for name, bathroom in bathrooms.items():
    airflow_m3h, airflow_cfm = calculate_airflow(bathroom['volume'], AIR_CHANGES_PER_HOUR)
    bathroom['airflow_m3h'] = airflow_m3h
    bathroom['airflow_cfm'] = airflow_cfm
    print(f"{name}: Volume {bathroom['volume']} m³, Required airflow: {airflow_m3h:.2f} m³/h ({airflow_cfm:.2f} CFM)")

# Calculate total airflow requirements for each fan
fan1_bathrooms = ['PB_Bano_C', 'PB_Bano_D', 'Cafeteria_Bano_C', 'Cafeteria_Bano_D']
fan2_bathrooms = ['Piso1_Bano_C', 'Piso1_Bano_D', 'Piso2_Bano_C', 'Piso2_Bano_D']

fan1_total_m3h = sum(bathrooms[name]['airflow_m3h'] for name in fan1_bathrooms)
fan1_total_cfm = sum(bathrooms[name]['airflow_cfm'] for name in fan1_bathrooms)

fan2_total_m3h = sum(bathrooms[name]['airflow_m3h'] for name in fan2_bathrooms)
fan2_total_cfm = sum(bathrooms[name]['airflow_cfm'] for name in fan2_bathrooms)

print(f"\nFan 1 (Ground Floor & Cafeteria): Total airflow: {fan1_total_m3h:.2f} m³/h ({fan1_total_cfm:.2f} CFM)")
print(f"Fan 2 (Floors 1 & 2): Total airflow: {fan2_total_m3h:.2f} m³/h ({fan2_total_cfm:.2f} CFM)")

# Define duct systems for pressure loss calculation
# This is a simplified calculation for the critical path of each fan system

# Fan 1 System - Ground Floor and Cafeteria
# Critical path through PB_Bano_C to fan
fan1_critical_path = [
    # Format: [duct_width_inch, duct_height_inch, length_m, duct_type]
    {'section': 'PB_Bano_C Grille', 'dim': (4, 11), 'length': 0, 'type': 'entry'},
    {'section': 'Duct 6x11', 'dim': (6, 11), 'length': 4, 'type': 'straight'},
    {'section': 'Expansion to 11x11', 'dim': (11, 11), 'length': 0, 'type': 'sudden_expansion'},
    {'section': 'Duct 11x11', 'dim': (11, 11), 'length': 2, 'type': 'straight'},
    {'section': '90° Elbow 11x11', 'dim': (11, 11), 'length': 0, 'type': 'elbow_90'},
    {'section': 'Transition to 12.5" circular', 'dim': (12.5, 12.5), 'length': 0, 'type': 'transition'},
    {'section': '12.5" Vertical Duct', 'dim': (12.5, 12.5), 'length': 6, 'type': 'straight'},
    {'section': '90° Elbow to Fan', 'dim': (12.5, 12.5), 'length': 0, 'type': 'elbow_90'},
    {'section': 'Fan Connection', 'dim': (12.5, 12.5), 'length': 0, 'type': 'exit'}
]

# Fan 2 System - Floors 1 & 2
# Critical path through Piso2 to fan
fan2_critical_path = [
    # Format: [duct_width_inch, duct_height_inch, length_m, duct_type]
    {'section': 'Piso2_Bano_C Grille', 'dim': (4, 11), 'length': 0, 'type': 'entry'},
    {'section': 'Duct 11x6', 'dim': (11, 6), 'length': 3.2, 'type': 'straight'},
    {'section': 'T-Connection', 'dim': (11, 6), 'length': 0, 'type': 'tee'},
    {'section': 'Main Duct 11x11', 'dim': (11, 11), 'length': 12, 'type': 'straight'},
    {'section': 'Transition to 12" circular', 'dim': (12, 12), 'length': 0, 'type': 'transition'},
    {'section': 'Fan Connection', 'dim': (12, 12), 'length': 0, 'type': 'exit'}
]

# Calculate pressure losses for Fan 1
def calculate_path_pressure_loss(path, total_airflow_m3h):
    straight_losses = []
    fitting_losses = []
    
    air_density = calculate_air_density(ALTITUDE, STANDARD_TEMP, STANDARD_RH)
    
    print("\nSection-by-section pressure loss calculation:")
    print("-" * 70)
    print(f"{'Section':<25} {'Velocity (m/s)':<15} {'Pressure Loss (Pa)':<20}")
    print("-" * 70)
    
    for section in path:
        # For circular ducts, both dimensions are the same (diameter)
        if section['dim'][0] == section['dim'][1]:  # Circular duct
            area = math.pi * (inch_to_meter(section['dim'][0])/2)**2
            equiv_diam = inch_to_meter(section['dim'][0])
        else:  # Rectangular duct
            area = calculate_area(section['dim'][0], section['dim'][1])
            equiv_diam = equivalent_diameter(section['dim'][0], section['dim'][1])
        
        velocity = calculate_velocity(total_airflow_m3h, area)
        
        if section['type'] == 'straight':
            loss = calculate_friction_loss(velocity, equiv_diam, section['length'])
            straight_losses.append(loss)
            loss_type = "Friction"
        else:
            k = get_k_factor(section['type'])
            loss = calculate_fitting_loss(velocity, k)
            fitting_losses.append(loss)
            loss_type = f"Fitting (k={k})"
        
        print(f"{section['section']:<25} {velocity:.2f} m/s      {loss:.2f} Pa ({loss_type})")
    
    # Convert losses to Pascals and apply air density correction
    total_loss_pa = (sum(straight_losses) + sum(fitting_losses)) * air_density
    total_loss_inwg = total_loss_pa * 0.00402  # Convert Pa to inches of water gauge
    
    print("-" * 70)
    print(f"Total pressure loss: {total_loss_pa:.2f} Pa ({total_loss_inwg:.4f} inWG)")
    
    return total_loss_pa, total_loss_inwg, straight_losses, fitting_losses

# Calculate pressure losses for both fans
print("\n=== Fan 1 (Ground Floor & Cafeteria) Pressure Loss Calculation ===")
fan1_pressure_pa, fan1_pressure_inwg, fan1_straight, fan1_fittings = calculate_path_pressure_loss(fan1_critical_path, fan1_total_m3h)

print("\n=== Fan 2 (Floors 1 & 2) Pressure Loss Calculation ===")
fan2_pressure_pa, fan2_pressure_inwg, fan2_straight, fan2_fittings = calculate_path_pressure_loss(fan2_critical_path, fan2_total_m3h)

# Fan selection recommendation
def recommend_fan(airflow_cfm, pressure_inwg, system_name):
    safety_factor = 1.3  # Add 30% safety factor
    design_airflow = airflow_cfm * safety_factor
    design_pressure = pressure_inwg * safety_factor
    
    print(f"\n=== Fan Recommendation for {system_name} ===")
    print(f"Required airflow: {airflow_cfm:.2f} CFM")
    print(f"Required static pressure: {pressure_inwg:.4f} inWG")
    print(f"Design airflow (with safety factor): {design_airflow:.2f} CFM")
    print(f"Design static pressure (with safety factor): {design_pressure:.4f} inWG")
    
    # Fan model recommendations based on performance requirements
    # This is a simplified selection - in reality, you'd look at fan curves
    if design_airflow < 500:
        fan_size = "Small"
    elif design_airflow < 1000:
        fan_size = "Medium"
    else:
        fan_size = "Large"
        
    if design_pressure < 0.5:
        pressure_class = "Low Pressure"
    elif design_pressure < 1.0:
        pressure_class = "Medium Pressure"
    else:
        pressure_class = "High Pressure"
    
    print(f"Recommended fan type: Centrifugal {fan_size}, {pressure_class} class")
    
    # Example commercial models (these would need to be replaced with actual models)
    if system_name == "Fan 1 (Ground Floor & Cafeteria)":
        if fan_size == "Small" and pressure_class == "Low Pressure":
            return "Greenheck SQ-85-VG, FANTECH FKD 10XL, or equivalent"
        elif fan_size == "Small" and pressure_class == "Medium Pressure":
            return "Greenheck SQ-100-VG, FANTECH FKD 12XL, or equivalent"
        elif fan_size == "Medium" and pressure_class == "Low Pressure":
            return "Greenheck SQ-120-VG, FANTECH FKD 14, or equivalent"
        elif fan_size == "Medium" and pressure_class == "Medium Pressure":
            return "Greenheck SQ-140-VG, FANTECH FKD 16XL, or equivalent"
        else:
            return "Greenheck SQ-160-VG, FANTECH FKD 18XL, or equivalent"
    else:  # Fan 2
        if fan_size == "Small" and pressure_class == "Low Pressure":
            return "Greenheck SQ-85-VG, FANTECH FKD 10XL, or equivalent"
        elif fan_size == "Small" and pressure_class == "Medium Pressure":
            return "Greenheck SQ-100-VG, FANTECH FKD 12XL, or equivalent"
        elif fan_size == "Medium" and pressure_class == "Low Pressure":
            return "Greenheck SQ-120-VG, FANTECH FKD 14, or equivalent"
        elif fan_size == "Medium" and pressure_class == "Medium Pressure":
            return "Greenheck SQ-140-VG, FANTECH FKD 16XL, or equivalent"
        else:
            return "Greenheck SQ-160-VG, FANTECH FKD 18XL, or equivalent"

# Check if airflow velocities are within recommended ranges
def check_airflow_velocities(path, total_airflow_m3h):
    print("\n=== Airflow Velocity Check ===")
    print("-" * 70)
    print(f"{'Section':<25} {'Size':<15} {'Velocity (m/s)':<15} {'Status':<15}")
    print("-" * 70)
    
    all_ok = True
    
    # Recommended velocity ranges based on ASHRAE standards
    # For bathroom ventilation ducts - these are simplified ranges
    min_velocity = 2.5  # m/s (to prevent settling)
    max_velocity = 10.0  # m/s (to prevent noise)
    
    # For grilles, the velocity should be lower
    min_grille_velocity = 1.0  # m/s
    max_grille_velocity = 3.0  # m/s
    
    for section in path:
        # For circular ducts, both dimensions are the same (diameter)
        if section['dim'][0] == section['dim'][1]:  # Circular duct
            area = math.pi * (inch_to_meter(section['dim'][0])/2)**2
            size_str = f"{section['dim'][0]}\" diameter"
        else:  # Rectangular duct
            area = calculate_area(section['dim'][0], section['dim'][1])
            size_str = f"{section['dim'][0]}\"x{section['dim'][1]}\""
        
        velocity = calculate_velocity(total_airflow_m3h, area)
        
        # Check if it's a grille or entry point
        if section['type'] == 'entry':
            status = "OK" if min_grille_velocity <= velocity <= max_grille_velocity else "REVIEW"
            if status == "REVIEW":
                all_ok = False
        else:
            status = "OK" if min_velocity <= velocity <= max_velocity else "REVIEW"
            if status == "REVIEW":
                all_ok = False
        
        print(f"{section['section']:<25} {size_str:<15} {velocity:.2f} m/s      {status:<15}")
    
    print("-" * 70)
    if all_ok:
        print("All velocities are within recommended ranges.")
    else:
        print("Some sections need review. Velocities should be:")
        print(f"- Ducts: {min_velocity} - {max_velocity} m/s")
        print(f"- Grilles: {min_grille_velocity} - {max_grille_velocity} m/s")
    
    return all_ok

# Get fan recommendations
fan1_recommendation = recommend_fan(fan1_total_cfm, fan1_pressure_inwg, "Fan 1 (Ground Floor & Cafeteria)")
fan2_recommendation = recommend_fan(fan2_total_cfm, fan2_pressure_inwg, "Fan 2 (Floors 1 & 2)")

print(f"\nRecommended model for Fan 1: {fan1_recommendation}")
print(f"Recommended model for Fan 2: {fan2_recommendation}")

# Check velocities
print("\nChecking velocities for Fan 1 system:")
fan1_velocities_ok = check_airflow_velocities(fan1_critical_path, fan1_total_m3h)

print("\nChecking velocities for Fan 2 system:")
fan2_velocities_ok = check_airflow_velocities(fan2_critical_path, fan2_total_m3h)

# Calculate proper load correction for high altitude and temperature/humidity
def altitude_correction_factor(altitude, temp, rh):
    # Standard air density at sea level (kg/m³)
    rho_0 = 1.225
    
    # Calculate actual air density
    rho_actual = calculate_air_density(altitude, temp, rh)
    
    # Correction factor is the ratio of standard to actual density
    return rho_0 / rho_actual

# Calculate correction factors for different conditions
std_correction = altitude_correction_factor(ALTITUDE, STANDARD_TEMP, STANDARD_RH)
hot_humid_correction = altitude_correction_factor(ALTITUDE, 28, 0.93)
cold_dry_correction = altitude_correction_factor(ALTITUDE, 16, 0.30)

print("\n=== Altitude and Climate Correction Factors ===")
print(f"Standard conditions (25°C, 60% RH): {std_correction:.3f}")
print(f"Hot humid conditions (28°C, 93% RH): {hot_humid_correction:.3f}")
print(f"Cold dry conditions (16°C, 30% RH): {cold_dry_correction:.3f}")

print("\nFinal Fan Requirements (with standard condition correction):")
print(f"Fan 1: {fan1_total_cfm * std_correction:.2f} CFM at {fan1_pressure_inwg:.4f} inWG")
print(f"Fan 2: {fan2_total_cfm * std_correction:.2f} CFM at {fan2_pressure_inwg:.4f} inWG")

# Visualize the pressure loss breakdown for both fans
def plot_pressure_losses(straight_losses, fitting_losses, fan_name):
    total_loss = sum(straight_losses) + sum(fitting_losses)
    
    # Calculate percentages
    straight_pct = sum(straight_losses) / total_loss * 100
    fitting_pct = sum(fitting_losses) / total_loss * 100
    
    # Create pie chart
    labels = ['Straight Duct Losses', 'Fitting Losses']
    sizes = [straight_pct, fitting_pct]
    colors = ['lightblue', 'lightgreen']
    
    plt.figure(figsize=(8, 6))
    plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    plt.axis('equal')
    plt.title(f'Pressure Loss Breakdown for {fan_name}')
    plt.tight_layout()
    
    # Save the chart
    filename = f"{fan_name.replace(' ', '_').replace('(', '').replace(')', '')}_pressure_loss.png"
    plt.savefig(filename)
    print(f"\nPressure loss breakdown chart saved as {filename}")
    
    # Close the figure to free memory
    plt.close()

# Generate pressure loss charts
plot_pressure_losses(fan1_straight, fan1_fittings, "Fan 1 (Ground Floor & Cafeteria)")
plot_pressure_losses(fan2_straight, fan2_fittings, "Fan 2 (Floors 1 & 2)")

print("\n=== Final Recommendations ===")
print("Based on the calculations for the bathroom ventilation system:")

if fan1_velocities_ok:
    print(f"1. For Fan 1 (Ground Floor & Cafeteria): {fan1_recommendation}")
    print(f"   - Required capacity: {fan1_total_cfm * std_correction:.2f} CFM at {fan1_pressure_inwg:.4f} inWG")
else:
    print("1. Fan 1 system needs duct size revisions to maintain proper velocities.")

if fan2_velocities_ok:
    print(f"2. For Fan 2 (Floors 1 & 2): {fan2_recommendation}")
    print(f"   - Required capacity: {fan2_total_cfm * std_correction:.2f} CFM at {fan2_pressure_inwg:.4f} inWG")
else:
    print("2. Fan 2 system needs duct size revisions to maintain proper velocities.")

print("\nNotes:")
print("- All fans should be centrifugal type with backward-curved blades for efficient operation.")
print("- Fans should be installed with vibration isolators to minimize noise transmission.")
print("- Fan motors should be TEFC (Totally Enclosed Fan Cooled) rated for continuous operation.")
print("- Consider variable speed drives for energy savings during periods of lower occupancy.")
print("- Ensure all duct connections are properly sealed to prevent air leakage.")
print("- Install backdraft dampers to prevent reverse airflow when fans are not operating.")

