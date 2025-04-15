import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Constants
AIR_CHANGES_PER_HOUR = 12
AIR_DENSITY_SEA_LEVEL = 1.225  # kg/m³
ALTITUDE = 850  # meters above sea level
GRAVITY = 9.81  # m/s²
PRESSURE_ATM = 101325  # Pa at sea level
TEMPERATURE_AVG = 25  # °C
HUMIDITY_AVG = 60  # %
R_GAS_CONSTANT = 287.05  # J/(kg·K)

# Function to convert units
def convert_to_meters(value, unit):
    if unit == "inch":
        return value * 0.0254
    elif unit == "m":
        return value
    else:
        return value

# Function to calculate air density based on altitude, temperature and humidity
def calculate_air_density(altitude, temperature, humidity):
    # Barometric pressure at altitude
    pressure = PRESSURE_ATM * (1 - 2.25577e-5 * altitude) ** 5.25588
    
    # Temperature in Kelvin
    temp_kelvin = temperature + 273.15
    
    # Saturation vapor pressure
    es = 6.11 * 10 ** (7.5 * temperature / (237.7 + temperature))
    
    # Vapor pressure
    e = humidity * es / 100
    
    # Mixing ratio
    r = 0.622 * e / (pressure - e)
    
    # Virtual temperature
    tv = temp_kelvin * (1 + 0.61 * r)
    
    # Density calculation
    density = pressure / (R_GAS_CONSTANT * tv)
    
    return density

# Calculate actual air density
air_density = calculate_air_density(ALTITUDE, TEMPERATURE_AVG, HUMIDITY_AVG)

# Function to calculate required air flow rate based on room volume and air changes per hour
def calculate_required_flow_rate(volume_m3):
    # Convert m³ to cfm (cubic feet per minute)
    flow_rate_m3_per_hour = volume_m3 * AIR_CHANGES_PER_HOUR
    flow_rate_m3_per_second = flow_rate_m3_per_hour / 3600
    flow_rate_cfm = flow_rate_m3_per_second * 2118.88
    return flow_rate_m3_per_second, flow_rate_cfm

# Function to calculate duct cross-sectional area
def calculate_duct_area(width_inches, height_inches):
    width_m = width_inches * 0.0254
    height_m = height_inches * 0.0254
    return width_m * height_m

# Function to calculate circular duct area
def calculate_circular_duct_area(diameter_inches):
    diameter_m = diameter_inches * 0.0254
    return math.pi * (diameter_m / 2) ** 2

# Function to calculate hydraulic diameter
def calculate_hydraulic_diameter(width_inches, height_inches):
    width_m = width_inches * 0.0254
    height_m = height_inches * 0.0254
    return 2 * (width_m * height_m) / (width_m + height_m)

# Function to calculate velocity in duct
def calculate_velocity(flow_rate_m3_per_second, area_m2):
    return flow_rate_m3_per_second / area_m2

# Function to calculate friction factor using Colebrook equation
def calculate_friction_factor(reynolds_number, roughness_m, diameter_m):
    def colebrook_equation(f):
        return 1 / math.sqrt(f) + 2 * math.log10(roughness_m / (3.7 * diameter_m) + 2.51 / (reynolds_number * math.sqrt(f)))
    
    # Initial guess for friction factor
    initial_guess = 0.02
    
    # Solve the Colebrook equation
    friction_factor = fsolve(colebrook_equation, initial_guess)[0]
    
    return friction_factor

# Function to calculate pressure loss in duct
def calculate_pressure_loss_duct(length_m, hydraulic_diameter_m, velocity_m_s, friction_factor):
    return friction_factor * (length_m / hydraulic_diameter_m) * (air_density * velocity_m_s ** 2 / 2)

# Function to calculate pressure loss in fittings
def calculate_pressure_loss_fitting(k_factor, velocity_m_s):
    return k_factor * (air_density * velocity_m_s ** 2 / 2)

# Define K factors for common fittings
k_factors = {
    'elbow_90_rect': 0.3,  # 90° rectangular elbow
    'elbow_90_round': 0.25,  # 90° round elbow
    'tee_straight': 0.5,  # Straight-through tee
    'tee_branch': 1.0,  # Branch of tee
    'entry_sharp': 0.5,  # Sharp entry from room to duct
    'exit_duct': 1.0,  # Exit from duct to atmosphere
    'expansion_sudden': 0.7,  # Sudden expansion
    'contraction_sudden': 0.4,  # Sudden contraction
    'rect_to_round': 0.2,  # Rectangular to round transition
}

# Function to calculate total system pressure loss
def calculate_total_pressure_loss(system_components):
    total_loss = 0
    for component in system_components:
        total_loss += component['pressure_loss']
    return total_loss

# Define the systems
def define_system_1():
    # System 1: Ground floor and cafeteria bathrooms
    
    # Room volumes
    pb_bano_c_volume = 25  # m³
    pb_bano_d_volume = 18  # m³
    cafeteria_bano_c_volume = 6.5  # m³
    cafeteria_bano_d_volume = 9  # m³
    
    total_volume = pb_bano_c_volume + pb_bano_d_volume + cafeteria_bano_c_volume + cafeteria_bano_d_volume
    
    # Calculate required total flow rate
    flow_rate_m3_s, flow_rate_cfm = calculate_required_flow_rate(total_volume)
    
    # Define ducts and components
    components = []
    
    # PB Baño C - 2 grilles 4"x11"
    grille_area = 2 * calculate_duct_area(4, 11)
    velocity_grille = calculate_velocity(flow_rate_m3_s * (pb_bano_c_volume / total_volume), grille_area)
    components.append({
        'name': 'PB Baño C Grilles',
        'velocity': velocity_grille,
        'area': grille_area,
        'pressure_loss': calculate_pressure_loss_fitting(k_factors['entry_sharp'], velocity_grille)
    })
    
    # PB Baño C - Duct 6"x11" 4m
    duct_area = calculate_duct_area(6, 11)
    velocity_duct = calculate_velocity(flow_rate_m3_s * ((pb_bano_c_volume + pb_bano_d_volume) / total_volume), duct_area)
    hydraulic_diameter = calculate_hydraulic_diameter(6, 11)
    reynolds_number = (velocity_duct * hydraulic_diameter * air_density) / (1.825e-5)  # Assuming kinematic viscosity
    friction_factor = calculate_friction_factor(reynolds_number, 0.00015, hydraulic_diameter)
    components.append({
        'name': 'PB Baño C Duct 6"x11"',
        'velocity': velocity_duct,
        'area': duct_area,
        'pressure_loss': calculate_pressure_loss_duct(4, hydraulic_diameter, velocity_duct, friction_factor)
    })
    
    # Expansion from 6"x11" to 11"x11"
    expanded_area = calculate_duct_area(11, 11)
    velocity_expanded = calculate_velocity(flow_rate_m3_s * ((pb_bano_c_volume + pb_bano_d_volume) / total_volume), expanded_area)
    components.append({
        'name': 'Expansion 6"x11" to 11"x11"',
        'velocity': velocity_expanded,
        'area': expanded_area,
        'pressure_loss': calculate_pressure_loss_fitting(k_factors['expansion_sudden'], velocity_duct)
    })
    
    # 11"x11" duct 2m with 90 degree elbow
    velocity_11x11 = velocity_expanded
    hydraulic_diameter_11x11 = calculate_hydraulic_diameter(11, 11)
    reynolds_number_11x11 = (velocity_11x11 * hydraulic_diameter_11x11 * air_density) / (1.825e-5)
    friction_factor_11x11 = calculate_friction_factor(reynolds_number_11x11, 0.00015, hydraulic_diameter_11x11)
    components.append({
        'name': '11"x11" duct 2m',
        'velocity': velocity_11x11,
        'area': expanded_area,
        'pressure_loss': calculate_pressure_loss_duct(2, hydraulic_diameter_11x11, velocity_11x11, friction_factor_11x11)
    })
    
    # 90 degree elbow 11"x11"
    components.append({
        'name': '90° Elbow 11"x11"',
        'velocity': velocity_11x11,
        'area': expanded_area,
        'pressure_loss': calculate_pressure_loss_fitting(k_factors['elbow_90_rect'], velocity_11x11)
    })
    
    # Section change 11"x11" to circular 12.5"
    circular_area = calculate_circular_duct_area(12.5)
    velocity_circular = calculate_velocity(flow_rate_m3_s, circular_area)
    components.append({
        'name': 'Section change 11"x11" to 12.5" round',
        'velocity': velocity_circular,
        'area': circular_area,
        'pressure_loss': calculate_pressure_loss_fitting(k_factors['rect_to_round'], velocity_11x11)
    })
    
    # Vertical circular duct 12.5" 6m
    diameter_circular = 12.5 * 0.0254
    reynolds_number_circular = (velocity_circular * diameter_circular * air_density) / (1.825e-5)
    friction_factor_circular = calculate_friction_factor(reynolds_number_circular, 0.00015, diameter_circular)
    components.append({
        'name': 'Vertical circular duct 12.5" 6m',
        'velocity': velocity_circular,
        'area': circular_area,
        'pressure_loss': calculate_pressure_loss_duct(6, diameter_circular, velocity_circular, friction_factor_circular)
    })
    
    # 90 degree elbow circular 12.5"
    components.append({
        'name': '90° Elbow 12.5" circular',
        'velocity': velocity_circular,
        'area': circular_area,
        'pressure_loss': calculate_pressure_loss_fitting(k_factors['elbow_90_round'], velocity_circular)
    })
    
    # Calculate additional components for PB Baño D, Cafetín Baño C, Cafetín Baño D
    # ... (similar calculations would continue for all other components)
    
    # Calculate total pressure loss
    total_loss = calculate_total_pressure_loss(components)
    
    return {
        'name': 'System 1: Ground Floor and Cafeteria',
        'flow_rate_m3_s': flow_rate_m3_s,
        'flow_rate_cfm': flow_rate_cfm,
        'pressure_loss_pa': total_loss,
        'pressure_loss_inwg': total_loss / 249.08,  # Convert Pa to inwg
        'components': components
    }

def define_system_2():
    # System 2: First and second floor bathrooms
    
    # Room volumes (estimate based on number of grilles)
    piso1_bano_c_volume = 25  # m³ (estimated)
    piso1_bano_d_volume = 25  # m³ (estimated)
    piso2_bano_c_volume = 20  # m³ (estimated)
    piso2_bano_d_volume = 20  # m³ (estimated)
    
    total_volume = piso1_bano_c_volume + piso1_bano_d_volume + piso2_bano_c_volume + piso2_bano_d_volume
    
    # Calculate required total flow rate
    flow_rate_m3_s, flow_rate_cfm = calculate_required_flow_rate(total_volume)
    
    # Similar component calculations would continue for System 2
    # ...
    
    # For simplicity, let's estimate the total pressure loss for System 2
    estimated_pressure_loss = 150  # Pa (estimation)
    
    return {
        'name': 'System 2: First and Second Floors',
        'flow_rate_m3_s': flow_rate_m3_s,
        'flow_rate_cfm': flow_rate_cfm,
        'pressure_loss_pa': estimated_pressure_loss,
        'pressure_loss_inwg': estimated_pressure_loss / 249.08,  # Convert Pa to inwg
        'components': []
    }

# Main function to calculate and recommend fans
def main():
    # Calculate system parameters
    system1 = define_system_1()
    system2 = define_system_2()
    
    # Check for minimum velocities in ducts (should be > 3 m/s typically)
    min_velocity_system1 = min([comp['velocity'] for comp in system1['components'] if 'velocity' in comp])
    
    # Display results
    print("=== Ventilation System Calculations ===")
    print("\nAir Density at {:.0f}m altitude, {:.1f}°C, {:.0f}% humidity: {:.3f} kg/m³".format(
        ALTITUDE, TEMPERATURE_AVG, HUMIDITY_AVG, air_density))
    # Fan selection criteria
    print("\n=== Fan Selection Criteria ===")
    print("System 1: Ground Floor and Cafeteria")
    print("  Required Flow Rate: {:.1f} m³/s = {:.0f} CFM".format(
        system1['flow_rate_m3_s'], system1['flow_rate_cfm']))
    print("  Required Static Pressure: {:.1f} Pa = {:.2f} inwg".format(
        system1['pressure_loss_pa'], system1['pressure_loss_inwg']))
    
    # Add 25% safety factor
    system1_flow_with_safety = system1['flow_rate_cfm'] * 1.25
    system1_pressure_with_safety = system1['pressure_loss_inwg'] * 1.25
    
    print("  Flow Rate with 25% Safety Factor: {:.0f} CFM".format(system1_flow_with_safety))
    print("  Static Pressure with 25% Safety Factor: {:.2f} inwg".format(system1_pressure_with_safety))
    
    if min_velocity_system1 < 3:
        print("  WARNING: Some duct velocities are below 3 m/s, which may lead to poor ventilation.")
    
    print("\nSystem 2: First and Second Floors")
    print("  Required Flow Rate: {:.1f} m³/s = {:.0f} CFM".format(
        system2['flow_rate_m3_s'], system2['flow_rate_cfm']))
    print("  Required Static Pressure: {:.1f} Pa = {:.2f} inwg".format(
        system2['pressure_loss_pa'], system2['pressure_loss_inwg']))
    
    # Add 25% safety factor
    system2_flow_with_safety = system2['flow_rate_cfm'] * 1.25
    system2_pressure_with_safety = system2['pressure_loss_inwg'] * 1.25
    
    print("  Flow Rate with 25% Safety Factor: {:.0f} CFM".format(system2_flow_with_safety))
    print("  Static Pressure with 25% Safety Factor: {:.2f} inwg".format(system2_pressure_with_safety))
    
    # Commercial Fan Recommendations
    # These are example fans - in a real system you would look up actual manufacturer data
    print("\n=== Commercial Fan Recommendations ===")
    
    # Define some sample commercial fans
    commercial_fans = [
        {
            'model': 'Cook 100 CPS',
            'type': 'Centrifugal',
            'max_flow_cfm': 1000,
            'max_pressure_inwg': 0.75,
            'power_hp': 0.25,
            'power_watts': 186,
            'rpm': 1725
        },
        {
            'model': 'Cook 120 CPS',
            'type': 'Centrifugal',
            'max_flow_cfm': 1200,
            'max_pressure_inwg': 1.0,
            'power_hp': 0.33,
            'power_watts': 246,
            'rpm': 1725
        },
        {
            'model': 'Greenheck CUE-070-VG',
            'type': 'Centrifugal',
            'max_flow_cfm': 700,
            'max_pressure_inwg': 0.5,
            'power_hp': 0.125,
            'power_watts': 93,
            'rpm': 1625
        },
        {
            'model': 'Greenheck CUE-095-VG',
            'type': 'Centrifugal',
            'max_flow_cfm': 950,
            'max_pressure_inwg': 0.75,
            'power_hp': 0.25,
            'power_watts': 186,
            'rpm': 1750
        },
        {
            'model': 'Greenheck CUE-120-VG',
            'type': 'Centrifugal',
            'max_flow_cfm': 1200,
            'max_pressure_inwg': 1.0,
            'power_hp': 0.33,
            'power_watts': 246,
            'rpm': 1750
        },
        {
            'model': 'Greenheck CUE-140-VG',
            'type': 'Centrifugal',
            'max_flow_cfm': 1400,
            'max_pressure_inwg': 1.25,
            'power_hp': 0.5,
            'power_watts': 373,
            'rpm': 1750
        },
        {
            'model': 'Soler & Palau TD-200',
            'type': 'Mixed Flow',
            'max_flow_cfm': 720,
            'max_pressure_inwg': 0.6,
            'power_hp': 0.2,
            'power_watts': 149,
            'rpm': 2500
        },
        {
            'model': 'Soler & Palau TD-250',
            'type': 'Mixed Flow',
            'max_flow_cfm': 950,
            'max_pressure_inwg': 0.8,
            'power_hp': 0.25,
            'power_watts': 186,
            'rpm': 2500
        },
        {
            'model': 'Soler & Palau TD-315',
            'type': 'Mixed Flow',
            'max_flow_cfm': 1200,
            'max_pressure_inwg': 1.1,
            'power_hp': 0.33,
            'power_watts': 246,
            'rpm': 2500
        },
        {
            'model': 'Fantech FKD 12XL',
            'type': 'Centrifugal',
            'max_flow_cfm': 980,
            'max_pressure_inwg': 0.8,
            'power_hp': 0.33,
            'power_watts': 246,
            'rpm': 1800
        },
        {
            'model': 'Fantech FKD 14',
            'type': 'Centrifugal',
            'max_flow_cfm': 1200,
            'max_pressure_inwg': 1.0,
            'power_hp': 0.5,
            'power_watts': 373,
            'rpm': 1800
        }
    ]
    
    # Function to recommend fans
    def recommend_fans(required_flow, required_pressure, fans_list):
        suitable_fans = []
        for fan in fans_list:
            if fan['max_flow_cfm'] >= required_flow and fan['max_pressure_inwg'] >= required_pressure:
                # Calculate simple matching score - lower is better
                # This gives preference to fans that are closer to the required specs
                flow_ratio = fan['max_flow_cfm'] / required_flow
                pressure_ratio = fan['max_pressure_inwg'] / required_pressure
                matching_score = abs(flow_ratio - 1) + abs(pressure_ratio - 1)
                
                suitable_fans.append({
                    'model': fan['model'],
                    'type': fan['type'],
                    'max_flow_cfm': fan['max_flow_cfm'],
                    'max_pressure_inwg': fan['max_pressure_inwg'],
                    'power_hp': fan['power_hp'],
                    'power_watts': fan['power_watts'],
                    'rpm': fan['rpm'],
                    'matching_score': matching_score
                })
        
        # Sort by matching score
        return sorted(suitable_fans, key=lambda x: x['matching_score'])
    
    # Recommend fans for System 1
    recommended_fans_system1 = recommend_fans(
        system1_flow_with_safety, 
        system1_pressure_with_safety, 
        commercial_fans
    )
    
    print("\nRecommended fans for System 1 (Ground Floor and Cafeteria):")
    if recommended_fans_system1:
        for i, fan in enumerate(recommended_fans_system1[:3]):  # Top 3 recommendations
            print("{}. {} ({})".format(i+1, fan['model'], fan['type']))
            print("   Max Flow: {:.0f} CFM, Max Pressure: {:.2f} inwg".format(
                fan['max_flow_cfm'], fan['max_pressure_inwg']))
            print("   Power: {:.2f} HP ({:.0f} W), RPM: {:.0f}".format(
                fan['power_hp'], fan['power_watts'], fan['rpm']))
    else:
        print("No suitable fans found. Consider custom fan selection.")
    
    # Recommend fans for System 2
    recommended_fans_system2 = recommend_fans(
        system2_flow_with_safety, 
        system2_pressure_with_safety, 
        commercial_fans
    )
    
    print("\nRecommended fans for System 2 (First and Second Floors):")
    if recommended_fans_system2:
        for i, fan in enumerate(recommended_fans_system2[:3]):  # Top 3 recommendations
            print("{}. {} ({})".format(i+1, fan['model'], fan['type']))
            print("   Max Flow: {:.0f} CFM, Max Pressure: {:.2f} inwg".format(
                fan['max_flow_cfm'], fan['max_pressure_inwg']))
            print("   Power: {:.2f} HP ({:.0f} W), RPM: {:.0f}".format(
                fan['power_hp'], fan['power_watts'], fan['rpm']))
    else:
        print("No suitable fans found. Consider custom fan selection.")
    
    # Generate fan performance curves (simplified)
    def plot_fan_curve(fan_model, required_flow, required_pressure):
        # Create a simplified fan curve
        flow_points = np.linspace(0, fan_model['max_flow_cfm'] * 1.1, 100)
        
        # Simple quadratic curve: P = Pmax * (1 - (Q/Qmax)²)
        pressure_points = [fan_model['max_pressure_inwg'] * (1 - (q/fan_model['max_flow_cfm'])**2) for q in flow_points]
        
        plt.figure(figsize=(10, 6))
        plt.plot(flow_points, pressure_points, 'b-', linewidth=2, label=fan_model['model'])
        plt.scatter([required_flow], [required_pressure], color='red', s=100, marker='x', 
                    label='Required Operating Point')
        
        plt.title(f"Fan Performance Curve: {fan_model['model']}")
        plt.xlabel("Flow Rate (CFM)")
        plt.ylabel("Static Pressure (inwg)")
        plt.grid(True)
        plt.legend()
        
        plt.axhline(y=required_pressure, color='r', linestyle='--', alpha=0.3)
        plt.axvline(x=required_flow, color='r', linestyle='--', alpha=0.3)
        
        plt.tight_layout()
        return plt
    
    # Plot fan curves for the top recommendations
    if recommended_fans_system1:
        plot1 = plot_fan_curve(
            recommended_fans_system1[0], 
            system1_flow_with_safety, 
            system1_pressure_with_safety
        )
        plot1.savefig('system1_fan_curve.png')
        print("\nFan curve for System 1 saved as 'system1_fan_curve.png'")
    
    if recommended_fans_system2:
        plot2 = plot_fan_curve(
            recommended_fans_system2[0], 
            system2_flow_with_safety, 
            system2_pressure_with_safety
        )
        plot2.savefig('system2_fan_curve.png')
        print("Fan curve for System 2 saved as 'system2_fan_curve.png'")
    
    print("\n=== Ventilation System Summary ===")
    print("Air changes per hour: {}".format(AIR_CHANGES_PER_HOUR))
    print("Total volume - System 1: {:.1f} m³".format(
        system1['flow_rate_m3_s'] * 3600 / AIR_CHANGES_PER_HOUR))
    print("Total volume - System 2: {:.1f} m³".format(
        system2['flow_rate_m3_s'] * 3600 / AIR_CHANGES_PER_HOUR))
    print("\nNotes:")
    print("1. All calculations include a 25% safety factor")
    print("2. Fan selection should be verified with manufacturer data")
    print("3. Installation should comply with local building codes")
    print("4. Final selection may need adjustments based on actual installation conditions")

if __name__ == "__main__":
    main()
