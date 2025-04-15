import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math

# Constants and conversion factors
INCH_TO_M = 0.0254
FT_TO_M = 0.3048
CFM_TO_M3S = 0.000471947
M3S_TO_CFM = 1/CFM_TO_M3S
PA_TO_INWG = 0.00402
INWG_TO_PA = 248.84

# Climate conditions
climate_conditions = {
    'average': {'temperature': 25, 'humidity': 60},
    'extreme_hot_humid': {'temperature': 28, 'humidity': 93},
    'extreme_cold_dry': {'temperature': 16, 'humidity': 30}
}

# Altitude
altitude = 850  # meters above sea level

# Air density calculation based on temperature, humidity and altitude
def air_density(temp_c, rel_humidity, alt_m):
    # Barometric pressure at altitude
    p0 = 101325  # sea level standard pressure, Pa
    g = 9.80665  # gravitational acceleration, m/s²
    M = 0.0289644  # molar mass of dry air, kg/mol
    R = 8.31447  # universal gas constant, J/(mol·K)
    T0 = 288.15  # sea level standard temperature, K
    L = 0.0065  # temperature lapse rate, K/m
    
    # Barometric formula
    temp_k = temp_c + 273.15
    p = p0 * (1 - L * alt_m / T0) ** (g * M / (R * L))
    
    # Saturation vapor pressure
    es = 611.2 * np.exp(17.67 * temp_c / (temp_c + 243.5))
    
    # Actual vapor pressure
    e = rel_humidity / 100 * es
    
    # Density of moist air
    Rd = 287.05  # specific gas constant for dry air, J/(kg·K)
    Rv = 461.495  # specific gas constant for water vapor, J/(kg·K)
    
    density = (p - e) / (Rd * temp_k) + e / (Rv * temp_k)
    
    return density

# Function to calculate duct area
def duct_area(width_inch, height_inch):
    width_m = width_inch * INCH_TO_M
    height_m = height_inch * INCH_TO_M
    return width_m * height_m

# Function to calculate circular duct area
def circular_duct_area(diameter_inch):
    radius_m = (diameter_inch / 2) * INCH_TO_M
    return np.pi * radius_m**2

# Function to calculate hydraulic diameter
def hydraulic_diameter(width_inch, height_inch):
    width_m = width_inch * INCH_TO_M
    height_m = height_inch * INCH_TO_M
    return 2 * (width_m * height_m) / (width_m + height_m)

# Function to calculate air velocity in a duct
def air_velocity(flow_rate_m3s, area_m2):
    return flow_rate_m3s / area_m2

# Function to calculate pressure loss in a duct
def pressure_loss(velocity, length, hydraulic_diam, k_factors, density):
    # Friction factor (assuming turbulent flow in commercial ducts)
    f = 0.02  
    
    # Pressure loss due to friction
    dp_friction = f * (length / hydraulic_diam) * (density * velocity**2 / 2)
    
    # Pressure loss due to fittings
    dp_fittings = sum(k_factors) * (density * velocity**2 / 2)
    
    return dp_friction + dp_fittings

# Function to calculate required air changes per hour based on room volume
def required_ach(room_volume_m3, min_velocity_m_s, grille_area_m2):
    # Calculate minimum flow rate to maintain minimum velocity at grilles
    min_flow_rate_m3s = min_velocity_m_s * grille_area_m2
    
    # Convert to air changes per hour
    ach = min_flow_rate_m3s * 3600 / room_volume_m3
    
    # Ensure minimum ACH for bathrooms (typically 8-15 ACH)
    return max(ach, 10)

# Define bathroom spaces
bathrooms = {
    'PB_Bano_C': {
        'volume_m3': 25,
        'grilles': [{'width': 4, 'height': 11, 'count': 2}],
        'ducts': [
            {'width': 6, 'height': 11, 'length': 4, 'k_factors': [0.5]},  # Main duct
            {'width': 11, 'height': 11, 'length': 2, 'k_factors': [0.5]},  # Expansion
            {'diameter': 12.5, 'length': 6, 'k_factors': [0.9, 0.9]}  # Vertical duct with 2 elbows
        ],
        'system': 'system1'
    },
    'PB_Bano_D': {
        'volume_m3': 18,
        'grilles': [{'width': 6, 'height': 11, 'count': 3}],
        'ducts': [
            {'width': 6, 'height': 11, 'length': 4, 'k_factors': [0.9, 0.5]}  # Main duct with elbow
        ],
        'connects_to': 'PB_Bano_C',
        'system': 'system1'
    },
    'Cafetin_Bano_C': {
        'volume_m3': 6.5,
        'grilles': [{'width': 7.5, 'height': 10, 'count': 1}],
        'ducts': [
            {'diameter': 10, 'length': 1.44, 'k_factors': [0.9]},  # Duct with elbow
        ],
        'connects_to': 'PB_Bano_C',
        'system': 'system1'
    },
    'Cafetin_Bano_D': {
        'volume_m3': 9,
        'grilles': [{'width': 7.5, 'height': 11, 'count': 2}],
        'ducts': [
            {'diameter': 10, 'length': 2, 'k_factors': [0.5]}  # Connecting duct
        ],
        'connects_to': 'Cafetin_Bano_C',
        'system': 'system1'
    },
    'Piso1_Bano_C': {
        'volume_m3': 10,
        'grilles': [{'width': 6, 'height': 14, 'count': 3}],
        'ducts': [
            {'width': 11, 'height': 6, 'length': 3.2, 'k_factors': [0.5]},
            {'width': 11, 'height': 6, 'length': 3.5, 'k_factors': [0.9, 1.0]}  # Vertical duct with elbow and T-junction
        ],
        'connects_to': 'Ducto_Principal_Piso2',
        'system': 'system2'
    },
    'Piso1_Bano_D': {
        'volume_m3': 10,
        'grilles': [{'width': 6, 'height': 14, 'count': 3}],
        'ducts': [
            {'width': 11, 'height': 6, 'length': 3.2, 'k_factors': [0.5]},
            {'width': 11, 'height': 6, 'length': 3.5, 'k_factors': [0.9, 1.0]}  # Vertical duct with elbow and T-junction
        ],
        'connects_to': 'Ducto_Principal_Piso2',
        'system': 'system2'
    },
    'Piso2_Bano_C': {
        'volume_m3': 10,
        'grilles': [{'width': 4, 'height': 11, 'count': 3}],
        'ducts': [
            {'width': 11, 'height': 6, 'length': 3.2, 'k_factors': [1.0]}  # Duct with T-junction
        ],
        'connects_to': 'Ducto_Principal_Piso2',
        'system': 'system2'
    },
    'Piso2_Bano_D': {
        'volume_m3': 10,
        'grilles': [{'width': 4, 'height': 11, 'count': 3}],
        'ducts': [
            {'width': 11, 'height': 6, 'length': 3.2, 'k_factors': [1.0]}  # Duct with T-junction
        ],
        'connects_to': 'Ducto_Principal_Piso2',
        'system': 'system2'
    },
    'Ducto_Principal_Piso2': {
        'volume_m3': 0,  # Not a room, just a duct
        'grilles': [],
        'ducts': [
            {'width': 11, 'height': 11, 'length': 12, 'k_factors': [0.5]},
            {'diameter': 12, 'length': 1, 'k_factors': [0.5]}  # Transition to circular
        ],
        'system': 'system2'
    }
}

# Minimum and maximum recommended velocities (m/s)
min_velocity_grille = 2.5  # m/s
max_velocity_grille = 4.0  # m/s
min_velocity_duct = 3.0    # m/s
max_velocity_duct = 10.0   # m/s

def analyze_ventilation_system(climate_condition='extreme_hot_humid'):
    # Get climate parameters
    temp = climate_conditions[climate_condition]['temperature']
    humidity = climate_conditions[climate_condition]['humidity']
    
    # Calculate air density
    density = air_density(temp, humidity, altitude)
    print(f"Air density at {temp}°C, {humidity}% RH, {altitude}m altitude: {density:.4f} kg/m³")
    
    # Calculate total grille area and required flow rate for each bathroom
    results = {}
    system_flow_rates = {'system1': 0, 'system2': 0}
    system_pressure_losses = {'system1': 0, 'system2': 0}
    
    for bathroom_name, bathroom in bathrooms.items():
        if bathroom_name == 'Ducto_Principal_Piso2':
            continue  # Skip the main duct for now
            
        # Calculate total grille area
        total_grille_area = 0
        for grille in bathroom['grilles']:
            grille_area = duct_area(grille['width'], grille['height']) * grille['count']
            total_grille_area += grille_area
        
        # Calculate required air changes per hour
        ach = required_ach(bathroom['volume_m3'], min_velocity_grille, total_grille_area)
        
        # Calculate required flow rate
        flow_rate_m3h = ach * bathroom['volume_m3']
        flow_rate_m3s = flow_rate_m3h / 3600
        
        # Calculate grille velocity
        grille_velocity = flow_rate_m3s / total_grille_area if total_grille_area > 0 else 0
        
        # Adjust flow rate if grille velocity is too low
        if grille_velocity < min_velocity_grille and total_grille_area > 0:
            flow_rate_m3s = min_velocity_grille * total_grille_area
            flow_rate_m3h = flow_rate_m3s * 3600
            ach = flow_rate_m3h / bathroom['volume_m3']
            grille_velocity = min_velocity_grille
        
        # Calculate duct velocities and pressure losses
        duct_velocities = []
        pressure_losses = []
        
        for duct in bathroom['ducts']:
            if 'diameter' in duct:
                area = circular_duct_area(duct['diameter'])
                hyd_diam = duct['diameter'] * INCH_TO_M
            else:
                area = duct_area(duct['width'], duct['height'])
                hyd_diam = hydraulic_diameter(duct['width'], duct['height'])
            
            velocity = air_velocity(flow_rate_m3s, area)
            duct_velocities.append(velocity)
            
            loss = pressure_loss(velocity, duct['length'], hyd_diam, duct['k_factors'], density)
            pressure_losses.append(loss)
        
        # Store results
        results[bathroom_name] = {
            'volume_m3': bathroom['volume_m3'],
            'total_grille_area_m2': total_grille_area,
            'air_changes_per_hour': ach,
            'flow_rate_m3h': flow_rate_m3h,
            'flow_rate_m3s': flow_rate_m3s,
            'flow_rate_cfm': flow_rate_m3s * M3S_TO_CFM,
            'grille_velocity_m_s': grille_velocity,
            'duct_velocities_m_s': duct_velocities,
            'pressure_losses_pa': pressure_losses,
            'total_pressure_loss_pa': sum(pressure_losses),
            'system': bathroom['system']
        }
        
        # Add to system totals
        system_flow_rates[bathroom['system']] += flow_rate_m3s
        
        # For pressure loss, we need to consider the path
        if 'connects_to' not in bathroom:
            system_pressure_losses[bathroom['system']] = max(
                system_pressure_losses[bathroom['system']],
                sum(pressure_losses)
            )
    
    # Now handle the main duct for system2
    main_duct = bathrooms['Ducto_Principal_Piso2']
    flow_rate_m3s = system_flow_rates['system2']
    
    duct_velocities = []
    pressure_losses = []
    
    for duct in main_duct['ducts']:
        if 'diameter' in duct:
            area = circular_duct_area(duct['diameter'])
            hyd_diam = duct['diameter'] * INCH_TO_M
        else:
            area = duct_area(duct['width'], duct['height'])
            hyd_diam = hydraulic_diameter(duct['width'], duct['height'])
        
        velocity = air_velocity(flow_rate_m3s, area)
        duct_velocities.append(velocity)
        
        loss = pressure_loss(velocity, duct['length'], hyd_diam, duct['k_factors'], density)
        pressure_losses.append(loss)
    
    # Store results for main duct
    results['Ducto_Principal_Piso2'] = {
        'flow_rate_m3s': flow_rate_m3s,
        'flow_rate_m3h': flow_rate_m3s * 3600,  # Make sure to include this key
        'flow_rate_cfm': flow_rate_m3s * M3S_TO_CFM,
        'duct_velocities_m_s': duct_velocities,
        'pressure_losses_pa': pressure_losses,
        'total_pressure_loss_pa': sum(pressure_losses),
        'system': 'system2'
    }
    
    # Add main duct pressure loss to system2
    system_pressure_losses['system2'] += sum(pressure_losses)
    
    # Calculate fan requirements
    fan_requirements = {}
    for system, flow_rate in system_flow_rates.items():
        fan_requirements[system] = {
            'flow_rate_m3s': flow_rate,
            'flow_rate_m3h': flow_rate * 3600,
            'flow_rate_cfm': flow_rate * M3S_TO_CFM,
            'pressure_pa': system_pressure_losses[system],
            'pressure_inwg': system_pressure_losses[system] * PA_TO_INWG
        }
    
    return results, fan_requirements

def suggest_improvements(results):
    suggestions = {}
    
    for bathroom_name, data in results.items():
        if bathroom_name == 'Ducto_Principal_Piso2':
            continue
            
        suggestions[bathroom_name] = []
        
        # Check grille velocity
        if 'grille_velocity_m_s' in data:
            if data['grille_velocity_m_s'] < min_velocity_grille:
                suggestions[bathroom_name].append(
                    f"Grille velocity ({data['grille_velocity_m_s']:.2f} m/s) is below minimum ({min_velocity_grille} m/s). "
                    f"Consider reducing grille size or increasing air changes per hour."
                )
            elif data['grille_velocity_m_s'] > max_velocity_grille:
                suggestions[bathroom_name].append(
                    f"Grille velocity ({data['grille_velocity_m_s']:.2f} m/s) is above maximum ({max_velocity_grille} m/s). "
                    f"Consider increasing grille size or reducing air changes per hour."
                )
        
        # Check duct velocities
        if 'duct_velocities_m_s' in data:
            for i, velocity in enumerate(data['duct_velocities_m_s']):
                if velocity < min_velocity_duct:
                    suggestions[bathroom_name].append(
                        f"Duct {i+1} velocity ({velocity:.2f} m/s) is below minimum ({min_velocity_duct} m/s). "
                        f"Consider reducing duct size or increasing flow rate."
                    )
                elif velocity > max_velocity_duct:
                    suggestions[bathroom_name].append(
                        f"Duct {i+1} velocity ({velocity:.2f} m/s) is above maximum ({max_velocity_duct} m/s). "
                        f"Consider increasing duct size or reducing flow rate."
                    )
    
    return suggestions

def recommend_fans(fan_requirements):
    fan_recommendations = {}
    
    # Define some commercial centrifugal fan models with their specifications
    commercial_fans = {
        'CF-100': {'max_flow_cfm': 100, 'max_pressure_inwg': 0.5, 'power_watts': 35},
        'CF-200': {'max_flow_cfm': 200, 'max_pressure_inwg': 0.75, 'power_watts': 70},
        'CF-300': {'max_flow_cfm': 300, 'max_pressure_inwg': 1.0, 'power_watts': 120},
        'CF-500': {'max_flow_cfm': 500, 'max_pressure_inwg': 1.25, 'power_watts': 180},
        'CF-750': {'max_flow_cfm': 750, 'max_pressure_inwg': 1.5, 'power_watts': 250},
        'CF-1000': {'max_flow_cfm': 1000, 'max_pressure_inwg': 2.0, 'power_watts': 370},
        'CF-1500': {'max_flow_cfm': 1500, 'max_pressure_inwg': 2.5, 'power_watts': 550},
        'CF-2000': {'max_flow_cfm': 2000, 'max_pressure_inwg': 3.0, 'power_watts': 750}
    }
    
    for system, requirements in fan_requirements.items():
        suitable_fans = []
        
        for fan_model, specs in commercial_fans.items():
            # Add 20% safety factor
            required_flow = requirements['flow_rate_cfm'] * 1.2
            required_pressure = requirements['pressure_inwg'] * 1.2
            
            if specs['max_flow_cfm'] >= required_flow and specs['max_pressure_inwg'] >= required_pressure:
                suitable_fans.append({
                    'model': fan_model,
                    'max_flow_cfm': specs['max_flow_cfm'],
                    'max_pressure_inwg': specs['max_pressure_inwg'],
                    'power_watts': specs['power_watts'],
                    'flow_margin': (specs['max_flow_cfm'] - required_flow) / required_flow * 100,
                    'pressure_margin': (specs['max_pressure_inwg'] - required_pressure) / required_pressure * 100
                })
        
        # Sort by closest match (lowest combined margin)
        if suitable_fans:
            suitable_fans.sort(key=lambda x: x['flow_margin'] + x['pressure_margin'])
            fan_recommendations[system] = suitable_fans[0]
        else:
            fan_recommendations[system] = "No suitable fan found. Consider custom solutions."
    
    return fan_recommendations

def generate_report(results, fan_requirements, suggestions, fan_recommendations, climate_condition):
    print("\n===== VENTILATION SYSTEM ANALYSIS REPORT =====")
    print(f"Climate condition: {climate_condition}")
    print("\n--- BATHROOM DETAILS ---")
    
    for bathroom_name, data in results.items():
        print(f"\n{bathroom_name}:")
        if bathroom_name != 'Ducto_Principal_Piso2' and 'volume_m3' in data:
            print(f"  Volume: {data['volume_m3']:.1f} m³")
            if 'air_changes_per_hour' in data:
                print(f"  Air changes per hour: {data['air_changes_per_hour']:.1f}")
        
        if 'flow_rate_m3h' in data and 'flow_rate_cfm' in data:
            print(f"  Flow rate: {data['flow_rate_m3h']:.1f} m³/h ({data['flow_rate_cfm']:.1f} CFM)")
        
        if bathroom_name != 'Ducto_Principal_Piso2' and 'grille_velocity_m_s' in data:
            print(f"  Grille velocity: {data['grille_velocity_m_s']:.2f} m/s")
        
        if 'duct_velocities_m_s' in data:
            print(f"  Duct velocities:")
            for i, velocity in enumerate(data['duct_velocities_m_s']):
                print(f"    Duct {i+1}: {velocity:.2f} m/s")
        
        if 'total_pressure_loss_pa' in data:
            print(f"  Total pressure loss: {data['total_pressure_loss_pa']:.1f} Pa ({data['total_pressure_loss_pa'] * PA_TO_INWG:.4f} inWG)")
    
    print("\n--- SYSTEM REQUIREMENTS ---")
    for system, requirements in fan_requirements.items():
        print(f"\n{system.upper()}:")
        print(f"  Total flow rate: {requirements['flow_rate_m3h']:.1f} m³/h ({requirements['flow_rate_cfm']:.1f} CFM)")
        print(f"  Total pressure: {requirements['pressure_pa']:.1f} Pa ({requirements['pressure_inwg']:.4f} inWG)")
    
    print("\n--- IMPROVEMENT SUGGESTIONS ---")
    for bathroom_name, bathroom_suggestions in suggestions.items():
        if bathroom_suggestions:
            print(f"\n{bathroom_name}:")
            for suggestion in bathroom_suggestions:
                print(f"  - {suggestion}")
    
    print("\n--- FAN RECOMMENDATIONS ---")
    for system, recommendation in fan_recommendations.items():
        print(f"\n{system.upper()}:")
        if isinstance(recommendation, dict):
            print(f"  Recommended model: {recommendation['model']}")
            print(f"  Maximum flow: {recommendation['max_flow_cfm']:.1f} CFM")
            print(f"  Maximum pressure: {recommendation['max_pressure_inwg']:.4f} inWG")
            print(f"  Power consumption: {recommendation['power_watts']} W")
            print(f"  Flow margin: {recommendation['flow_margin']:.1f}%")
            print(f"  Pressure margin: {recommendation['pressure_margin']:.1f}%")
        else:
            print(f"  {recommendation}")

def plot_system_diagram(results, fan_recommendations):
    # Create a new figure with a larger size
    plt.figure(figsize=(15, 10))
    
    # Create two subplots for the two systems
    plt.subplot(1, 2, 1)
    plt.title('System 1: Ground Floor & Cafeteria')
    
    # Plot system 1 components
    system1_bathrooms = [name for name, data in results.items() if data['system'] == 'system1']
    
    # Simple representation with boxes and lines
    y_pos = 0
    for bathroom in system1_bathrooms:
        plt.text(0.1, y_pos, bathroom, fontsize=10, bbox=dict(facecolor='lightblue', alpha=0.5))
        y_pos += 1
    
    # Add fan at the top
    if isinstance(fan_recommendations['system1'], dict):
        fan_text = f"Fan: {fan_recommendations['system1']['model']}\n{fan_recommendations['system1']['max_flow_cfm']} CFM, {fan_recommendations['system1']['max_pressure_inwg']:.2f} inWG"
    else:
        fan_text = "Fan: Not found"
    
    plt.text(0.5, y_pos + 1, fan_text, fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
    
    # Connect bathrooms to fan
    for i in range(len(system1_bathrooms)):
        plt.plot([0.3, 0.5], [i, y_pos + 1], 'k-')
    
    plt.axis('off')
    
    # System 2
    plt.subplot(1, 2, 2)
    plt.title('System 2: Floor 1 & Floor 2')
    
    # Plot system 2 components
    system2_bathrooms = [name for name, data in results.items() if data['system'] == 'system2' and name != 'Ducto_Principal_Piso2']
    
    # Simple representation with boxes and lines
    y_pos = 0
    for bathroom in system2_bathrooms:
        plt.text(0.1, y_pos, bathroom, fontsize=10, bbox=dict(facecolor='lightblue', alpha=0.5))
        y_pos += 1
    
    # Add main duct
    plt.text(0.5, y_pos/2, 'Ducto_Principal_Piso2', fontsize=10, bbox=dict(facecolor='lightyellow', alpha=0.5))
    
    # Add fan at the top
    if isinstance(fan_recommendations['system2'], dict):
        fan_text = f"Fan: {fan_recommendations['system2']['model']}\n{fan_recommendations['system2']['max_flow_cfm']} CFM, {fan_recommendations['system2']['max_pressure_inwg']:.2f} inWG"
    else:
        fan_text = "Fan: Not found"
    
    plt.text(0.5, y_pos + 1, fan_text, fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
    
    # Connect bathrooms to main duct
    for i in range(len(system2_bathrooms)):
        plt.plot([0.3, 0.5], [i, y_pos/2], 'k-')
    
    # Connect main duct to fan
    plt.plot([0.5, 0.5], [y_pos/2, y_pos + 1], 'k-')
    
    plt.axis('off')
    
    plt.tight_layout()
    plt.savefig('ventilation_system_diagram.png', dpi=300)
    plt.close()

def main():
    # Analyze for worst climate case
    climate_condition = 'extreme_hot_humid'
    results, fan_requirements = analyze_ventilation_system(climate_condition)
    
    # Generate suggestions for improvements
    suggestions = suggest_improvements(results)
    
    # Recommend fans
    fan_recommendations = recommend_fans(fan_requirements)
    
    # Generate report
    generate_report(results, fan_requirements, suggestions, fan_recommendations, climate_condition)
    
    # Plot system diagram
    plot_system_diagram(results, fan_recommendations)
    
    # Save results to CSV
    results_df = pd.DataFrame()
    for bathroom_name, data in results.items():
        row = {
            'Bathroom': bathroom_name,
            'System': data['system']
        }
        
        if bathroom_name != 'Ducto_Principal_Piso2':
            row.update({
                'Volume (m³)': data['volume_m3'],
                'Air Changes per Hour': data['air_changes_per_hour'],
                'Grille Velocity (m/s)': data['grille_velocity_m_s']
            })
        
        row.update({
            'Flow Rate (m³/h)': data['flow_rate_m3h'],
            'Flow Rate (CFM)': data['flow_rate_cfm'],
            'Total Pressure Loss (Pa)': data['total_pressure_loss_pa']
        })
        
        results_df = pd.concat([results_df, pd.DataFrame([row])], ignore_index=True)
    
    results_df.to_csv('ventilation_results.csv', index=False)
    
    # Save fan recommendations to CSV
    fan_df = pd.DataFrame()
    for system, recommendation in fan_recommendations.items():
        if isinstance(recommendation, dict):
            row = {
                'System': system,
                'Model': recommendation['model'],
                'Max Flow (CFM)': recommendation['max_flow_cfm'],
                'Max Pressure (inWG)': recommendation['max_pressure_inwg'],
                'Power (W)': recommendation['power_watts'],
                'Flow Margin (%)': recommendation['flow_margin'],
                'Pressure Margin (%)': recommendation['pressure_margin']
            }
        else:
            row = {
                'System': system,
                'Model': 'Not found',
                'Max Flow (CFM)': None,
                'Max Pressure (inWG)': None,
                'Power (W)': None,
                'Flow Margin (%)': None,
                'Pressure Margin (%)': None
            }
        
        fan_df = pd.concat([fan_df, pd.DataFrame([row])], ignore_index=True)
    
    fan_df.to_csv('fan_recommendations.csv', index=False)
    
    print("\nResults saved to 'ventilation_results.csv'")
    print("Fan recommendations saved to 'fan_recommendations.csv'")
    print("System diagram saved to 'ventilation_system_diagram.png'")

if __name__ == "__main__":
    main()
