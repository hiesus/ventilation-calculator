import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Constants
AIR_DENSITY_STD = 1.2  # kg/m³ at sea level, 20°C
GRAVITY = 9.81  # m/s²
AIR_VISCOSITY = 1.81e-5  # Pa·s at 25°C
ROUGHNESS = 0.09e-3  # m (galvanized steel)
CHANGES_PER_HOUR = 12  # Air changes per hour

# Altitude correction
ALTITUDE = 850  # meters above sea level
TEMP_AVG = 25  # °C
RH_AVG = 60  # %

# Convert imperial to metric
def inch_to_m(inch):
    return inch * 0.0254

def cfm_to_m3s(cfm):
    return cfm * 0.000471947

def m3s_to_cfm(m3s):
    return m3s / 0.000471947

def pa_to_inwg(pa):
    return pa * 0.00402

def inwg_to_pa(inwg):
    return inwg * 248.84

# Air density correction for altitude and temperature
def air_density(altitude, temperature, rh):
    # Barometric pressure at altitude (Pa)
    p0 = 101325  # sea level standard pressure (Pa)
    T0 = 288.15  # sea level standard temperature (K)
    g = 9.80665  # gravitational acceleration (m/s²)
    L = 0.0065  # temperature lapse rate (K/m)
    R = 8.31447  # universal gas constant (J/(mol·K))
    M = 0.0289644  # molar mass of dry air (kg/mol)
    
    # Barometric pressure at altitude
    p = p0 * (1 - (L * altitude / T0)) ** (g * M / (R * L))
    
    # Saturation vapor pressure (Pa)
    T = temperature + 273.15  # Convert to Kelvin
    es = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
    
    # Actual vapor pressure (Pa)
    e = rh / 100 * es
    
    # Density of moist air (kg/m³)
    Rd = 287.058  # specific gas constant for dry air (J/(kg·K))
    Rv = 461.495  # specific gas constant for water vapor (J/(kg·K))
    
    density = (p - e) / (Rd * T) + e / (Rv * T)
    
    return density

# Calculate required airflow based on room volume and air changes per hour
def calculate_airflow(volume, changes_per_hour):
    # Volume in m³, result in m³/s
    return volume * changes_per_hour / 3600

# Calculate hydraulic diameter for rectangular ducts
def hydraulic_diameter(width, height):
    # width and height in meters
    return 2 * (width * height) / (width + height)

# Calculate friction factor using Colebrook equation
def friction_factor(reynolds, roughness, diameter):
    def colebrook(f):
        return 1/np.sqrt(f) + 2*np.log10(roughness/(3.7*diameter) + 2.51/(reynolds*np.sqrt(f)))
    
    # Initial guess using Swamee-Jain equation
    f_initial = 0.25 / (np.log10(roughness/(3.7*diameter) + 5.74/reynolds**0.9))**2
    
    try:
        f = fsolve(colebrook, f_initial)[0]
        return f
    except:
        # Fallback to Swamee-Jain if fsolve fails
        return f_initial

# Calculate pressure drop in straight ducts
def pressure_drop_straight(flow_rate, length, diameter, density, viscosity, roughness):
    # flow_rate in m³/s, length in m, diameter in m
    velocity = flow_rate / (np.pi * (diameter/2)**2)
    reynolds = density * velocity * diameter / viscosity
    f = friction_factor(reynolds, roughness, diameter)
    
    return f * density * length * velocity**2 / (2 * diameter)

# Calculate pressure drop in rectangular ducts
def pressure_drop_rect(flow_rate, length, width, height, density, viscosity, roughness):
    # flow_rate in m³/s, length in m, width and height in m
    area = width * height
    velocity = flow_rate / area
    dh = hydraulic_diameter(width, height)
    reynolds = density * velocity * dh / viscosity
    f = friction_factor(reynolds, roughness, dh)
    
    return f * density * length * velocity**2 / (2 * dh)

# Calculate pressure drop in fittings (K-method)
def pressure_drop_fitting(flow_rate, k_factor, area, density):
    # flow_rate in m³/s, area in m²
    velocity = flow_rate / area
    return k_factor * 0.5 * density * velocity**2

# K factors for common fittings
k_factors = {
    'elbow_90': 0.3,  # 90° smooth bend
    'elbow_90_sharp': 1.2,  # 90° sharp bend
    'tee_straight': 0.5,  # Flow through straight section of tee
    'tee_branch': 1.0,  # Flow through branch of tee
    'entry': 0.5,  # Duct entry
    'exit': 1.0,  # Duct exit
    'sudden_expansion': lambda A1, A2: (1 - A1/A2)**2,  # A1=smaller area, A2=larger area
    'sudden_contraction': lambda A1, A2: 0.4 * (1 - A2/A1)  # A1=larger area, A2=smaller area
}

# Define bathroom systems
def define_bathroom_systems():
    systems = {
        'system1': {  # PB Baño C, PB Baño D, Cafetín Baño C, Cafetín Baño D
            'rooms': [
                {
                    'name': 'PB Baño C',
                    'volume': 25,  # m³
                    'grilles': 2,
                    'grille_size': (inch_to_m(4), inch_to_m(11))
                },
                {
                    'name': 'PB Baño D',
                    'volume': 18,  # m³
                    'grilles': 3,
                    'grille_size': (inch_to_m(6), inch_to_m(11))
                },
                {
                    'name': 'Cafetín Baño C',
                    'volume': 6.5,  # m³
                    'grilles': 1,
                    'grille_size': (inch_to_m(7.5), inch_to_m(10))
                },
                {
                    'name': 'Cafetín Baño D',
                    'volume': 9,  # m³
                    'grilles': 2,
                    'grille_size': (inch_to_m(7.5), inch_to_m(11))
                }
            ],
            'ducts': [
                {
                    'name': 'PB Baño C main duct',
                    'type': 'rectangular',
                    'dimensions': (inch_to_m(6), inch_to_m(11)),
                    'length': 4  # m
                },
                {
                    'name': 'PB expansion',
                    'type': 'expansion',
                    'from_dimensions': (inch_to_m(6), inch_to_m(11)),
                    'to_dimensions': (inch_to_m(11), inch_to_m(11))
                },
                {
                    'name': 'PB larger duct',
                    'type': 'rectangular',
                    'dimensions': (inch_to_m(11), inch_to_m(11)),
                    'length': 2  # m
                },
                {
                    'name': 'PB 90° elbow',
                    'type': 'elbow',
                    'dimensions': (inch_to_m(11), inch_to_m(11)),
                    'angle': 90
                },
                {
                    'name': 'PB rect to round',
                    'type': 'transition',
                    'from_dimensions': (inch_to_m(11), inch_to_m(11)),
                    'to_dimensions': inch_to_m(12.5)  # circular diameter
                },
                {
                    'name': 'PB vertical duct',
                    'type': 'circular',
                    'dimensions': inch_to_m(12.5),
                    'length': 6  # m
                },
                {
                    'name': 'PB final elbow',
                    'type': 'elbow',
                    'dimensions': inch_to_m(12.5),
                    'angle': 90
                },
                {
                    'name': 'PB Baño D duct',
                    'type': 'rectangular',
                    'dimensions': (inch_to_m(6), inch_to_m(11)),
                    'length': 4  # m
                },
                {
                    'name': 'PB Baño D elbows',
                    'type': 'elbow',
                    'dimensions': (inch_to_m(6), inch_to_m(11)),
                    'angle': 90,
                    'count': 3
                },
                {
                    'name': 'Cafetín Baño C duct',
                    'type': 'circular',
                    'dimensions': inch_to_m(10),
                    'length': 1.44  # m
                },
                {
                    'name': 'Cafetín Baño C elbow',
                    'type': 'elbow',
                    'dimensions': inch_to_m(10),
                    'angle': 90
                },
                {
                    'name': 'Cafetín transition',
                    'type': 'transition',
                    'from_dimensions': inch_to_m(10),
                    'to_dimensions': inch_to_m(12.5)
                },
                {
                    'name': 'Cafetín Baño D duct',
                    'type': 'circular',
                    'dimensions': inch_to_m(10),
                    'length': 2  # m
                }
            ]
        },
        'system2': {  # Piso 1 Baño C, Piso 1 Baño D, Piso 2 Baño C, Piso 2 Baño D
            'rooms': [
                {
                    'name': 'Piso 1 Baño C',
                    'volume': 25,  # Estimated volume
                    'grilles': 3,
                    'grille_size': (inch_to_m(6), inch_to_m(14))
                },
                {
                    'name': 'Piso 1 Baño D',
                    'volume': 25,  # Estimated volume
                    'grilles': 3,
                    'grille_size': (inch_to_m(6), inch_to_m(14))
                },
                {
                    'name': 'Piso 2 Baño C',
                    'volume': 20,  # Estimated volume
                    'grilles': 3,
                    'grille_size': (inch_to_m(4), inch_to_m(11))
                },
                {
                    'name': 'Piso 2 Baño D',
                    'volume': 20,  # Estimated volume
                    'grilles': 3,
                    'grille_size': (inch_to_m(4), inch_to_m(11))
                }
            ],
            'ducts': [
                {
                    'name': 'Piso 1 Baño C duct',
                    'type': 'rectangular',
                    'dimensions': (inch_to_m(11), inch_to_m(6)),
                    'length': 3.2  # m
                },
                {
                    'name': 'Piso 1 Baño C elbow',
                    'type': 'elbow',
                    'dimensions': (inch_to_m(11), inch_to_m(6)),
                    'angle': 90
                },
                {
                    'name': 'Piso 1 Baño C vertical',
                    'type': 'rectangular',
                    'dimensions': (inch_to_m(11), inch_to_m(6)),
                    'length': 3.5  # m
                },
                {
                    'name': 'Piso 1 Baño C tee',
                    'type': 'tee',
                    'dimensions': (inch_to_m(11), inch_to_m(6)),
                    'main_dimensions': (inch_to_m(11), inch_to_m(11))
                },
                {
                    'name': 'Piso 1 Baño D duct',
                    'type': 'rectangular',
                    'dimensions': (inch_to_m(11), inch_to_m(6)),
                    'length': 3.2  # m
                },
                {
                    'name': 'Piso 1 Baño D elbow',
                    'type': 'elbow',
                    'dimensions': (inch_to_m(11), inch_to_m(6)),
                    'angle': 90
                },
                {
                    'name': 'Piso 1 Baño D vertical',
                    'type': 'rectangular',
                    'dimensions': (inch_to_m(11), inch_to_m(6)),
                    'length': 3.5  # m
                },
                {
                    'name': 'Piso 1 Baño D tee',
                    'type': 'tee',
                    'dimensions': (inch_to_m(11), inch_to_m(6)),
                    'main_dimensions': (inch_to_m(11), inch_to_m(11))
                },
                {
                    'name': 'Piso 2 Baño C duct',
                    'type': 'rectangular',
                    'dimensions': (inch_to_m(11), inch_to_m(6)),
                    'length': 3.2  # m
                },
                {
                    'name': 'Piso 2 Baño C tee',
                    'type': 'tee',
                    'dimensions': (inch_to_m(11), inch_to_m(6)),
                    'main_dimensions': (inch_to_m(11), inch_to_m(11))
                },
                {
                    'name': 'Piso 2 Baño D duct',
                    'type': 'rectangular',
                    'dimensions': (inch_to_m(11), inch_to_m(6)),
                    'length': 3.2  # m
                },
                {
                    'name': 'Piso 2 Baño D tee',
                    'type': 'tee',
                    'dimensions': (inch_to_m(11), inch_to_m(6)),
                    'main_dimensions': (inch_to_m(11), inch_to_m(11))
                },
                {
                    'name': 'Ducto Principal Piso 2',
                    'type': 'rectangular',
                    'dimensions': (inch_to_m(11), inch_to_m(11)),
                    'length': 12  # m
                },
                {
                    'name': 'Transition to fan',
                    'type': 'transition',
                    'from_dimensions': (inch_to_m(11), inch_to_m(11)),
                    'to_dimensions': inch_to_m(12)  # circular diameter
                }
            ]
        }
    }
    return systems

# Calculate total airflow requirements for each system
def calculate_system_airflows(systems):
    results = {}
    
    # Calculate corrected air density
    density = air_density(ALTITUDE, TEMP_AVG, RH_AVG)
    
    for system_name, system in systems.items():
        total_flow = 0
        room_flows = {}
        
        for room in system['rooms']:
            flow = calculate_airflow(room['volume'], CHANGES_PER_HOUR)
            room_flows[room['name']] = flow
            total_flow += flow
        
        results[system_name] = {
            'total_flow': total_flow,
            'total_flow_cfm': m3s_to_cfm(total_flow),
            'room_flows': room_flows,
            'density': density
        }
    
    return results

# Calculate pressure drops in the system
def calculate_pressure_drops(systems, airflows):
    results = {}
    
    for system_name, system in systems.items():
        total_pressure_drop = 0
        component_drops = {}
        density = airflows[system_name]['density']
        total_flow = airflows[system_name]['total_flow']
        
        # Calculate pressure drops for each duct component
        for duct in system['ducts']:
            if duct['type'] == 'rectangular':
                width, height = duct['dimensions']
                area = width * height
                dh = hydraulic_diameter(width, height)
                
                dp = pressure_drop_rect(
                    total_flow, 
                    duct['length'], 
                    width, 
                    height, 
                    density, 
                    AIR_VISCOSITY, 
                    ROUGHNESS
                )
                
                component_drops[duct['name']] = dp
                total_pressure_drop += dp
                
            elif duct['type'] == 'circular':
                diameter = duct['dimensions']
                area = np.pi * (diameter/2)**2
                
                dp = pressure_drop_straight(
                    total_flow, 
                    duct['length'], 
                    diameter, 
                    density, 
                    AIR_VISCOSITY, 
                    ROUGHNESS
                )
                
                component_drops[duct['name']] = dp
                total_pressure_drop += dp
                
            elif duct['type'] == 'elbow':
                if isinstance(duct['dimensions'], tuple):
                    # Rectangular elbow
                    width, height = duct['dimensions']
                    area = width * height
                else:
                    # Circular elbow
                    diameter = duct['dimensions']
                    area = np.pi * (diameter/2)**2
                
                k = k_factors['elbow_90']
                count = duct.get('count', 1)
                
                dp = pressure_drop_fitting(total_flow, k * count, area, density)
                
                component_drops[duct['name']] = dp
                total_pressure_drop += dp
                
            elif duct['type'] == 'tee':
                width, height = duct['dimensions']
                area = width * height
                
                k = k_factors['tee_branch']
                
                dp = pressure_drop_fitting(total_flow, k, area, density)
                
                component_drops[duct['name']] = dp
                total_pressure_drop += dp
                
            elif duct['type'] == 'expansion':
                width1, height1 = duct['from_dimensions']
                width2, height2 = duct['to_dimensions']
                
                area1 = width1 * height1
                area2 = width2 * height2
                
                k = k_factors['sudden_expansion'](area1, area2)
                
                dp = pressure_drop_fitting(total_flow, k, area1, density)
                
                component_drops[duct['name']] = dp
                total_pressure_drop += dp
                
            elif duct['type'] == 'transition':
                # Handle different types of transitions
                if isinstance(duct['from_dimensions'], tuple) and not isinstance(duct['to_dimensions'], tuple):
                    # Rectangular to circular
                    width, height = duct['from_dimensions']
                    area1 = width * height
                    
                    diameter = duct['to_dimensions']
                    area2 = np.pi * (diameter/2)**2
                elif not isinstance(duct['from_dimensions'], tuple) and isinstance(duct['to_dimensions'], tuple):
                    # Circular to rectangular
                    diameter = duct['from_dimensions']
                    area1 = np.pi * (diameter/2)**2
                    
                    width, height = duct['to_dimensions']
                    area2 = width * height
                elif isinstance(duct['from_dimensions'], tuple) and isinstance(duct['to_dimensions'], tuple):
                    # Rectangular to rectangular
                    width1, height1 = duct['from_dimensions']
                    width2, height2 = duct['to_dimensions']
                    
                    area1 = width1 * height1
                    area2 = width2 * height2
                else:
                    # Circular to circular
                    diameter1 = duct['from_dimensions']
                    diameter2 = duct['to_dimensions']
                    
                    area1 = np.pi * (diameter1/2)**2
                    area2 = np.pi * (diameter2/2)**2
                
                if area1 < area2:
                    k = k_factors['sudden_expansion'](area1, area2)
                else:
                    k = k_factors['sudden_contraction'](area1, area2)
                
                dp = pressure_drop_fitting(total_flow, k, min(area1, area2), density)
                
                component_drops[duct['name']] = dp
                total_pressure_drop += dp
        
        # Add entry and exit losses
        # Entry loss (grilles)
        total_grille_area = 0
        for room in system['rooms']:
            width, height = room['grille_size']
            grille_area = width * height * room['grilles']
            total_grille_area += grille_area
        
        entry_dp = pressure_drop_fitting(total_flow, k_factors['entry'], total_grille_area, density)
        component_drops['Grille entry losses'] = entry_dp
        total_pressure_drop += entry_dp
        
        # Exit loss (fan outlet)
        if system_name == 'system1':
            exit_area = np.pi * (inch_to_m(12.5)/2)**2
        else:  # system2
            exit_area = np.pi * (inch_to_m(12)/2)**2
            
        exit_dp = pressure_drop_fitting(total_flow, k_factors['exit'], exit_area, density)
        component_drops['Fan exit losses'] = exit_dp
        total_pressure_drop += exit_dp
        
        # Add safety factor (20%)
        total_pressure_drop *= 1.2
        
        results[system_name] = {
            'total_pressure': total_pressure_drop,
            'total_pressure_inwg': pa_to_inwg(total_pressure_drop),
            'component_drops': component_drops
        }
    
    return results

# Select appropriate fans based on requirements
def select_fans(airflows, pressure_drops):
    fan_options = {
        'system1': [],
        'system2': []
    }
    
    # Define some commercial centrifugal fans with their performance data
    # Format: [Model, Max Flow (CFM), Max Static Pressure (inWG), Power (W), Diameter (inches)]
    commercial_fans = [
        ['CF-100', 200, 0.5, 35, 4],
        ['CF-150', 350, 0.8, 60, 6],
        ['CF-200', 500, 1.0, 90, 8],
        ['CF-250', 700, 1.2, 120, 10],
        ['CF-300', 900, 1.5, 180, 12],
        ['CF-350', 1200, 1.8, 250, 14],
        ['CF-400', 1500, 2.0, 350, 16],
        ['CF-500', 2000, 2.5, 500, 20],
        ['CF-600', 2500, 3.0, 750, 24],
        ['CF-700', 3000, 3.5, 1000, 28]
    ]
    
    for system_name in airflows:
        required_flow_cfm = airflows[system_name]['total_flow_cfm']
        required_pressure_inwg = pressure_drops[system_name]['total_pressure_inwg']
        
        # Filter fans that meet requirements with 20% margin
        suitable_fans = []
        for fan in commercial_fans:
            model, max_flow, max_pressure, power, diameter = fan
            
            # Check if fan meets requirements with margin
            if max_flow >= required_flow_cfm * 1.2 and max_pressure >= required_pressure_inwg * 1.2:
                suitable_fans.append({
                    'model': model,
                    'max_flow_cfm': max_flow,
                    'max_pressure_inwg': max_pressure,
                    'power_w': power,
                    'diameter_inches': diameter,
                    'flow_margin': (max_flow / required_flow_cfm - 1) * 100,
                    'pressure_margin': (max_pressure / required_pressure_inwg - 1) * 100
                })
        
        # Sort by closest match (minimize oversizing)
        if suitable_fans:
            suitable_fans.sort(key=lambda x: x['flow_margin'] + x['pressure_margin'])
            fan_options[system_name] = suitable_fans
    
    return fan_options

# Check if duct velocities are within standards
def check_duct_velocities(systems, airflows):
    results = {}
    
    # Recommended velocity ranges (m/s)
    velocity_ranges = {
        'main_ducts': (5, 10),  # 5-10 m/s for main ducts
        'branch_ducts': (3, 8),  # 3-8 m/s for branch ducts
        'grilles': (2, 4)  # 2-4 m/s for grilles
    }
    
    for system_name, system in systems.items():
        total_flow = airflows[system_name]['total_flow']
        velocity_issues = []
        
        # Check duct velocities
        for duct in system['ducts']:
            if duct['type'] in ['rectangular', 'circular']:
                if duct['type'] == 'rectangular':
                    width, height = duct['dimensions']
                    area = width * height
                else:  # circular
                    diameter = duct['dimensions']
                    area = np.pi * (diameter/2)**2
                
                velocity = total_flow / area
                
                # Determine if main or branch duct
                is_main = False
                if 'main' in duct['name'].lower() or 'principal' in duct['name'].lower():
                    is_main = True
                elif any(dim >= inch_to_m(10) for dim in duct['dimensions']) if duct['type'] == 'rectangular' else duct['dimensions'] >= inch_to_m(10):
                    is_main = True
                
                range_key = 'main_ducts' if is_main else 'branch_ducts'
                min_vel, max_vel = velocity_ranges[range_key]
                
                if velocity < min_vel or velocity > max_vel:
                    velocity_issues.append({
                        'component': duct['name'],
                        'velocity': velocity,
                        'recommended_range': (min_vel, max_vel),
                        'status': 'Too low' if velocity < min_vel else 'Too high'
                    })
        
        # Check grille velocities
        for room in system['rooms']:
            width, height = room['grille_size']
            grille_area = width * height * room['grilles']
            room_flow = airflows[system_name]['room_flows'][room['name']]
            
            velocity = room_flow / grille_area
            min_vel, max_vel = velocity_ranges['grilles']
            
            if velocity < min_vel or velocity > max_vel:
                velocity_issues.append({
                    'component': f"{room['name']} grilles",
                    'velocity': velocity,
                    'recommended_range': (min_vel, max_vel),
                    'status': 'Too low' if velocity < min_vel else 'Too high'
                })
        
        results[system_name] = velocity_issues
    
    return results

# Generate a comprehensive report
def generate_report(systems, airflows, pressure_drops, fan_options, velocity_issues):
    report = {
        'systems': {},
        'recommendations': {}
    }
    
    for system_name in systems:
        system_report = {
            'airflow': {
                'total_m3s': airflows[system_name]['total_flow'],
                'total_cfm': airflows[system_name]['total_flow_cfm'],
                'rooms': {}
            },
            'pressure_drop': {
                'total_pa': pressure_drops[system_name]['total_pressure'],
                'total_inwg': pressure_drops[system_name]['total_pressure_inwg'],
                'components': {}
            },
            'velocity_issues': velocity_issues[system_name]
        }
        
        # Add room airflows
        for room_name, flow in airflows[system_name]['room_flows'].items():
            system_report['airflow']['rooms'][room_name] = {
                'm3s': flow,
                'cfm': m3s_to_cfm(flow)
            }
        
        # Add component pressure drops
        for component, dp in pressure_drops[system_name]['component_drops'].items():
            system_report['pressure_drop']['components'][component] = {
                'pa': dp,
                'inwg': pa_to_inwg(dp)
            }
        
        report['systems'][system_name] = system_report
        
        # Add fan recommendations
        if fan_options[system_name]:
            best_fan = fan_options[system_name][0]
            report['recommendations'][system_name] = {
                'fan_model': best_fan['model'],
                'max_flow_cfm': best_fan['max_flow_cfm'],
                'max_pressure_inwg': best_fan['max_pressure_inwg'],
                'power_w': best_fan['power_w'],
                'diameter_inches': best_fan['diameter_inches'],
                'flow_margin_percent': best_fan['flow_margin'],
                'pressure_margin_percent': best_fan['pressure_margin']
            }
        else:
            report['recommendations'][system_name] = {
                'message': 'No suitable fan found. Consider custom solutions or redesigning the duct system.'
            }
    
    return report

# Main function to run the calculations
def calculate_ventilation_system():
    # Define the systems
    systems = define_bathroom_systems()
    
    # Calculate airflow requirements
    airflows = calculate_system_airflows(systems)
    
    # Calculate pressure drops
    pressure_drops = calculate_pressure_drops(systems, airflows)
    
    # Check duct velocities
    velocity_issues = check_duct_velocities(systems, airflows)
    
    # Select appropriate fans
    fan_options = select_fans(airflows, pressure_drops)
    
    # Generate report
    report = generate_report(systems, airflows, pressure_drops, fan_options, velocity_issues)
    
    return report

# Function to print the report in a readable format
def print_report(report):
    print("=== VENTILATION SYSTEM REPORT ===\n")
    
    for system_name, system_data in report['systems'].items():
        print(f"SYSTEM: {system_name}")
        print("-" * 40)
        
        # Airflow information
        print("AIRFLOW REQUIREMENTS:")
        print(f"Total: {system_data['airflow']['total_m3s']:.4f} m³/s ({system_data['airflow']['total_cfm']:.1f} CFM)")
        print("Room breakdown:")
        for room, flow in system_data['airflow']['rooms'].items():
            print(f"  - {room}: {flow['m3s']:.4f} m³/s ({flow['cfm']:.1f} CFM)")
        
        # Pressure drop information
        print("\nPRESSURE DROP ANALYSIS:")
        print(f"Total: {system_data['pressure_drop']['total_pa']:.2f} Pa ({system_data['pressure_drop']['total_inwg']:.4f} inWG)")
        print("Component breakdown:")
        for component, dp in system_data['pressure_drop']['components'].items():
            print(f"  - {component}: {dp['pa']:.2f} Pa ({dp['inwg']:.4f} inWG)")
        
        # Velocity issues
        print("\nVELOCITY ISSUES:")
        if system_data['velocity_issues']:
            for issue in system_data['velocity_issues']:
                print(f"  - {issue['component']}: {issue['velocity']:.2f} m/s")
                print(f"    Recommended: {issue['recommended_range'][0]}-{issue['recommended_range'][1]} m/s")
                print(f"    Status: {issue['status']}")
        else:
            print("  No velocity issues detected.")
        
        # Fan recommendation
        print("\nFAN RECOMMENDATION:")
        if system_name in report['recommendations']:
            rec = report['recommendations'][system_name]
            if 'fan_model' in rec:
                print(f"  Model: {rec['fan_model']}")
                print(f"  Max Flow: {rec['max_flow_cfm']} CFM")
                print(f"  Max Static Pressure: {rec['max_pressure_inwg']} inWG")
                print(f"  Power: {rec['power_w']} W")
                print(f"  Diameter: {rec['diameter_inches']} inches")
                print(f"  Flow Margin: {rec['flow_margin_percent']:.1f}%")
                print(f"  Pressure Margin: {rec['pressure_margin_percent']:.1f}%")
            else:
                print(f"  {rec['message']}")
        
        print("\n" + "=" * 40 + "\n")

# Function to visualize the system
def visualize_system(report):
    # Create a figure with subplots for each system
    fig, axs = plt.subplots(len(report['systems']), 2, figsize=(15, 8 * len(report['systems'])))
    
    for i, (system_name, system_data) in enumerate(report['systems'].items()):
        # Airflow distribution pie chart
        ax1 = axs[i, 0] if len(report['systems']) > 1 else axs[0]
        room_names = list(system_data['airflow']['rooms'].keys())
        room_flows = [flow['cfm'] for flow in system_data['airflow']['rooms'].values()]
        
        ax1.pie(room_flows, labels=room_names, autopct='%1.1f%%', startangle=90)
        ax1.set_title(f'{system_name} - Airflow Distribution (CFM)')
        
        # Pressure drop bar chart
        ax2 = axs[i, 1] if len(report['systems']) > 1 else axs[1]
        component_names = list(system_data['pressure_drop']['components'].keys())
        component_pressures = [dp['pa'] for dp in system_data['pressure_drop']['components'].values()]
        
        # Sort by pressure drop
        sorted_indices = np.argsort(component_pressures)[::-1]
        sorted_names = [component_names[j] for j in sorted_indices]
        sorted_pressures = [component_pressures[j] for j in sorted_indices]
        
        # Limit to top 10 components for readability
        if len(sorted_names) > 10:
            sorted_names = sorted_names[:10]
            sorted_pressures = sorted_pressures[:10]
            sorted_names.append('Others')
            sorted_pressures.append(sum(component_pressures) - sum(sorted_pressures))
        
        bars = ax2.barh(sorted_names, sorted_pressures)
        ax2.set_title(f'{system_name} - Pressure Drop by Component (Pa)')
        ax2.set_xlabel('Pressure Drop (Pa)')
        
        # Add values to bars
        for bar in bars:
            width = bar.get_width()
            ax2.text(width + 1, bar.get_y() + bar.get_height()/2, f'{width:.1f}', 
                    ha='left', va='center')
    
    plt.tight_layout()
    plt.savefig('ventilation_analysis.png')
    plt.close()

# Function to export report to CSV
def export_to_csv(report):
    # Export airflow data
    airflow_data = []
    for system_name, system_data in report['systems'].items():
        for room, flow in system_data['airflow']['rooms'].items():
            airflow_data.append({
                'System': system_name,
                'Room': room,
                'Flow (m³/s)': flow['m3s'],
                'Flow (CFM)': flow['cfm']
            })
    
    airflow_df = pd.DataFrame(airflow_data)
    airflow_df.to_csv('airflow_requirements.csv', index=False)
    
    # Export pressure drop data
    pressure_data = []
    for system_name, system_data in report['systems'].items():
        for component, dp in system_data['pressure_drop']['components'].items():
            pressure_data.append({
                'System': system_name,
                'Component': component,
                'Pressure Drop (Pa)': dp['pa'],
                'Pressure Drop (inWG)': dp['inwg']
            })
    
    pressure_df = pd.DataFrame(pressure_data)
    pressure_df.to_csv('pressure_drops.csv', index=False)
    
    # Export fan recommendations
    fan_data = []
    for system_name, rec in report['recommendations'].items():
        if 'fan_model' in rec:
            fan_data.append({
                'System': system_name,
                'Model': rec['fan_model'],
                'Max Flow (CFM)': rec['max_flow_cfm'],
                'Max Pressure (inWG)': rec['max_pressure_inwg'],
                'Power (W)': rec['power_w'],
                'Diameter (inches)': rec['diameter_inches'],
                'Flow Margin (%)': rec['flow_margin_percent'],
                'Pressure Margin (%)': rec['pressure_margin_percent']
            })
    
    if fan_data:
        fan_df = pd.DataFrame(fan_data)
        fan_df.to_csv('fan_recommendations.csv', index=False)
    
    # Export velocity issues
    velocity_data = []
    for system_name, system_data in report['systems'].items():
        for issue in system_data['velocity_issues']:
            velocity_data.append({
                'System': system_name,
                'Component': issue['component'],
                'Velocity (m/s)': issue['velocity'],
                'Min Recommended (m/s)': issue['recommended_range'][0],
                'Max Recommended (m/s)': issue['recommended_range'][1],
                'Status': issue['status']
            })
    
    if velocity_data:
        velocity_df = pd.DataFrame(velocity_data)
        velocity_df.to_csv('velocity_issues.csv', index=False)

if __name__ == "__main__":
    report = calculate_ventilation_system()
    print_report(report)
    visualize_system(report)
    export_to_csv(report)
