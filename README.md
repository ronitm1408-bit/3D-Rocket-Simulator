6-DOF Rocket Trajectory Simulation and 3D Visualisation

Project Overview
This project simulates the flight of a small solid-fuel rocket using a
12-state, six-degree-of-freedom (6-DOF) model implemented in Python. It
calculates the rocket's translational and rotational state over time,
exports the simulation data to a CSV file, and creates an interactive 3D
Plotly animation of the trajectory.

The simulation includes: - A time-dependent rocket mass model. - A
constant-thrust motor model. - An International Standard Atmosphere
(ISA) density model. - Altitude-dependent gravitational acceleration. -
Quadratic aerodynamic drag. - Numerical integration using SciPy's
solve_ivp solver. - CSV data export using pandas. - Interactive 3D
visualisation with a rocket model, trajectory markers, play/pause
controls, and a time slider.

Important model-verification note: The project report identifies a
frame-consistency issue in the translational equations. The velocity
states are treated as body-frame components for position updates,
while their derivatives are assigned inertial-frame acceleration.
Consequently, the reported horizontal trajectory, apogee timing,
range, and flight time are not yet validated physical predictions. See
Limitations and Verification below before interpreting the
results.

Objectives
1. Represent the rocket's mass, geometry, motor, aerodynamic
   properties, and moments of inertia using a Python class.
2. Implement atmospheric density, gravity, mass variation, and thrust
   as separate functions.
3. Formulate and numerically integrate the 12-state 6-DOF equations of
   motion.
4. Calculate flight quantities such as apogee, time to apogee, speed,
   acceleration, downrange distance, and flight time.
5. Export the simulation history to Rocket_6DOF_Trajectory.csv.
6. Generate an interactive 3D trajectory animation in
   Rocket_6DOF_3D_Trajectory.html.
7. Compare the simulation with an independent inertial-frame check
   model and document the model's limitations.
   
Rocket and Simulation Parameters
  Parameter                                                   Value
  Initial mass                                               5.0 kg
Dry mass                                                   3.0 kg
Propellant mass                                            2.0 kg
Rocket length                                              1.50 m
Rocket diameter                                            0.08 m
Frontal area                                           0.00503 m²
Drag coefficient, CD                                       0.40
Motor thrust                                                150 N
Burn time                                                   3.0 s
Mass flow rate                                         0.667 kg/s
Launch elevation                                              80°
Initial speed                                             0.1 m/s
Moments of inertia         Ixx = 0.02, Iyy = Izz = 0.80 kg·m²

The integration is configured for 0--120 seconds, with 6,000 requested
output samples (nominal interval 0.02 seconds), relative tolerance
1e-7, absolute tolerance 1e-9, and maximum solver step 0.02
seconds.
Requirements
- Python 3.9 or later is recommended.
- The following Python packages are required:
  - NumPy
  - SciPy
  - pandas
  - Plotly
The Python standard library modules os and webbrowser are also used;
they do not need to be installed separately.
Project Structure
A suggested packaged project layout is:
rocket-6dof-project/
├── README.md
├── pyproject.toml
├── src/
│   └── ...                 # Rocket simulation package source
├── examples/
│   └── ...                 # Example run scripts
└── tests/
    └── ...                 # Verification or unit tests
If you are using the standalone source-code listing supplied with the
project report, place the Python file (for example,
OpenElective_Project_original.py) in the project folder. The script
writes its output files to the current working directory.
Installation
1. Create or activate a Python environment (recommended)
python -m venv .venv
Activate it:

Windows PowerShell
.venv\Scripts\Activate.ps1
macOS / Linux
source .venv/bin/activate
2. Install the dependencies
For the standalone script:
python -m pip install numpy scipy pandas plotly
If the packaged project includes a valid pyproject.toml, install it
from the project root instead:
python -m pip install .
Running the Program
Standalone Python script
1. Open a terminal or command prompt.
2. Change directory to the folder containing the Python source file.
3. Run:
python OpenElective_Project_original.py
Replace the filename with the actual name of your Python script if it
differs.

Installed package
After installing the package with python -m pip install ., run the
example or entry-point command provided by that package. The exact
command depends on the package's pyproject.toml configuration; the
report does not specify a console-script command.

Output Files
After a successful run, the program generates:
  File                                Description
  Rocket_6DOF_Trajectory.csv        Tabular simulation history
containing time, position, speed,
acceleration, mass, and thrust.
Rocket_6DOF_3D_Trajectory.html    Interactive 3D Plotly visualisation
of the flight.

The CSV contains 6,000 requested samples and 10 columns:
- Time_s
- X_m
- Y_m
- Altitude_m
- Velocity_mps
- Velocity_kmph
- Acceleration_mps2
- Acceleration_g
- Mass_kg
- Thrust_N
Open the HTML file in a modern web browser. Use the PLAY and
PAUSE controls to animate the rocket, the time slider to move
through the flight, and mouse controls to rotate or zoom the 3D scene.

The solver in the submitted implementation continues beyond ground
impact. The report therefore evaluates the flight only through the first
sample at or below ground level after apogee. Later CSV rows represent
below-ground mathematical continuation, not actual flight.
Model Methodology
The program follows this sequence:
1. Create the rocket object and load its physical and motor parameters.
2. Calculate atmospheric temperature, pressure, and density using the
   ISA model.
3. Calculate gravity as a function of altitude.
4. Calculate rocket mass and motor thrust as functions of time.
5. Construct the body-to-inertial rotation matrix from Euler angles.
6. Evaluate the 12-state equations in the dynamics() function.
7. Integrate the equations using scipy.integrate.solve_ivp.
8. Extract position, velocity, and acceleration and identify apogee and
   impact.
9. Save the results to CSV.
10. Construct the rocket geometry and interactive 3D Plotly animation.
11. Write the HTML visualisation and attempt to open it in a browser.
Reported Results
The submitted program reports the following values:

  Quantity                                                 Program result
  Maximum altitude (apogee)                                      474.78 m
Time to apogee                                                  20.18 s
Maximum acceleration                                39.63 m/s² (4.04 g)
Downrange distance at first below-ground sample              1,583.51 m
Lateral distance                                                 0.00 m
Time to first below-ground sample                               31.25 s

The report notes that the program's printed maximum speed over the
entire 120-second integration includes below-ground samples. Within the
flight rows up to the first below-ground sample, the maximum speed is
approximately 186.0 m/s.

These are reported program outputs, not fully validated flight
predictions.

Limitations and Verification
The report's independent inertial-frame check model gives an apogee of
409.8 m, downrange distance at impact of 313.4 m, and flight time of
20.1 s. The differences reveal that the submitted translational
equations require correction before the trajectory results can be
trusted.

The report identifies the main issues and simplifications:
- Translational velocity components are used as body-frame quantities,
  but inertial-frame acceleration is assigned to their derivatives. A
  consistent formulation should either retain velocity in the inertial
  frame or include the required body-frame rotational and gravity
  terms.
- No aerodynamic or control moments are included. The body rates
  remain zero, so the rotational dynamics do not become active and the
  attitude does not respond to the flight path.
- The drag coefficient is constant; Mach-number and angle-of-attack
  effects, lift, and fin forces are not modelled.
- Thrust and mass flow are simplified as constant during the burn.
- Wind, Earth rotation, launch-rail dynamics, and parachute recovery
  are not included.
- The solver is not terminated at ground impact; the impact time is
  taken from the first sampled point at or below ground, so it has a
  small sampling overshoot.
- The graphical rocket is enlarged for visualisation and is not drawn
  to its real physical scale.
  
Suggested Future Improvements
1. Reformulate the translational equations consistently in one
   reference frame.
2. Add aerodynamic forces and moments, including centre-of-pressure and
   stability effects.
3. Use a measured motor thrust curve and a Mach-dependent drag
   coefficient.
4. Include wind, launch-rail dynamics, and parachute recovery.
5. Stop integration at impact using a solver event function.
6. Validate the updated model against flight-test data or an
   established simulation tool such as OpenRocket.
   
Main Python Functions
  Function / class                     Purpose
  Rocket                             Stores rocket, motor, geometry, and
inertia parameters.

  atmosphere(h)                      Returns ISA temperature, pressure,
and density.

  gravity(h)                         Returns altitude-dependent
gravitational acceleration.

  rocket_mass(t)                     Calculates rocket mass during and
after the burn.

  rocket_thrust(t)                   Returns motor thrust as a function
of time.

  rotation_matrix(phi, theta, psi)   Computes the body-to-inertial
rotation matrix.

  dynamics(t, state)                 Computes the 12 state derivatives.

  fix_rocket_distortion(vertices)    Corrects the graphical rocket mesh
for scene aspect-ratio distortion.

  rocket_model(index)                Builds the rocket's graphical
components at a trajectory sample.
exhaust_trace(index)               Creates the exhaust visual while the motor is burning.

Troubleshooting
- ModuleNotFoundError: Install the required packages using
  python -m pip install numpy scipy pandas plotly.
- HTML file does not open automatically: Open
  Rocket_6DOF_3D_Trajectory.html manually in a web browser.
- Output files are difficult to find: They are saved in the
  program's current working directory. Run the script from the folder
  where you want the outputs.
- Package installation fails: Confirm that the project root
  contains a valid pyproject.toml. If you only have the standalone
  .py file, install the dependencies and run the script directly.
  
Acknowledgement
This README describes the implementation and findings documented in the
project report, 6-DOF Rocket Trajectory Simulation and 3D
Visualisation.
