# ================================================================
#       ROCKET 6-DOF TRAJECTORY SIMULATION
#       3D TRAJECTORY + PROPER UPRIGHT ROCKET MODEL
#       (FIXED: rocket model no longer distorted by the scene's
#        non-uniform manual aspect ratio)
# ================================================================
#
# Installing required libraries:
#
# pip install numpy scipy pandas plotly
#
# ================================================================

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.integrate import solve_ivp
import webbrowser
import os


# ================================================================
# 1. ROCKET PARAMETERS
# ================================================================

class Rocket:

    def __init__(self):

        # MASS
        self.dry_mass = 3.0
        self.propellant_mass = 2.0
        self.initial_mass = self.dry_mass + self.propellant_mass

        # PHYSICAL ROCKET
        self.length = 1.50
        self.diameter = 0.08
        self.radius = self.diameter / 2.0
        self.area = np.pi * self.radius ** 2

        # AERODYNAMICS
        self.CD = 0.40

        # MOTOR
        self.thrust = 150.0
        self.burn_time = 3.0
        self.mass_flow_rate = self.propellant_mass / self.burn_time

        # MOMENTS OF INERTIA
        self.Ixx = 0.02
        self.Iyy = 0.80
        self.Izz = 0.80


rocket = Rocket()


# ================================================================
# 2. ATMOSPHERIC MODEL
# ================================================================

def atmosphere(h):

    h = max(float(h), 0.0)

    T0 = 288.15
    P0 = 101325.0
    R = 287.05
    g0 = 9.80665
    lapse = 0.0065

    if h <= 11000.0:
        T = T0 - lapse * h
        P = P0 * (T / T0) ** (g0 / (R * lapse))
    else:
        T = 216.65
        P = 22632.06 * np.exp(-g0 * (h - 11000.0) / (R * T))

    rho = P / (R * T)

    return T, P, rho


# ================================================================
# 3. GRAVITY
# ================================================================

def gravity(h):

    Re = 6371000.0
    g0 = 9.80665

    return g0 * (Re / (Re + max(h, 0.0))) ** 2


# ================================================================
# 4. ROCKET MASS
# ================================================================

def rocket_mass(t):

    if t < rocket.burn_time:
        m = rocket.initial_mass - rocket.mass_flow_rate * t
        return max(m, rocket.dry_mass)

    return rocket.dry_mass


# ================================================================
# 5. THRUST
# ================================================================

def rocket_thrust(t):

    if 0.0 <= t <= rocket.burn_time:
        return rocket.thrust

    return 0.0


# ================================================================
# 6. EULER ROTATION MATRIX
# ================================================================

def rotation_matrix(phi, theta, psi):

    cp, sp = np.cos(phi), np.sin(phi)
    ct, st = np.cos(theta), np.sin(theta)
    cs, ss = np.cos(psi), np.sin(psi)

    R = np.array([
        [ct * cs, sp * st * cs - cp * ss, cp * st * cs + sp * ss],
        [ct * ss, sp * st * ss + cp * cs, cp * st * ss - sp * cs],
        [-st,     sp * ct,                cp * ct]
    ])

    return R


# ================================================================
# 7. 6-DOF EQUATIONS OF MOTION
# ================================================================

def dynamics(t, state):

    X, Y, Z = state[0], state[1], state[2]
    u, v, w = state[3], state[4], state[5]
    phi, theta, psi = state[6], state[7], state[8]
    p, q, r = state[9], state[10], state[11]

    _, _, rho = atmosphere(Z)

    m = rocket_mass(t)

    R = rotation_matrix(phi, theta, psi)

    V_body = np.array([u, v, w])
    V_inertial = R @ V_body

    V = np.linalg.norm(V_body)
    V_safe = max(V, 1e-8)

    q_dynamic = 0.5 * rho * V_safe ** 2

    drag = q_dynamic * rocket.area * rocket.CD
    drag_body = -drag * V_body / V_safe

    thrust_body = np.array([rocket_thrust(t), 0.0, 0.0])

    F_body = thrust_body + drag_body
    F_inertial = R @ F_body

    F_gravity = np.array([0.0, 0.0, -m * gravity(Z)])

    acceleration = (F_inertial + F_gravity) / m

    cos_theta = np.cos(theta)
    if abs(cos_theta) < 1e-6:
        cos_theta = 1e-6

    phi_dot = p + q * np.sin(phi) * np.tan(theta) + r * np.cos(phi) * np.tan(theta)
    theta_dot = q * np.cos(phi) - r * np.sin(phi)
    psi_dot = q * np.sin(phi) / cos_theta + r * np.cos(phi) / cos_theta

    p_dot = (rocket.Iyy - rocket.Izz) * q * r / rocket.Ixx
    q_dot = (rocket.Izz - rocket.Ixx) * p * r / rocket.Iyy
    r_dot = (rocket.Ixx - rocket.Iyy) * p * q / rocket.Izz

    return np.array([
        V_inertial[0], V_inertial[1], V_inertial[2],
        acceleration[0], acceleration[1], acceleration[2],
        phi_dot, theta_dot, psi_dot,
        p_dot, q_dot, r_dot
    ])


# ================================================================
# 8. INITIAL CONDITIONS
# ================================================================

launch_angle = np.radians(80.0)
initial_speed = 0.1

initial_state = np.array([
    0.0, 0.0, 0.0,
    initial_speed, 0.0, 0.0,
    0.0, -launch_angle, 0.0,
    0.0, 0.0, 0.0
])


# ================================================================
# 9. RUN SIMULATION
# ================================================================

print()
print("=" * 70)
print("6-DOF ROCKET TRAJECTORY SIMULATION")
print("=" * 70)
print()

t_final = 120.0
t_eval = np.linspace(0.0, t_final, 6000)

solution = solve_ivp(
    dynamics,
    (0.0, t_final),
    initial_state,
    t_eval=t_eval,
    rtol=1e-7,
    atol=1e-9,
    max_step=0.02
)

if not solution.success:
    raise RuntimeError("Simulation failed: " + solution.message)


# ================================================================
# 10. EXTRACT POSITION
# ================================================================

time = solution.t
x = solution.y[0]
y = solution.y[1]
z = solution.y[2]


# ================================================================
# 11. VELOCITY
# ================================================================

u = solution.y[3]
v = solution.y[4]
w = solution.y[5]
velocity = np.sqrt(u ** 2 + v ** 2 + w ** 2)


# ================================================================
# 12. ACCELERATION
# ================================================================

acceleration = np.zeros(len(time))

for i in range(len(time)):
    derivative = dynamics(time[i], solution.y[:, i])
    acceleration[i] = np.linalg.norm(derivative[3:6])


# ================================================================
# 13. APOGEE
# ================================================================

apogee_index = np.argmax(z)
maximum_altitude = z[apogee_index]
time_to_apogee = time[apogee_index]


# ================================================================
# 14. IMPACT
# ================================================================

impact_index = len(z) - 1

for i in range(apogee_index, len(z)):
    if z[i] <= 0.0:
        impact_index = i
        break

flight_time = time[impact_index]
downrange_distance = x[impact_index]
lateral_distance = y[impact_index]


# ================================================================
# 15. MAXIMUM VELOCITY
# ================================================================

maximum_velocity = np.max(velocity)


# ================================================================
# 16. MAXIMUM ACCELERATION
# ================================================================

maximum_acceleration = np.max(acceleration)
maximum_g = maximum_acceleration / 9.80665


# ================================================================
# 17. DISPLAY RESULTS
# ================================================================

print()
print("=" * 70)
print("SIMULATION RESULTS")
print("=" * 70)
print()
print(f"Maximum altitude       : {maximum_altitude:.2f} m")
print(f"Time to apogee         : {time_to_apogee:.2f} s")
print(f"Maximum velocity       : {maximum_velocity:.2f} m/s")
print(f"Maximum velocity       : {maximum_velocity * 3.6:.2f} km/h")
print(f"Maximum acceleration   : {maximum_acceleration:.2f} m/s\u00b2")
print(f"Maximum acceleration   : {maximum_g:.2f} g")
print(f"Downrange distance     : {downrange_distance:.2f} m")
print(f"Lateral distance       : {lateral_distance:.2f} m")
print(f"Flight time            : {flight_time:.2f} s")
print()


# ================================================================
# 18. SAVE SIMULATION DATA
# ================================================================

data = pd.DataFrame({
    "Time_s": time,
    "X_m": x,
    "Y_m": y,
    "Altitude_m": z,
    "Velocity_mps": velocity,
    "Velocity_kmph": velocity * 3.6,
    "Acceleration_mps2": acceleration,
    "Acceleration_g": acceleration / 9.80665,
    "Mass_kg": [rocket_mass(ti) for ti in time],
    "Thrust_N": [rocket_thrust(ti) for ti in time]
})

csv_file = os.path.abspath("Rocket_6DOF_Trajectory.csv")
data.to_csv(csv_file, index=False)


# ================================================================
# 19. TRAJECTORY DATA
# ================================================================

plot_end = impact_index + 1

px = x[:plot_end]
py = y[:plot_end]
pz = np.maximum(z[:plot_end], 0.0)


# ================================================================
# 19b. AXIS RANGES + ASPECT-DISTORTION CORRECTION  ***THE FIX***
# ================================================================
#
# WHY THIS IS NEEDED:
#
# The 3D scene uses a MANUAL aspect ratio (x:1.35, y:0.75, z:1.10)
# so the trajectory stays readable even when it barely moves
# sideways. But Plotly stretches each axis by a DIFFERENT amount
# to hit that ratio, because the real data ranges (downrange
# distance, lateral distance, altitude) are very different sizes
# -- e.g. lateral motion is nearly zero and gets floored to a
# 50 m span, while downrange/altitude can be hundreds of metres.
#
# The rocket mesh is built in plain x/y/z metres. Dropped into a
# scene where "1 metre" does NOT render the same size on every
# axis, the rocket gets sheared into a flat, twisted blob -- this
# is exactly what you saw in the screenshot.
#
# THE FIX: work out how many render-units one metre is worth on
# each axis, then pre-scale the rocket geometry (using X as the
# reference axis) so one metre looks the same size on every axis
# again, before the rocket is ever placed into the scene.
#
# ================================================================

x_min = np.min(px)
x_max = np.max(px)

y_min = np.min(py)
y_max = np.max(py)

z_max = np.max(pz)

x_span = max(x_max - x_min, 100.0)
y_span = max(y_max - y_min, 50.0)
z_span = max(z_max * 1.10, 100.0)

ASPECT_X = 1.35
ASPECT_Y = 0.75
ASPECT_Z = 1.10

render_units_per_metre_x = ASPECT_X / x_span
render_units_per_metre_y = ASPECT_Y / y_span
render_units_per_metre_z = ASPECT_Z / z_span

# X is the reference axis; Y and Z get corrected to match it so
# the rocket keeps true, undistorted proportions.

ROCKET_CORRECTION_X = 1.0
ROCKET_CORRECTION_Y = render_units_per_metre_x / render_units_per_metre_y
ROCKET_CORRECTION_Z = render_units_per_metre_x / render_units_per_metre_z


def fix_rocket_distortion(vertices):
    """
    Rescales rocket-part vertices so the rocket renders as a
    correctly-proportioned rocket despite the scene's non-uniform
    manual aspect ratio.
    """

    v = np.array(vertices, dtype=float).copy()

    v[:, 0] *= ROCKET_CORRECTION_X
    v[:, 1] *= ROCKET_CORRECTION_Y
    v[:, 2] *= ROCKET_CORRECTION_Z

    return v


# ================================================================
# 20. GRAPHICAL ROCKET SIZE
# ================================================================
#
# The real rocket is only 1.5 m long. A trajectory may be
# hundreds of metres long, so the graphical model is enlarged
# ONLY for visualization.
#
# ================================================================

trajectory_scale = max(np.ptp(px), np.ptp(py), np.ptp(pz), 100.0)

ROCKET_LENGTH = np.clip(trajectory_scale * 0.055, 20.0, 50.0)


# ================================================================
# 21. ROCKET GRAPHICAL DIMENSIONS
# ================================================================

BODY_LENGTH = ROCKET_LENGTH * 0.58
NOSE_LENGTH = ROCKET_LENGTH * 0.24
BODY_RADIUS = ROCKET_LENGTH * 0.055
FIN_ROOT_LENGTH = ROCKET_LENGTH * 0.27
FIN_HEIGHT = ROCKET_LENGTH * 0.25
FIN_THICKNESS = ROCKET_LENGTH * 0.035
NOZZLE_LENGTH = ROCKET_LENGTH * 0.08
NOZZLE_RADIUS = BODY_RADIUS * 0.72


# ================================================================
# 22. CYLINDER ALONG Z AXIS
# ================================================================
#
# The rocket is built along +Z. +Z = UP. Therefore the nose
# ALWAYS points upward.
#
# ================================================================

def cylinder_z(radius, z0, z1, n=32):

    vertices = []
    faces = []

    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)

    for a in angles:
        vertices.append([radius * np.cos(a), radius * np.sin(a), z0])

    for a in angles:
        vertices.append([radius * np.cos(a), radius * np.sin(a), z1])

    for i in range(n):
        j = (i + 1) % n
        faces.append([i, j, n + j])
        faces.append([i, n + j, n + i])

    bottom_center = len(vertices)
    vertices.append([0, 0, z0])

    for i in range(n):
        j = (i + 1) % n
        faces.append([bottom_center, j, i])

    top_center = len(vertices)
    vertices.append([0, 0, z1])

    for i in range(n):
        j = (i + 1) % n
        faces.append([top_center, n + i, n + j])

    return np.array(vertices), np.array(faces)


# ================================================================
# 23. NOSE CONE
# ================================================================

def nose_cone_z(radius, z_base, length, n=32):

    vertices = []
    faces = []

    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)

    for a in angles:
        vertices.append([radius * np.cos(a), radius * np.sin(a), z_base])

    tip = len(vertices)
    vertices.append([0.0, 0.0, z_base + length])

    for i in range(n):
        j = (i + 1) % n
        faces.append([i, j, tip])

    return np.array(vertices), np.array(faces)


# ================================================================
# 24. FIN
# ================================================================

def create_fin(angle):
    """Creates one aerodynamic fin, positioned radially around the
    rocket body. Four copies are used."""

    z_front = BODY_LENGTH * 0.12
    z_rear = BODY_LENGTH * 0.52

    z_outer_front = z_front - FIN_ROOT_LENGTH * 0.18
    z_outer_rear = z_rear + FIN_ROOT_LENGTH * 0.05

    r_root = BODY_RADIUS * 0.85
    r_outer = BODY_RADIUS + FIN_HEIGHT

    profile = [
        [r_root, 0, z_front],
        [r_root, 0, z_rear],
        [r_outer, 0, z_outer_rear],
        [r_outer, 0, z_outer_front]
    ]

    vertices = []

    for p in profile:
        vertices.append(p)

    for p in profile:
        vertices.append([p[0], FIN_THICKNESS, p[2]])

    c = np.cos(angle)
    s = np.sin(angle)

    for i in range(len(vertices)):
        X, Y, Z = vertices[i]
        vertices[i] = [X * c - Y * s, X * s + Y * c, Z]

    faces = [
        [0, 1, 2], [0, 2, 3],
        [4, 6, 5], [4, 7, 6],
        [0, 4, 5], [0, 5, 1],
        [1, 5, 6], [1, 6, 2],
        [2, 6, 7], [2, 7, 3],
        [3, 7, 4], [3, 4, 0]
    ]

    return np.array(vertices), np.array(faces)


# ================================================================
# 25. CREATE ROCKET COMPONENTS  (corrected geometry applied here)
# ================================================================

body_v, body_f = cylinder_z(BODY_RADIUS, 0.0, BODY_LENGTH)
body_v = fix_rocket_distortion(body_v)

nose_v, nose_f = nose_cone_z(BODY_RADIUS, BODY_LENGTH, NOSE_LENGTH)
nose_v = fix_rocket_distortion(nose_v)

fin_components = []

for angle in [0.0, np.pi / 2.0, np.pi, 3.0 * np.pi / 2.0]:
    fv, ff = create_fin(angle)
    fv = fix_rocket_distortion(fv)
    fin_components.append((fv, ff))

nozzle_v, nozzle_f = cylinder_z(NOZZLE_RADIUS, -NOZZLE_LENGTH, 0.0)
nozzle_v = fix_rocket_distortion(nozzle_v)


# ================================================================
# 26. TRANSFORM COMPONENT TO ROCKET POSITION
# ================================================================

def move_mesh(vertices, position):
    return vertices + position


# ================================================================
# 27. CREATE MESH TRACE
# ================================================================

def mesh_trace(vertices, faces, position, name, color):

    V = move_mesh(vertices, position)

    return go.Mesh3d(
        x=V[:, 0], y=V[:, 1], z=V[:, 2],
        i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
        color=color,
        opacity=1.0,
        flatshading=False,
        lighting=dict(ambient=0.30, diffuse=0.85, specular=0.80, roughness=0.25),
        name=name,
        hoverinfo="skip"
    )


# ================================================================
# 28. ROCKET MODEL
# ================================================================

def rocket_model(index):

    position = np.array([px[index], py[index], pz[index]])

    traces = [
        mesh_trace(body_v, body_f, position, "Rocket Body", "lightgray"),
        mesh_trace(nose_v, nose_f, position, "Pointed Nose", "red"),
    ]

    for number, (fv, ff) in enumerate(fin_components):
        traces.append(mesh_trace(fv, ff, position, f"Fin {number + 1}", "red"))

    traces.append(mesh_trace(nozzle_v, nozzle_f, position, "Engine Nozzle", "black"))

    return traces


# ================================================================
# 29. EXHAUST
# ================================================================

def exhaust_trace(index):

    position = np.array([px[index], py[index], pz[index]])

    if time[index] > rocket.burn_time:
        return go.Scatter3d(x=[], y=[], z=[], mode="lines", showlegend=False)

    # Corrected the same way as the nozzle so it stays attached
    start_z = -NOZZLE_LENGTH * ROCKET_CORRECTION_Z
    end_z = (-NOZZLE_LENGTH - ROCKET_LENGTH * 0.35) * ROCKET_CORRECTION_Z

    return go.Scatter3d(
        x=[position[0], position[0]],
        y=[position[1], position[1]],
        z=[position[2] + start_z, position[2] + end_z],
        mode="lines",
        line=dict(color="orange", width=14),
        showlegend=False,
        hoverinfo="skip"
    )


# ================================================================
# 30. TRAJECTORY
# ================================================================

trajectory_trace = go.Scatter3d(
    x=px, y=py, z=pz, mode="lines",
    line=dict(color="blue", width=7),
    name="6-DOF Trajectory"
)


# ================================================================
# 31. GROUND TRACK
# ================================================================

ground_trace = go.Scatter3d(
    x=px, y=py, z=np.zeros(len(px)), mode="lines",
    line=dict(color="gray", width=2, dash="dash"),
    name="Ground Track"
)


# ================================================================
# 32. LAUNCH MARKER
# ================================================================

launch_marker = go.Scatter3d(
    x=[px[0]], y=[py[0]], z=[pz[0]],
    mode="markers+text",
    marker=dict(size=10, color="green"),
    text=["LAUNCH"], textposition="top center",
    name="Launch"
)


# ================================================================
# 33. APOGEE MARKER
# ================================================================

apogee_marker = go.Scatter3d(
    x=[px[apogee_index]], y=[py[apogee_index]], z=[pz[apogee_index]],
    mode="markers+text",
    marker=dict(size=12, color="yellow", line=dict(color="black", width=2)),
    text=["APOGEE<br>" + f"{maximum_altitude:.1f} m"],
    textposition="top center",
    name="Apogee"
)


# ================================================================
# 34. IMPACT MARKER
# ================================================================

impact_marker = go.Scatter3d(
    x=[px[impact_index]], y=[py[impact_index]], z=[0],
    mode="markers+text",
    marker=dict(size=10, color="red"),
    text=["IMPACT"], textposition="top center",
    name="Impact"
)


# ================================================================
# 35. INITIAL ROCKET
# ================================================================

rocket_traces = rocket_model(0)
initial_exhaust = exhaust_trace(0)


# ================================================================
# 36. FIGURE
# ================================================================

fig = go.Figure(
    data=[trajectory_trace, ground_trace, launch_marker, apogee_marker, impact_marker]
    + rocket_traces
    + [initial_exhaust]
)


# ================================================================
# 37. TRACE INDEX
# ================================================================
#
# 0 = trajectory, 1 = ground, 2 = launch, 3 = apogee, 4 = impact
# 5 = body, 6 = nose, 7-10 = fins, 11 = nozzle, 12 = exhaust
#
# ================================================================

rocket_trace_indices = [5, 6, 7, 8, 9, 10, 11, 12]


# ================================================================
# 38. ANIMATION FRAMES
# ================================================================

frame_indices = np.linspace(0, impact_index, 160, dtype=int)
frame_indices = np.unique(frame_indices)

frames = []

for frame_number, index in enumerate(frame_indices):

    components = rocket_model(index)
    exhaust = exhaust_trace(index)

    frames.append(go.Frame(
        name=str(frame_number),
        data=components + [exhaust],
        traces=rocket_trace_indices
    ))

fig.frames = frames


# ================================================================
# 39. SLIDER
# ================================================================

slider_steps = []

for frame_number, index in enumerate(frame_indices):

    slider_steps.append({
        "label": f"{time[index]:.1f} s",
        "method": "animate",
        "args": [
            [str(frame_number)],
            {
                "mode": "immediate",
                "frame": {"duration": 0, "redraw": True},
                "transition": {"duration": 0}
            }
        ]
    })


# ================================================================
# 40. PLAY / PAUSE
# ================================================================

fig.update_layout(
    updatemenus=[{
        "type": "buttons",
        "direction": "left",
        "x": 0.02,
        "y": 1.08,
        "buttons": [
            {
                "label": "\u25b6 PLAY",
                "method": "animate",
                "args": [None, {
                    "frame": {"duration": 60, "redraw": True},
                    "transition": {"duration": 0},
                    "fromcurrent": True,
                    "mode": "immediate"
                }]
            },
            {
                "label": "\u23f8 PAUSE",
                "method": "animate",
                "args": [[None], {
                    "frame": {"duration": 0, "redraw": False},
                    "mode": "immediate"
                }]
            }
        ]
    }]
)


# ================================================================
# 41. SLIDER
# ================================================================

fig.update_layout(
    sliders=[{
        "active": 0,
        "currentvalue": {"prefix": "Time: "},
        "steps": slider_steps
    }]
)


# ================================================================
# 42. FINAL 3D VISUALIZATION
# ================================================================
# (axis ranges were already computed in section 19b, reused here)
# ================================================================

fig.update_layout(
    title={
        "text": "\U0001f680 6-DOF Rocket Trajectory Simulation",
        "x": 0.5,
        "font": {"size": 24}
    },

    scene={
        "xaxis": {
            "title": "Downrange X (m)",
            "range": [x_min - 0.05 * x_span, x_max + 0.05 * x_span]
        },
        "yaxis": {
            "title": "Lateral Y (m)",
            "range": [y_min - 0.15 * y_span, y_max + 0.15 * y_span]
        },
        "zaxis": {
            "title": "Altitude Z (m)",
            "range": [0, z_span]
        },

        # Z is vertical. Rocket is always constructed along +Z,
        # so the nose points UP. Its geometry has already been
        # corrected (section 19b) to render undistorted despite
        # this non-uniform manual aspect ratio.

        "aspectmode": "manual",
        "aspectratio": {"x": ASPECT_X, "y": ASPECT_Y, "z": ASPECT_Z},

        "camera": {
            "eye": {"x": 1.5, "y": 1.5, "z": 1.2},
            "center": {"x": 0.0, "y": 0.0, "z": 0.2}
        }
    },

    width=1350,
    height=850,
    showlegend=True,
    margin=dict(l=0, r=0, t=80, b=0)
)


# ================================================================
# 43. SAVE HTML
# ================================================================

html_file = os.path.abspath("Rocket_6DOF_3D_Trajectory.html")
fig.write_html(html_file, include_plotlyjs=True, auto_open=False)


# ================================================================
# 44. OPEN BROWSER
# ================================================================

print()
print("=" * 70)
print("3D ROCKET VISUALIZATION READY")
print("=" * 70)
print()
print("Rocket model:")
print("  \u2713 Cylindrical fuselage")
print("  \u2713 Pointed nose cone")
print("  \u2713 Four aerodynamic fins")
print("  \u2713 Engine nozzle")
print("  \u2713 Exhaust")
print("  \u2713 Nose permanently pointing upward")
print("  \u2713 Corrected for scene aspect distortion")
print("  \u2713 Animated along 6-DOF trajectory")
print()
print("HTML file:")
print(html_file)
print()
print("CSV file:")
print(csv_file)
print()
print("=" * 70)

try:
    webbrowser.open("file:///" + html_file.replace("\\", "/"))
except Exception:
    print("Open this file manually:")
    print(html_file)