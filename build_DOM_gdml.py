#!/usr/bin/env python3
"""
Generate a GDML geometry that instantiates an OuterTube_logical (shell)
*and* a bare Tube_logical (PMT) at every (x,y,z) coordinate given below.

The Tube_logical copy is shifted 0.5 mm farther along its look‑direction
so that its front face sits flush with the outer shell’s open end.
"""

import math
from pathlib import Path
import numpy as np
# --------------------------------------------------------------------
# 1.  User‑supplied coordinates and tube axes
# --------------------------------------------------------------------
x= [-1.3515787011974438e-14, 0.014519848978303414, -93.98700187061975, -94.00152171959806, -0.014519848978317153, 93.98700187061976, 94.00152171959806, 0.01876423690348331, -166.9753290817535, -164.79502539694286, -83.47895639584895, -164.81378963384634, 83.49637268590448, -0.018764236903487032, 166.97532908175347, 164.79502539694286, 83.47895639584895, 164.81378963384634, -83.49637268590449, 0.007588696125238299, -190.30999558119555, -144.60122657909912, -95.18599955861427, -144.60881527522434, 95.12399602258122, -0.007588696125244882, 190.30999558119552, 144.60122657909912, 95.18599955861428, 144.60881527522434, -95.12399602258124]
y= [-157.59110560319414, -127.2078144162383, -127.2078144162383, -127.20781441623785, -127.2078144162383, -127.2078144162383, -127.2078144162383, -16.21282012919346, -67.98421361231476, -16.21282012919346, -67.98421361231476, -16.21282012919346, -67.98421361231476, -16.21282012919346, -67.98421361231476, -16.21282012919346, -67.98421361231476, -16.21282012919346, -67.98421361231476, 144.72952571926407, 92.9405067056814, 144.72952571926407, 92.9405067056814, 144.72952571926362, 92.9405067056814, 144.72952571926407, 92.9405067056814, 144.72952571926407, 92.9405067056814, 144.72952571926407, 92.9405067056814]
z= [-1.2835610710505035e-14, 108.53522469937256, 54.28018690776063, -54.255037791612004, -108.53522469937259, -54.280186907760616, 54.25503779161196, 190.29973808529974, -0.010055299751782388, 95.16611934849095, -144.609904439941, -95.13361873680888, -144.59984914018918, -190.29973808529974, 0.010055299751841983, -95.1661193484909, 144.609904439941, 95.13361873680884, 144.59984914018918, 166.9754955169366, 0.03579775821946067, 83.4943197620944, -164.79539188830992, -83.4811757548423, -164.8311896465293, -166.9754955169366, -0.03579775821937266, -83.49431976209435, 164.79539188830995, 83.48117575484225, 164.8311896465293]
dx= [-6.838714497526159e-17, 7.334204760100329e-05, -0.4747431722858066, -0.4748165143334076, -7.334204760107293e-05, 0.47474317228580665, 0.47481651433340766, 9.478265156253965e-05, -0.8434531507819575, -0.8324191147117075, -0.4216825873814833, -0.83251389736327, 0.4217705634004739, -9.478265156255871e-05, 0.8434531507819575, 0.8324191147117076, 0.4216825873814833, 0.83251389736327, -0.421770563400474, 3.8331405102609264e-05, -0.9612534487736069, -0.7303979633475706, -0.48078331393604345, -0.7304362947526732, 0.480470134837563, -3.8331405102642745e-05, 0.9612534487736066, 0.7303979633475707, 0.4807833139360434, 0.7304362947526731, -0.48047013483756307]
dy= [-0.8058814181034215, -0.6422099790192162, -0.6422099790192162, -0.642209979019216, -0.6422099790192162, -0.6422099790192162, -0.6422099790192162, -0.08156120720162435, -0.3430841543735022, -0.08156120720162435, -0.3430841543735021, -0.0815612072016243, -0.3430841543735022, -0.08156120720162435, -0.3430841543735022, -0.08156120720162435, -0.3430841543735022, -0.0815612072016243, -0.3430841543735022, 0.7313842722065715, 0.4697846380716747, 0.7313842722065715, 0.4697846380716747, 0.7313842722065715, 0.46978463807167464, 0.7313842722065715, 0.46978463807167464, 0.7313842722065715, 0.46978463807167464, 0.7313842722065715, 0.46978463807167464]
dz= [-5.032556088088008e-17, 0.5482285406812201, 0.2741777864169982, -0.2740507542642222, -0.5482285406812201, -0.2741777864169981, 0.274050754264222, 0.9612495227041828, -5.0792978246163326e-05, 0.48070684553618287, -0.730477251968325, -0.48054267716800053, -0.7304264589900785, -0.9612495227041828, 5.079297824649338e-05, -0.48070684553618254, 0.730477251968325, 0.48054267716800025, 0.7304264589900785, 0.8434130521556031, 0.0001808140368125949, 0.42173972204838345, -0.832379499094941, -0.42167333010722025, -0.832560313131753, -0.8434130521556031, -0.00018081403681212128, -0.4217397220483831, 0.8323794990949409, 0.42167333010721997, 0.832560313131753]


xyz   = np.column_stack([x, y, z]).astype(float)
udir  = np.column_stack([dx, dy, dz]).astype(float)

# normalise the direction vectors (in case they aren’t already)
udir /= np.linalg.norm(udir, axis=1)[:, None]

# angle between position vector and axis
dot   = np.sum(xyz * udir, axis=1) / np.linalg.norm(xyz, axis=1)
angle = np.degrees(np.arccos(np.clip(dot, -1.0, 1.0)))   # 0° = perfect

tol = 1e-3             # angular tolerance in degrees (adjust to taste)
bad = np.where(angle > tol)[0]

if bad.size == 0:
    print(f"✔  All {len(x)} tubes point outward (|θ| ≤ {tol}°).")
else:
    print(f"✘  {bad.size} / {len(x)} tubes outside tolerance {tol}°:")
    for i in bad:
        print(f"   #{i:2d}: angle = {angle[i]:.4f}°  "
              f"pos=({x[i]:.2f},{y[i]:.2f},{z[i]:.2f})  "
              f"axis=({dx[i]:.3g},{dy[i]:.3g},{dz[i]:.3g})")

# Ensure all four arrays are the same length
assert len({len(x), len(y), len(z), len(dx)}) == 1, "Array length mismatch!"
n = len(x)

# --------------------------------------------------------------------
# 2.  Helper:  direction‑vector ➜ Euler angles (deg)
#     We rotate the tube’s local +Z axis to the desired unit vector v.
#     A convenient intrinsic sequence is  RotY(−θ) · RotZ(φ)
#     where φ = atan2(dy,dx) and θ = acos(dz).
# --------------------------------------------------------------------
from scipy.spatial.transform import Rotation as R

def direction_to_euler(dx, dy, dz):
    direction = np.array([dx, dy, dz])
    direction /= np.linalg.norm(direction)

    # Original axis is +Z
    original = np.array([0, 0, 1])

    # Check if vectors are already aligned or opposite
    if np.allclose(direction, original):
        return (0.0, 0.0, 0.0)
    elif np.allclose(direction, -original):
        # Rotate 180 deg about X or Y axis (arbitrary choice)
        return (180.0, 0.0, 0.0)

    # Axis-angle rotation from original to direction
    rot_axis = np.cross(original, direction)
    rot_axis /= np.linalg.norm(rot_axis)
    rot_angle = np.arccos(np.clip(np.dot(original, direction), -1.0, 1.0))

    # Convert axis-angle to rotation matrix
    rotvec = rot_axis * rot_angle
    rotation = R.from_rotvec(rotvec)

    # GDML rotations are extrinsic XYZ rotations
    euler_angles = rotation.as_euler('xyz', degrees=True)

    return tuple(euler_angles)

# --------------------------------------------------------------------
# 3.  Fixed parts of the GDML (header & footer taken from your template)
# --------------------------------------------------------------------
HEAD = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
      xsi:noNamespaceSchemaLocation="">

  <!-- 1) Define any matrices, etc. (optional) -->
  <define>
    <matrix coldim="2" name="REFLECTIVITY_MIRROR" values="1.55e-06 1.0  1.55e-05 1.0" />
    <matrix coldim="2" name="REFLECTIVITYZero" values="1.55e-06 0.0  1.55e-05 0.0" />
    <matrix coldim="2" name="REFLECTIVITYMAX" values="1.55e-06 1  1.55e-05 1" />
    <matrix coldim="2" name="RINDEX_MIRROR" values="1.55e-06 1.0 1.55e-05 1.0" />
    <matrix coldim="2" name="RINDEX_WATER" values="1.55e-06 1.333 1.55e-05 1.333" />
    <matrix coldim="2" name="ABSLENGTH_WATER" values="1.55e-06 100.0*m 1.55e-05 100.0*m" />
    <matrix coldim="2" name="SCATTERLENGTH_WATER" values="1.55e-06 50.0*m 1.55e-05 50.0*m" />
    <matrix coldim="2" name="RINDEX_GLASS" values="1.55e-06 1.5  1.55e-05 1.5" />
    <matrix coldim="2" name="EFFICIENCYDET" values="2.034*eV 1 2.068*eV 1 2.103*eV 1 2.139*eV 1 2.177*eV 1 2.216*eV 1 2.256*eV 1 2.298*eV 1 2.341*eV 1 2.386*eV 1 2.433*eV 1 2.481*eV 1 2.532*eV 1 2.585*eV 1 2.64*eV 1 2.697*eV 1 2.757*eV 1 2.82*eV 1 2.885*eV 1 2.954*eV 1 3.026*eV 1 3.102*eV 1 3.181*eV  1 3.265*eV 1 3.353*eV 1 3.446*eV 1 3.545*eV 1 3.649*eV 1 3.76*eV 1 3.877*eV 1 4.002*eV 1 4.136*eV 1 5.0*eV 1 6.0*eV 1 7.0*eV 1 8.0*eV 1 9.0*eV 1 9.0*eV 1 10.0*eV 1 11.0*eV 1 11.6*eV 1" />
  </define>

  <!-- 2) Materials -->
  <materials>
    <element name="H" formula="H" Z="1">
      <atom value="1.0079" unit="g/mole"/>
    </element>
    <element name="O" formula="O" Z="8">
      <atom value="15.999" unit="g/mole"/>
    </element>
    <material name="WaterMaterial" state="liquid">
      <D value="1.0" unit="g/cm3"/>
      <fraction n="0.1119" ref="H"/>
      <fraction n="0.8881" ref="O"/>
      <property name="RINDEX" ref="RINDEX_WATER"/>
      <property name="ABSLENGTH" ref="ABSLENGTH_WATER"/>
      <property name="SCATTERLENGTH" ref="SCATTERLENGTH_WATER"/>
    </material>
    <material name="MirrorMaterial" state="solid">
      <D value="2.7" unit="g/cm3"/>
      <fraction n="1.0" ref="O"/>
      <property name="REFLECTIVITY" ref="REFLECTIVITY_MIRROR"/>
      <property name="RINDEX" ref="RINDEX_MIRROR"/>
    </material>
    <material name="GlassMaterial" state="solid">
      <D value="2.5" unit="g/cm3"/>
      <fraction n="1.0" ref="O"/>
      <property name="RINDEX" ref="RINDEX_GLASS"/>
    </material>
    <material name="TubeMaterial" state="solid">
      <D value="1.2" unit="g/cm3"/>
      <fraction n="1.0" ref="O"/>
    </material>
    <material name="OuterTubeMaterial" state="solid">
      <D value="1.4" unit="g/cm3"/>
      <fraction n="1.0" ref="O"/>
    </material>
  </materials>

  <!-- 3) Solids: geometry -->
  <solids>
    <box name="WorldBox" x="100000" y="100000" z="100000" lunit="mm" />
    <sphere name="MirrorSphere" rmin="300" rmax="301" deltaphi="6.28318530718" deltatheta="1.5707963268" aunit="rad" lunit="mm"/>
    <box name="GlassSphere" x="30" y="30" z="30" lunit="mm"/>

 <!-- Tube primitives -->
    <tube name="Tube1" rmin="0" rmax="23"   z="120"   deltaphi="360" aunit="deg" lunit="mm"/>
    <tube name="Tube2" rmin="0" rmax="23.1" z="121"  deltaphi="360" aunit="deg" lunit="mm"/>
    <tube name="Tube3" rmin="0" rmax="23"   z="120.1" deltaphi="360" aunit="deg" lunit="mm"/>

    <!-- Outer tube shell created by Boolean subtraction (inner tube shifted +0.5 mm) -->
    <subtraction name="OuterTube">
      <first  ref="Tube2"/>
      <second ref="Tube3">
        <positionref ref="Tube1Offset"/>
      </second>
    </subtraction>

  <!-- Optical Surfaces -->
    <opticalsurface name="MirrorSurface" type="dielectric_metal" model="glisur" finish="polished" value="1">
      <property name="REFLECTIVITY" ref="REFLECTIVITY_MIRROR" />
      <property name="RINDEX" ref="RINDEX_MIRROR" />
    </opticalsurface>
    <opticalsurface finish="0" model="0" name="Det_optical" type="0" value="1">
      <property name="EFFICIENCY" ref="EFFICIENCYDET"/>
      <property name="REFLECTIVITY" ref="REFLECTIVITYMAX"/>
    </opticalsurface>
    <opticalsurface finish="0" model="0" name="Absorb_optical" type="0" value="1">
      <property name="REFLECTIVITY" ref="REFLECTIVITYZero"/>
    </opticalsurface>
  </solids>

  <!-- 5) Structure -->
  <structure>
    <volume name="Mirror_logical">
      <materialref ref="MirrorMaterial"/>
      <solidref ref="MirrorSphere"/>
    </volume>
    <volume name="Tube_logical">
      <materialref ref="TubeMaterial"/>
      <auxiliary auxtype="SensDet" auxvalue="PhotonDetector" />
      <solidref ref="Tube1"/>
    </volume>
<volume name="OuterTube_logical">
  <materialref ref="OuterTubeMaterial"/>
  <solidref ref="OuterTube"/>
  
  <!-- Now nest the inner tube inside the outer one -->
  <physvol>
    <volumeref ref="Tube_logical"/>
    <position name="InnerTubePos" unit="mm" x="0" y="0" z="0.5"/>
    <rotation name="InnerTubeRot" unit="deg" x="0" y="0" z="0"/>
  </physvol>
</volume>


    <volume name="World_logical">
      <materialref ref="WaterMaterial"/>
      <solidref ref="WorldBox"/>
      <physvol>
        <volumeref ref="Mirror_logical"/>
        <position name="MirrorPos" unit="mm" x="0" y="0" z="0"/>
      </physvol>
"""

FOOT = """
    </volume>   <!-- World_logical -->

    <!-- Skin‑surfaces (unchanged) -->
    <skinsurface name="MirrorSkinSurface" surfaceproperty="MirrorSurface">
      <volumeref ref="Mirror_logical"/>
    </skinsurface>
    <skinsurface name="GlassSkinSurface" surfaceproperty="Det_optical">
      <volumeref ref="Tube_logical"/>
    </skinsurface>
    <skinsurface name="GlassSkinSurface" surfaceproperty="Absorb_optical">
      <volumeref ref="OuterTube_logical"/>
    </skinsurface>
  </structure>

  <!-- 6) Setup -->
  <setup name="Default" version="1.0">
    <world ref="World_logical"/>
  </setup>
</gdml>
"""

# --------------------------------------------------------------------
# 4.  Build placement strings
# --------------------------------------------------------------------
placements = []

import numpy as np

def rotation_matrix_from_vectors(vec1, vec2):
    """Returns the rotation matrix that aligns vec1 to vec2"""
    a, b = (vec1 / np.linalg.norm(vec1)), (vec2 / np.linalg.norm(vec2))
    v = np.cross(a, b)
    c = np.dot(a, b)
    if np.isclose(c, -1.0):
        # 180 degree rotation special case
        orthogonal = np.array([1, 0, 0])
        if (np.abs(a[0]) > 0.9):
            orthogonal = np.array([0, 1, 0])
        v = np.cross(a, orthogonal)
        v /= np.linalg.norm(v)
        return R.from_rotvec(np.pi * v).as_matrix()
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

placements = []

ROT = {
    0:  (-90.00,   0.00, 0.0),
    1:  (-56.79,   0.00, 0.0),
    2:  (-74.09,  29.20, 0.0),
    3:  (-106.18, 29.20, 0.0),
    4:  (-146.79,  0.00, 0.0),
    5:  (-106.18,-29.54, 0.0),
    6:  (-74.09, -29.54, 0.0),
    7:  (-16.00,   0.00, 0.0),
    8:  (-90.00,  57.30, 0.0),
    9:  (-30.30,  71.40, 0.0),
    10: (-136.96, 38.12, 0.0),
    11: (-118.72, 71.67, 0.0),
    12: (-136.95,-38.12, 0.0),
    13: (-164.00,  0.00, 0.0),
    14: (-90.00, -57.49, 0.0),
    15: (-118.74,-71.66, 0.0),
    16: (-43.04, -38.12, 0.0),
    17: (-61.26, -71.67, 0.0),
    18: (-43.05,  38.12, 0.0),
    19: ( 57.50,   0.00, 0.0),
    20: ( 89.99,  74.00, 0.0),
    21: ( 65.05,  53.64, 0.0),
    22: (146.33,  60.00, 0.0),
    23: (114.95,  53.65, 0.0),
    24: (146.36, -60.00, 0.0),
    25: (147.50,   0.00, 0.0),
    26: ( 89.99, -74.00, 0.0),
    27: (114.95, -53.64, 0.0),
    28: ( 33.67, -60.00, 0.0),
    29: ( 65.05, -53.65, 0.0),
    30: ( 33.64,  60.00, 0.0),
}

# --------------------------------------------------------------------
# 5.  Build <physvol> blocks with your fixed angles
# --------------------------------------------------------------------
placements = []
shift_mm = 0.5

for i in range(n):
    # --- position vectors ---
    px, py, pz = x[i], y[i], z[i]
    vx, vy, vz = dx[i], dy[i], dz[i]

    # --- Euler angles for this index ---
    rx, ry, rz = ROT[i]          # OuterTubeRot_i  (and TubeRot_i)

    # ---------- outer shell ----------
    placements.append(f"""
      <!-- PMT #{i} – outer shell -->
      <physvol>
        <volumeref ref="OuterTube_logical"/>
        <position name="OuterTubePos_{i}" unit="mm"
                  x="{px:.6f}" y="{py:.6f}" z="{pz:.6f}"/>
        <rotation name="OuterTubeRot_{i}" unit="deg"
                  x="{rx:.2f}" y="{ry:.2f}" z="{rz:.2f}"/>
      </physvol>""")

    # ---------- inner tube (shifted 0.5 mm) ----------
    norm = math.sqrt(vx*vx + vy*vy + vz*vz)
    p_inner = (px + shift_mm*vx/norm,
               py + shift_mm*vy/norm,
               pz + shift_mm*vz/norm)

    placements.append(f"""
      <!-- PMT #{i} – sensitive inner tube -->
      <physvol>
        <volumeref ref="Tube_logical"/>
        <position name="TubePos_{i}" unit="mm"
                  x="{p_inner[0]:.6f}" y="{p_inner[1]:.6f}" z="{p_inner[2]:.6f}"/>
        <rotation name="TubeRot_{i}" unit="deg"
                  x="{rx:.2f}" y="{ry:.2f}" z="{rz:.2f}"/>
      </physvol>""")


# --------------------------------------------------------------------
# 5.  Write out full GDML
# --------------------------------------------------------------------
gdml_text = HEAD + "\n".join(placements) + FOOT
Path("auto_pmt_array.gdml").write_text(gdml_text)
print(f"Wrote {n} PMT pairs to auto_pmt_array.gdml")
