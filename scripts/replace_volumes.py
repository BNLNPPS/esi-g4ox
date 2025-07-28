#!/usr/bin/env python3
"""
Replace all
    <volumeref ref="MirrorPyramid0x56008bfb13e0"/>
with
    <volumeref ref="MirrorPyramidAbsorbing"/>
*only* inside <physvol> blocks whose name attribute is one of
the names listed in `TARGET_NAMES`.
Usage:
    python3 replace_mirror_pyramid.py in.gdml out.gdml
"""

import argparse
import xml.etree.ElementTree as ET
from pathlib import Path

# ----------------------------------------------------------------------
# 1.  Configuration
# ----------------------------------------------------------------------
OLD_REF = "MirrorPyramid0x56008bfb13e0"
NEW_REF = "MirrorPyramidAbsorbing"

TARGET_NAMES = {
    "MirrorPyramidEdge_516", "MirrorPyramidEdge_515", "MirrorPyramidEdge_514",
    "MirrorPyramidEdge_513", "MirrorPyramidEdge_512",
    "MirrorPyramidEdge_616", "MirrorPyramidEdge_615", "MirrorPyramidEdge_614",
    "MirrorPyramidEdge_613", "MirrorPyramidEdge_612",
    "MirrorPyramidEdge_312", "MirrorPyramidEdge_311", "MirrorPyramidEdge_310",
    "MirrorPyramidEdge_309", "MirrorPyramidEdge_308",
    "MirrorPyramidEdge_412", "MirrorPyramidEdge_411", "MirrorPyramidEdge_410",
    "MirrorPyramidEdge_409", "MirrorPyramidEdge_408",
    "MirrorPyramidEdge_217", "MirrorPyramidEdge_216", "MirrorPyramidEdge_215",
    "MirrorPyramidEdge_214", "MirrorPyramidEdge_213",
    "MirrorPyramidEdge_207", "MirrorPyramidEdge_206", "MirrorPyramidEdge_200",
    "MirrorPyramid_negX", "MirrorPyramid_posX",
    "MirrorPyramid_negY",  "MirrorPyramid_posY",
}

# ----------------------------------------------------------------------
# 2.  Command‑line arguments
# ----------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Replace MirrorPyramid references in selected physvols.")
parser.add_argument("input_file",  type=Path, help="Path to the original GDML")
parser.add_argument("output_file", type=Path, help="Path to write the patched GDML")
args = parser.parse_args()

# ----------------------------------------------------------------------
# 3.  Parse, update, write
# ----------------------------------------------------------------------
tree = ET.parse(args.input_file)
root = tree.getroot()

for physvol in root.iter("physvol"):
    if physvol.get("name") in TARGET_NAMES:
        for child in physvol:
            if child.tag == "volumeref" and child.get("ref") == OLD_REF:
                child.set("ref", NEW_REF)        # <-- perform the substitution
                break                            # (at most one per physvol)

# Pretty‑print output (optional)
ET.indent(tree, space="  ", level=0)  # Python ≥3.9

tree.write(args.output_file, encoding="utf-8", xml_declaration=True)
print(f"Patched file written to {args.output_file}")
