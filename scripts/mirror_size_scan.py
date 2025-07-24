import subprocess
import re
import shutil

GDML_FILE = "lo3.gdml"
BACKUP_FILE = GDML_FILE + ".bak"
RUN_CMD = ["./build/src/simg4ox", "-g", GDML_FILE, "-m", "esi-g4ox/run.mac"]
RMIN_RANGE = range(300, 2001, 50)  # 90, 100, ..., 300

def update_gdml(rmin, rmax):
    with open(BACKUP_FILE, "r") as f:
        lines = f.readlines()
    # Replace MirrorSphere rmin and rmax
    new_lines = []
    for line in lines:
        if '<sphere name="MirrorSphere"' in line:
            newline = re.sub(r'rmin="[^"]+"', f'rmin="{rmin}"', line)
            newline = re.sub(r'rmax="[^"]+"', f'rmax="{rmax}"', newline)
            new_lines.append(newline)
        else:
            new_lines.append(line)
    with open(GDML_FILE, "w") as f:
        f.writelines(new_lines)

def run_and_get_hits():
    result = subprocess.run(RUN_CMD, capture_output=True, text=True)
    for line in result.stdout.splitlines():
        m = re.search(r'PhotonSD::EndOfEvent Number of PhotonHits:\s*(\d+)', line)
        if m:
            return int(m.group(1))
    return None

def main():
    # Backup GDML first
    shutil.copyfile(GDML_FILE, BACKUP_FILE)
    results = []
    try:
        for rmin in RMIN_RANGE:
            rmax = rmin + 1
            print(f"Running for rmin={rmin}, rmax={rmax}")
            update_gdml(rmin, rmax)
            hits = run_and_get_hits()
            if hits is not None:
                print(f"  PhotonHits: {hits}")
                results.append((rmin, hits))
            else:
                print("  ERROR: PhotonHits not found!")
                results.append((rmin, "ERROR"))
    finally:
        # Restore original GDML
        shutil.move(BACKUP_FILE, GDML_FILE)
    # Write results
    with open("rmin_vs_photonhits.csv", "w") as f:
        f.write("rmin,PhotonHits\n")
        for rmin, hits in results:
            f.write(f"{rmin},{hits}\n")
    print("Done. Results in rmin_vs_photonhits.csv")

if __name__ == "__main__":
    main()

