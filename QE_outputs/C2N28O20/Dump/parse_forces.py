import sys
import json
out = sys.argv[1]
txt = open(out).read()
if "JOB DONE" not in txt:
    sys.exit(1)
blocks = txt.split("Forces acting on atoms")
if len(blocks) < 2:
    sys.exit(1)
forces = []
for line in blocks[-1].splitlines():
    s = line.strip()
    if s.startswith("atom") and "force" in s:
        p = s.split()
        forces.append([float(p[-3]), float(p[-2]), float(p[-1])])
    elif forces and not s:
        break
print(json.dumps({"name": out, "forces": forces}))
