import json, sys

path = sys.argv[1]
cell_idx = int(sys.argv[2])
d = json.load(open(path, 'r', encoding='utf-8'))
c = d['cells'][cell_idx]

print("=== SOURCE ===")
for line in c['source']:
    print(line, end='')

print("\n\n=== OUTPUTS ===")
for out in c.get('outputs', []):
    print(f"[{out['output_type']}] name={out.get('name','N/A')}")
    for t in out.get('text', []):
        print(t, end='')
    print()
