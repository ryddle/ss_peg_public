import json, sys

path = sys.argv[1]
d = json.load(open(path, 'r', encoding='utf-8'))
cells = d['cells']
for i, c in enumerate(cells):
    ct = c['cell_type']
    cid = c['id']
    ex = c.get('execution_count', 'N/A')
    outs = len(c.get('outputs', []))
    src = c['source']
    first = src[0][:80].strip() if src else '(empty)'
    print(f"{i}: [{ct}] id={cid} exec={ex} outputs={outs} | {first}")
