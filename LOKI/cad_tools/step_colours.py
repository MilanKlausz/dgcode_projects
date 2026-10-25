"""Map STEP products to their surface colours (a proxy for the material), without FreeCAD.

Follows STYLED_ITEM -> styles -> COLOUR_RGB, and the styled geometric item -> its shape representation ->
PRODUCT_DEFINITION_SHAPE -> PRODUCT_DEFINITION -> PRODUCT. Prints, per product name, the colours found on its solids/faces.
usage: python step_colours.py FILE.stp [OUT.json]"""
import re
import sys
import json
import collections

ENT = re.compile(r'#(\d+)\s*=\s*([A-Z0-9_]+)\s*\((.*)\)\s*;', re.S)
REF = re.compile(r'#(\d+)')


def parse(fn):
    txt = open(fn, 'r', errors='replace').read()
    data = txt.split('DATA;', 1)[1]
    ents = {}
    for stmt in data.split(';\n'):
        stmt = stmt.strip()
        if not stmt.startswith('#'):
            continue
        m = ENT.match(stmt + ';')
        if m:
            ents[int(m.group(1))] = (m.group(2), m.group(3))
        else:                                    # complex entity: #n = ( A() B() ... );
            mm = re.match(r'#(\d+)\s*=\s*\((.*)\)\s*$', stmt, re.S)
            if mm:
                ents[int(mm.group(1))] = ('COMPLEX', mm.group(2))
    return ents


def refs(s):
    return [int(x) for x in REF.findall(s)]


def main(fn, out=None):
    E = parse(fn)
    # colour of each STYLED_ITEM's item
    def colour_of(eid, depth=0):
        t, a = E.get(eid, ('', ''))
        if t == 'COLOUR_RGB':
            v = [float(x) for x in re.findall(r'[-+]?\d*\.?\d+(?:E[-+]?\d+)?', a.split(',', 1)[1])]
            return tuple(round(x, 3) for x in v[:3])
        if depth > 8:
            return None
        for r in refs(a):
            c = colour_of(r, depth + 1)
            if c:
                return c
        return None
    item_col = {}
    for eid, (t, a) in E.items():
        if t in ('STYLED_ITEM', 'OVER_RIDING_STYLED_ITEM'):
            rr = refs(a)
            if not rr:
                continue
            item = rr[-1]
            c = None
            for r in rr[:-1]:
                c = colour_of(r)
                if c:
                    break
            if c:
                item_col.setdefault(item, set()).add(c)
    # which representation contains which item (one level, plus faces of solids via the brep)
    rep_items = {}
    for eid, (t, a) in E.items():
        if t.endswith('SHAPE_REPRESENTATION') or t == 'SHAPE_REPRESENTATION':
            rep_items[eid] = refs(a)[:-1]
    # product of each representation
    pd_name = {}
    for eid, (t, a) in E.items():
        if t == 'PRODUCT':
            pd_name[eid] = re.findall(r"'([^']*)'", a)[0]
    pdf_prod = {eid: refs(a)[-1] for eid, (t, a) in E.items() if t in ('PRODUCT_DEFINITION_FORMATION', 'PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE')}
    pdef_pdf = {eid: refs(a)[0] for eid, (t, a) in E.items() if t == 'PRODUCT_DEFINITION'}
    pds_pdef = {eid: refs(a)[0] for eid, (t, a) in E.items() if t == 'PRODUCT_DEFINITION_SHAPE'}
    rep_prod = {}
    for eid, (t, a) in E.items():
        if t == 'SHAPE_DEFINITION_REPRESENTATION':
            pds, rep = refs(a)[:2]
            p = pd_name.get(pdf_prod.get(pdef_pdf.get(pds_pdef.get(pds))))
            if p:
                rep_prod[rep] = p
    # brep representations linked to the product's shape representation
    for eid, (t, a) in E.items():
        if t in ('SHAPE_REPRESENTATION_RELATIONSHIP',) or ('SHAPE_REPRESENTATION_RELATIONSHIP' in a and t == 'COMPLEX'):
            rr = refs(a)
            if len(rr) >= 2:
                r1, r2 = rr[0], rr[1]
                if r1 in rep_prod and r2 not in rep_prod:
                    rep_prod[r2] = rep_prod[r1]
                elif r2 in rep_prod and r1 not in rep_prod:
                    rep_prod[r1] = rep_prod[r2]
    # map solids' faces to their solid: CLOSED_SHELL -> faces, MANIFOLD_SOLID_BREP -> shell
    face_solid = {}
    for eid, (t, a) in E.items():
        if t == 'MANIFOLD_SOLID_BREP':
            sh = refs(a)[0]
            for f in refs(E.get(sh, ('', ''))[1]):
                face_solid[f] = eid
    prod_col = collections.defaultdict(collections.Counter)
    for rep, items in rep_items.items():
        p = rep_prod.get(rep)
        if not p:
            continue
        for it in items:
            for c in item_col.get(it, ()):
                prod_col[p][c] += 1
    for it, cs in item_col.items():
        s = face_solid.get(it)
        if s is None:
            continue
        for rep, items in rep_items.items():
            if s in items and rep in rep_prod:
                for c in cs:
                    prod_col[rep_prod[rep]][c] += 1
    res = {p: {','.join(map(str, c)): n for c, n in cnt.items()} for p, cnt in prod_col.items()}
    allc = collections.Counter()
    for p, cnt in prod_col.items():
        for c in cnt:
            allc[c] += 1
    print(f'{len(res)} products with colours; colours (number of products):')
    for c, n in allc.most_common():
        ex = [p for p, cnt in prod_col.items() if c in cnt][:8]
        print(f'  {c}: {n}  e.g. {ex}')
    if out:
        json.dump(res, open(out, 'w'), indent=0)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
