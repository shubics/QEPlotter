"""
Client-side 3D structure viewer (3Dmol.js).

We only build an HTML/JS snippet here; the actual WebGL rendering happens in the
user's browser via 3Dmol.js loaded from a CDN. The server does *not* render
anything, which keeps its load near zero (a hard project requirement).

The snippet is meant to be embedded with ``streamlit.components.v1.html``.
"""
import json
import os
import tempfile
import colorsys
import hashlib

from ase.io import write as ase_write

# Matplotlib-free style presets mapped to 3Dmol style dicts.
_STYLES = {
    "ball-stick": {"sphere": {"scale": 0.25}},
    "stick": {"sphere": {"scale": 0.13}},
    "spacefill": {"sphere": {}},
    "wireframe": {"sphere": {"scale": 0.10, "opacity": 0.72}},
}

_BOND_RADII = {"ball-stick": 0.11, "stick": 0.14,
               "wireframe": 0.025, "spacefill": 0.0}

_CDN = "https://3Dmol.org/build/3Dmol-min.js"

# Deliberately not the conventional CPK colours: several CPK entries are very
# similar on a dark canvas.  This palette optimises category separation so a
# three-element structure always has three clearly different colours.
_ELEMENT_PALETTE = (
    "#38BDF8", "#FB7185", "#A3E635", "#FBBF24", "#C084FC",
    "#2DD4BF", "#F97316", "#818CF8", "#F472B6", "#4ADE80",
    "#FACC15", "#22D3EE", "#E879F9", "#60A5FA", "#BEF264",
)


def element_colors(atoms):
    """Return a stable, high-contrast colour for every element in *atoms*."""
    elements = list(dict.fromkeys(atoms.get_chemical_symbols()))
    colors = list(_ELEMENT_PALETTE)
    for i in range(len(colors), len(elements)):
        # Golden-angle hue stepping extends the palette without ever reusing a
        # colour in unusually composition-rich structures.
        hue = (i * 0.61803398875) % 1.0
        red, green, blue = colorsys.hsv_to_rgb(hue, 0.68, 0.96)
        colors.append(f"#{round(red * 255):02X}{round(green * 255):02X}{round(blue * 255):02X}")
    return dict(zip(elements, colors))


def _atoms_to_cif(atoms):
    # ASE's CIF writer expects a binary file handle, so go via a temp file.
    fd, path = tempfile.mkstemp(suffix=".cif")
    os.close(fd)
    try:
        ase_write(path, atoms, format="cif")
        with open(path, "r") as f:
            return f.read()
    finally:
        os.remove(path)


def build_3dmol_html(
    atoms,
    style="ball-stick",
    show_cell=True,
    supercell=(1, 1, 1),
    height=480,
    background="#15191E",
    highlight=None,
    bond_tol=1.15,
    periodic_bonds=False,
):
    """
    Build a self-contained 3Dmol.js viewer HTML string for an ASE ``Atoms``.

    Parameters
    ----------
    atoms : ase.Atoms
    style : str
        One of ``ball-stick``, ``stick``, ``spacefill``, ``wireframe``.
    show_cell : bool
        Draw the unit-cell box.
    supercell : tuple of int
        (nx, ny, nz) repetitions, applied server-side before export.
    height : int
        Viewer height in pixels.
    background : str
        CSS/hex background colour.

    Returns
    -------
    str
        HTML to pass to ``st.components.v1.html(..., height=height)``.
    """
    nx, ny, nz = (max(1, int(n)) for n in supercell)
    disp = atoms.repeat((nx, ny, nz)) if (nx, ny, nz) != (1, 1, 1) else atoms

    base_cif = _atoms_to_cif(atoms)
    cif = base_cif if disp is atoms else _atoms_to_cif(disp)
    style_dict = _STYLES.get(style, _STYLES["ball-stick"])
    colors = element_colors(disp)

    # 3Dmol's automatic CIF bonding ignores the UI tolerance. Replace it with
    # the exact same ASE/PBC geometry used by the analysis table.
    from qeplotter.structure.bonds import analyse_bonds, atom_labels
    _, visible_bonds = analyse_bonds(
        disp, tol=bond_tol, include_periodic=periodic_bonds)
    bond_segments = [{"points": record["points"],
                      "elements": record["point_elements"]}
                     for record in visible_bonds]
    labels = atom_labels(disp)

    highlight = highlight or {}

    cif_js = json.dumps(cif)
    style_js = json.dumps(style_dict)
    colors_js = json.dumps(colors)
    bonds_js = json.dumps(bond_segments)
    labels_js = json.dumps(labels)
    highlight_js = json.dumps(highlight)
    bond_radius_js = json.dumps(_BOND_RADII.get(style, 0.11))
    bg_js = json.dumps(background)
    show_cell_js = "true" if show_cell else "false"
    camera_key_js = json.dumps(
        "qeplotter-camera-" + hashlib.sha1(base_cif.encode("utf-8")).hexdigest()[:16])

    # Unique id so multiple embeds / Streamlit reruns don't collide.
    viewer_state = (cif, style, show_cell, supercell, float(bond_tol),
                    bool(periodic_bonds),
                    json.dumps(highlight, sort_keys=True))
    div_id = f"viewer_{abs(hash(viewer_state)) % (10 ** 8)}"

    return f"""
<div style="position:relative">
  <div id="{div_id}" style="width:100%;height:{height}px;position:relative;
       border-radius:5px;overflow:hidden;border:1px solid #3A424C;"></div>
  <div id="{div_id}_legend" style="position:absolute;left:12px;bottom:12px;
       display:flex;flex-wrap:wrap;gap:6px;padding:7px 9px;border-radius:3px;
       background:rgba(17,20,24,.90);color:#E7EAED;font:12px sans-serif"></div>
  <div id="{div_id}_views" style="position:absolute;right:12px;top:12px;
       display:flex;gap:4px;padding:5px;border-radius:3px;
       background:rgba(17,20,24,.92);font:12px sans-serif">
    <button data-view="top">Top</button><button data-view="front">Front</button>
    <button data-view="right">Right</button><button data-view="left">Left</button>
  </div>
  <div style="position:absolute;right:12px;bottom:12px;padding:6px 9px;border-radius:3px;
       background:rgba(17,20,24,.90);color:#C6CDD5;font:12px sans-serif">
       Click atom → pin identity</div>
  {('<div style="position:absolute;left:12px;top:12px;padding:7px 10px;border-radius:3px;'
     'background:rgba(17,20,24,.92);color:#E7EAED;font:12px sans-serif">'
     + ('<span style="color:#FDE047">●</span> selected bond'
        if highlight.get("kind") == "bond" else
        '<span style="color:#22D3EE">◎</span> selected symmetry orbit'
        if highlight.get("kind") == "orbit" else
        '<span style="color:#FDE047">━</span> angle arms &nbsp; '
        '<span style="color:#FB7185">●</span> vertex &nbsp; '
        '<span style="color:#22D3EE">⌒</span> angle') + '</div>')
    if highlight.get("kind") else ''}
</div>
<script src="{_CDN}"></script>
<script>
(function() {{
  function draw() {{
    if (typeof $3Dmol === "undefined") {{ setTimeout(draw, 50); return; }}
    var el = document.getElementById("{div_id}");
    if (!el) return;
    var viewer = $3Dmol.createViewer(el, {{ backgroundColor: {bg_js} }});
    var model = viewer.addModel({cif_js}, "cif");
    var colors = {colors_js};
    var atoms = model.selectedAtoms({{}});
    // Remove CIF-inferred bonds. All sticks below come from the tolerance-aware
    // neighbour list, so the slider and the viewer can no longer disagree.
    atoms.forEach(function(atom) {{ atom.bonds = []; atom.bondOrder = []; }});
    Object.keys(colors).forEach(function(elem) {{
      var atomStyle = JSON.parse(JSON.stringify({style_js}));
      Object.keys(atomStyle).forEach(function(part) {{
        atomStyle[part].color = colors[elem];
      }});
      viewer.setStyle({{elem: elem}}, atomStyle);
    }});
    function point(values) {{ return {{x:values[0], y:values[1], z:values[2]}}; }}
    function midpoint(a, b) {{ return {{x:(a.x+b.x)/2, y:(a.y+b.y)/2, z:(a.z+b.z)/2}}; }}
    function cylinder(a, b, radius, color, opacity) {{
      viewer.addCylinder({{start:a, end:b, radius:radius, color:color,
                          opacity:opacity === undefined ? 1 : opacity,
                          fromCap:1, toCap:1}});
    }}
    function selectedAtomLabel(position, text, borderColor) {{
      viewer.addLabel(text, {{position:position, inFront:true,
        backgroundColor:"#020617", backgroundOpacity:.90,
        borderColor:borderColor, borderThickness:1,
        fontColor:"#FFFFFF", fontSize:11, padding:3}});
    }}

    var bondRadius = {bond_radius_js};
    if (bondRadius > 0) {{
      {bonds_js}.forEach(function(bond) {{
        var start = point(bond.points[0]), end = point(bond.points[1]);
        var middle = midpoint(start, end);
        cylinder(start, middle, bondRadius, colors[bond.elements[0]], 0.92);
        cylinder(middle, end, bondRadius, colors[bond.elements[1]], 0.92);
      }});
    }}

    var selected = {highlight_js};
    if (selected.kind && selected.points && selected.points.length) {{
      var selectedPoints = selected.points.map(point);
      selected.point_labels = selected.point_labels || selectedPoints.map(function(_, index) {{
        return "Atom " + (index + 1);
      }});
      if (selected.kind === "bond" && selectedPoints.length === 2) {{
        cylinder(selectedPoints[0], selectedPoints[1], 0.19, "#FDE047", 1);
        selectedPoints.forEach(function(p, index) {{
          selectedAtomLabel(p, selected.point_labels[index], "#FDE047");
        }});
        viewer.addLabel(selected.value_label, {{position:midpoint(selectedPoints[0], selectedPoints[1]),
          inFront:true, backgroundColor:"#0F172A", backgroundOpacity:.94,
          borderColor:"#FDE047", borderThickness:2, fontColor:"#FFFFFF",
          fontSize:15, padding:6}});
      }} else if (selected.kind === "angle" && selectedPoints.length === 3) {{
        var left = selectedPoints[0], vertex = selectedPoints[1], right = selectedPoints[2];
        cylinder(vertex, left, 0.16, "#FDE047", 1);
        cylinder(vertex, right, 0.16, "#FDE047", 1);
        selectedAtomLabel(left, selected.point_labels[0], "#FDE047");
        selectedAtomLabel(vertex, selected.point_labels[1] + " (vertex)", "#FB7185");
        selectedAtomLabel(right, selected.point_labels[2], "#FDE047");

        var v1={{x:left.x-vertex.x,y:left.y-vertex.y,z:left.z-vertex.z}};
        var v2={{x:right.x-vertex.x,y:right.y-vertex.y,z:right.z-vertex.z}};
        var n1=Math.hypot(v1.x,v1.y,v1.z)||1, n2=Math.hypot(v2.x,v2.y,v2.z)||1;
        v1={{x:v1.x/n1,y:v1.y/n1,z:v1.z/n1}};
        v2={{x:v2.x/n2,y:v2.y/n2,z:v2.z/n2}};
        var arcRadius=Math.min(n1,n2)*.38, previous=null, arcMiddle=null;
        for (var step=0; step<=24; step++) {{
          var t=step/24;
          var direction={{x:(1-t)*v1.x+t*v2.x,y:(1-t)*v1.y+t*v2.y,z:(1-t)*v1.z+t*v2.z}};
          var dn=Math.hypot(direction.x,direction.y,direction.z)||1;
          var current={{x:vertex.x+arcRadius*direction.x/dn,
                       y:vertex.y+arcRadius*direction.y/dn,
                       z:vertex.z+arcRadius*direction.z/dn}};
          if (previous) viewer.addLine({{start:previous,end:current,color:"#22D3EE",linewidth:4}});
          if (step===12) arcMiddle=current;
          previous=current;
        }}
        viewer.addLabel(selected.value_label, {{position:arcMiddle || vertex, inFront:true,
          backgroundColor:"#0F172A", backgroundOpacity:.95, borderColor:"#22D3EE",
          borderThickness:2, fontColor:"#FFFFFF", fontSize:15, padding:6}});
      }} else if (selected.kind === "orbit") {{
        selectedPoints.forEach(function(p, index) {{
          viewer.addSphere({{center:p, radius:.52, color:"#22D3EE",
                            opacity:.75, wireframe:true}});
          selectedAtomLabel(p, selected.point_labels[index], "#22D3EE");
        }});
      }}
    }}
    var legend = document.getElementById("{div_id}_legend");
    legend.innerHTML = Object.keys(colors).map(function(elem) {{
      return '<span style="display:flex;align-items:center;gap:5px">' +
        '<i style="width:9px;height:9px;border-radius:50%;background:' +
        colors[elem] + '"></i>' + elem + '</span>';
    }}).join('');
    if ({show_cell_js}) {{
      viewer.addUnitCell(model, {{ box: {{ color: "#94A3B8" }} }});
    }}
    var atomLabels = {labels_js};
    atoms.forEach(function(atom, index) {{ atom.qepLabel = atomLabels[index] || atom.elem; }});
    viewer.setHoverable({{}}, true, function(atom, vw) {{
      if (!atom.hoverLabel) {{
        atom.hoverLabel = vw.addLabel(atom.qepLabel,
          {{ position: atom, backgroundColor: "#1E293B", fontColor: "white",
             fontSize: 11 }});
      }}
    }}, function(atom, vw) {{
      if (atom.hoverLabel) {{ vw.removeLabel(atom.hoverLabel); delete atom.hoverLabel; }}
    }});
    var pinnedLabel = null, pinnedHalo = null;
    viewer.setClickable({{}}, true, function(atom, vw) {{
      if (pinnedLabel) vw.removeLabel(pinnedLabel);
      if (pinnedHalo) vw.removeShape(pinnedHalo);
      pinnedHalo = vw.addSphere({{center:atom, radius:.56, color:"#FFFFFF",
                                  opacity:.58, wireframe:true}});
      pinnedLabel = vw.addLabel(atom.qepLabel + "  •  atom " + (atom.index + 1), {{
        position:atom, inFront:true, backgroundColor:"#020617", backgroundOpacity:.96,
        borderColor:"#FFFFFF", borderThickness:2, fontColor:"#FFFFFF",
        fontSize:14, padding:6
      }});
      vw.render();
    }});
    var cameraKey = {camera_key_js};
    function saveCamera() {{
      try {{ localStorage.setItem(cameraKey, JSON.stringify(viewer.getView())); }} catch (error) {{}}
    }}
    viewer.zoomTo();
    viewer.zoom(1.15);
    try {{
      var savedCamera = localStorage.getItem(cameraKey);
      if (savedCamera) viewer.setView(JSON.parse(savedCamera));
    }} catch (error) {{}}

    var orientations = {{
      top:[0.7071068,0,0,0.7071068], front:[0,0,0,1],
      right:[0,0.7071068,0,0.7071068], left:[0,-0.7071068,0,0.7071068]
    }};
    document.querySelectorAll("#{div_id}_views button").forEach(function(button) {{
      button.style.cssText = "border:1px solid #4B5561;border-radius:3px;background:#20252B;" +
        "color:#E7EAED;padding:4px 7px;cursor:pointer";
      button.addEventListener("click", function() {{
        viewer.zoomTo(); viewer.zoom(1.15);
        var view = viewer.getView(), quaternion = orientations[button.dataset.view];
        view[4]=quaternion[0]; view[5]=quaternion[1];
        view[6]=quaternion[2]; view[7]=quaternion[3];
        viewer.setView(view); viewer.render(); saveCamera();
      }});
    }});
    el.addEventListener("pointerup", function() {{ setTimeout(saveCamera, 0); }});
    el.addEventListener("wheel", function() {{ setTimeout(saveCamera, 80); }}, {{passive:true}});
    viewer.render();
  }}
  draw();
}})();
</script>
"""
