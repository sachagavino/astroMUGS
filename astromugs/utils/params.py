from dataclasses import dataclass, fields, is_dataclass, field
import html
from typing import Literal, Any


# =======================================================
# Utility: scientific number formatting
# =======================================================
def fmt_value(v, sci_threshold=1e3):
    if isinstance(v, (float, int)):
        if v != 0 and (abs(v) >= sci_threshold or abs(v) < 1e10):
            return f"{v:.3e}"
    return repr(v)


# =======================================================
# HTML REPR FOR JUPYTER (collapsible panel)
# =======================================================

def _html_repr(obj):
    """Return a collapsible HTML panel with MathJax-enabled table."""
    if not is_dataclass(obj):
        return None

    cls_name = obj.__class__.__name__
    flds = fields(obj)

    rows = ""
    for f in flds:
        desc_raw = f.metadata.get("desc", "")
        # HTML escape text except math (MathJax handles $...$)
        desc = html.escape(desc_raw)
        desc = desc.replace("&dollar;", "$")  # Allow $ in math mode

        val = fmt_value(getattr(obj, f.name))
        typ = f.type.__name__ if hasattr(f.type, "__name__") else str(f.type)

        rows += f"""
        <tr>
            <td><b>{f.name}</b></td>
            <td>{typ}</td>
            <td><code>{val}</code></td>
            <td style="text-align:left; white-space:normal;">{desc}</td>
        </tr>
        """

    # note: MathJax automatically renders inline $...$ or block $$...$$
    html_block = f"""
    <details style="margin:8px 0; padding:6px; border:1px solid #ccc; border-radius:6px;">
      <summary style="font-size:16px; font-weight:bold; cursor:pointer;">{cls_name}</summary>

      <div style="padding:10px;">
        <table style="border-collapse: collapse; margin-top:10px;">
          <thead>
            <tr style="text-align:left; border-bottom:1px solid #888;">
              <th style="padding-right:20px;">Field</th>
              <th style="padding-right:20px;">Type</th>
              <th style="padding-right:20px;">Value</th>
              <th style="text-align:left;">Description</th>
            </tr>
          </thead>
          <tbody>
            {rows}
          </tbody>
        </table>
      </div>

    </details>
    """

    return html_block

# =======================================================
# TEXT FALLBACK (for terminal printing)
# =======================================================
def fancy_repr(obj) -> str:
    """Non-HTML pretty repr for terminal."""
    if not is_dataclass(obj):
        return repr(obj)

    cls_name = obj.__class__.__name__
    flds = fields(obj)

    max_name = max(len(f.name) for f in flds)
    max_type = max(len(f.type.__name__) if hasattr(f.type, "__name__") else len(str(f.type))
                   for f in flds)
    max_desc = max(len(f.metadata.get("desc", "")) for f in flds)

    VALUE_COL_WIDTH = 14

    lines = [f"{cls_name}:"]
    lines.append(
        f"{'Field'.ljust(max_name)}  {'Type'.ljust(max_type)}  "
        f"{'Value'.ljust(VALUE_COL_WIDTH)}  {'Description'.ljust(max_desc)}"
    )
    lines.append(
        f"{'-'*max_name}  {'-'*max_type}  {'-'*VALUE_COL_WIDTH}  {'-'*max_desc}"
    )

    for f in flds:
        name = f.name.ljust(max_name)
        typ = (
            f.type.__name__.ljust(max_type)
            if hasattr(f.type, "__name__")
            else str(f.type).ljust(max_type)
        )
        desc = f.metadata.get("desc", "").ljust(max_desc)
        value = fmt_value(getattr(obj, f.name)).ljust(VALUE_COL_WIDTH)
        lines.append(f"{name}  {typ}  {value}  {desc}")

    return "\n".join(lines)

# =======================================================
# PARAMETER DATACLASSES
# =======================================================

@dataclass
class EnvelopeParams:
    acc_rate: float = field(default=1.000e-05, 
        metadata={'desc': r'[Msun yr-1] Accretion rate '})
    cav_fact: float = field(default=1.000e+01, 
        metadata={'desc': r'Factor of decrease of the density in the cavity'})
    cavpl: float = field(default=1.000e+00, 
        metadata={'desc': r'Opening angle and shape of the outflow'})
    cavz0: float = field(default=1.000e+01, 
        metadata={'desc': r'[AU]'})
    coordsystem: Literal['spherical', 'cylindrical'] = field(default='spherical', 
        metadata={'desc': 'Coordinate system for radiative transfer'})
    dtogas: float = field(default=1.000e-02, 
        metadata={'desc': r'Dust to gas mass ratio'})
    dust_env_mass: float = field(default=1.000e-05, 
        metadata={'desc': r'[Msun] Total dust mass of the envelope'})
    r_centri: float = field(default=1.000e+02, 
        metadata={'desc': r'[AU] Centrifugal radius (critical radius inside of which the envelope flattens)'})
    rmax: float = field(default=5.000e+03, 
        metadata={'desc': r'[AU] Outer radius of the envelope)'})
    rmin: float = field(default=1.000e+00, 
        metadata={'desc': r'[AU] inner radius of the envelope)'})

    def __repr__(self):  # terminal
        return fancy_repr(self)

    def _repr_html_(self):  # Jupyter
        return _html_repr(self)

@dataclass
class DiskParams:
    acc_rate: float = field(default=1.000e-05, 
        metadata={'desc': r'[Msun yr-1] Accretion rate '})
    alpha: float = field(default=1.000e-02, 
        metadata={'desc': 'Viscosity coefficient'})
    coordsystem: Literal['spherical', 'cylindrical'] = field(default='spherical', 
        metadata={'desc': 'Coordinate system for radiative transfer'})
    d_exp: float = field(default=3.500e+00, 
        metadata={'desc': "Size distribution exponent. WARNING: don't type 4"})
    disk_mass: float = field(default=0.15, 
        metadata={'desc': r'[Msun] Total mass of the disk'})
    dtogas: float = field(default=1.000e-02, 
        metadata={'desc': "Dust to gas mass ratio"})
    dust_mass: float = field(default=1.000e-04, 
        metadata={'desc': r"[Msun] Total dust mass in the disk"})
    h0: float = field(default=1.000e+01, 
        metadata={'desc': "[AU] Scale height at reference"})
    isothermal: bool = field(default=False, 
        metadata={'desc': 'Is vertically isothermal. Default is non-isothermal'})
    lim_h: float = field(default=1.000e+00, 
        metadata={'desc': "Value in fraction of scale height until which accretion heating is computed"})
    max_H: float = field(default=4.000e+00, 
        metadata={'desc': "Maximum altitude in fraction of scale height e.g. n=4 means that the largest computed z is 4*Hg"})
    nz_chem: int = field(default=5.000e+01, 
        metadata={'desc': "Number of vertical spatial points for chemistry"})
    p_exp: float = field(default=1.500e+00, 
        metadata={'desc': "Surface density exponent"})
    q_c: float = field(default=6.000e+00, 
        metadata={'desc': "Extinction efficiency at resonance"})
    q_exp: float = field(default=4.000e-01, 
        metadata={'desc': "Exponent for the radial variation of temperature and gas scale height"})
    ref_radius: float = field(default=1.000e+02, 
        metadata={'desc': "[AU] Reference radius for parametric laws"}) 
    rho_m: float = field(default=3.000e+00, 
        metadata={'desc': "[g.cm-3] Material density of the grains"}) 
    rin: float = field(default=1.000e+00,
        metadata={'desc': "[AU] inner radius"})
    rout: float = field(default=3.000e+02, 
        metadata={'desc': "[AU] Outer radius"})
    schmidtnumber: float = field(default=1.000e+00, 
        metadata={'desc': "Schmidt number"})
    settfact: float = field(default=9.000e-01,
        metadata={'desc': "Settling factor."})
    settling: bool = field(default=True,
        metadata={'desc': 'Is vertical settling. Default is True'})       
    sigma_gas_ref: float = field(default=3.350e-01, 
        metadata={'desc': "[g.cm-2] Surface density of the gas at reference radius"})
    sigma_t: float = field(default=2.000e+00, 
        metadata={'desc': "Stiffness of the vertical temperature profile"})
    star_mass: float = field(default=5.799e-01, 
        metadata={'desc': r"[Msun] Mass of the central star"})
    tatmos_ref: float = field(default=5.000e+01,
        metadata={'desc': "[K] Atmospheric temperature at the reference radius (where z=n*H)"})
    tmidplan_ref: float = field(default=1.000e+01, 
        metadata={'desc': "[K] Midplane temperature at the reference radius"})

    def __repr__(self):
        # add a blank line before DiskParams for readability
        return "\n" + fancy_repr(self)

    def _repr_html_(self):
        return "<br/>" + _html_repr(self)   # blank line in Jupyter
    

# ---------------- Structure Params ----------------
@dataclass
class StructureParams:
    envelope: EnvelopeParams = field(default_factory=EnvelopeParams)
    disk: DiskParams = field(default_factory=DiskParams)

    def __repr__(self):
        # blank line between dataclasses
        return fancy_repr(self.envelope) + "\n\n" + fancy_repr(self.disk)

    def _repr_html_(self):
        # Jupyter collapsible panels with spacing
        return _html_repr(self.envelope) + "<br/><br/>" + _html_repr(self.disk)