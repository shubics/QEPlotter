"""Setyawan--Curtarolo (AFLOW) k-path recipes for all 14 Bravais lattices.

Coordinates are fractional in the reciprocal basis of the AFLOW-canonical
primitive cell.  Formulae follow Comput. Mater. Sci. 49, 299--312 (2010).
This module contains data/formulae only; cell identification and BZ geometry
live in separate modules.
"""
import numpy as np


G = "GAMMA"


def _segments(*branches):
    result = []
    for branch in branches:
        labels = [G if label == "G" else label for label in branch.split("-")]
        result.extend(zip(labels, labels[1:]))
    return list(result)


def _recipe(variant, labels, coords, branches, parameters=None):
    labels = [G if label == "G" else label for label in labels]
    return {"variant": variant,
            "points": {label: np.asarray(coord, dtype=float)
                       for label, coord in zip(labels, coords)},
            "path": _segments(*branches),
            "parameters": parameters or {}}


def _common(name):
    data = {
        "CUB": (["G","M","R","X"], [[0,0,0],[.5,.5,0],[.5,.5,.5],[0,.5,0]],
                ["G-X-M-G-R-X", "M-R"]),
        "FCC": (["G","K","L","U","W","X"],
                [[0,0,0],[3/8,3/8,3/4],[.5,.5,.5],[5/8,.25,5/8],[.5,.25,.75],[.5,0,.5]],
                ["G-X-W-K-G-L-U-W-L-K", "U-X"]),
        "BCC": (["G","H","P","N"], [[0,0,0],[.5,-.5,.5],[.25,.25,.25],[0,0,.5]],
                ["G-H-N-G-P-H", "P-N"]),
        "TET": (["G","A","M","R","X","Z"],
                [[0,0,0],[.5,.5,.5],[.5,.5,0],[0,.5,.5],[0,.5,0],[0,0,.5]],
                ["G-X-M-G-Z-R-A-Z", "X-R", "M-A"]),
        "ORC": (["G","R","S","T","U","X","Y","Z"],
                [[0,0,0],[.5,.5,.5],[.5,.5,0],[0,.5,.5],[.5,0,.5],[.5,0,0],[0,.5,0],[0,0,.5]],
                ["G-X-S-Y-G-Z-U-R-T-Z", "Y-T", "U-X", "S-R"]),
        "HEX": (["G","M","K","A","L","H"],
                [[0,0,0],[.5,0,0],[1/3,1/3,0],[0,0,.5],[.5,0,.5],[1/3,1/3,.5]],
                ["G-M-K-G-A-L-H-A", "L-M", "K-H"]),
    }
    labels, coords, branches = data[name]
    return _recipe(name, labels, coords, branches)


def get_sc_recipe(lattice):
    """Return recipe for an AFLOW-canonical ASE Bravais lattice object."""
    name, variant, p = lattice.name, lattice.variant, lattice.vars()
    if name in {"CUB", "FCC", "BCC", "TET", "ORC", "HEX"}:
        return _common(name)
    if name == "BCT":
        return _bct(variant, p["a"], p["c"])
    if name == "ORCF":
        return _orcf(variant, p["a"], p["b"], p["c"])
    if name == "ORCI":
        return _orci(p["a"], p["b"], p["c"])
    if name == "ORCC":
        return _orcc(p["a"], p["b"], p["c"])
    if name == "RHL":
        return _rhl(variant, p["alpha"])
    if name == "MCL":
        return _mcl(p["a"], p["b"], p["c"], p["alpha"])
    if name == "MCLC":
        return _mclc(variant, p["a"], p["b"], p["c"], p["alpha"])
    if name == "TRI":
        return _tri(variant)
    raise ValueError(f"No Setyawan–Curtarolo recipe for lattice {name}")


def _bct(variant, a, c):
    a2, c2 = a*a, c*c
    if variant == "BCT1":
        eta = .25 * (1 + c2/a2)
        return _recipe(variant, ["G","M","N","P","X","Z","Z1"],
                       [[0,0,0],[-.5,.5,.5],[0,.5,0],[.25,.25,.25],[0,0,.5],
                        [eta,eta,-eta],[-eta,1-eta,eta]],
                       ["G-X-M-G-Z-P-N-Z1-M", "X-P"], {"eta":eta})
    eta = .25 * (1 + a2/c2)
    zeta = .5 * a2/c2
    return _recipe(variant, ["G","N","P","S","S1","X","Y","Y1","Z"],
                   [[0,0,0],[0,.5,0],[.25,.25,.25],[-eta,eta,eta],
                    [eta,1-eta,-eta],[0,0,.5],[-zeta,zeta,.5],
                    [.5,.5,-zeta],[.5,.5,-.5]],
                   ["G-X-Y-S-G-Z-S1-N-P-Y1-Z", "X-P"],
                   {"eta":eta, "zeta":zeta})


def _orcf(variant, a, b, c):
    a2,b2,c2=a*a,b*b,c*c
    xm=.25*(1+a2/b2-a2/c2); xp=.25*(1+a2/b2+a2/c2)
    if variant in {"ORCF1","ORCF3"}:
        zeta,eta=xm,xp
        coords=[[0,0,0],[.5,.5+zeta,zeta],[.5,.5-zeta,1-zeta],[.5,.5,.5],
                [1,.5,.5],[0,eta,eta],[1,1-eta,1-eta],[.5,0,.5],[.5,.5,0]]
        branches = (["G-Y-T-Z-G-X-A1-Y", "T-X1", "X-A-Z", "L-G"]
                    if variant=="ORCF1" else ["G-Y-T-Z-G-X-A1-Y", "X-A-Z", "L-G"])
        return _recipe(variant,["G","A","A1","L","T","X","X1","Y","Z"],
                       coords,branches,{"eta":eta,"zeta":zeta})
    phi=.25*(1+c2/b2-c2/a2); delta=.25*(1+b2/a2-b2/c2); eta=xm
    coords=[[0,0,0],[.5,.5-eta,1-eta],[.5,.5+eta,eta],[.5-delta,.5,1-delta],
            [.5+delta,.5,delta],[.5,.5,.5],[1-phi,.5-phi,.5],[phi,.5+phi,.5],
            [0,.5,.5],[.5,0,.5],[.5,.5,0]]
    return _recipe(variant,["G","C","C1","D","D1","L","H","H1","X","Y","Z"],
                   coords,["G-Y-C-D-X-G-Z-D1-H-C", "C1-Z", "X-H1", "H-Y", "L-G"],
                   {"eta":eta,"phi":phi,"delta":delta})


def _orci(a,b,c):
    a2,b2,c2=a*a,b*b,c*c
    zeta=.25*(1+a2/c2); eta=.25*(1+b2/c2)
    delta=.25*(b2-a2)/c2; mu=.25*(a2+b2)/c2
    coords=[[0,0,0],[-mu,mu,.5-delta],[mu,-mu,.5+delta],[.5-delta,.5+delta,-mu],
            [0,.5,0],[.5,0,0],[0,0,.5],[.25,.25,.25],[-zeta,zeta,zeta],
            [zeta,1-zeta,-zeta],[eta,-eta,eta],[1-eta,eta,-eta],[.5,.5,-.5]]
    return _recipe("ORCI",["G","L","L1","L2","R","S","T","W","X","X1","Y","Y1","Z"],
                   coords,["G-X-L-T-W-R-X1-Z-G-Y-S-W", "L1-Y", "Y1-Z"],
                   {"zeta":zeta,"eta":eta,"delta":delta,"mu":mu})


def _orcc(a,b,c):
    zeta=.25*(1+a*a/(b*b))
    coords=[[0,0,0],[zeta,zeta,.5],[-zeta,1-zeta,.5],[0,.5,.5],[0,.5,0],
            [-.5,.5,.5],[zeta,zeta,0],[-zeta,1-zeta,0],[-.5,.5,0],[0,0,.5]]
    return _recipe("ORCC",["G","A","A1","R","S","T","X","X1","Y","Z"],coords,
                   ["G-X-S-R-A-Z-G-Y-X1-A1-T-Y", "Z-T"],{"zeta":zeta})


def _rhl(variant, alpha):
    ar=np.radians(alpha)
    if variant=="RHL1":
        ca=np.cos(ar); eta=(1+4*ca)/(2+4*ca); nu=.75-.5*eta
        coords=[[0,0,0],[eta,.5,1-eta],[.5,1-eta,eta-1],[.5,.5,0],[.5,0,0],
                [0,0,-.5],[eta,nu,nu],[1-nu,1-nu,1-eta],[nu,nu,eta-1],
                [1-nu,nu,0],[nu,0,-nu],[.5,.5,.5]]
        return _recipe(variant,["G","B","B1","F","L","L1","P","P1","P2","Q","X","Z"],
                       coords,["G-L-B1", "B-Z-G-X", "Q-F-P1-Z", "L-P"],
                       {"eta":eta,"nu":nu})
    eta=1/(2*np.tan(ar/2)**2); nu=.75-.5*eta
    coords=[[0,0,0],[.5,-.5,0],[.5,0,0],[1-nu,-nu,1-nu],[nu,nu-1,nu-1],
            [eta,eta,eta],[1-eta,-eta,-eta],[.5,-.5,.5]]
    return _recipe(variant,["G","F","L","P","P1","Q","Q1","Z"],coords,
                   ["G-P-Z-Q-G-F-P1-Q1-L-Z"],{"eta":eta,"nu":nu})


def _mcl(a,b,c,alpha):
    ar=np.radians(alpha); ca=np.cos(ar); sa=np.sin(ar)
    eta=(1-b*ca/c)/(2*sa**2); nu=.5-eta*c*ca/b
    coords=[[0,0,0],[.5,.5,0],[0,.5,.5],[.5,0,.5],[.5,0,-.5],[.5,.5,.5],
            [0,eta,1-nu],[0,1-eta,nu],[0,eta,-nu],[.5,eta,1-nu],
            [.5,1-eta,nu],[.5,eta,-nu],[0,.5,0],[0,0,.5],[0,0,-.5],[.5,0,0]]
    return _recipe("MCL",["G","A","C","D","D1","E","H","H1","H2","M","M1","M2","X","Y","Y1","Z"],
                   coords,["G-Y-H-C-E-M1-A-X-H1", "M-D-Z", "Y-D"],
                   {"eta":eta,"nu":nu})


def _mclc(variant,a,b,c,alpha):
    v=int(variant[-1]); ar=np.radians(alpha); ca=np.cos(ar); sa=np.sin(ar); s2=sa*sa
    a2,b2=a*a,b*b
    if v in (1,2):
        zeta=(2-b*ca/c)/(4*s2); eta=.5+2*zeta*c*ca/b
        psi=.75-a2/(4*b2*s2); phi=psi+(.75-psi)*b*ca/c
        coords=[[0,0,0],[.5,0,0],[0,-.5,0],[1-zeta,1-zeta,1-eta],[zeta,zeta,eta],
                [-zeta,-zeta,1-eta],[1-zeta,-zeta,1-eta],[phi,1-phi,.5],
                [1-phi,phi-1,.5],[.5,.5,.5],[.5,0,.5],[1-psi,psi-1,0],
                [psi,1-psi,0],[psi-1,-psi,0],[.5,.5,0],[-.5,-.5,0],[0,0,.5]]
        branches=(["G-Y-F-L-I", "I1-Z-F1", "Y-X1", "X-G-N", "M-G"] if v==1
                  else ["G-Y-F-L-I", "I1-Z-F1", "N-G-M"])
        params={"zeta":zeta,"eta":eta,"psi":psi,"phi":phi}
        labels=["G","N","N1","F","F1","F2","F3","I","I1","L","M","X","X1","X2","Y","Y1","Z"]
    elif v in (3,4):
        mu=.25*(1+b2/a2); delta=b*c*ca/(2*a2)
        zeta=mu-.25+(1-b*ca/c)/(4*s2); eta=.5+2*zeta*c*ca/b
        phi=1+zeta-2*mu; psi=eta-2*delta
        coords=[[0,0,0],[1-phi,1-phi,1-psi],[phi,phi-1,psi],[1-phi,-phi,1-psi],
                [zeta,zeta,eta],[1-zeta,-zeta,1-eta],[-zeta,-zeta,1-eta],
                [.5,-.5,.5],[.5,0,.5],[.5,0,0],[0,-.5,0],[.5,-.5,0],
                [mu,mu,delta],[1-mu,-mu,-delta],[-mu,-mu,-delta],[mu,mu-1,delta],[0,0,.5]]
        branches=(["G-Y-F-H-Z-I-F1", "H1-Y1-X-G-N", "M-G"] if v==3
                  else ["G-Y-F-H-Z-I", "H1-Y1-X-G-N", "M-G"])
        params={"mu":mu,"delta":delta,"zeta":zeta,"eta":eta,"phi":phi,"psi":psi}
        labels=["G","F","F1","F2","H","H1","H2","I","M","N","N1","X","Y","Y1","Y2","Y3","Z"]
    else:
        zeta=.25*(b2/a2+(1-b*ca/c)/s2); eta=.5+2*zeta*c*ca/b
        mu=.5*eta+b2/(4*a2)-b*c*ca/(2*a2); nu=2*mu-zeta
        omega=(4*nu-1-b2*s2/a2)*c/(2*b*ca); delta=zeta*c*ca/b+omega/2-.25
        rho=1-zeta*a2/b2
        coords=[[0,0,0],[nu,nu,omega],[1-nu,1-nu,1-omega],[nu,nu-1,omega],
                [zeta,zeta,eta],[1-zeta,-zeta,1-eta],[-zeta,-zeta,1-eta],
                [rho,1-rho,.5],[1-rho,rho-1,.5],[.5,.5,.5],[.5,0,.5],
                [.5,0,0],[0,-.5,0],[.5,-.5,0],[mu,mu,delta],[1-mu,-mu,-delta],
                [-mu,-mu,-delta],[mu,mu-1,delta],[0,0,.5]]
        branches=["G-Y-F-L-I", "I1-Z-H-F1", "H1-Y1-X-G-N", "M-G"]
        params={"zeta":zeta,"eta":eta,"mu":mu,"nu":nu,"omega":omega,"delta":delta,"rho":rho}
        labels=["G","F","F1","F2","H","H1","H2","I","I1","L","M","N","N1","X","Y","Y1","Y2","Y3","Z"]
    return _recipe(variant,labels,coords,branches,params)


def _tri(variant):
    labels=["G","L","M","N","R","X","Y","Z"]
    if variant in {"TRI1a","TRI2a"}:
        coords=[[0,0,0],[.5,.5,0],[0,.5,.5],[.5,0,.5],[.5,.5,.5],[.5,0,0],[0,.5,0],[0,0,.5]]
    else:
        coords=[[0,0,0],[.5,-.5,0],[0,0,.5],[-.5,-.5,.5],[0,-.5,.5],[0,-.5,0],[.5,0,0],[-.5,0,.5]]
    return _recipe(variant,labels,coords,["X-G-Y", "L-G-Z", "N-G-M", "R-G"])
