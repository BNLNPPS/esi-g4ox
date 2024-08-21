import numpy as np
import plotly.graph_objects as go

from opticks.CSG.CSGFoundry import CSGFoundry

mesh = None

def set_csg(path):
    global mesh
    cf = CSGFoundry.Load(path)
    mesh = cf.sim.stree.mesh


def create_mesh(mesh_subdir):
    vtx = mesh_subdir.vtx
    tri = mesh_subdir.tri

    # close triangles and create edge coordinates
    tri_close = np.hstack((tri, tri[:,0][:,np.newaxis]))
    tri_sides = np.array([[None]*3 if i is None else vtx[i].tolist() for t in tri_close for i in list(t)+[None]], dtype=float)

    faces = go.Mesh3d(x=vtx.T[0], y=vtx.T[1], z=vtx.T[2], i=tri.T[0], j=tri.T[1], k=tri.T[2], color='green', opacity=0.2)
    edges = go.Scatter3d(x=tri_sides.T[0], y=tri_sides.T[1], z=tri_sides.T[2], mode='lines', line_width=5)

    return faces, edges


def fig_show(data):
    fig = go.Figure(data=data)
    fig.show()
    return fig


def fig_add(fig, data):
    fig.add_traces(data)
    return fig
