#!/usr/bin/env python3

import argparse
import logging
import pathlib as pl
import plotly.graph_objects as go
import numpy as np

from g4ox import csg

logging.basicConfig(level=logging.INFO)


def read_mesh(path_csg, path_out):
    tri = np.load(path_csg/'tri.npy')
    vtx = np.load(path_csg/'vtx.npy')

    faces, edges = csg.mesh(tri, vtx)

    logging.info(path_csg.name)

    fig = go.Figure(data=[faces, edges])
    html_file = path_out / (str(path_csg).replace("/", "_").replace("\\", "_").replace(".", "") + '.html')
    fig.write_html(html_file, include_plotlyjs='cdn')

    logging.info(f'Saved image to: {html_file}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('path', help="A path to CSG dir")
    parser.add_argument('-o', '--outdir', default='/home/web/g4ox', help="A path to output dir")

    args = parser.parse_args()

    path_csg = pl.Path(args.path)
    path_out = pl.Path(args.outdir)
    logging.info(f'Reading data from dir: {path_csg}')

    read_mesh(path_csg, path_out)
