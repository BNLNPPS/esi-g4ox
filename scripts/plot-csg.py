#!/usr/bin/env python3

import argparse
import logging
import pathlib as pl
import plotly.graph_objects as go
import numpy as np

from g4ox import csg

logging.basicConfig(level=logging.INFO)


def read_mesh(path_csg, path_rec, path_out):
    tri = np.load(path_csg/'tri.npy')
    vtx = np.load(path_csg/'vtx.npy')

    rays = []
    if path_rec:
        rec = np.load(path_rec/'record.npy')
        rays = csg.rays(rec)

    faces, edges = csg.mesh(tri, vtx)


    logging.info(path_csg.name)

    fig = go.Figure(data=[faces, edges, *rays])
    html_file = path_out / (str(path_csg).replace("/", "_").replace("\\", "_").replace(".", "") + '.html')
    fig.write_html(html_file, include_plotlyjs='cdn')

    logging.info(f'Saved image to: {html_file}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('path', help="A path to CSG dir with tri.npy and vtx.npy files")
    parser.add_argument('-r', '--rays', help="A path to record.npy file")
    parser.add_argument('-o', '--outdir', default='/home/web/g4ox', help="A path to output dir")

    args = parser.parse_args()

    path_csg = pl.Path(args.path)
    path_rec = pl.Path(args.rays) if args.rays else None
    path_out = pl.Path(args.outdir)
    logging.info(f'Reading data from dir: {path_csg}')

    read_mesh(path_csg, path_rec, path_out)
