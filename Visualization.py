#!/bin/env python

import plotly
import plotly.graph_objs as go
import plotly.express as px

def create_figure(objs, name):
    fig = go.Figure(objs)
    size = 800
    fig.update_layout(
            template="plotly_dark",
            paper_bgcolor="rgb(0,0,0,1)",
            scene = dict(
                xaxis = dict(nticks=10, range=[-300,300],),
                yaxis = dict(nticks=10, range=[-300,300],),
                zaxis = dict(nticks=10, range=[-300,300],),
                aspectratio=dict(x=1, y=1, z=1),
                ),
            width=1000,
            height=1000,
            margin=dict(r=20, l=10, b=10, t=10));
    fig.write_html(name)

def get_modules_goScatter3D(list_of_detids, geom):
    zs=[]
    xs=[]
    ys=[]
    j = geom.geom_data
    for detid in list_of_detids:
        z_=[j[detid][0][0], j[detid][1][0], j[detid][2][0], j[detid][3][0], j[detid][0][0], None]
        x_=[j[detid][0][1], j[detid][1][1], j[detid][2][1], j[detid][3][1], j[detid][0][1], None]
        y_=[j[detid][0][2], j[detid][1][2], j[detid][2][2], j[detid][3][2], j[detid][0][2], None]
        zs += z_
        xs += x_
        ys += y_
    mods = go.Scatter3d(x=xs,
                        y=ys,
                        z=zs,
                        mode="lines",
                        line=dict(
                            color="rgb(0,255,255,0.5)",
                            ),
                        hoverinfo='text',
                        )
    return mods
