import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import json


#read json file
with open('it7.json','r') as file:
	json_obj = json.load(file)
data = json_obj['data']
clusters = json_obj['clusters']
nClusters = len(clusters)

#update axis names
x = data['x']
y = data['y']
z = data['z']

gr_arr = []

#add data
gr_arr.append(go.Scatter3d(x=x,y=y,z=z,mode='markers',marker=dict(
		size = 4, color = 'rgba(132,242,201,0.0)', line=dict(
			color = 'rgba(132, 242, 201, 1.)', width = 30))))

for i in range(nClusters):
	idx = str(i)
	a = clusters[idx]['eigenVal_0']
	b = clusters[idx]['eigenVal_1']
	c = clusters[idx]['eigenVal_2']
	
	op = clusters[idx]['mixing_coeff_norm']
	
	# compute ellipsoid coordinates on standard basis
	# e1=(1, 0, 0), e2=(0, 1, 0), e3=(0, 0, 1)
	u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
	x0 = clusters[idx]['mu_x'] 
	y0 = clusters[idx]['mu_y']
	z0 = clusters[idx]['mu_z']
	
	x1 = a * np.cos(u) * np.sin(v) 
	y1 = b * np.sin(u) * np.sin(v) 
	z1 = c * np.cos(v)
	# points on the ellipsoid
	points = np.stack([t.flatten() for t in [x1, y1, z1]])
	
	
	#eigenvectors
	v1 = clusters[idx]['eigenVec_0']
	v2 = clusters[idx]['eigenVec_1']
	v3 = clusters[idx]['eigenVec_2']
	# 3x3 transformation matrix
	T = np.array([v1, v2, v3])
	
	# transform coordinates to the new orthonormal basis
	new_points = T @ points
	x2 = new_points[0, :] + x0
	y2 = new_points[1, :] + y0
	z2 = new_points[2, :] + z0
	x2, y2, z2 = [t.reshape(x1.shape) for t in [x2, y2, z2]]

	#add ellipsoids
	gr_arr.append(go.Surface(x=x2, y=y2, z=z2, opacity=op, colorscale="aggrnyl", surfacecolor=y1, cmin=y1.min(), cmax=y1.max(), colorbar=dict(len=0.6, yanchor="bottom", y=0, x=0.95))),


# scale vector for better visualization
#scale = 5
#v1, v2, v3 = [scale * t for t in [v1, v2, v3]]

'''
fig = go.Figure([
    # axis on the standard base
    go.Scatter3d(x=[0, 5], y=[0, 0], z=[0, 0], mode="lines", name="x1", line=dict(width=5, color="red")),
    go.Scatter3d(x=[0, 0], y=[0, 5], z=[0, 0], mode="lines", name="y1", line=dict(width=5, color="green")),
    go.Scatter3d(x=[0, 0], y=[0, 0], z=[0, 5], mode="lines", name="z1", line=dict(width=5, color="blue")),
    # axis on the new orthonormal base
    go.Scatter3d(x=[0, v1[0]], y=[0, v1[1]], z=[0, v1[2]], mode="lines", name="x2", line=dict(width=2, color="red")),
    go.Scatter3d(x=[0, v2[0]], y=[0, v2[1]], z=[0, v2[2]], mode="lines", name="y2", line=dict(width=2, color="green")),
    go.Scatter3d(x=[0, v3[0]], y=[0, v3[1]], z=[0, v3[2]], mode="lines", name="z2", line=dict(width=2, color="blue")),
    # original ellipsoid aligned to the standard base
    #go.Surface(x=x1, y=y1, z=z1, opacity=0.35, colorscale="plotly3", surfacecolor=y1, cmin=y1.min(), cmax=y1.max(), colorbar=dict(len=0.6, yanchor="bottom", y=0)),
    # final ellipsoid aligned to the new orthonormal base
    go.Surface(x=x2, y=y2, z=z2, opacity=0.5, colorscale="aggrnyl", surfacecolor=y1, cmin=y1.min(), cmax=y1.max(), colorbar=dict(len=0.6, yanchor="bottom", y=0, x=0.95)),
    #scatter
    go.Scatter3d(x=x,y=y,z=z,mode='markers',marker=dict(size = 8, color = z, colorscale = 'Viridis'))
    #slices
    #go.Surface(x=x2, y=y2, z=z3, opacity=1, colorscale="agsunset", surfacecolor=y1, cmin=y1.min(), cmax=y1.max(), colorbar=dict(len=0.6, yanchor="bottom", y=0, x=0.95)),
    #go.Surface(x=x2, y=y2, z=z4, opacity=1, colorscale="agsunset", surfacecolor=y1, cmin=y1.min(), cmax=y1.max(), colorbar=dict(len=0.6, yanchor="bottom", y=0, x=0.95))
])
'''

fig = go.Figure(gr_arr)

fig.update_layout({"scene": {"aspectmode": "auto"}})
fig.show()
