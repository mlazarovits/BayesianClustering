import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import json
import argparse
import os


def plot_json(jsonfile, dataonly = False):
	#read json file
	with open(jsonfile,'r') as file:
		json_obj = json.load(file)
	data = json_obj['data']
	clusters = json_obj['clusters']
	nClusters = len(clusters)
	
	plotname = jsonfile[:jsonfile.find(".json")]


	#update axis names
	x = data['x']
	y = data['y']
	z = data['z']
	
	gr_arr = []

	
	#add data
	gr_arr.append(go.Scatter3d(x=x,y=y,z=z,mode='markers',marker=dict(
			size = 4, color = 'rgba(132,242,201,1.)', line=dict(
				color = 'rgba(132, 242, 201, 1.)', width = 30)), showlegend = False))
	if dataonly is True:
		fig = go.Figure(gr_arr)
		fig.update_layout({"scene": {"aspectmode": "auto"}},title=plotname, template=None)
		return fig
	
	
	for i in range(nClusters):
		idx = str(i)
		a = clusters[idx]['eigenVal_0']
		b = clusters[idx]['eigenVal_1']
		c = clusters[idx]['eigenVal_2']
		
		op = clusters[idx]['mixing_coeff_norm']
		
		# compute ellipsoid coordinates on standard basis
		# e1=(1, 0, 0), e2=(0, 1, 0), e3=(0, 0, 1)
		u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:40j]
		x0 = clusters[idx]['mu_x'] 
		y0 = clusters[idx]['mu_y']
		z0 = clusters[idx]['mu_z']
		
		x1 = np.sqrt(a) * np.cos(u) * np.sin(v) 
		y1 = np.sqrt(b) * np.sin(u) * np.sin(v) 
		z1 = np.sqrt(c) * np.cos(v)
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
		gr_arr.append(go.Surface(x=x2, y=y2, z=z2, opacity=op, colorscale=['rgba(132, 242, 201, 1.)','rgba(132, 242, 201, 1.)'], surfacecolor=y1, cmin=y1.min(), cmax=y1.max(), showscale = False, showlegend = False)),

		#draw means
		gr_arr.append(go.Scatter3d(x=[x0],y=[y0],z=[z0],mode='markers',marker=dict(
			size = 4, color = 'rgba(0,0,0,1.)', symbol='x', line=dict(
				color = 'rgba(0, 0, 0, 1.)', width = 30)), showlegend = False))
	
	fig = go.Figure(gr_arr)
	fig.update_layout({"scene": {"aspectmode": "auto"}},title=plotname, template=None)
	return fig


parser = argparse.ArgumentParser()
parser.add_argument('--dir','-d',help='directory with json files')
parser.add_argument('--json','-j',help='json file to plot')
parser.add_argument('--data',help='plot data only',action='store_true')
args = parser.parse_args()

if args.dir is None and args.json is None:
	print("Error: please provide either directory with jsons or single json file.")
	exit()

files = []
outname = ''
if args.dir is not None:
	#sort files by iteration
	it_to_file = {}
	for j in os.listdir(args.dir):
		if ".json" not in j:
			continue
		it = int(j[j.find("it")+2:j.find(".json")])
		it_to_file[it] = j
	it_to_file = dict(sorted(it_to_file.items()))
	for j in it_to_file:
		files.append(it_to_file[j])
	outname += args.dir+'/'

if len(args.dir) < 1:
	exit()

if args.json is not None:
	files.append(args.json)


for f in files:
	if ".json" not in f:
		continue
	fig = plot_json(outname+f,args.data)
	if args.dir is not None:
		name = f[:f.find(".json")]
		print("Writing to",args.dir+"/"+name+".pdf")
		fig.write_image(args.dir+"/"+name+".pdf")
	if args.data:
		break
fig.show()
if args.dir is not None:
	os.system("convert -delay 50 -loop 1 "+args.dir+"/*.pdf "+args.dir+"/total.gif");