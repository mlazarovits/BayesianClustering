import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import json
import argparse
import os


class JsonPlotter:
	def __init__(self, jsonfile):
		#read json file
		with open(jsonfile,'r') as file:
			self.json_obj = json.load(file)
		self.jsonfile = jsonfile

	def buildColorDict(self):
		colors = {}
		#leaf color - hollow points
		colors[-1] = 'rgba(0, 0, 0, 0.1)'
		#light green
		colors[0] = 'rgba(173, 226, 209, 1.)'
		#light red
		colors[1] = 'rgba(248, 147, 186, 1.)'
		#light blue
		colors[2] = 'rgba(132, 215, 227, 1.)'
		#light purple
		colors[3] = 'rgba(156, 180, 252, 1.)' 
		#light orange
		colors[4] = 'rgba(255, 205, 162, 1.)'
		#light pink
		colors[5] = 'rgba(221, 178, 214, 1.)'
		#dark green
		colors[6] = 'rgba(51, 99, 84, 1.)'
		#dark red
		colors[7] = 'rgba(135, 61, 89, 1.)'
		#dark blue
		colors[8] = 'rgba(39, 129, 143, 1.)'
		#dark purple
		colors[9] = 'rgba(39, 55, 105, 1.)' 
		#dark orange
		colors[10] = 'rgba(156, 92, 36, 1.)'
		#dark pink
		colors[11] = 'rgba(138, 48, 123, 1.)'
		#dark magenta
		colors[12] = 'rgba(136, 48, 138, 1.)'
		#light magenta
		colors[13] = 'rgba(188, 119, 189, 1.)'
		#dark yellow
		colors[14] = 'rgba(148, 129, 34, 1.)'
		#light yellow-green
		colors[15] = 'rgba(187, 201, 123, 1.)'
		#light blue-green
		colors[16] = 'rgba(123, 201, 174, 1.)'
		#dark blue-green
		colors[17] = 'rgba(33, 99, 76, 1.)'
		return colors
	
	
	
	
	def plot_json(self, dataonly = False, numlevels = 0):
		levels = self.json_obj['levels']
			
		plotname = self.jsonfile[:self.jsonfile.find(".json")]
		
		figs = []
		if numlevels > 0:
			nLevels = numlevels
			print("plotting",nLevels,"levels")
		else:
			nLevels = len(levels)
			print("json has",nLevels,"levels")
		for l in range(nLevels):
			fig = self.plot_level(l, plotname+'_level'+str(l))
			figs.append(fig)
		return figs
			
	
	
	def plot_level(self, l, filename):
		#get level l for each tree -> will return clusters_level_l for each tree
		#if tree has less levels than l, plot all data in tree as leaves
		level = self.json_obj["levels"]["level_"+str(l)]
		nTrees = len(level)
	
		gr_arr = []
		minLevel = 0
		for t in range(nTrees):
			gr_arr.append(self.plot_tree(level, t))	
	
		
		#make sure arr is flat
		gr_arr = [gr for i in gr_arr for gr in i]
		#this level is one plot
		fig = go.Figure(gr_arr)
		fig.update_layout({"scene": {"aspectmode": "auto"}},title=filename, template=None)
		return fig
	
		
	
	def plot_tree(self, level, t):
		tree = level["tree_"+str(t)]
		nClusters = len(tree)
		#print("Tree",t,"has",nClusters,"clusters")
		gr_arr = []
		for c in range(nClusters):
			cluster = tree["cluster_"+str(c)]
			gr_arr.append(self.plot_cluster(cluster, t))
		
		#make sure arr is flat
		gr_arr = [gr for i in gr_arr for gr in i]
		return gr_arr
	
	
	def plot_cluster(self, cluster, c, dataonly = False):
		#update axis names
		data = cluster['data']
		x = data['x']
		y = data['y']
		z = data['z']
		
		gr_arr = []
	
		colors = self.buildColorDict()
		
		#leaf
		if(len(x) == 1):
			cl = -1
		else:
			cl = c % 16	
	
		#add data
		gr_arr.append(go.Scatter3d(x=x,y=y,z=z,mode='markers',marker=dict(
				size = 4, color = colors[cl], line=dict(
					color = colors[cl], width = 30)), showlegend=False))
		if dataonly is True:
			return gr_arr
		
		#don't plot subclusters for clusters with points less than minPoints
		minPoints = 2
		if(len(x) < minPoints):
			return gr_arr
		
		nSubClusters = len(cluster['subclusters'])
		for i in range(nSubClusters):
			idx = str(i)
			subcluster = cluster['subclusters']['subcluster_'+idx]
			#should be subcluster_i	
			
			a = subcluster['eigenVal_0']
			b = subcluster['eigenVal_1']
			c = subcluster['eigenVal_2']
			
			op = subcluster['mixing_coeff_norm']
			
			# compute ellipsoid coordinates on standard basis
			# e1=(1, 0, 0), e2=(0, 1, 0), e3=(0, 0, 1)
			u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
			x0 = subcluster['mu_x'] 
			y0 = subcluster['mu_y']
			z0 = subcluster['mu_z']
			
			x1 = a * np.cos(u) * np.sin(v) 
			y1 = b * np.sin(u) * np.sin(v) 
			z1 = c * np.cos(v)
			# points on the ellipsoid
			points = np.stack([t.flatten() for t in [x1, y1, z1]])
			
			
			#eigenvectors
			v1 = subcluster['eigenVec_0']
			v2 = subcluster['eigenVec_1']
			v3 = subcluster['eigenVec_2']
			# 3x3 transformation matrix
			T = np.array([v1, v2, v3])
			
			# transform coordinates to the new orthonormal basis
			new_points = T @ points
			x2 = new_points[0, :] + x0
			y2 = new_points[1, :] + y0
			z2 = new_points[2, :] + z0
			x2, y2, z2 = [t.reshape(x1.shape) for t in [x2, y2, z2]]
		
			#add ellipsoids
			gr_arr.append(go.Surface(x=x2, y=y2, z=z2, opacity=op, colorscale=[colors[cl],colors[cl]], showscale = False)),
		
		return gr_arr

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--json','-j',help='json file to plot')
	parser.add_argument('--data',help='plot data only',action='store_true')
	parser.add_argument('--nlevels',help='number of levels to plot (from top)',default=0)
	args = parser.parse_args()
	
	if args.json is None:
		print("Error: please provide either directory with jsons or single json file.")
		exit()
	
	if ".json" not in args.json:
		print("Error: please provide either directory with jsons or single json file.")
		exit()

	f = args.json	
	jp = JsonPlotter(f)	
	figs = jp.plot_json(args.data, int(args.nlevels))
	name = f[:f.find(".json")]
	print("Writing to directory",name)
	if(os.path.exists(name)):
		#remake files
		os.system("rm -rf "+name)	
	os.mkdir(name)
	for f, fig in enumerate(figs):
		fig.write_image(name+"/level_"+str(f)+".pdf")
		if f < 10:	
			fig.show()

if __name__ == "__main__":
	main()
