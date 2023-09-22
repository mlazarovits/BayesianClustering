import numpy as np
import matplotlib.pyplot as plt
from plotly.colors import sample_colorscale
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
		self.dirname = jsonfile[:jsonfile.find(".json")]
		self.cl_gifs = []
		self.minPoints = 3

	def setVerb(self, v):
		self._v = v

	def SetWraparound(self, w):
		self.wraparound = w

	def makeOpacities(self, w):
		wmin = min(w);
		wmax = max(w);

		op = [(i - wmin)/(wmax - wmin) for i in w]
		return op

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
		#dark blue-green
		colors[11] = 'rgba(33, 99, 76, 1.)'
		#dark pink
		colors[12] = 'rgba(138, 48, 123, 1.)'
		#dark yellow
		colors[13] = 'rgba(148, 129, 34, 1.)'
		#light magenta
		colors[14] = 'rgba(188, 119, 189, 1.)'
		#light yellow-green
		colors[15] = 'rgba(187, 201, 123, 1.)'
		#light blue-green
		colors[16] = 'rgba(123, 201, 174, 1.)'
		#brown
		colors[17] = 'rgba(97, 58, 18, 1.)'
		return colors
	
	
	
	
	def plotDataset(self, numlevels = 0, onlydata = False):
		levels = self.json_obj['levels']
			
		plotname = self.jsonfile[:self.jsonfile.find(".json")]
		
		figs = []
		if numlevels > 0:
			nLevels = numlevels
			print("Plotting",nLevels,"levels")
		else:
			nLevels = len(levels)
			print("json has",nLevels,"levels")
		for l in range(nLevels):
			fig = self.plot_level(l, plotname+'_level'+str(l), onlydata)
			figs.append(fig)
		return figs
			
	
	
	def plot_level(self, l, filename, onlydata = False):
		#get level l for each tree -> will return clusters_level_l for each tree
		#if tree has less levels than l, plot all data in tree as leaves
		level = self.json_obj["levels"]["level_"+str(l)]
		nTrees = len(level)
		if self._v > 0:
			print("level",l,"has",nTrees,"trees")
		gr_arr = []
		minLevel = 0
		for t in range(nTrees):
			gr_arr.append(self.plot_tree(l, t, onlydata))	
	
		
		#make sure arr is flat
		gr_arr = [gr for i in gr_arr for gr in i]
		#this level is one plot
		fig = go.Figure(gr_arr)
		fig.update_layout({"scene": {"aspectmode": "auto"}},title=filename, template=None)
		return fig
	
		
	
	def plot_tree(self, l, t, onlydata = False):
		level = self.json_obj["levels"]["level_"+str(l)]
		tree = level["tree_"+str(t)]
		nClusters = len(tree)
		#print("Tree",t,"has",nClusters,"clusters")
		gr_arr = []
	
		if self._v > 0:
			print("	tree",t,"has",nClusters,"clusters")

		for c in range(nClusters):
			cluster = tree["cluster_"+str(c)]
			#plot whole data
			gr_arr.append(self.plot_cluster(cluster, t, True))
			if onlydata:
				continue	
		
			#plot individual clusters
			gr_cl = self.plot_cluster(cluster,t,False)
			#if npts in data of cluster c is less than minPoints (defined in plot_cluster) subclusters won't be drawn, only want to draw clusters with subclusters so skip these
			if(len(cluster["data"]["x"]) < self.minPoints):
				continue
			#this level is one plot
			fig = go.Figure(gr_cl)
			filename = self.jsonfile[:self.jsonfile.find(".json")]+"_cluster"+str(t)+"_level"+str(l)
			fig.update_layout({"scene": {"aspectmode": "auto"}},title=filename, template=None)
			#write to file
			#make sure cluster directory exists
			if(os.path.exists(self.dirname+"/cluster"+str(t))):
				#remake files
				os.system("rm -rf "+self.dirname+"/cluster"+str(t))	
			os.mkdir(self.dirname+"/cluster"+str(t))
			fig.write_image(self.dirname+"/cluster"+str(t)+"/level_"+str(l)+".pdf")
			if l == 0:
				fig.show()
				gifcmd = "convert -delay 50 -loop 1 -reverse "
				self.cl_gifs.append(gifcmd)
			self.cl_gifs[c] += self.dirname+"/cluster"+str(t)+"/level_"+str(l)+" "
				
		return gr_arr
	
	
	def plot_cluster(self, cluster, c, alldata = False):
		#update axis names
		data = cluster['data']
		x = data['x']
		y = data['y']
		z = data['z']
		w = data['w']
		
		
		nSubClusters = len(cluster['subclusters'])

		if self._v > 0:
			print("		cluster",c,"has",nSubClusters,"subclusters")

		for k in range(nSubClusters):
			if self._v > 0:
				print("			subcluster",k,"has",cluster['subclusters']['subcluster_'+str(k)]['color'],'weight and',len(x),'points')
			#w.append(cluster['subclusters']['subcluster_'+str(k)]["color"])
	
		colors = self.buildColorDict()


		#leaf
		if(len(x) == 1):
			cl = -1
		else:
			cl = c
		if cl > 17:
			cl = cl % 16	
	
		cols = [colors[cl] for i in range(len(w))]
		name = "cluster "+str(c)
		if(len(x) >= self.minPoints or nSubClusters > 1):
			name += " has "+str(len(x))+" points ("+str(round(sum(w),2))+")"
			name += " in "+str(nSubClusters)+" subclusters" 
		
		if(min(w) != max(w)):
			opacities = self.makeOpacities(w)
			cols = [i.replace("1.",str(opacities[idx])) for idx, i in enumerate(cols)]
		

		if alldata is True:
			return go.Scatter3d(x=x,y=y,z=z,mode='markers',marker=dict(
				size = 4, color = cols, line=dict(
					color = colors[cl], width = 30)), showlegend=True, name = name)
		
		#don't plot subclusters for clusters with points less than minPoints
		if(len(x) < self.minPoints):
			return go.Scatter3d(x=x,y=y,z=z,mode='markers',marker=dict(
				size = 4, color = cols, line=dict(
					color = colors[cl], width = 30)), showlegend=True, name = name)

		gr_arr = []
		#transform the data points into local coordinates + have color be energy scale
		if(self.wraparound is True):
			for i in y:
				if i < 0:
					i += 2*np.pi

		avg_x = np.mean(x)
		avg_y = np.mean(y)
		avg_z = np.mean(z)
		for i, pt in enumerate(zip(x,y,z)):
			x[i] = pt[0] - avg_x
			#transform phi coord
			y[i] = pt[1] - avg_y
			z[i] = pt[2] - avg_z
		
		#data in "local"/delta coordinates
		gr_arr.append(go.Scatter3d(x=x,y=y,z=z,mode='markers',marker=dict(
			size = 4, cmax = max(w), cmin = min(w), color = w, colorscale = 'Plotly3', colorbar=dict(title="Energy (GeV)", x=-0.2)),
			showlegend=True, name = name))
		
		for i in range(nSubClusters):
			idx = str(i)
			subcluster = cluster['subclusters']['subcluster_'+idx]
			#should be subcluster_i	
		
	
			a = round(subcluster['eigenVal_0'], 10)
			b = round(subcluster['eigenVal_1'], 10)
			c = round(subcluster['eigenVal_2'], 10)
			
			op = subcluster['mixing_coeff_norm']
			
			# compute ellipsoid coordinates on standard basis
			# e1=(1, 0, 0), e2=(0, 1, 0), e3=(0, 0, 1)
			u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
			x0 = subcluster['mu_x'] - avg_x 
			y0 = subcluster['mu_y'] - avg_y
			z0 = subcluster['mu_z'] - avg_z
		
			x1 = np.sqrt(a) * np.cos(u) * np.sin(v) 
			y1 = np.sqrt(b) * np.sin(u) * np.sin(v) 
			z1 = np.sqrt(c) * np.cos(v)
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

			#make ellipsoid color of average energy across points (responsibilities)
			if max(w) == min(w):
				cl_w = (subcluster["color"])# - min(w))/(max(w) - min(w))
			else:
				cl_w = (subcluster["color"] - min(w))/(max(w) - min(w))
				scale = True
			cl = sample_colorscale("Plotly3",cl_w)
			cl = np.array([cl,cl]).flatten()
		
			ell_name = "subcluster "+str(i)+" with weight "+str(round(subcluster["color"],2)) 
	
			#add ellipsoids
			#gr_arr.append(go.Surface(x=x2, y=y2, z=z2, opacity=op, colorscale=cl, surfacecolor=y1, cmin=y1.min(), cmax=y1.max(), showscale = False, showlegend = False)),
			gr_arr.append(go.Surface(x=x2, y=y2, z=z2, opacity=op, colorscale=cl, showscale = False, name = ell_name,showlegend=True)),
		return gr_arr

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--json','-j',help='json file to plot')
	parser.add_argument('--data',help='plot data only',action='store_true')
	parser.add_argument('--nlevels',help='number of levels to plot (from top)',default=0,type=int)
	parser.add_argument('--verbosity','-v',help='verbosity',default=0,type=int)
	parser.add_argument('--noViz',action='store_true',help='do not make plots',default=False)
	parser.add_argument('--noWrap',action='store_false',help='no (phi) wraparound',default=True)
	args = parser.parse_args()
	
	if args.json is None:
		print("Error: please provide either directory with jsons or single json file.")
		exit()
	
	if ".json" not in args.json:
		print("Error: please provide either directory with jsons or single json file.")
		exit()
	

	f = args.json	
	jp = JsonPlotter(f)
	name = jp.dirname
	if args.data and os.path.exists(name):
		#remake files
		os.system("rm -rf "+name)	
	os.mkdir(name)
	jp.SetWraparound(args.noWrap)	
	jp.setVerb(args.verbosity)
	#draw all data - also plots individual clusters with GMM components
	figs = jp.plotDataset(int(args.nlevels),args.data)
	if args.noViz:
		exit()
	print("Writing to directory",name)
	files = []
	for f, fig in enumerate(figs):
		fig.write_image(name+"/level_"+str(f)+".pdf")
		files.append(name+"/level_"+str(f)+".pdf")
		if f < 11:	
			fig.show()
	gifcmd = "convert -delay 50 -loop 1 -reverse "
	for f in files:
		gifcmd += f+" "


	gifcmd += name+"/total.gif"

	if args.nlevels == 0:
		os.system(gifcmd)
if __name__ == "__main__":
	main()
