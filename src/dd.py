#!/usr/bin/python

import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import plotly.plotly as py
import numpy as np
import sys
import base64

import gzip
import StringIO

STEP=2
YMAX=1.1
MAXSNP=100000	# max snp allowed, otherwize thin
MAXFSIZE=10000000
F_fai='/media/Setaria_1/Share/Sviridis_vphytozome1.1/Sviridis_311_v1.0.fa.fai'

# genome size dic
f=open(F_fai)
chcolordic={0:'rgba(200,255,0,.2)',1:'rgba(50,150,0,.15)'}
startdic={}
pos,n,chs=0,0,[]
for l in f:
        if l.startswith('scaffold'):
		continue
	ls=l.split()
	chs.append(ls[0])
        startdic[ls[0]]=(pos,pos+int(ls[1]),chcolordic[n%2])
        pos+=int(ls[1])
	n+=1
f.close()


# data entry from vcf.gz file
def readvcf(infile, LOW):
	d=[]
	ch,n,thinning=0,0,1
	if len(infile)>MAXSNP:
		thinning=len(infile)/MAXSNP
	for l in infile:
		if l.startswith("#"):
			continue
		if l.startswith('scaffold'):
			continue
		ls=l.split()
		if len(ls[4])!=1:   # multiple alternative
			continue
		if "AD" not in ls[8]:
			sys.exit(0) # no AD term, could not be analyzed
		ch,relpos=ls[0],int(ls[1])
		abspos=startdic[ch][0]+relpos
		wt,mut=map(int,ls[9].split(':')[1].split(','))
		frequency=mut/(wt+mut+0.0000001)
		if frequency<LOW:
			continue
		disruptive=0+('MODERATE' in ls[7] or 'HIGH' in ls[7])
		if n%thinning!=0:
			n+=1
			continue
		for item in ls[7].split(';'):
			if 'EFF' in item:
				eff=item.split('=')[1]
		annotation=str(ch)+":"+str(relpos)+'<br>'+str(mut)+'/'+str(mut+wt)+' mut/total allels<br>'+'<br>'.join(eff.split(','))
		d.append([ch,relpos,relpos,abspos,frequency,disruptive,annotation,wt,mut])
		n+=1
	d=map(list,zip(*d))
	colordic={0:'rgba(0,70,190,.7)',1:'rgba(255,30,30,1)'}
	sizedic={0:7,1:10}
	snpcolors=[colordic[i%2] for i in d[5]]
	snpsizes=[sizedic[i] for i in d[5]]
	return d, snpcolors, snpsizes

def smoother(d,WIN,STEP):
	n,w_pre=0,0
	smoothpos,smoothfreq,w_depth=[],[],[]
	for n in xrange(len(d[0])): 
		if n-w_pre>=STEP and n>=WIN:
			w_pos=sum(d[3][n-WIN:n])/WIN
			w_wt=sum(d[7][n-WIN:n])
			w_mut=sum(d[8][n-WIN:n])
			w_depth.append(w_wt+w_mut)
			smoothpos.append(w_pos)
			smoothfreq.append(w_mut/(w_wt+w_mut+0.0000001))
			w_pre=n
		n+=1
	return smoothpos, smoothfreq


#figure2=go.Figure(
#	data=[go.Histogram(x=w_depth,histnorm='probability')],
#	layout=go.Layout(title='Histogram for smoothing window read coverage',
#		xaxis=dict(title='Window read coverage')
#		)
#	)


### call server
app = dash.Dash()

# page layout
title_text = '''
### BSA plot for *Setaria viridis*

Bulked Segregant Analysis by deep sequencing.

Please cite [Huang et al. *Sparse panicle1* is required for inflorescence development in *Setaria viridis* and maize. Nature plants 3.5 (2017): 17054.](https://www.nature.com/articles/nplants201754)
'''

app.layout = html.Div(children=[
	dcc.Markdown(children=title_text),
	dcc.Upload(id='file_upload',children=html.Button('Upload .vcf File')),
	dcc.Graph(	id='BSA plot',
			#figure=go.Figure(   
			#	data=[go.Scatter(x=[1,2,3],y=[4,5,6])], # default figure
			#	layout=go.Layout(title='sample'))
		),
	html.Div([
			html.Label('Low frequency cutoff'),
			dcc.Slider(id='freq-slider',min=0,max=100,value=25,step=None,marks={str(WIN): str(WIN) for WIN in range(0,100,5)}),
			dcc.Markdown('_'),
			html.Label('Sliding window size'),
			dcc.Slider(id='win-slider',min=5,max=50,value=10,step=None,marks={str(WIN): str(WIN) for WIN in range(5,50)}),
		],style={'columnCount': 1}),
	#dcc.Graph(figure=figure2,id='test hist'),
])


@app.callback(
	dash.dependencies.Output('BSA plot', 'figure'),
	[
		dash.dependencies.Input('win-slider', 'value'),
		dash.dependencies.Input('freq-slider', 'value'),
		dash.dependencies.Input('file_upload', 'contents'),
		dash.dependencies.Input('file_upload', 'filename'),
	])
def update_figure(WIN,LOW,contents,fname):
	# file reading
	if fname.endswith('.gz'):
		tmp=base64.b64decode(contents.split(',')[1])
		if len(tmp)>MAXFSIZE:
			return go.Figure(layout=go.Layout(title='Warning, too large file for analysis, consider thinning data'))
		infile=gzip.GzipFile(fileobj=StringIO.StringIO(tmp)).read().splitlines()
	else:
		infile=base64.b64decode(contents.split(',')[1]).splitlines()
	if len(infile)>MAXSNP:
		return go.Figure(layout=go.Layout(title='Warning, too large file for analysis, consider thinning data'))
	
	
	# base snp plot
	d, snpcolors, snpsizes = readvcf(infile, LOW/100.)

	# calculate smoothcurve
	smoothpos, smoothfreq = smoother(d,WIN,STEP)
	
	return go.Figure(
		data=[
			go.Scatter(x=d[3], y=d[4], name=u'SNPs',mode = 'markers', marker=dict(color=snpcolors,size=snpsizes),text=d[6]),
			go.Scatter(x=smoothpos, y=smoothfreq, name=u'Smooth line',mode='line', line=dict(color='rgb(255,50,50)'),hoverinfo='none'),
		],
		layout=go.Layout(title='BSA plot for "'+fname+'"',
			yaxis=dict(title='Derived Allele Frequency', range=(0,YMAX)),
			xaxis=dict(title='Genomic position'),
			shapes=[dict(type='rect',x0=startdic[ch][0],x1=startdic[ch][1],y0=0,y1=YMAX,fillcolor=startdic[ch][2], \
				line=dict(width=0)) for ch in chs]),
		)
		




if __name__ == '__main__':
	app.run_server(debug=True, host='0.0.0.0')











