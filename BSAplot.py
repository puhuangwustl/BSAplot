#!/usr/bin/python

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table_experiments as dt
import plotly.graph_objs as go
import plotly.plotly as py
import numpy as np
import sys
import base64
#import json
import re
import urllib

import commands
import gzip
import StringIO

STEP=2
YMAX=1.1
MAXSNP=100000	# max snp allowed, otherwize thin
MAXFSIZE=10000000
F_fai='/media/Setaria_1/Share/Sviridis_vphytozome1.1/Sviridis_311_v1.0.fa.fai'
F_anno='/media/Setaria_1/Share/Sviridis_vphytozome1.1/Sviridis_311_v1.1.annotation_info.txt'
ColumnOrder=['Pos','Mut/Total','Annotation','Impact','Codon','Allele','Gene','AtLog','AtName','AtAnno','Transcript']
#F_fai='Sviridis_311_v1.0.fa.fai'
#F_fai='Sviridis_311_v1.0.fa.fai'
#F_anno='Sviridis_311_v1.1.annotation_info.txt'

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
	colordic={0:'rgba(0,70,190,.7)',1:'rgba(255,30,30,1)'}
	sizedic={0:7,1:10}
	d,dd=[],[]
	ch,n=0,0
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
		for item in ls[7].split(';'):
			if 'EFF' in item:
				eff='OLD,'+item.split('=')[1] # snpEff version old
			if 'ANN' in item:
				eff='NEW,'+item.split('=')[1] # snpEff version new
		annotation=str(ch)+":"+str(relpos)+'<br>'+str(mut)+'/'+str(mut+wt)+' mut/total allels<br>'+'<br>'.join(eff.split(','))
		chartcolor=colordic[disruptive%2]
		chartsize=sizedic[disruptive%2]
		d.append([ch,relpos,relpos,abspos,frequency,disruptive,annotation,wt,mut,chartcolor,chartsize])
		if disruptive:
			dd.append([ch,relpos,relpos,abspos,frequency,disruptive,annotation,wt,mut,chartcolor,chartsize])
		n+=1
	d=map(list,zip(*d))
	dd=map(list,zip(*dd))
	return d,dd

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

#encoded_image = base64.b64encode(open('testimg.png', 'rb').read())

app.layout = html.Div(children=[
	dcc.Markdown(children=title_text),
	dcc.Upload(id='file_upload',children=html.Button('Upload .vcf File')),
	dcc.Graph(id='BSA_plot'),
	html.Div([
			html.Label('Low frequency cutoff'),
			dcc.Slider(id='freq-slider',min=0,max=100,value=25,step=None,marks={str(WIN): str(WIN) for WIN in range(0,100,5)}),
			dcc.Markdown('_'),
			html.Label('Sliding window size'),
			dcc.Slider(id='win-slider',min=5,max=50,value=10,step=None,marks={str(WIN): str(WIN) for WIN in range(5,50)}),
		],style={'columnCount': 1}),
		#dcc.Graph(figure=figure2,id='test hist'),
	
	dcc.Markdown((""" ### SNP annotation """)),
	html.A(children='Download Table' , id='download',download="BSAdata.tsv", href="",target="_blank"),
		#html.Label('SNP annotation'),
	html.Div([
		dt.DataTable(
			id='table',
			# Initialise the rows
			rows=[{}],
			columns=ColumnOrder,
			column_widths=[100,100,100,100,100,200,100,100,100,400,100],
			row_selectable=True,
			#filterable=True,
			sortable=True,
			selected_row_indices=[],
		),
	],style={'fontSize': 12}),
	#html.Div([html.Img(src='data:image/png;base64,{}'.format(encoded_image))]),
])


@app.callback(
	dash.dependencies.Output('BSA_plot', 'figure'),
	[
		dash.dependencies.Input('win-slider', 'value'),
		dash.dependencies.Input('freq-slider', 'value'),
		dash.dependencies.Input('file_upload', 'contents'),
		dash.dependencies.Input('file_upload', 'filename'),
	])
def update_figure(WIN,LOW,contents,fname):
	# file reading
	if not fname:
		return go.Figure(   
				data=[go.Scatter(x=[1],y=[1])], # default figure
				layout=go.Layout(title='BSA plot server by Pu Huang',
					yaxis=dict(zeroline=False,title='Derived Allele Frequency', range=(0,YMAX)),
					xaxis=dict(zeroline=False,title='Genomic position',range=(0,pos)),
					shapes=[dict(type='rect',x0=startdic[ch][0],x1=startdic[ch][1],y0=0,y1=YMAX,fillcolor=startdic[ch][2], \
						line=dict(width=0)) for ch in chs],
					hovermode='closest',
					)
		)
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
	d, dd = readvcf(infile, LOW/100.)

	# calculate smoothcurve
	smoothpos, smoothfreq = smoother(d,WIN,STEP)
	
	return go.Figure(
		data=[
			go.Scatter(x=dd[3], y=dd[4], name=u'Disruptive SNPs',mode = 'markers', marker=dict(color=dd[9],size=dd[10]),hoverinfo='none'),
			go.Scatter(x=d[3], y=d[4], name=u'SNPs',mode = 'markers', marker=dict(color=d[9],size=d[10]),text=d[6]),
			go.Scatter(x=smoothpos, y=smoothfreq, name=u'Smooth line',mode='line', line=dict(color='rgb(255,50,50)'),hoverinfo='none'),
		],
		layout=go.Layout(title='BSA plot for "'+fname+'"',
			yaxis=dict(zeroline=False,title='Derived Allele Frequency', range=(0,YMAX)),
			xaxis=dict(zeroline=False,title='Genomic position',range=(0,pos)),
			shapes=[dict(type='rect',x0=startdic[ch][0],x1=startdic[ch][1],y0=0,y1=YMAX,fillcolor=startdic[ch][2], \
				line=dict(width=0)) for ch in chs],
			hovermode='closest',
		)
	)


@app.callback(
	dash.dependencies.Output('table', 'rows'),
	[
		dash.dependencies.Input('BSA_plot', 'selectedData'),
		dash.dependencies.Input('file_upload', 'filename'),
		#dash.dependencies.Input('BSA_plot', 'selectedData'),
	])
def display_selected_data(selectedData,fname):
	if selectedData==None:
		return [{}]
	points=selectedData[u'points']
	out=[]
	text=[]
	for point in points:
		if point[u'curveNumber']==1:
			text.append(point[u'text'])
	if text==[]:
		return [{}]
	fulltexts=[re.split('<br>',t) for t in text]
	out=[]
	effdic={'HIGH':'1 ','MODERATE':'2 ','LOW':'3 ','MODIFIER':'4 '}
	for snp_entry in fulltexts:
		for item in snp_entry[3:]:
			'''ColumnOrder=['Pos','Mut/Total','Annotation','Impact','Codon','Allele','Gene','AtLog','AtName','AtAnno','Transcript'] '''
			dic={}
			anno=item.replace('(','|').replace(')','|').split('|')
			# SNPeff new version got different format, be cautious
			if snp_entry[2]=='OLD':
				if anno[6]:
					cmd='grep '+anno[6]+' '+F_anno+' |  head -1 | cut -f11-'
					geneanno=commands.getstatusoutput(cmd)[1].split('\t')
				else:
					geneanno=['','','','','']
				values=[snp_entry[0],snp_entry[1].split()[0],anno[0],effdic[anno[1]]+anno[1],anno[3],anno[4],anno[6],geneanno[0],geneanno[1],geneanno[2],anno[9]]
				dic={ColumnOrder[i]:values[i] for i in xrange(len(ColumnOrder))}
					
			elif snp_entry[2]=='NEW':
				if anno[3]:
					genes=anno[3].split('-')
					for gene in genes:
						cmd='grep '+gene+' '+F_anno+' |  head -1 | cut -f11-'
						geneanno=commands.getstatusoutput(cmd)[1].split('\t')
						values=[snp_entry[0],snp_entry[1].split()[0],anno[1],effdic[anno[2]]+anno[2],anno[9],anno[10],gene,geneanno[0],geneanno[1],geneanno[2],anno[6]]
						dic={ColumnOrder[i]:values[i] for i in xrange(len(ColumnOrder))}
						out.append(dic)
				else:
					geneanno=['','','','','']
					values=[snp_entry[0],snp_entry[1].split()[0],anno[1],effdic[anno[2]]+anno[2],anno[9],anno[10],anno[3],geneanno[0],geneanno[1],geneanno[2],anno[6]]
					out.append(dic)
	return out

@app.callback(
	dash.dependencies.Output('download','href'),
	[dash.dependencies.Input('table','rows'),])
def update_downloader(contents):
	if (contents==[{}]):
		return 'data:text/tsv;charset=utf-8,'+urllib.quote(' ')
	tsvstring='\n'.join(['\t'.join(ColumnOrder)]+['\t'.join(i.values()) for i in contents])
	tsvstring='\n'.join(['\t'.join(ColumnOrder)]+['\t'.join([i[c] for c in ColumnOrder]) for i in contents])
	tsvstring = "data:text/tsv;charset=utf-8," + urllib.quote(tsvstring)
	return tsvstring

if __name__ == '__main__':
	app.run_server(debug=True, host='0.0.0.0')




