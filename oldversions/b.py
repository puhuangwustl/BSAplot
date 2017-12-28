#!/usr/bin/python

import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import plotly.plotly as py
import numpy as np
import sys

import gzip


LOW=.25 # low frequency cut
WIN=10
STEP=2
YMAX=1.2
F_fai='/media/Setaria_1/Share/Sviridis_vphytozome1.1/Sviridis_311_v1.0.fa.fai'
F_vcf='mut.eff.vcf.gz'

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
d=[]
ch=0
with gzip.open(F_vcf, 'rb') as f:
	for l in f:
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
				eff=item.split('=')[1]
		annotation=str(ch)+":"+str(relpos)+'<br>'+str(mut)+'/'+str(mut+wt)+' mut/total allels<br>'+'<br>'.join(eff.split(','))
		d.append([ch,relpos,relpos,abspos,frequency,disruptive,annotation,wt,mut])

d=map(list,zip(*d))

# calculate smoothing line
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


# call plot

colordic={0:'rgba(0,70,190,.7)',1:'rgba(255,30,30,1)'}
sizedic={0:7,1:10}
snpcolors=[colordic[i%2] for i in d[5]]
snpsizes=[sizedic[i] for i in d[5]]


figure=go.Figure(
	data=[
		go.Scatter(x=d[3], y=d[4], name=u'SNPs',mode = 'markers', marker=dict(color=snpcolors,size=snpsizes),text=d[6]),
		go.Scatter(x=smoothpos, y=smoothfreq, name=u'Smooth line',mode='line', line=dict(color='rgb(255,50,50)'),hoverinfo='none'),
		],
	layout=go.Layout(title='BSA plot for "'+F_vcf.split('/')[-1]+'"',
			yaxis=dict(title='Derived Allele Frequency', range=(0,YMAX)),
			xaxis=dict(title='Genomic position'),
			shapes=[dict(type='rect',x0=startdic[ch][0],x1=startdic[ch][1],y0=0,y1=YMAX,fillcolor=startdic[ch][2], \
				line=dict(width=0)) for ch in chs]),
)
#py.image.save_as(figure, filename='a-simple-plot.png')



hisdata = [go.Histogram(x=w_depth,
                     histnorm='probability')]

figure2=go.Figure(
	data=hisdata,
	layout=go.Layout(title='Histogram for smoothing window read coverage',
		xaxis=dict(title='Window read coverage')
		)
	)


### call server
app = dash.Dash()

# page layout
markdown_text = '''
### BSA plot for *Setaria viridis*

Bulked Segregant Analysis by deep sequencing.

Please cite 'Huang et al. *Sparse panicle1* is required for inflorescence development in *Setaria viridis* and maize. Nature plants 3.5 (2017): 17054.' [link](https://www.nature.com/articles/nplants201754)
'''

app.layout = html.Div(children=[
	dcc.Markdown(children=markdown_text),
	dcc.Graph(figure=figure,id='BSA plot'),
	dcc.Graph(figure=figure2,id='test hist'),
])










if __name__ == '__main__':
	app.run_server(debug=True, host='0.0.0.0')











