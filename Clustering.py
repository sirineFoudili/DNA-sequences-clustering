
'''
Foudili Sirine
M1
Bio-Info/USTHB
Fouille de données
==================
Clustering des séquences d'ADN
==================
'''

import random
from random import randint
import sys
from re import *
from PyQt5 import *
from PyQt5.Qt import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtGui import QImage
from PyQt5.QtGui import QIcon
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from time import *


'''
global a
a=None
'''
class App(QMainWindow):
	global f 
	global k
	global it
	global pts
	global minpts
	global cg
	f=None
	k=0
	it=None
	pts=0
	minpts=0
	cg=None

	'''	
################################################################################################################################
 ################################################################################################################################
 
                                                  K-MEANS
 
 ################################################################################################################################
 ################################################################################################################################
	'''
	def tri_ADN(self,l):
		c=''
		for i in range(39,len(l)):
			if(l[i] in {'A','C','G','T','-'}):
				c=c+l[i]
		return c

	def trans_fichier(self,fichier):
		f=open(fichier,'r')
		i=0
		seq={}
		for ligne in f:
			ligne=f.readline()
			s=self.tri_ADN(ligne)
			seq[i]=s
			i=i+1
		f.close()
		return seq

	def longueur(self,f):
		i=0
		for ligne in f:
			i=i+1
		return i

	def comparaison(self,c1,c2):
		if c1==c2:
			return 1
		else:
			return -1

	def score(self,match,mismatch,c1,c2):
		if self.comparaison(c1,c2)==1:
			return match
		else:
			return mismatch

	def levenshtein(self,chaine1,chaine2,match,mismatch,indel):
		mat=[ [0] * (len(chaine1)+2) for i in range(len(chaine2)+2)] 
		for i in range(len(chaine2)+2):
			for j in range(len(chaine1)+2):
				if(i==0):
					if(j>=2):
						mat[0][j]=chaine1[j-2] 
				if(j==0):
					if(i>=2):
						mat[i][0]=chaine2[i-2] 
				if(i==1):
					if(j>=2):
						mat[i][j]=j-1
				if( j==1):
					if(i>=2):
						mat[i][j]=i-1
				if((i>=2) and (j>=2)):
					mat[i][j]=min( mat[i-1][j]+indel , mat[i][j-1]+indel , mat[i-1][j-1]+self.score(match,mismatch,mat[0][j],mat[i][0]) )
		return mat[len(chaine2)+1][len(chaine1)+1]


	def existe(self,dic,seq):
		for valeur in dic.values():
			if valeur==seq:
				return 1
			else:
				return 0

	def centre(self,f,k):
		centre={}
		for i in range(k):
			key = random.choice(list(f))
			j=0
			while(self.existe(centre,f[key])==1):
				key = random.choice(list(f))
				j=j+1
			centre[i]=f[key]
		return centre

	def longest(self,l):
		m=''
		for i in range(len(l)-1):
			if len(l[i]) >= len(m):
				m=l[i]
		return m

	def car_frequent(self,l):
		a=l.count('A')
		c=l.count('C')
		g=l.count('G')
		t=l.count('T')
		gap=l.count('-')
		m=max(a,c,g,t,gap)
		if(m==a):
			return 'A'
		if(m==c):
			return 'C'
		if(m==g):
			return 'G'
		if(m==t):
			return 'T'
		if(m==gap):
			return '-'

	def mean(self,l):
		mn=''
		for p in range(len(self.longest(l))):
			c=''
			for i in range(len(l)):
				if(p<=len(l[i])-1):
					c=c+(l[i])[p]
			mn=mn+self.car_frequent(c)
		m=mn
		return m
	def liste_adn(self,l,f):
		n=[]
		if (isinstance(l,int)==True):
			n.append(f[l])
		else:
			for i in range(len(l)):
				n.append(f[l[i]])
		return n

	def centroid(self,c,f):
		centre={}
		for i in range(len(c)):
			centre[i]=self.mean(self.liste_adn(c[i],f))
		return centre

	def comp(self,centres1,centres2,k):
		c=0
		for i in range(k):
			c=c+self.levenshtein(centres1[i],centres2[i],0,1,1)
			if c==0:
				return 1
			else:
				return 0
# Indice du centre le plus proche à l'adn (l'indice du cluster)
	def plus_proche(self,adn,centres):
		p=0
		m=self.levenshtein(adn,centres[0],0,1,1)
		for i in range(1,len(centres)):
			if self.levenshtein(adn,centres[i],0,1,1)<=m:
				p=i
				m=self.levenshtein(adn,centres[i],0,1,1)
		return p

	def clust(self,f,centre,k):
		c={} 
		for i in range(k):
			c[i]=[]
		for i in range(len(f)):
			p=self.plus_proche(f[i],centre)
			c[p].append(i)
		return c

	def k_means(self,f,k):
		print("="*50,'\n',' '*20,'K-Means',' '*20,'\n',"="*50)
		global it
		print('\n itération: 0 \n')
		ctr=self.centre(f,k)
		print('centre: ',ctr,'\n')
		clust1=self.clust(f,ctr,k)
		print('cluster: ',clust1,'\n\n')
		ctrd=self.centroid(clust1,f)
		it=0
		while(self.comp(ctr,ctrd,k)==0):
			print('iteration: ',it+1,'\n')
			ctr=ctrd
			print('centre: ',ctr,'\n')
			a=self.clust(f,ctrd,k)
			print('cluster: ',a,'\n\n')
			ctrd=self.centroid(a,f)
			it=it+1
		print(it)
		return self.clust(f,ctrd,k)
	'''
 ################################################################################################################################
 ################################################################################################################################
 
                                                    K-MEDOID
 
 ################################################################################################################################
 ################################################################################################################################
 	'''
	def medo(self,l,f):
		distance={}
		d=0
		c=self.liste_adn(l,f)
		if (isinstance(l,int)==True):
			d=d+self.levenshtein(c[0],c[0],0,1,1)
			distance[l]=d
		else:
			for i in range(len(l)):
				for j in range(len(l)):
					d=d+self.levenshtein(c[i],c[j],0,1,1)
				distance[l[i]]=d
		if len(distance.values())==0:
			return ' '
		else:
			m=min(distance.values())
			for c,v in distance.items():
				if v==m:
					medo=c
			return f[medo]

	def medoid(self,c,f):
		med={}
		for i in range(len(c)):
			med[i]=self.medo(c[i],f)
		return med



	def k_medoid(self,f,k):
		print("="*50,'\n',' '*20,'K-Medoïd',' '*20,'\n',"="*50)
		global it
		ctr=self.centre(f,k)
		clust1=self.clust(f,ctr,k)
		print('centre: ',ctr,'\n')
		clust1=self.clust(f,ctr,k)
		print('cluster: ',clust1,'\n\n')
		ctrd=self.medoid(clust1,f)
		it=0
		while(self.comp(ctr,ctrd,k)==0):
			print('itération:',it)
			ctr=ctrd
			print('centre: ',ctr,'\n')
			a=self.clust(f,ctrd,k)
			print('cluster: ',a,'\n\n')
			ctrd=self.medoid(a,f)
			it=it+1
		print(it)
		return self.clust(f,ctrd,k)

	'''
################################################################################################################################
 ################################################################################################################################
 
                                                 DBSCAN
 
 ################################################################################################################################
 ################################################################################################################################
	'''

 #Le point fera partie du voisinage  
	def voisinage(self,f,point,eps):
		voisins=[]
		for i in range(len(f)):
			if self.levenshtein(point,f[i],0,1,1)<=eps:
				voisins.append(i)
		return voisins

	def coeur(self,f,point,eps,minpts):
		if len(self.voisinage(f,point,eps))>=minpts:
			return 1
		else:
			return 0


	def classifie(self,dic,seq):
		existe=0
		for valeur in dic.values():
			for i in range(len(valeur)):
				if valeur[i]==seq:
					existe=1
		return existe

	def dbScan(self,f,eps,minpts):
		print("="*50,'\n',' '*20,'DBSCAN',' '*20,'\n',"="*50)
		global it
		it=0
		clus={}
		noise=[]
		for point in f.keys():
			print('itération: ',it,'\n')
			if self.classifie(clus,point)==0:
				if self.coeur(f,f[point],eps,minpts-1)==1:
					print('Coeur: ',f[point],'\n')
					c=self.voisinage(f,f[point],eps)
					clus[it]=c
					print('Voisinage: ',c,'\n')
					it=it+1
				else:
					noise.append(point)
		print('noise: ',noise,'\n')
		return clus

	'''
################################################################################################################################
 ################################################################################################################################
 
                                                 AGNES
 
 ################################################################################################################################
 ################################################################################################################################
	'''


	def distance(self,tmp,f):
		dist=[ [0] * (len(f)) for i in range(len(f))]
		for i in tmp.keys():
			for j in tmp.keys():
				dist[i][j]=self.levenshtein(f[i],f[j],0,1,1) 
		return dist

	def minimum(self,mat):
		x=0
		y=0
		while(mat[x][y]==0 and x<=len(mat)-1):
			y=y+1
			if y==len(mat)-1:
				y=0
				x=x+1
		m=mat[x][y]
		l=[x,y]
		for i in range(len(mat)):
			for j in range(len(mat)):
				if(min(mat[i][j],m)==mat[i][j] and mat[i][j]!=0):
					m=min(mat[i][j],m)
					l=[i,j]
		return l

	def remplacer(self,liste,tmp,f):
		tmp[liste[0]]=self.mean(self.liste_adn(liste,f))
		del tmp[liste[1]]
		return tmp


	def mis_a_jour(self,c,liste):
		if(isinstance(c[liste[0]],int)==True) and ( isinstance(c[liste[1]],int)==True ):
			l=[]
			l.append(c[liste[0]])
			l.append(c[liste[1]])
		else:
			if(isinstance(c[liste[0]],int)==False) and ( isinstance(c[liste[1]],int)==False):
				l=c[liste[0]]+c[liste[1]]
			if(isinstance(c[liste[0]],int)==True):
				t=[]
				t.append(c[liste[0]])
				l=c[liste[1]]+t
			if(isinstance(c[liste[1]],int)==True):
				t=[]
				t.append(c[liste[1]])
				l=c[liste[0]]+t
		c[liste[0]]=l
		del c[liste[1]]
		return c


	def agnes(self,f):
		print("="*50,'\n',' '*20,'Agnes',' '*20,'\n',"="*50)
		global it
		print('\n itération: 0 \n')
		it=0
		c={}
		clust={}
		tmp=dict(f)
		for i in range(len(tmp)):
			c[i]=i
			clust[it]=dict(c)
		it=it+1
		print('\n itération: ',it,' \n')
		i=0
		while len(c)>1 :
			m=self.minimum(self.distance(tmp,f))
			self.remplacer(m,tmp,f)
			self.mis_a_jour(c,m)
			clust[it]=dict(c)
			print('cluster :',i,'\n',c,'\n')
			i=i+1
			it=it+1
		return clust


	def affichage_agnes(self,CLUS):
		message=''
		p=0
		for p in CLUS.keys():
			message=message+'\n\nItération :'+str(p)+'\n\n'
			for i in CLUS[p].keys():
				message=message+'\nCluster '+str(i)+': '+str((CLUS[p])[i])+'\n'
		return message


	def affichage_detaillé_agnes(self,CLUS,f):
		message=''
		p=0
		#print('f :',f,'\n')
		for p in CLUS.keys():
			message=message+'\n\n Itération :'+str(p)+'\n\n'
			for i in CLUS[p].keys():
				message=message+'\nCluster '+str(i)+': \n'
				if (isinstance((CLUS[p])[i],int)==True):
					a=int((CLUS[p])[i])
					message=message+str(f[a])+'\n'
				else:
					for j in range(len((CLUS[p])[i])):
						b=int((((CLUS[p])[i])[j])) 
						message=message+str(f[b])  +'\n'
			#if p>0:
			print('====it',p)
			message=message+"\n\n\nPreformances :\n inertie intra-classe = "+str(self.inertie_intra(CLUS[p],cg)/10)+"\ninertie inter-classe = "+str(self.inertie_inter(CLUS[p],cg)/10)+'\n'   
		return message

	'''
################################################################################################################################
 ################################################################################################################################
 
                                                 PERFORMANCES
 
 ################################################################################################################################
 ################################################################################################################################
	'''
	def cg(self,c,f):
		moyenne=self.centroid(c,f)
		return moyenne

	def inertie_intra(self,cluster,cg):
		print('====cluster  :',cluster,'len',len(cluster))
		global f
		
		d=0
		for k in cluster.keys():
			print('cluster[k]',cluster[k])
			if (isinstance(cluster[k],int)==True):
				d=d+self.levenshtein(cg,f[cluster[k]],0,1,1)
			else:
				for j in range(len(cluster[k])):
					d=d+self.levenshtein(cg,f[cluster[k][j]],0,1,1)
		print("intra :",d)
		return d

	def inertie_inter(self,cluster,cg):
		print('===cluster :',cluster,'len',len(cluster))
		d=0
		for k in cluster.keys():
			d=d+self.levenshtein(cg,self.medo(cluster[k],f),0,1,1)
		print('inter',d)
		return d


	def affichage(self,CLUS):
		message=''
		for i in range(len(CLUS)):
			message=message+'\n\nCluster '+str(i+1)+': '+str(CLUS[i])+'\n'
		return message
	def affichage_detaillé(self,CLUS,f):
		global cg
		global it
		it=it+1
		message=''
		for i in range(len(CLUS)):
			message=message+'\n\nCluster '+str(i+1)+': \n '
			for j in range(len(CLUS[i])):
				message=message+str(CLUS[i][j])+' : '+str(f[CLUS[i][j]])+'.\n'
		message=message+"\n\n\nLe nombre d'itération :"+str(it)+'\n'
		message=message+"\n\n\nPreformances :\n inertie intra-classe = "+str(self.inertie_intra(CLUS,cg)/10)+"\ninertie inter-classe = "+str(self.inertie_inter(CLUS,cg)/10)+'\n'
		return message

	def affichage_detaillé_db(self,CLUS,f):
		global cg
		global it
		it=it+1
		message=''
		for i in range(len(CLUS)):
			message=message+'\n\nCluster '+str(i+1)+': \n '
			for j in range(len(CLUS[i])):
				message=message+str(CLUS[i][j])+' : '+str(f[CLUS[i][j]])+'.\n'
		message=message+"\n\n\nLe nombre d'itération :"+str(it)+'\n'
		message=message+"\n\n\nPreformances :\n inertie intra-classe = "+str((self.inertie_intra(CLUS,cg)/10)-self.inertie_inter(CLUS,cg)/10)+"\ninertie inter-classe = "+str(self.inertie_inter(CLUS,cg)/10)+'\n'
		return message

	def Temps(self,algo,t1):
		algo
		t2=clock()
		return (t2 - t1)

	def fermer(self,f):
		open('dna_examples.txt','r').close()
	'''
################################################################################################################################
 ################################################################################################################################
                                           GRAPHES
 
 ################################################################################################################################
 ################################################################################################################################
	'''

	def graph(self,cluster):
		plt.style.use('ggplot')
		fig,ax =plt.subplots()
		color=['b','g','r','c','m','y','k']
		for k in range(len(cluster)):
			
			if k<=len(color)-1:
				coul=color[k]
			else:
				coul=color[random.randint(0,len(color)-1)]
			y=[]
			df={'X':[],'Y':[]}
			for j in range(len(cluster[k])):
				y.append(k)
			df['Y'].extend(cluster[k])
			df['X'].extend(y)
			ax.scatter(df['Y'],df['X'], c=coul,s=30, alpha=0.5 )
		ax.set_ylabel("Clusters")
		ax.set_xlabel("Séquences")
		ax.set_title("Classification des séquences")
		plt.show()


	'''
 ################################################################################################################################
 ################################################################################################################################
                                          INTREFACE GRAPHIQUE
 
 ################################################################################################################################
 ################################################################################################################################
	'''
	def __init__(self):
		"""Iitialisation de l'interface graphique"""
		super(QMainWindow,self).__init__()
		#Initialisation du titre
		self.title = "Classification de séquences d'ADN"
		#Initialisation de la taille de la fenetre
		self.left = 200
		self.top = 100
		self.width = 1000
		self.height = 600
		self.initUI()

	def initUI(self):
		#self.setWindowIcon(QtGui.QIcon("icon.png"))
		self.setWindowTitle(self.title)
		self.setGeometry(self.left, self.top, self.width, self.height)
		 # Set window background color
		self.setAutoFillBackground(True)
		p = self.palette()
		p.setColor(self.backgroundRole(), Qt.white)
		self.setPalette(p)
		
		oImage = QImage("bg.png")
		sImage = oImage.scaled(QSize(1400,1000))                   # resize Image to widgets size
		palette = QPalette()
		palette.setBrush(10, QBrush(sImage))                     # 10 = Windowrole
		self.setPalette(palette)
		
		'''
		self.label = QLabel('Test', self)                        # test, if it's really backgroundimage
		self.label.setGeometry(50,50,200,50)
		'''
		mainMenu = self.menuBar() 
		'''
		helpMenu = mainMenu.addMenu('A Propos ?')
		helpButton = QAction('?', self)
		helpButton.setShortcut('Ctrl+H')
		helpButton.setStatusTip('A Propos de L\'application')
		helpButton.triggered.connect(self.menu_click0)
		helpMenu.addAction(helpButton)
		'''
		#Ajout, positionnement et configuration des boutons
		button1 = QPushButton('Nombre de clusters', self)
		button1.move(700,80)
		button1.resize(200,40)
		button1.clicked.connect(self.on_click1)

		button2 = QPushButton("Importer des séquences d'ADN", self)
		button2.move(500,80)
		button2.resize(200,40) 
		button2.clicked.connect(self.on_click2)

		button3 = QPushButton('K-means', self)
		button3.move(150,410)
		button3.resize(200,40) 
		button3.clicked.connect(self.on_click3)

		button4 = QPushButton('K-medoïd', self)
		button4.move(150,460)
		button4.resize(200,40) 
		button4.clicked.connect(self.on_click4)

		button5 = QPushButton('Dbscan', self)
		button5.move(150,510)
		button5.resize(200,40) 
		button5.clicked.connect(self.on_click5)

		button6 = QPushButton('Agnes', self)
		button6.move(150,570)
		button6.resize(200,40) 
		button6.clicked.connect(self.on_click6)

		button7 = QPushButton('Rayon de voisinage', self)
		button7.move(900,80)
		button7.resize(200,40)
		button7.clicked.connect(self.on_click7)

		button8 = QPushButton('Nombre minimal de voisins', self)
		button8.move(1100,80)
		button8.resize(200,40)
		button8.clicked.connect(self.on_click8)

		self.textArea = QPlainTextEdit(self)
		self.textArea.move(500,350)
		self.textArea.resize(700,400)
		self.textArea.insertPlainText("Clustering...\n\n")
		self.show()

	def openFileNameDialog(self):
		"""La fonction openFileNameDialog demande à l'utilisateur de sélectionner un fichier parmi 
		   le fichiers existants dans la machine et return l'emlacement du fichier"""
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()",
								 "","All Files (*);;Python Files (*.py)", options=options)
		if fileName:
			return fileName




	def getInteger(self):
		i, okPressed = QInputDialog.getInt(self, "Nombre de clusters","Insérer le nombre de clusters k:", 1, 0, 500, 1)
		if okPressed:
			return (i)

	def getInteger_minpts(self):
		"""La fonction getInteger demande à l'utilisateur de saisir un entier et renvoi l'entier saisi"""
		i, okPressed = QInputDialog.getInt(self, "MinPts","Insérer le nombre minimal de voisins MinPts :", 1, 0, 500, 1)
		if okPressed:
			return (i)

	def getInteger_eps(self):
		"""La fonction getInteger demande à l'utilisateur de saisir un entier et renvoi l'entier saisi"""
		i, okPressed = QInputDialog.getInt(self, "EPS","Insérer le rayon de voisinage eps :", 1, 0, 500, 1)
		if okPressed:
			return (i)






###########################

#Configuration des cliques sur les boutons
	@pyqtSlot()
	def on_click1(self):
		global k
		k=self.getInteger()
		self.textArea.insertPlainText("> Nombre de clusters: "+str(k)+" \n\n")


	@pyqtSlot()
	def on_click2(self):
		global f
		global cg 
		name = self.openFileNameDialog()
		f= self.trans_fichier(name)
		cg=self.mean(f)
		self.textArea.insertPlainText("\n\n> Fichier :\n"+ name +" \n\n")

	@pyqtSlot()
	def on_click3(self):
		r=self.k_means(f,k)
		self.textArea.insertPlainText("\n\n> Méthode K-means :\n Démarage . . . \n\nRésultat: \n"+self.affichage(r)+"\n\n\nRésultat Détaillé :\n"+self.affichage_detaillé(r,f)+"\n")
		self.graph(r)
		t1=clock()
		t=self.Temps(self.k_means(f,k),t1)
		self.textArea.insertPlainText("> Temps d'éxecution : "+str(t)+"\n\n")
	@pyqtSlot()
	def on_click4(self):
		r=self.k_medoid(f,k)
		self.textArea.insertPlainText("\n\n> Méthode K-medoïd :\n Démarage . . . \n\nRésultat: \n"+self.affichage(r)+"\n\n\nRésultat Détaillé :\n"+self.affichage_detaillé(r,f)+"\n\n")
		self.graph(r)
		t1=clock()
		t=self.Temps(self.k_medoid(f,k),t1)
		self.textArea.insertPlainText("> Temps d'éxecution : "+str(t)+"\n\n")
	@pyqtSlot()
	def on_click5(self):
		r=self.dbScan(f,pts,minpts)
		self.textArea.insertPlainText("\n\n> Méthode DBSCAN :\n Démarage . . . \n\nRésultat: \n"+self.affichage(r)+"\n\n\nRésultat Détaillé :\n"+self.affichage_detaillé_db(r,f)+"\n\n")
		self.graph(r)
		t1=clock()
		t=self.Temps(self.dbScan(f,pts,minpts),t1)
		self.textArea.insertPlainText("> Temps d'éxecution : "+str(t)+"\n\n")
	@pyqtSlot()
	def on_click6(self):
		r=self.agnes(f)
		self.textArea.insertPlainText("\n\n> Méthode Agnes :\n Démarage . . . \n\nRésultat: \n"+self.affichage_agnes(r)+"\n\n\nRésultat Détaillé :\n"+self.affichage_detaillé_agnes(r,f)+"\n\n")
		t1=clock()
		t=self.Temps(self.agnes(f),t1)
		self.textArea.insertPlainText("> Temps d'éxecution : "+str(t)+"\n\n")
	@pyqtSlot()
	def on_click7(self):
		global pts
		pts=self.getInteger_eps()
		self.textArea.insertPlainText("> Rayon de voisinage EPS: "+str(pts)+" \n\n")

	@pyqtSlot()
	def on_click8(self):
		global minpts
		minpts=self.getInteger_minpts()
		self.textArea.insertPlainText("> Nombre minimal de voisins MinPts: "+str(minpts)+" \n\n")

	
	@pyqtSlot()
	def menu_click0(self):
		QMessageBox.information(self,"A Propos"," \n\n \n-")




###########################

if __name__ == '__main__':

	app = QApplication(sys.argv)
	ex = App()
	sys.exit(app.exec_())



