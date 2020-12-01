# -*- coding: utf-8 -*-
"""
Created on Thu May 13 18:38:02 2020

@author: FOUDILI sirine
"""
#Lire le fichier fasta
from Bio import SeqIO
#générer le vecteur de comaraison
from itertools import combinations 
from itertools import product
import numpy
from numpy import array
#Needelman wunch alignement
from Bio import pairwise2
from Bio import Align
from Bio.pairwise2 import format_alignment

from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO

#Phylogénie
#from Bio.Cluster import treecluster
#from scipy.cluster.hierarchy import average, dendrogram
import biotite.sequence.phylo as phylo
import matplotlib.pyplot as plt
import biotite.sequence.graphics as graphics

from Bio import Phylo

#Affichage
from flask import render_template
from bokeh.plotting import figure
from bokeh.embed import components

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot
from bokeh.io import output_file, show
from bokeh.resources import CDN
import time
import os

#filtre de bloom
import math
import hashlib

class Alignement:
    #entrer methode,dictionnaire vide, fichier (sous forme de dict)
    def __init__(self,m,k,sequences={}): 
        self.adn = sequences.copy()
        self.methode=int(m)
        self.vecteurs_kmers={}
        self.seqs=[seq for seq  in self.adn.keys()]
        self.k=k
        self.kmers=self.k_mer()
        self.scores_matrice=self.score()
        self.f= None
        self.temps=time.time()
        self.details=None
        self.taille=len(self.adn)
        self.it=None
        
    '''
            K-mers
        
    '''
    #un dictionnaire {kmer:{id_seq:nbr kmers}}
    def k_mer(self):
        kfreq={}       
        for id,seq in self.adn.items():
            for i in range(0 ,len(seq) - self.k + 1) :
                kmer=seq[i:i+ self.k]
                if kmer in kfreq:
                    if id in kfreq[kmer]:
                        kfreq[kmer][id]+=1
                    else:
                        kfreq[kmer].update({id:1})                    
                else:
                    kseq={}
                    kseq[id]=1
                    kfreq[kmer]=kseq.copy()
        return kfreq
    '''
            FILTRAGE
    '''
    #nbr(id,'aac') --> nbr de kmers aac dans la seq id
    def nbr(self,seq,kmer):
        return self.kmers[kmer][seq]
    
    #vérifier si 2 seq ont des kmers communs
    def commun(self,kmer,seq1,seq2):
        if seq1 in self.kmers[kmer].keys() and seq2 in self.kmers[kmer].keys():
            
            return True
        else:
            return False
        
    # nbr de kmers communs
    def comptage(self,kmer,seq1,seq2):
        if(self.commun(kmer,seq1,seq2)):        
            cmpt=min(self.nbr(self.kmers, seq1, kmer),self.nbr(self.kmers, seq2, kmer))
            return cmpt
        else:
            return 0

    #vérifier si le nbr de kmers communs>= hit
    def filtrage(self,seq1,seq2,hit):
        c=0
        for k in self.kmers.keys():
            c=c+ self.comptage(k,seq1,seq2)
        if(c>= hit):
            return True
        else:
            return False  
        
    '''
       COMPARAISON AVEC K-MERS VECTEURS
    '''
    #vecteur global, contient tous les kmers possibles selon la taille k (taille_alphabet^k)
    def vecteur_kmer(self):    
        v=[''.join(c) for c in product('ACGT', repeat=self.k)]
        self.vecteurs_kmers['ref']=v
        return v
    #convertir une sequen en vecteur a partir de son id
    def seq_to_vecteur(self,seq,v):
        vecteur=[]
        for kmer in v:
            if kmer in self.kmers.keys() and seq in self.kmers[kmer].keys():
                vecteur.append(self.nbr(seq,kmer))
            else:
                vecteur.append(0)
        return numpy.array(vecteur)

    #distance euclidienne entre 2 seq (vecteurs)
    def distance_kmers(self,a,b):
        v=self.vecteur_kmer()
        v1=self.seq_to_vecteur(a,v)
        print("v(",a,")",v1)
        v2=self.seq_to_vecteur(b,v)
        print("v(",b,")",v2)
        dist = numpy.linalg.norm(v1-v2)
        print("d (",a,"-",b,") :",dist)
        return round(dist, 4)
    
        '''
           COMPARAISON AVEC K-MERS VECTEURS & FILTRE DE BLOOM
        '''
    
        #Appel par ID   
    def sequence_vers_vecteur(self,s,v):
        vecteur=[0]*len(v)
        seq=self.adn[s]
        taille=min(len(v),len(seq)-self.k +1)
        print("n:",taille)
        bloom = BloomFilter(taille,0.1)
        print("k:",bloom.k)
        print("m:",bloom.N)
        for i in range(0 ,len(seq) - self.k + 1) :
            kmer=seq[i:i+ self.k]
            position=v.index(kmer)
            if bloom.contains(kmer)==False:
                #ajouter le kmer au filtre
                #incrémenter 1
                bloom.add(kmer)
                vecteur[position]+=1
            else:
                #cas de repetition
                vecteur[position]+=1
        print("S:",s,"\nbloom:\n",bloom.bitfield)
        return numpy.array(vecteur)
        
    #Appel par ID
    def distance_kmers_bloom(self,a,b):
        v=self.vecteur_kmer()
        if a in self.vecteurs_kmers.keys():
            v1=self.vecteurs_kmers[a]
            
        else:
            v1=self.sequence_vers_vecteur(a, v)
            self.vecteurs_kmers[a]=v1
            print("v(",a,"):",v1)
            
        if b in self.vecteurs_kmers.keys():
            v2=self.vecteurs_kmers[b]
           
        else:
            v2=self.sequence_vers_vecteur(b, v)
            self.vecteurs_kmers[b]=v2
            print("v(",b,"):",v2)
        dist = numpy.linalg.norm(v1-v2)
        
        return round(dist, 4)
        
        '''
            COMPARAISON AVEC NEEDELMAN & WUNCH
       '''
    def needelman_wunch_align(self,S1,S2):
        a=self.adn[S1]
        b=self.adn[S2]
        if a==b:
            return 0
        else:  
           
            alignments = pairwise2.align.globalmx(a,b,0,1)
            self.f.write(format_alignment(*alignments[0]))
            self.details.write('>')
            self.details.write(str(S1)+'-'+self.it)
            self.details.write('\n')
            self.details.write(str(alignments[0][0]))
            self.details.write('\n')
            self.details.write('>')
            self.details.write(str(S2)+'-'+str(self.it))
            self.details.write('\n')
            self.details.write(str(alignments[0][1]))
            self.details.write('\n')
            self.details.write('\n')
            score=alignments[0].score
            return score
      
    #score d'alignement entre 2 seq selon la méthode
    def alignement(self,S1,S2):
        if self.methode==1:
            print("choix: NW")
            return self.needelman_wunch_align(S1,S2)
        if self.methode==3:
            print("choix: kmer bloom")
            return self.distance_kmers_bloom(S1,S2)
        if self.methode==2:
            print("choix: kmer simple")
            return self.distance_kmers(S1,S2)
   
    #Matrice des scores     
    def score(self):
        
        sequences=self.adn.copy()
        scores=numpy.zeros((len(self.adn),len(self.adn)))
        if self.methode==1:
            self.f=open("static/files/AlignementPaire.fasta", "w")
            self.details=open("static/files/AlignementNWk.fasta", "w")
        i=0
        for S1 in self.adn.keys():
            
            j=0
            for S2 in self.adn.keys():
                if S2 in sequences:
                    if self.methode==1:
                        self.it=str(i)+'|'+str(j)
                        self.f.write('---------------------\n séquences : ')
                        self.f.write(str(S1))
                        self.f.write(' et ')
                        self.f.write(str(S2))
                        self.f.write('\n')
                    scores[i][j]=self.alignement(S1,S2)
                    scores[j][i]=scores[i][j]
                    if j==0:
                        del sequences[S1]
                    j=j+1
                else:
                    j=j+1
                
            i=i+1
        if self.methode==1:
            self.f.close()
            self.details.close()
        return scores
    

class Phylogenie:
    def __init__(self,matrice,seqs): 
        self.arbre = None
        tree=phylo.upgma(matrice)
        self.arbre=tree.to_newick(include_distance=True,labels=seqs)
        print("---------\n-----------\nFIGURE")
        fig, ax = plt.subplots(figsize=(6.0, 6.0))
        graphics.plot_dendrogram(ax, tree, labels=seqs)
        fig.tight_layout()
        
        
        image= "static/img/phylogenie/arbre.png"
        if os.path.isfile(image):
            print("is file")
            os.remove(image) 
        plt.savefig("static/img/phylogenie/arbre.png")
        

        
class Affichage:
    
   
    #transforme une liste de cars ['a','c','-',...] à une liste de couleurs
    def get_colors(self,seqs):
        """make colors for bases in sequence"""
        text = [i for s in list(seqs) for i in s]
        clrs =  {'A':'red','T':'green','G':'orange','C':'blue','-':'white','N':'black'}
        colors=[]
        for i in text:
            colors.append(clrs[i])
        return colors
    def k_mers_colors(self,seqs):
        text = [i for s in list(seqs) for i in s]
        #vecteur de ref
        v=[''.join(c) for c in product('ACGT', repeat=self.k)]
        #vecteur de 16 couleurs
        couleurs=['peachpuff', 'peru', 'pink', 'plum', 'powderblue', 'purple', 'red', 'rosybrown', 'royalblue', 'saddlebrown', 'salmon', 'sandybrown', 'seagreen', 'sienna', 'silver', 'skyblue', 'slateblue', 'slategray', 'slategrey', 'snow', 'springgreen', 'steelblue', 'magenta', 'maroon', 'mediumaquamarine', 'mediumblue', 'mediumorchid', 'mediumpurple', 'mediumseagreen', 'mediumslateblue', 'mediumspringgreen', 'mediumturquoise', 'mediumvioletred', 'midnightblue', 'tan', 'teal', 'thistle', 'tomato', 'turquoise', 'violet', 'wheat', 'yellow', 'yellowgreen', 'mistyrose', 'moccasin', 'navajowhite', 'navy', 'olive', 'olivedrab', 'orange', 'orangered', 'orchid', 'palegoldenrod', 'palegreen', 'paleturquoise', 'palevioletred', 'lightgray', 'lightgreen', 'lightgrey', 'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue', 'lightslategray', 'lightslategrey', 'lightsteelblue', 'lightblue', 'lightcoral', 'lawngreen', 'hotpink', 'indianred', 'indigo', 'gold', 'goldenrod', 'gray', 'green', 'greenyellow', 'grey', 'forestgreen', 'fuchsia', 'gainsboro', 'crimson', 'cyan', 'darkblue', 'darkcyan', 'darkgoldenrod', 'darkgray', 'darkgreen', 'darkgrey', 'darkkhaki', 'darkmagenta', 'darkolivegreen', 'darkorange', 'darkorchid', 'darkred', 'darksalmon', 'darkseagreen', 'darkslateblue', 'darkslategray', 'darkslategrey', 'darkturquoise', 'darkviolet', 'deeppink', 'deepskyblue', 'dimgray', 'dimgrey', 'dodgerblue', 'firebrick', 'bisque', 'black', 'blanchedalmond', 'blue', 'blueviolet', 'brown', 'burlywood', 'cadetblue', 'chartreuse', 'chocolate', 'coral', 'cornflowerblue', 'aqua', 'aquamarine','peachpuff', 'peru', 'pink', 'plum', 'powderblue', 'purple', 'red', 'rosybrown', 'royalblue', 'saddlebrown', 'salmon', 'sandybrown', 'seagreen', 'sienna', 'silver', 'skyblue', 'slateblue', 'slategray', 'slategrey', 'snow', 'springgreen', 'steelblue', 'magenta', 'maroon', 'mediumaquamarine', 'mediumblue', 'mediumorchid', 'mediumpurple', 'mediumseagreen', 'mediumslateblue', 'mediumspringgreen', 'mediumturquoise', 'mediumvioletred', 'midnightblue', 'tan', 'teal', 'thistle', 'tomato', 'turquoise', 'violet', 'wheat', 'yellow', 'yellowgreen', 'mistyrose', 'moccasin', 'navajowhite', 'navy', 'olive', 'olivedrab', 'orange', 'orangered', 'orchid', 'palegoldenrod', 'palegreen', 'paleturquoise', 'palevioletred', 'lightgray', 'lightgreen', 'lightgrey', 'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue', 'lightslategray', 'lightslategrey', 'lightsteelblue', 'lightblue', 'lightcoral', 'lawngreen', 'hotpink', 'indianred', 'indigo', 'gold', 'goldenrod', 'gray', 'green', 'greenyellow', 'grey', 'forestgreen', 'fuchsia', 'gainsboro', 'crimson', 'cyan', 'darkblue', 'darkcyan', 'darkgoldenrod', 'darkgray', 'darkgreen', 'darkgrey', 'darkkhaki', 'darkmagenta', 'darkolivegreen', 'darkorange', 'darkorchid', 'darkred', 'darksalmon', 'darkseagreen', 'darkslateblue', 'darkslategray', 'darkslategrey', 'darkturquoise', 'darkviolet', 'deeppink', 'deepskyblue', 'dimgray', 'dimgrey', 'dodgerblue', 'firebrick', 'bisque', 'black', 'blanchedalmond', 'blue', 'blueviolet', 'brown', 'burlywood', 'cadetblue', 'chartreuse', 'chocolate', 'coral', 'cornflowerblue', 'aqua', 'aquamarine','peachpuff', 'peru', 'pink', 'plum', 'powderblue', 'purple', 'red', 'rosybrown', 'royalblue', 'saddlebrown', 'salmon', 'sandybrown', 'seagreen', 'sienna', 'silver', 'skyblue', 'slateblue', 'slategray', 'slategrey', 'snow', 'springgreen', 'steelblue', 'magenta', 'maroon', 'mediumaquamarine', 'mediumblue', 'mediumorchid', 'mediumpurple', 'mediumseagreen', 'mediumslateblue', 'mediumspringgreen', 'mediumturquoise', 'mediumvioletred', 'midnightblue', 'tan', 'teal', 'thistle', 'tomato', 'turquoise', 'violet', 'wheat', 'yellow', 'yellowgreen', 'mistyrose', 'moccasin', 'navajowhite', 'navy', 'olive', 'olivedrab', 'orange', 'orangered', 'orchid', 'palegoldenrod', 'palegreen', 'paleturquoise', 'palevioletred', 'lightgray', 'lightgreen', 'lightgrey', 'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue', 'lightslategray', 'lightslategrey', 'lightsteelblue', 'lightblue', 'lightcoral', 'lawngreen', 'hotpink', 'indianred', 'indigo', 'gold', 'goldenrod', 'gray', 'green', 'greenyellow', 'grey', 'forestgreen', 'fuchsia', 'gainsboro', 'crimson', 'cyan', 'darkblue', 'darkcyan', 'darkgoldenrod', 'darkgray', 'darkgreen', 'darkgrey', 'darkkhaki', 'darkmagenta', 'darkolivegreen', 'darkorange', 'darkorchid', 'darkred', 'darksalmon', 'darkseagreen', 'darkslateblue', 'darkslategray', 'darkslategrey', 'darkturquoise', 'darkviolet', 'deeppink', 'deepskyblue', 'dimgray', 'dimgrey', 'dodgerblue', 'firebrick', 'bisque', 'black', 'blanchedalmond', 'blue', 'blueviolet', 'brown', 'burlywood', 'cadetblue', 'chartreuse', 'chocolate', 'coral', 'cornflowerblue', 'aqua', 'aquamarine','peachpuff', 'peru', 'pink', 'plum', 'powderblue', 'purple', 'red', 'rosybrown', 'royalblue', 'saddlebrown', 'salmon', 'sandybrown', 'seagreen', 'sienna', 'silver', 'skyblue', 'slateblue', 'slategray', 'slategrey', 'snow', 'springgreen', 'steelblue', 'magenta', 'maroon', 'mediumaquamarine', 'mediumblue', 'mediumorchid', 'mediumpurple', 'mediumseagreen', 'mediumslateblue', 'mediumspringgreen', 'mediumturquoise', 'mediumvioletred', 'midnightblue', 'tan', 'teal', 'thistle', 'tomato', 'turquoise', 'violet', 'wheat', 'yellow', 'yellowgreen', 'mistyrose', 'moccasin', 'navajowhite', 'navy', 'olive', 'olivedrab', 'orange', 'orangered', 'orchid', 'palegoldenrod', 'palegreen', 'paleturquoise', 'palevioletred', 'lightgray', 'lightgreen', 'lightgrey', 'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue', 'lightslategray', 'lightslategrey', 'lightsteelblue', 'lightblue', 'lightcoral', 'lawngreen', 'hotpink', 'indianred', 'indigo', 'gold', 'goldenrod', 'gray', 'green', 'greenyellow', 'grey', 'forestgreen', 'fuchsia', 'gainsboro', 'crimson', 'cyan', 'darkblue', 'darkcyan', 'darkgoldenrod', 'darkgray', 'darkgreen', 'darkgrey', 'darkkhaki', 'darkmagenta', 'darkolivegreen', 'darkorange', 'darkorchid', 'darkred', 'darksalmon', 'darkseagreen', 'darkslateblue', 'darkslategray', 'darkslategrey', 'darkturquoise', 'darkviolet', 'deeppink', 'deepskyblue', 'dimgray', 'dimgrey', 'dodgerblue', 'firebrick', 'bisque', 'black', 'blanchedalmond', 'blue', 'blueviolet', 'brown', 'burlywood', 'cadetblue', 'chartreuse', 'chocolate', 'coral', 'cornflowerblue', 'aqua', 'aquamarine']
        clrs={}
        for kmer,i in zip(v,range(len(v))):
            clrs[kmer]=couleurs[i]
        #colors = [clrs[k] for k in self.kmers]
        colors=[]
        for k in self.kmers:
            colors.append(clrs[k])
        return colors
    def view_alignment(self,aln, fontsize="9pt", plot_width=800):
    
        """Bokeh sequence alignment view"""
    
        #make sequence and id lists from the aln object
        seqs = [rec.seq for rec in (aln)]
        ids = [rec.id for rec in aln]
        #plus longue seq
        sizes = [len(rec) for rec in seqs]
        #taille=max(sizes)
        taille=min(sizes)
        #résultat :['a','c','-',...]
        text = [i for s in list(seqs) for i in s]
        #résultat :['green','red','white',...]
        colors = self.get_colors(seqs)
        #longueur des chaines (fixe)
        #N = len(seqs[0])
        N=taille
        #Nombre de séquences
        S = len(seqs)    
        width = .4
    
        x = numpy.arange(1,N+1)
        y = numpy.arange(0,S,1)
        #creates a 2D grid of coords from the 1D arrays
        xx, yy = numpy.meshgrid(x, y)
        #flattens the arrays/maillage (Affichage de longueur)
        gx = xx.ravel()
        gy = yy.flatten()
        #use recty for rect coords with an offset
        recty = gy+.5
        h= 1/S
    
        #now we can create the ColumnDataSource with all the arrays
        source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
        plot_height = len(seqs)*15+50
        x_range = Range1d(0,N+1, bounds='auto')
        if N>100:
            viewlen=100
        else:
            viewlen=N
        #view_range is for the close up view
        view_range = (0,viewlen)
        tools="xpan, xwheel_zoom, reset, save"
    
    
        #entire sequence view (no text, with zoom)
        p = figure(title=None, plot_width= plot_width, plot_height=50,
                   x_range=x_range, y_range=(0,S), tools=tools,
                   min_border=0, toolbar_location='below')
        rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                     line_color=None, fill_alpha=0.6)
        p.add_glyph(source, rects)
        p.yaxis.visible = False
        p.grid.visible = False  
    
        #sequence text view with ability to scroll along x axis
        p1 = figure(title=None, plot_width=plot_width, plot_height=plot_height,
                    x_range=view_range, y_range=ids, tools="xpan,reset",
                    min_border=0, toolbar_location='below')#, lod_factor=1)          
        glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                    text_font="monospace",text_font_size=fontsize)
        rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                    line_color=None, fill_alpha=0.4)
        p1.add_glyph(source, glyph)
        p1.add_glyph(source, rects)
    
        p1.grid.visible = False
        p1.xaxis.major_label_text_font_style = "bold"
        p1.yaxis.minor_tick_line_width = 0
        p1.yaxis.major_tick_line_width = 0
    
        p = gridplot([[p],[p1]], toolbar_location='below')
        script, div = components(p)
        cdn_js=CDN.js_files[0]
        return p


class Node:
    # Added parsing of the "na" variable to name and value.
    # Added a parent field
    def __init__(self, name_val):
        name, val_str = name_val[::-1].split(":")
        self.name = name
        self.value = float(val_str)
        self.children = []
        self.parent = None

    # Method to get the depth of the node (for printing)
    def get_depth(self):
        current_node = self
        depth = 0
        while current_node.parent:
            current_node = current_node.parent
            depth += 1
        return depth

    # String representation
    def __str__(self):
        #return "{}:{}".format(self.name, self.value)
        return self.name

class Arbre:
    def __init__(self, tree):
        self.newick = tree
        self.liste=self.arbre_vers_liste()
    
    def arbre_vers_liste(self):
        root = None
        # na was not defined before.
        na = ""
        stack = []
        for i in list(reversed(self.newick)):
            if i == ')':
                if na != "":
                    node = Node(na)
                    na = ""
                    if len(stack):
                        stack[-1].children.append(node)
                        node.parent = stack[-1]
                    else:
                        root = node
                    stack.append(node)
        
            elif i == '(':
                if (na != ""):
                    node = Node(na)
                    na = ""
                    stack[-1].children.append(node)
                    node.parent = stack[-1]
                stack.pop()
            elif i == ',':
                if (na != ""):
                    node = Node(na)
                    na = ""
                    stack[-1].children.append(node)
                    node.parent = stack[-1]
            else:
                # n was not defined before, changed to i.
                na += i
        
        # Just to print the parsed tree.
        print_stack = [root]
        l=[]
        while len(print_stack):
            node = print_stack.pop()
            l.append(node.name)
           # print(" " * node.get_depth(), node)
            print_stack.extend(node.children)
        return l

class AlignementMultiple:
       def __init__(self,tree,adn,methode):
           self.seqs= adn
           self.arbre=tree
           self.m=methode
           
            
       def tree_indexes(self):
           indexes=[]
           #colors = [clrs[i] for i in text]
           sequences={}
           for i,id in zip(range(len(self.seqs)),self.seqs.keys()):
               sequences[str(i)]=id
           for noeud in reversed(self.arbre): 
               if noeud!='':
                   indexes.append(sequences[noeud])
           return indexes
       
       def fichier(self):
           f = open("static/files/alignementKmers.fasta", "w")
           for indexe in self.tree_indexes():
               f.write('>')
               f.write(indexe)
               f.write('\n')
               f.write(str(self.seqs[indexe]))
               f.write('\n')
           f.close()



class BloomFilter:
   # Constructeur
    def __init__(self,capacity,errRate):
        # Dimensionnement du filtre
        self.k,self.N = self.bloom_size(capacity,errRate)
        # Bitfield de taille N
        self.bitfield = [0] * self.N                
    def add(self,key):
        print("key:",key)
        for n in self.makeHashs(key,self.k,self.N):
            #print("n:",n)
            self.bitfield[n] = 1
    def contains(self,key):
        hashs = self.makeHashs(key,self.k,self.N)
        return all([self.bitfield[n] == 1 for n in hashs])

    
    # Fabrique une liste de n hashs basés sur SHA1.
    # Chaque hash est un entier entre 0 et max-1
    # Cette fonction simule la création de n fonctions différentes (la base théorique du filtre # de Bloom) en renvoyant n valeurs "suffisamment indépendantes".
    def makeHashs(self,key,n,max):
        res = []
        m = hashlib.sha1()
        #hashlib.sha256("a".encode('utf-8')).hexdigest()
        
        # On sérialise la clef pour rendre la
        # fonction makeHashs polymorphique.
        s = str(key)
        for i in range(n):
            #m.update(s) 
            m.update(s.encode('utf-8'))
            # Et on re-hache la même valeur
            #print("s:",s)
            hash = int(m.hexdigest(),16) % max
            #print("hash:",hash)
            res.append(hash)
        print("res:",res)
        return res
    
    # Calcule un dimensionnement correct (nombre de hashs, taille du bitset) pour un filtre de Bloom
    def bloom_size(self,capacity,errRate):
        factor = math.log(1/(math.pow(2.0, math.log(2.0))))
        N = math.ceil((capacity * math.log(errRate))/factor)
        k = math.ceil(math.log(2.0) * N / capacity)
        return int(k),int(N)
