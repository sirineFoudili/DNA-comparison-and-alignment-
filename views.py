"""
Foudili Sirine
11.05.2020
Alignement de séquences d'ADN
"""

from flask import Flask, escape, request,render_template, redirect, url_for
#from flask_wtf import FlaskForm 
from flask import flash
# import tableau
#from flask import flask_table
from flask_table import Table, Col

import os
import numpy
#vérifier la sécurité des fichiers télécargés
from werkzeug.utils import secure_filename
#Lire le fichier fasta
from Bio import SeqIO
#générer le vecteur de comaraison
from itertools import combinations 
from itertools import product
#Needelman wunch alignement
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Align
#class des fonctions
from Alignement import Alignement
from Alignement import Phylogenie
from Alignement import Affichage
#Affichage
from bokeh.plotting import figure
from bokeh.embed import components






import os, io, random
import re
import string
import numpy as np

from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot
from bokeh.io import output_file, show
from bokeh.resources import CDN
import time






app = Flask(__name__)
app.secret_key = 'admin'
app.config.from_object('config')

#Télécharger fichier
app.config["FILE_UPLOADS"] = "static/files"
app.config["ALLOWED_FILE_EXTENSIONS"] = ["TXT", "FASTA"]




'''
Fichier fasta
'''
def fichier_autorise(filename):

    # vérifier la validité
    if not "." in filename:
        return False

    # extraire l'extention
    ext = filename.rsplit(".", 1)[1]

    # vérifier le type
    if ext.upper() in app.config["ALLOWED_FILE_EXTENSIONS"]:
        return True
    else:
        return False

def lecture_fichier(filename):
    dict={}
    print("lecture fichier")
    for seq_record in SeqIO.parse(filename, "fasta"):
        dict[seq_record.id]=str(seq_record.seq).replace('N','')
    return dict
def seqs_to_file(sequences):
    fich=open("static/files/Sequences.fasta", "w")
    fich.write(sequences)
    fich.close()
    return fich


#/nom de la page, pour une nouvelle page changer le nom de la fonction aussi

#@app.route('/', methods=["GET", "POST"])
@app.route('/')

def index():
    return render_template('index.html')



@app.route('/resultat', methods=["GET", "POST"])


def resultat():
    
    print("RESULTAT")
    '''          formulaire
    '''
    if request.method == "POST":
        
        
        #Récupérer méthode
        
        methode = request.form['methode']
        print("select:",str(methode))
        
        #Séquence saisies
        
        sequences=request.form['sequence']
        print("ADN:",sequences)
        
        #Récuperer fichier
        
        if request.files:
            fichier = request.files["fichier"]
            print('fichier :',fichier)
            #lecture_fichier(fichier)
            
        
        if (fichier.filename == "") and (sequences==""):
            print("No filename")
            flash("Fichier ou séquences introuvables", "danger")
            return render_template('index.html')
            
        if (fichier.filename != "") and not fichier_autorise(fichier.filename):
            flash("Le fichier doit etre sous format .fasta", "danger")
            print("fichier non autorisé")
            return render_template('index.html')
        else:
            if (fichier.filename == ""):
                s=seqs_to_file(sequences)
                f="static/files/Sequences.fasta"
                filename=''
            else:
                #par défaut 
                filename = secure_filename(fichier.filename)
                fichier.save(os.path.join(app.config["FILE_UPLOADS"],filename))
                f="static/files/"+filename 
            print('f:',f)
            d="static/files/AlignementNWk.fasta"
            a=lecture_fichier(f)
            if(a=={}):
                print("Non conforme")
                flash("le format des séquences n'est pas accepté", "danger")
                return render_template('index.html')
            if(len(a)<2):
                print("une seq")
                flash("Veuillez saisir au moins 2 séquences", "danger")
                return render_template('index.html')
            start_time = time.time()
            aligne=Alignement(int(methode),3,a)
            aligne.methode=methode
            matrice=aligne.scores_matrice.tolist()
            t=time.time() - start_time
            sequences=aligne.seqs
            print(sequences)
            phyl=Phylogenie(aligne.scores_matrice,aligne.seqs)
            noeuds=phyl.arbre
            #Alignement multiple
            if methode=='1':
                aln=list(SeqIO.parse(d,"fasta"))
            else:
                aln=list(SeqIO.parse(f,"fasta"))
            fig=Affichage()
            fig.methode=methode
            p = fig.view_alignment(aln,fontsize="9pt" ,plot_width=900)
            script, div = components(p)
            cdn_js=CDN.js_files[0]
            taille=len(sequences)
            if request.form['action'] == 'aligner':
                #test
                kwargs = {'script': script, 'div': div}
                return render_template('resultat.html', script=script,div=div,cdn_js=cdn_js,methode= aligne.methode,fichier=filename,matrice=matrice,sequences=sequences,taille=taille,noeuds=noeuds,temps=t)
            else:
                return render_template('resultat_arbre.html', noeuds=noeuds)
    return render_template('resultat.html', script=script,div=div)

@app.route('/resultat_details', methods=["GET", "POST"])
def resultat_details():
    print("details")
    methode = request.args.get('methode', None)
    matrice = request.args.getlist('matrice')
    seqs=request.args.getlist('sequences')
    if methode=='1':
        f=open("static/files/AlignementPaire.fasta","r")
        contenue=f.readlines()
        c=[]
        for line in contenue:
            line.replace('|', '')
            c.append(line.replace('\n',''))
    else:
        c=""
    return render_template('resultat_details.html',contenue=c,methode=methode,matrice=matrice,sequences=seqs)

@app.route('/resultat_arbre', methods=["GET", "POST"])
def resultat_arbre():
    noeuds = request.args.get('noeuds', None)
    return render_template('resultat_arbre.html',noeuds=noeuds)

@app.route('/scoreMatrice', methods=["GET", "POST"])

def scoreMatrice():
    print("matrice")
    methode = request.args.get('methode', None)
    mat = request.args.getlist('matrice')

    matrice=[]
    for vec in mat:
        v=list(re.split(",|[|]",vec))
        matrice.append(v)
    seqs=request.args.getlist('sequences')
    print("matrice:",matrice,"\nsequences:",seqs)
    for s, pos in zip(seqs, range(len(matrice))):
        matrice[pos].insert(0,s)
    print("scores:",matrice)

    return render_template('scoreMatrice.html',sequences=seqs,matrice=matrice)

    
if __name__ == '__main__':
    app.run() 
    


