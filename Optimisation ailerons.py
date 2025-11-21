"""
----- READ ME -----

V1.0 - 18/11/2025
Basé sur StabTraj v3-4_5

Auteurs : Clément Gissler (clement.gissler@gmail.com)
          Tommaso Bogoni (tom.bogoni@gmail.com)

/!\ Mise en garde /!\
Ce programme ne remplace en aucun cas le StabTraj et peut
devenir obsolète en cas de mise à jour majeure.
Il est donc recommandé de l'utiliser en parallèle de StabTraj
et de vérifier la cohérence des résultats.
À noter que dans cette version, seules les fusées mono-diamètre
et mono-empennage sont prises en charge.

Ce programme est la propriété de ses auteurs et ne doit en aucun
cas être diffusé en dehors de l'Ensma Space Project. Pour toute
remarque ou suggestion, merci de nous contacter par mail.
-------------------
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

#Paramètres à renseigner
masse = ...                 #en g SANS PROPULSEUR
centre_de_masse = ...       #en mm SANS PROPULSEUR
longueur_totale = ...       #en mm
diametre = ...              #en mm
hauteur_coiffe = ...        #en mm
nombre_ailerons = ...

type_fusee = ...            #parmi "exp", "mini"
type_moteur = ...           #parmi "Barasinga", "Pandora"
type_ogive = ...            #parmi "Parabolique", "Ogivale", "Conique"

MS_nominale = ...           #point de fonctionnement souhaité
Portance_nominale = ...     #cf StabTraj

depassement = ...           #dépassement minimal en mm des ailerons par rapport au bas de la fusée. 0 si sans importance 
plage_emplanture = ...      #plages de test
plage_saumon = ...          #[min, max]
plage_fleche = ...          #[x,x] pour fixer à x
plage_envergure = ...       #ne pas commencer à 0
pas_simulation = ...        #précision sur chaque paramètre géométrique

#Vous pouvez lancer le script

#------------------------------------------------
#Dictionnaire des caractéristiques moteurs [MpropuPlein, MpropuVide, XpropuPlein, XpropuVide, Longueur]
dic_moteurs = dict()
dic_moteurs["Barasinga"] = [1685, 652, 250, 240, 488]
dic_moteurs["Pandora"] = [159.9, 84.3, 114, 114, 288]

#Dictionnaire des caractéristiques des ogives
dic_ogives = dict()
dic_ogives["Parabolique"] = 1/2
dic_ogives["Ogivale"] = 7/15
dic_ogives["Conique"] = 2/3

#Dictionnaires des conditions de stabilité [Finessmin, Finessemax, Portancemin, Portancemax, MSmin, MSmax, Couplemin, Couplemax]
conditions = dict()
conditions["exp"] = [10, 35, 15, 40, 2, 6, 40, 100]
conditions["mini"] = [10, 20, 15, 30, 1.5, 6, 30, 100]
#------------------------------------------------

def calcul_stabilite(emplanture, saumon, fleche, envergure):
    """
    Fonction permettant d'accéder aux différents coefficients de Finesse, Portance, MargeStat. et Couple du fichier StrabTraj.
    Attention : valable uniquement pour une fusée mono-diamètre et mono-empennage.
    """
    
    #Calcul de la Finesse
    Finesse = longueur_totale/diametre
   
    #Calcul de la Portance
    Portance = (4*nombre_ailerons*(envergure/diametre)**2*(1+(diametre/(2*envergure+diametre))))/(1+np.sqrt(1+(2*np.sqrt((fleche+saumon/2-emplanture/2)**2+envergure**2)/(emplanture+saumon))**2))+2
   
    #Calcul de la Marge Statique
    Marge_Statique_Min = 1/diametre*(((Portance-2)*(longueur_totale-emplanture+(fleche*(emplanture+2*saumon)/(3*(emplanture+saumon)))+1/6*(emplanture+saumon-((emplanture*saumon)/(saumon+emplanture))))+2*dic_ogives[type_ogive]*hauteur_coiffe)/Portance-(centre_de_masse*masse+(longueur_totale-dic_moteurs[type_moteur][4]+dic_moteurs[type_moteur][2])*dic_moteurs[type_moteur][0])/(masse+dic_moteurs[type_moteur][0]))          
    Marge_Statique_Max = 1/diametre*(((Portance-2)*(longueur_totale-emplanture+(fleche*(emplanture+2*saumon)/(3*(emplanture+saumon)))+1/6*(emplanture+saumon-((emplanture*saumon)/(saumon+emplanture))))+2*dic_ogives[type_ogive]*hauteur_coiffe)/Portance-(centre_de_masse*masse+(longueur_totale-dic_moteurs[type_moteur][4]+dic_moteurs[type_moteur][3])*dic_moteurs[type_moteur][1])/(masse+dic_moteurs[type_moteur][1]))
   
    #Calcul du Couple
    Couple_Min = Marge_Statique_Min*Portance
    Couple_Max = Marge_Statique_Max*Portance
   
    return Finesse, Portance, Marge_Statique_Min, Marge_Statique_Max, Couple_Min, Couple_Max

def illustration_aileron(geometrie, coefficients):
    """
    Fonction permettant d'afficher la forme des ailerons par rapport à la fusée
    """
    
    #Géométrie
    emplanture, saumon, fleche, envergure = geometrie
    
    #Coefficients
    finesse = round(coefficients[0],2)
    portance = round(coefficients[1],2)
    MSmin = round(coefficients[2],2)
    MSmax = round(coefficients[3],2)
    Cmin = round(coefficients[4],2)
    Cmax = round(coefficients[5],2)
    erreurMS = round((MSmin+MSmax)/2 - MS_nominale,2)
    erreurPortance = round(portance - Portance_nominale,2)
    
    #Création de la figure et des axes
    fig, ax = plt.subplots(figsize=(6, 8))
    
    #Corps
    x_corps = -diametre / 2
    y_corps = 0
    corps = Rectangle((x_corps, y_corps), diametre, longueur_totale - hauteur_coiffe, color='lightgray', ec='black')
    
    #Ogive arrondie
    theta = np.linspace(0, np.pi, 50)
    r = diametre / 2
    x_ogive = r * np.cos(theta)
    y_ogive = r * np.sin(theta) * hauteur_coiffe / r + longueur_totale-hauteur_coiffe
    ogive = Polygon(
        np.column_stack([x_ogive, y_ogive]),
        closed=True, color='lightgray', ec='black'
    )
    
    #Ailerons
    aileron_gauche = Polygon([
        (x_corps, y_corps),
        (x_corps - envergure, y_corps - (fleche + saumon - emplanture)),
        (x_corps - envergure, y_corps + emplanture - fleche),
        (x_corps, y_corps + emplanture)
    ], closed=True, color='lightgray', ec='black')
    
    aileron_droit = Polygon([
        (x_corps + diametre, y_corps),
        (x_corps + diametre + envergure, y_corps - (fleche + saumon - emplanture)),
        (x_corps + diametre + envergure, y_corps + emplanture - fleche),
        (x_corps + diametre, y_corps + emplanture)
    ], closed=True, color='lightgray', ec='black')
    
    #Ajout des formes
    ax.add_patch(corps)
    ax.add_patch(ogive)
    ax.add_patch(aileron_gauche)
    ax.add_patch(aileron_droit)
    
    #Légende
    texte_leg = (
        f"Masse : {masse} g\n"
        f"Cdm : {centre_de_masse} mm\n"
        f"Hauteur : {longueur_totale} mm\n"
        f"Diamètre : {diametre} mm\n"
        f"Emplanture : {emplanture} mm\n"
        f"Saumon : {saumon} mm\n"
        f"Flèche : {fleche} mm\n"
        f"Envergure : {envergure} mm\n"
        "\n"
        f"Finesse : {finesse}\n"
        f"Portance : {portance}\n"
        f"MSmin : {MSmin}\n"
        f"MSmax : {MSmax}\n"
        f"Couplemin : {Cmin}\n"
        f"Couplemax : {Cmax}\n"
        f"ΔMS : {erreurMS}\n"
        f"ΔPortance : {erreurPortance}"
    )
    
    #Réglages d'affichage
    ax.set_ylim(-1.1 * (fleche + saumon - emplanture), (longueur_totale - hauteur_coiffe) * 1.15)
    ax.set_xlim(-(envergure + diametre) * 1.5, (envergure + diametre) * 1.5)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.text(0.7, 0.775, texte_leg, transform=ax.transAxes, fontsize=10, va='top', ha='left',
            bbox=dict(facecolor='white', edgecolor='None'))
    plt.show()
    print()
    print(texte_leg)

#Calcul de la distance euclidienne
def distance(X1, Y1, X2, Y2):
    return (X1-X2)**2 + (Y1-Y2)**2

def boucle_optimisation(plage_emplanture, plage_saumon, plage_fleche, plage_envergure):
    """
    Fonction qui teste toutes les géométries possibles des ailerons et renvoie celui qui se rapproche le plus du point nominal de stabilité
    """
    
    #Listes des géométries
    emplantures = np.arange(plage_emplanture[0], plage_emplanture[-1] + pas_simulation, pas_simulation)
    saumons = np.arange(plage_saumon[0], plage_saumon[-1] + pas_simulation, pas_simulation)
    fleches = np.arange(plage_fleche[0], plage_fleche[-1] + pas_simulation, pas_simulation)
    envergures = np.arange(plage_envergure[0], plage_envergure[-1] + pas_simulation, pas_simulation)
    
    #Initialisation de la recherche
    distance_min = float('inf')
    geometrie = []
    coefficients = []
    
    #Taux d'avancement de la recherche
    N_total = len(emplantures)*len(saumons)*len(fleches)*len(envergures)
    num_aileron = 0
    avancement = -1
    
    #Boucle sur la géométrie des ailerons
    for emplanture in emplantures:
        for saumon in saumons:
            for fleche in fleches:
                for envergure in envergures:
                    
                    #Calcul des coefficients de stabilité
                    Calculs = calcul_stabilite(emplanture, saumon, fleche, envergure)
                    MS_moy = (Calculs[2]+Calculs[3])/2
                    
                    #Conditions de stabilité
                    if (conditions[type_fusee][0] <= Calculs[0] <= conditions[type_fusee][1]) and (
                    conditions[type_fusee][2] <= Calculs[1] <= conditions[type_fusee][3]) and (
                    conditions[type_fusee][4] <= Calculs[2] <= conditions[type_fusee][5]) and (
                    conditions[type_fusee][4] <= Calculs[3] <= conditions[type_fusee][5]) and (
                    conditions[type_fusee][6] <= Calculs[4] <= conditions[type_fusee][7]) and (
                    conditions[type_fusee][6] <= Calculs[5] <= conditions[type_fusee][7]):
                        
                        dist = distance(MS_moy, Calculs[1], MS_nominale, Portance_nominale)
                        if dist < distance_min:                                              #aileron plus optimisé                            
                            if not depassement or (fleche+saumon>=emplanture+depassement):   #depassement (ou non)
                                distance_min = dist
                                geometrie = [emplanture, saumon, fleche, envergure]
                                coefficients = Calculs
                                
                    #Avancement
                    num_aileron += 1
                    taux = 100*num_aileron/N_total
                    if taux>=avancement+1:
                        avancement += 1
                        print(round(taux,2),"%")
    
    #Résultats              
    if distance_min == float('inf'):
        return 'Aucune configuration stable'
    return illustration_aileron(geometrie, coefficients)

boucle_optimisation(plage_emplanture, plage_saumon, plage_fleche, plage_envergure)
