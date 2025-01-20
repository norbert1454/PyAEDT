"""
Y-junction circulator with a Y access topology on a ferrite substrate

This script builds a 3-port Y-junction ferrite circulator with a Y access
topology on a ferrite substrate. Due to limitations from the pyaedt API,
the following user inputs are necessary to setup de magnetic properties and
magnetic bias of the ferrite once the project is saved and the aedt Desktop
has been released:
    Define the magnetic permeability as non-linear in the ferrite's material property menu
    Define the magnetic saturation as $ferrite_Mr
    Define the magnetic losses as $ferrite_delta_H
    Define the magnetic losses measurement frequency as $ferrite_freq_delta_H
    Add the magnetic Bias to the ferrite along the Z-axis whose value is Hint_apm 

Keywords: **HFSS**, **ferrite**, **circulator**.

Created on Mon Jan 20 15:08:17 2025

@author: parker
"""

import numpy as np
import ansys.aedt.core
from ansys.aedt.core.application.variables import Variable
import ansys.aedt.core.modeler
from pathlib import Path

###############################
# Paramètres de la simulation #
###############################

# Paramètres du setup
setup_frequency = "18GHz" # freq_unit
sweep_start = "16GHz" # freq_unit
sweep_stop = "20GHz" # freq_unit
sweep_step = "0.05GHz" # freq_unit

max_passes = 30
max_delta_S = 0.02
percent_refinement = 20

############################
# Propriétés des matériaux #
############################

# Propriétés du ferrite
ferrite_epsilon = "20"
ferrite_tand = "0.005"
Hk = "18000"      # en Oe
Mr = "3500"       # en Gauss
# Le coefficient démagnétisant est déterminé à partir du modèle de Aharoni pour une plaquette rectangulaire directement dans HFSS
delta_H = "200Oe"
freq_delta_H = "40GHz"

########################
# Paramètres du design #
########################

# Dimensions du substrat
hauteur_substrat = "100um"
longueur_substrat_arriere = "2000um"
longueur_substrat_avant = "2000um"
longueur_substrat = "longueur_substrat_avant + longueur_substrat_arriere"
largeur_substrat = "4000um"
largeur_taper = "100um"

# Dimensions de la jonction Y
epaisseur_metallisation = "4um"
rayon_jonction = "950um"
longueur_adaptation = "600um"
largeur_adaptation = "48um"
largeur_50_Ohm = "48um"
rayon_courbure = "500um"
ecartement_ports = "1.6mm"

# Dimensions du ferrite
rayon_ferrite = rayon_jonction
hauteur_ferrite = hauteur_substrat

######################
# Paramètres de HFSS #
######################

aedt_version = "2024.2"
non_graphical = False
new_desktop = True

#############
# Fonctions #
#############

def Polder_Mu_eff(freq,Hk,Nz,Mr,dH,f_dH):
    gyro_ratio = 2.8e6 #Mz/Oe aka Gamma

    Hint = Hk-Nz*Mr
    
    damping = (gyro_ratio*dH)/(2*f_dH)
    
    w = freq
    w_0 = gyro_ratio * Hint
    w_m = gyro_ratio*Mr

    polder_mu = 1 + ((w_0+1j*damping*w)*w_m)/((w_0+1j*damping*w)**2-w**2)
    polder_kappa = (w*w_m)/((w_0+1j*damping*w)**2-w**2)
    polder_mu_eff = (polder_mu**2 - polder_kappa**2)/polder_mu

    return polder_mu_eff

##########################
# Initialisation de HFSS #
##########################

# Création du nom du projet
    # Pour obtenir le chemin du script
script_path = Path(__file__).resolve()
    # Pour obtenir le dossier dans lequel est le script
script_dir = script_path.parent
    # Définition du nom du projet            
project_name = "Circulateur en Y"
project_dir = "Designs/" + project_name

project = script_dir / project_dir / (project_name+".aedt")
project.parent.parent.mkdir(parents=True,exist_ok=True)
project.parent.mkdir(parents=True,exist_ok=True)

# Création du projet HFSS
Circulateur = ansys.aedt.core.Hfss(project = str(project),
                                     version = aedt_version,
                                     design = "Circulateur",
                                     non_graphical = non_graphical,
                                     new_desktop = new_desktop,
                                     solution_type = "Modal")

################################
# Initialisation des variables #
################################

    # Propriétés des matériaux
        # Ferrite
Circulateur["$ferrite_epsilon"] = ferrite_epsilon
Circulateur["$ferrite_tand"] = ferrite_tand
Circulateur["$ferrite_Mr"] = Mr+"Gauss"
Circulateur["$ferrite_delta_H"] = delta_H
Circulateur["$ferrite_freq_delta_H"] = freq_delta_H
    # Dimensions du substrat
Circulateur["hauteur_substrat"] = hauteur_substrat
Circulateur["longueur_substrat_arriere"] = longueur_substrat_arriere
Circulateur["longueur_substrat_avant"] = longueur_substrat_avant
Circulateur["longueur_substrat"] = longueur_substrat
Circulateur["largeur_substrat"] = largeur_substrat
Circulateur["largeur_taper"] = largeur_taper
    # Dimensions de la jonction Y
Circulateur["rayon_jonction"] = rayon_jonction
Circulateur["longueur_adaptation"] = longueur_adaptation
Circulateur["longueur_adaptation_1"] = "longueur_adaptation"
Circulateur["longueur_adaptation_2"] = "longueur_adaptation"
Circulateur["longueur_adaptation_3"] = "longueur_adaptation"
Circulateur["largeur_adaptation"] = largeur_adaptation
Circulateur["largeur_adaptation_1"] = "largeur_adaptation"
Circulateur["largeur_adaptation_2"] = "largeur_adaptation"
Circulateur["largeur_adaptation_3"] = "largeur_adaptation"
    # Propriétés métallisation
Circulateur["epaisseur_metallisation"] = epaisseur_metallisation
Circulateur["largeur_50_Ohm"] = largeur_50_Ohm
Circulateur["ecartement_ports"] = ecartement_ports
Circulateur["rayon_courbure"]  = rayon_courbure
Circulateur["rayon_courbure_2"]  = "rayon_courbure"
Circulateur["rayon_courbure_3"]  = "rayon_courbure"
Circulateur["longueur_50_Ohm_avant_courbe_2"] = "(-ecartement_ports/2-(rayon_jonction+longueur_adaptation_2)*sin(60deg)-rayon_courbure_2*(1-cos(60deg)))/sin(60deg)"
Circulateur["longueur_50_Ohm_avant_courbe_3"] = "(-ecartement_ports/2-(rayon_jonction+longueur_adaptation_3)*sin(60deg)-rayon_courbure_3*(1-cos(60deg)))/sin(60deg)"
# Circulateur["longueur_50_Ohm_avant_courbe_2"] = "(-ecartement_ports-(rayon_jonction+longueur_adaptation_2)*sin(30deg)-rayon_courbure_2*(1-cos(30deg)))/sin(30deg)"
# Circulateur["longueur_50_Ohm_avant_courbe_3"] = "(-ecartement_ports-(rayon_jonction+longueur_adaptation_3)*sin(30deg)-rayon_courbure_3*(1-cos(30deg)))/sin(30deg)"
Circulateur["longueur_50_Ohm_2"] = "largeur_substrat/2-(rayon_jonction+longueur_adaptation_2+longueur_50_Ohm_avant_courbe_2)*cos(30deg)-rayon_courbure_2*sin(30deg)"
Circulateur["longueur_50_Ohm_3"] = "largeur_substrat/2-(rayon_jonction+longueur_adaptation_3+longueur_50_Ohm_avant_courbe_3)*cos(30deg)-rayon_courbure_3*sin(30deg)"
    # Calcul de la dimension des ports
        # Equations from https://emtalk.com/waveport_calc.htm
w = Circulateur.variable_manager.decompose("largeur_50_Ohm")[0]
h = Circulateur.variable_manager.decompose("hauteur_substrat")[0]
if w/h < 1:
    largeur_port = 20*np.round(w)
    hauteur_port = np.round(h/w)*np.ceil(h)
else:
    largeur_port = np.round(w*(w/h)*(w/h))
    hauteur_port = np.round(h*w/h)
        # Dimensions des ports    
Circulateur["hauteur_port"] = str(hauteur_port)+"um"
Circulateur["largeur_port"] = str(largeur_port)+"um"
Circulateur["epaisseur_pec"] = "10um"

    # Dimensions du ferrite
Circulateur["rayon_ferrite"] = rayon_ferrite
Circulateur["hauteur_ferrite"] = hauteur_ferrite
    # Propriétés du ferrite
Circulateur["Hk"] = Hk
Circulateur["Mr"] = Mr
        # Calcul de Nz
# Nz = Aharoni(longueur_substrat,largeur_substrat,hauteur_substrat)
Circulateur["N1"] = "(largeur_substrat^2 - hauteur_substrat^2)/(2*largeur_substrat*hauteur_substrat)*ln((sqrt(longueur_substrat^2 + largeur_substrat^2 + hauteur_substrat^2) - longueur_substrat)/(sqrt(longueur_substrat^2 + largeur_substrat^2 + hauteur_substrat^2) + longueur_substrat))"
Circulateur["N2"] = "(longueur_substrat^2 - hauteur_substrat^2)/(2*longueur_substrat*hauteur_substrat)*ln((sqrt(longueur_substrat^2 + largeur_substrat^2 + hauteur_substrat^2) - largeur_substrat)/(sqrt(longueur_substrat^2 + largeur_substrat^2 + hauteur_substrat^2) + largeur_substrat))"
Circulateur["N3"] = "largeur_substrat/(2*hauteur_substrat)*ln((sqrt(longueur_substrat^2 + largeur_substrat^2) + longueur_substrat)/(sqrt(longueur_substrat^2 + largeur_substrat^2) - longueur_substrat))"
Circulateur["N4"] = "longueur_substrat/(2*hauteur_substrat)*ln((sqrt(longueur_substrat^2 + largeur_substrat^2) + largeur_substrat)/(sqrt(longueur_substrat^2 + largeur_substrat^2) - largeur_substrat))"
Circulateur["N5"] = "hauteur_substrat/(2*longueur_substrat)*ln((sqrt(largeur_substrat^2 + hauteur_substrat^2) - largeur_substrat)/(sqrt(largeur_substrat^2 + hauteur_substrat^2) + largeur_substrat))"
Circulateur["N6"] = "hauteur_substrat/(2*largeur_substrat)*ln((sqrt(longueur_substrat^2 + hauteur_substrat^2) - longueur_substrat)/(sqrt(longueur_substrat^2 + hauteur_substrat^2) + longueur_substrat))"
Circulateur["N7"] = "2*atan((longueur_substrat*largeur_substrat)/(hauteur_substrat*sqrt(longueur_substrat^2+largeur_substrat^2+hauteur_substrat^2)))"
Circulateur["N8"] = "(longueur_substrat^3+largeur_substrat^3-2*hauteur_substrat^3)/(3*longueur_substrat*largeur_substrat*hauteur_substrat)"
Circulateur["N9"] = "(longueur_substrat^2+largeur_substrat^2-2*hauteur_substrat^2)/(3*longueur_substrat*largeur_substrat*hauteur_substrat)*sqrt(longueur_substrat^2+largeur_substrat^2+hauteur_substrat^2)"
Circulateur["N10"] = "hauteur_substrat/(longueur_substrat*largeur_substrat)*(sqrt(longueur_substrat^2+hauteur_substrat^2)+sqrt(largeur_substrat^2+hauteur_substrat^2))"
Circulateur["N11"] = "((longueur_substrat^2+largeur_substrat^2)^(3/2)+(largeur_substrat^2+hauteur_substrat^2)^(3/2)+(hauteur_substrat^2+longueur_substrat^2)^(3/2))/(3*longueur_substrat*largeur_substrat*hauteur_substrat)"
Circulateur["Nz"] = "(N1+N2+N3+N4+N5+N6+N7+N8+N9+N10-N11)/pi"
        # Champ interne
Circulateur["Hint"] = "Hk-Nz*Mr"
Circulateur["Hint_apm"] = "1000*Hint/(4*pi)"

        # Taper autour du ferrite
var_setup_frequency = Variable(setup_frequency)
var_setup_frequency.rescale_to("Hz")

var_ferrite_freq_delta_H = Variable(freq_delta_H)
var_ferrite_freq_delta_H.rescale_to("Hz")

polder_mu_eff = Polder_Mu_eff(var_setup_frequency.numeric_value,
                              Circulateur.variable_manager.decompose("Hk")[0],
                              Circulateur.variable_manager.decompose("Nz")[0],
                              Circulateur.variable_manager.decompose("Mr")[0],
                              Circulateur.variable_manager.decompose("$ferrite_delta_H")[0],
                              var_ferrite_freq_delta_H.numeric_value)

Circulateur["$ferrite_mu_effective"] = np.real(polder_mu_eff)

#######################
# Ajout des matériaux #
#######################

# Création du dielectrique
Taper = Circulateur.materials.add_material("mon_ferrite_taper")
Taper.permittivity.value = "$ferrite_epsilon"
Taper.dielectric_loss_tangent.value = "$ferrite_tand"
Taper.permeability = "$ferrite_mu_effective"


# Création du ferrite
Ferrite = Circulateur.materials.add_material("mon_ferrite")
Ferrite.permittivity = "$ferrite_epsilon"
Ferrite.dielectric_loss_tangent = "$ferrite_tand"
#Ferrite.permeability.type = "nonlinear"
#Ferrite.magnetic_saturation = "$ferrite_Mr"

###############################
# Modélisation du circulateur #
###############################

# Création du plan de masse
GND = Circulateur.modeler.create_box(origin = ["-longueur_substrat_arriere", "-largeur_substrat/2", "0mm"],
                                       sizes = ["longueur_substrat","largeur_substrat","-epaisseur_metallisation"],
                                       name = "GND",
                                       material = "gold")
GND.material_appearance = True

# Création du substrat
Substrat = Circulateur.modeler.create_box(origin = ["-longueur_substrat_arriere+largeur_taper", "-largeur_substrat/2+largeur_taper", "0mm"],
                                       sizes = ["longueur_substrat-2*largeur_taper","largeur_substrat-2*largeur_taper","hauteur_substrat"],
                                       name = "Ferrite",
                                       material = "mon_ferrite")
Substrat.color = (64,64,64)

# Création du Taper
Taper = Circulateur.modeler.create_box(origin = ["-longueur_substrat_arriere", "-largeur_substrat/2", "0mm"],
                                         sizes = ["longueur_substrat","largeur_substrat","hauteur_substrat"],
                                         name = "Taper",
                                         material = "mon_ferrite_taper")

Circulateur.modeler.subtract(blank_list = "Taper",
                               tool_list = "Ferrite",
                               keep_originals = True)

Taper.color = (160,160,160)

# Métallisation supérieure
    # Création du résonateur
Jonction = Circulateur.modeler.create_cylinder(orientation = 'Z',
                                                 origin = [0, 0, "hauteur_substrat"],
                                                 radius = "rayon_jonction",
                                                 height = "epaisseur_metallisation",
                                                 name = "Jonction",
                                                 material = "gold")

    # Création des lignes d'adaptation
Ligne_adaptation_1 = Circulateur.modeler.create_box(origin = [0, "-largeur_adaptation_1/2", "hauteur_substrat"],
                                                      sizes = ["rayon_jonction+longueur_adaptation_1","largeur_adaptation_1","epaisseur_metallisation"],
                                                      name = "Ligne_adaptation_1",
                                                      material = "gold")

Ligne_adaptation_2 = Circulateur.modeler.create_box(origin = [0, "-largeur_adaptation_2/2", "hauteur_substrat"],
                                                      sizes = ["rayon_jonction+longueur_adaptation_2","largeur_adaptation_2","epaisseur_metallisation"],
                                                      name = "Ligne_adaptation_2",
                                                      material = "gold")
Ligne_adaptation_2.rotate(axis = 'Z',
                          angle = "120deg")

Ligne_adaptation_3 = Circulateur.modeler.create_box(origin = [0, "-largeur_adaptation_3/2", "hauteur_substrat"],
                                                      sizes = ["rayon_jonction+longueur_adaptation_3","largeur_adaptation_3","epaisseur_metallisation"],
                                                      name = "Ligne_adaptation_3",
                                                      material = "gold")
Ligne_adaptation_3.rotate(axis = 'Z',
                          angle = "240deg")

    # Création des lignes d'accès
        # Acces 1
            # Ligne simple
Ligne_50_Ohm_1 = Circulateur.modeler.create_box(origin = ["rayon_jonction+longueur_adaptation_1", "-largeur_50_Ohm/2", "hauteur_substrat"],
                                                  sizes = ["longueur_substrat_avant-rayon_jonction-longueur_adaptation_1","largeur_50_Ohm","epaisseur_metallisation"],
                                                  name = "Ligne_50_Ohm_1",
                                                  material = "gold")

        # Acces 2 & 3
            # Lignes avant les courbes
Ligne_50_Ohm_avant_courbe_2 = Circulateur.modeler.create_box(origin = ["rayon_jonction+longueur_adaptation_2", "-largeur_50_Ohm/2", "hauteur_substrat"],
                                                               sizes = ["longueur_50_Ohm_avant_courbe_2","largeur_50_Ohm","epaisseur_metallisation"],
                                                               name = "Ligne_50_Ohm_avant_courbe_2",
                                                               material = "gold")
Ligne_50_Ohm_avant_courbe_2.rotate(axis = 'Z',
                                   angle = "120deg")

Ligne_50_Ohm_avant_courbe_3 = Circulateur.modeler.create_box(origin = ["rayon_jonction+longueur_adaptation_3", "-largeur_50_Ohm/2", "hauteur_substrat"],
                                                               sizes = ["longueur_50_Ohm_avant_courbe_3","largeur_50_Ohm","epaisseur_metallisation"],
                                                               name = "Ligne_50_Ohm_avant_courbe_3",
                                                               material = "gold")
Ligne_50_Ohm_avant_courbe_3.rotate(axis = 'Z',
                                   angle = "240deg")

            # Courbes
Ligne_50_Ohm_Courbe_2 = Circulateur.modeler.create_equationbased_surface(x_uv = "rayon_courbure_2*cos(_v)-sin(30deg)*(rayon_jonction+longueur_adaptation_2+longueur_50_Ohm_avant_courbe_2)-cos(30deg)*rayon_courbure_2",
                                                                            y_uv = "rayon_courbure_2*sin(_v)+cos(30deg)*(rayon_jonction+longueur_adaptation_2+longueur_50_Ohm_avant_courbe_2)-sin(30deg)*rayon_courbure_2",
                                                                            z_uv = "_u",
                                                                            u_start = "hauteur_substrat",
                                                                            u_end = "hauteur_substrat+epaisseur_metallisation",
                                                                            v_start = "30deg",
                                                                            v_end = "60deg",
                                                                            name = "Ligne_50_Ohm_Courbe_2")

Circulateur.modeler.thicken_sheet(assignment = "Ligne_50_Ohm_Courbe_2",
                                    thickness = "epaisseur_metallisation")

Ligne_50_Ohm_Courbe_2.material_name = "gold"

Ligne_50_Ohm_Courbe_3 = Circulateur.modeler.create_equationbased_surface(x_uv = "rayon_courbure_3*cos(_v)-sin(30deg)*(rayon_jonction+longueur_adaptation_3+longueur_50_Ohm_avant_courbe_3)-cos(30deg)*rayon_courbure_3",
                                                                            y_uv = "rayon_courbure_3*sin(_v)-cos(30deg)*(rayon_jonction+longueur_adaptation_3+longueur_50_Ohm_avant_courbe_3)+sin(30deg)*rayon_courbure_3",
                                                                            z_uv = "_u",
                                                                            u_start = "hauteur_substrat",
                                                                            u_end = "hauteur_substrat+epaisseur_metallisation",
                                                                            v_start = "270deg",
                                                                            v_end = "330deg",
                                                                            name = "Ligne_50_Ohm_Courbe_3")

Circulateur.modeler.thicken_sheet(assignment = "Ligne_50_Ohm_Courbe_3",
                                    thickness = "epaisseur_metallisation")

Ligne_50_Ohm_Courbe_3.material_name = "gold"

            # Lignes entre la courbe et le bord du substrat
Ligne_50_Ohm_2 = Circulateur.modeler.create_box(origin = ["-longueur_arriere", "ecartement_ports/2-largeur_50_Ohm/2", "hauteur_substrat"],
                                                  sizes = ["longueur_50_Ohm_2","largeur_50_Ohm","epaisseur_metallisation"],
                                                  name = "Ligne_50_Ohm_2",
                                                  material = "gold")

Ligne_50_Ohm_3 = Circulateur.modeler.create_box(origin = ["-longueur_arriere", "-ecartement_ports/2-largeur_50_Ohm/2", "hauteur_substrat"],
                                                  sizes = ["longueur_50_Ohm_3","largeur_50_Ohm","epaisseur_metallisation"],
                                                  name = "Ligne_50_Ohm_3",
                                                  material = "gold")

    # Union des différents éléments
Circulateur.modeler.unite(["Jonction",
                             "Ligne_adaptation_1","Ligne_50_Ohm_1",
                             "Ligne_adaptation_2","Ligne_50_Ohm_avant_courbe_2","Ligne_50_Ohm_Courbe_2","Ligne_50_Ohm_2",
                             "Ligne_adaptation_3","Ligne_50_Ohm_avant_courbe_3","Ligne_50_Ohm_Courbe_3","Ligne_50_Ohm_3"])

Jonction.material_appearance = True

###################
# Ajout des ports #
###################

# Ajout du port 1
PEC_Port_1 = Circulateur.modeler.create_box(origin = ["longueur_substrat_avant", "-largeur_port/2", 0],
                                              sizes = ["epaisseur_pec","largeur_port","hauteur_port"],
                                              name = "PEC_Port_1",
                                              material = "pec")

Circulateur.wave_port(PEC_Port_1.bottom_face_x,
                        name = "1")

# Ajout du port 2
PEC_Port_2 = Circulateur.modeler.create_box(origin = ["-longueur_arriere", "ecartement_ports/2-largeur_port/2", 0],
                                              sizes = ["largeur_port","epaisseur_pec","hauteur_port"],
                                              name = "PEC_Port_2",
                                              material = "pec")

Circulateur.wave_port(PEC_Port_2.bottom_face_y,
                        name = "2")

# Ajout du port 3
PEC_Port_3 = Circulateur.modeler.create_box(origin = ["-longueur_arriere", "-ecartement_ports/2-largeur_port/2", 0],
                                              sizes = ["largeur_port","-epaisseur_pec","hauteur_port"],
                                              name = "PEC_Port_3",
                                              material = "pec")

Circulateur.wave_port(PEC_Port_3.top_face_y,
                        name = "3")

# Couleur des ports
PEC_Port_1.color = (255,0,255)
PEC_Port_2.color = (255,0,255)
PEC_Port_3.color = (255,0,255)

########################
# Paramétrage du setup #
########################

# Création de la boite d'air et des conditions aux limites
Circulateur.set_auto_open(enable = True)

# Analysis setup setting
Setup = Circulateur.setups[0]
    # Sweep setup
var_sweep_start = Variable(sweep_start)
var_sweep_stop = Variable(sweep_stop)
var_sweep_step = Variable(sweep_step)

var_sweep_stop.rescale_to(var_sweep_start.units)
var_sweep_step.rescale_to(var_sweep_start.units)

Setup.create_linear_step_sweep(name = "Sweep",
                               unit = var_sweep_start.units,
                               start_frequency = var_sweep_start.numeric_value,
                               stop_frequency = var_sweep_stop.numeric_value,
                               step_size = var_sweep_step.numeric_value,
                               sweep_type = "Interpolating")
    # Setup setup
Setup.properties["Name"] = "Setup"
Setup.properties["Solution Freq"] = setup_frequency
Setup.properties["Delta S"] = max_delta_S
Setup.properties["Passes"] = max_passes
Setup.properties["Percent Refinement"] = percent_refinement

###################
# Post processing #
###################

plot_tous_les_ports = Circulateur.post.create_report(expressions = ["dB(S(1,1))", "dB(S(2,1))", "dB(S(3,1))",
                                                                      "dB(S(2,2))", "dB(S(3,2))", "dB(S(1,2))",
                                                                      "dB(S(3,3))", "dB(S(1,3))", "dB(S(2,3))"])
plot_tous_les_ports.plot_name = "Tous les ports"
plot_tous_les_ports.hide_legend(solution_name = False,
                                trace_name = True,
                                variation_key = False,
                                font_size = 10)
plot_tous_les_ports.edit_y_axis_scaling(min_scale = -30,
                                        max_scale = 0)

plot_Port_1 = Circulateur.post.create_report(expressions = ["db(S11)", "db(S21)", "db(S31)"])
plot_Port_1.plot_name = "Port 1"
plot_Port_1.hide_legend(solution_name = True,
                        trace_name = True,
                        variation_key = False,
                        font_size = 10)
plot_Port_1.edit_y_axis_scaling(min_scale = -30,
                                max_scale = 0)

plot_Port_2 = Circulateur.post.create_report(expressions = ["dB(S(2,2))", "dB(S(3,2))", "dB(S(1,2))"])
plot_Port_2.plot_name = "Port 2"
plot_Port_2.hide_legend(solution_name = True,
                        trace_name = True,
                        variation_key = False,
                        font_size = 10)
plot_Port_2.edit_y_axis_scaling(min_scale = -30,
                                max_scale = 0)

plot_Port_3 = Circulateur.post.create_report(expressions = ["dB(S(3,3))", "dB(S(1,3))", "dB(S(2,3))"])
plot_Port_3.plot_name = "Port 3"
plot_Port_3.hide_legend(solution_name = True,
                        trace_name = True,
                        variation_key = False,
                        font_size = 10)
plot_Port_3.edit_y_axis_scaling(min_scale = -30,
                                max_scale = 0)

plot_Adaptation = Circulateur.post.create_report(expressions = ["dB(S(1,1))", "dB(S(2,2))", "dB(S(3,3))"])
plot_Adaptation.plot_name = "Adaptation"
plot_Adaptation.hide_legend(solution_name = True,
                            trace_name = True,
                            variation_key = False,
                            font_size = 10)
plot_Adaptation.edit_y_axis_scaling(min_scale = -30,
                                    max_scale = 0)

plot_Isolation = Circulateur.post.create_report(expressions = ["dB(S(2,1))", "dB(S(3,2))", "dB(S(1,3))"])
plot_Isolation.plot_name = "Isolation"
plot_Isolation.hide_legend(solution_name = True,
                           trace_name = True,
                           variation_key = False,
                           font_size = 10)
plot_Isolation.edit_y_axis_scaling(min_scale = -30,
                                   max_scale = 0)


plot_Transmission = Circulateur.post.create_report(expressions = ["dB(S(3,1))", "dB(S(1,2))", "dB(S(2,3))"])
plot_Transmission.plot_name = "Transmission"
plot_Transmission.hide_legend(solution_name = True,
                              trace_name = True,
                              variation_key = False,
                              font_size = 10)
plot_Transmission.edit_y_axis_scaling(min_scale = -30,
                                      max_scale = 0)


#  Circulateur.release_desktop()

