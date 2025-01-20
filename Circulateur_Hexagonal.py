"""
Y-junction circulator with a Y access topology on a ferrite substrate

This script builds a 3-port Y-junction ferrite circulator with a Y access
topology on a ferrite substrate. Due to limitations from the pyaedt API,
the following user inputs are necessary to setup de magnetic properties and
magnetic bias of the ferrite:
    Define the magnetic permeability as non-linear in the ferrite's material property menu
    Define the magnetic saturation as $ferrite_Mr
    Add the magnetic Bias to the ferrite along the Z-axis whose value is Hint_apm     

Keywords: **HFSS**, **ferrite**, **circulator**.

Created on Mon Jan 13 10:31:37 2025

@author: parker
"""

import numpy as np
from math import pi
import ansys.aedt.core
from pathlib import Path
from ansys.aedt.core.application.variables import Variable

###############################
# Paramètres de la simulation #
###############################

# Paramètres du setup
setup_frequency = "24GHz"
sweep_start = "20GHz"
sweep_stop = "28GHz"
sweep_step = "0.05GHz"

max_passes = 30
max_delta_S = 0.02
percent_refinement = 20

############################
# Propriétés des matériaux #
############################

# Propriétés du diélectrique
dielectrique_epsilon = "21"
dielectrique_tand = "0.003"

# Propriétés du ferrite
ferrite_epsilon = "20"
ferrite_tand = "0.005"
Hk = "14000"      # en Oe
Mr = "3500"       # en Gauss
Nzm_Chen_filename = "Nzm_Chen.tab" # Le coefficient démagnétisant est obtenu directement depuis les valeurs déterminées par Chen
delta_H = "200Oe"
freq_delta_H = "40GHz"

########################
# Paramètres du design #
########################

# Dimensions du substrat
hauteur_substrat = "352um"
longueur_substrat = "4000um"

# Dimensions de la jonction Y
epaisseur_metallisation = "4um"
rayon_jonction = "1100um"
longueur_adaptation = "1600um"
largeur_adaptation = "720um"
largeur_50_Ohm = "115um"
longueur_50_Ohm_min = "1000um"

# Dimensions du ferrite
rayon_ferrite = rayon_jonction
hauteur_ferrite = hauteur_substrat

######################
# Paramètres de HFSS #
######################

aedt_version = "2024.2"
non_graphical = False
new_desktop = True

##########################
# Initialisation de HFSS #
##########################

# Création du nom du projet
    # Pour obtenir le chemin du script
script_path = Path(__file__).resolve()
    # Pour obtenir le dossier dans lequel est le script
script_dir = script_path.parent
    # Définition du nom du projet            
project_name = "Circulateur Hexagonal"
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
        # Dielectrique
Circulateur["$dielectrique_epsilon"] = dielectrique_epsilon
Circulateur["$dielectrique_tand"] = dielectrique_tand
        # Ferrite
Circulateur["$ferrite_epsilon"] = ferrite_epsilon
Circulateur["$ferrite_tand"] = ferrite_tand
Circulateur["$ferrite_Mr"] = Mr+"Gauss"
Circulateur["$ferrite_delta_H"] = delta_H
Circulateur["$ferrite_freq_delta_H"] = freq_delta_H
    # Propriétés métallisation
    # Dimensions de la jonction Y
Circulateur["epaisseur_metallisation"] = epaisseur_metallisation
Circulateur["rayon_jonction"] = rayon_jonction
Circulateur["longueur_adaptation"] = longueur_adaptation
Circulateur["largeur_adaptation"] = largeur_adaptation
Circulateur["largeur_50_Ohm"] = largeur_50_Ohm
Circulateur["longueur_50_Ohm_min"] = longueur_50_Ohm_min
    # Dimensions du substrat
Circulateur["hauteur_substrat"] = hauteur_substrat
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
Circulateur["longueur_substrat"] = str(np.max([1.5*largeur_port/(2*np.tan(30*pi/180)),Circulateur.variable_manager.decompose("rayon_jonction")[0]+Circulateur.variable_manager.decompose("longueur_adaptation")[0]+Circulateur.variable_manager.decompose("longueur_50_Ohm_min")[0]]))+"um"
    # Dimensions du ferrite
Circulateur["rayon_ferrite"] = rayon_ferrite
Circulateur["hauteur_ferrite"] = hauteur_ferrite
    # Propriétés du ferrite
Circulateur["Hk"] = Hk
Circulateur["Mr"] = Mr
        # Obtention du Nz à partir du dataset Nzm_Chen qui doit se trouver dans le même dossier que ce script
            # Pour obtenir le chemin du script
script_path = Path(__file__).resolve()
            # Pour obtenir le dossier dans lequel est le script
script_dir = script_path.parent
Nzm_Chen_path = script_dir / Nzm_Chen_filename
Circulateur["Nzm_Chen"] = Circulateur.import_dataset1d(input_file = str(Nzm_Chen_path), is_project_dataset = False)
Circulateur["gamma_Chen"] = "hauteur_ferrite/(2*rayon_ferrite)"
Circulateur["Nz"] = "pwl(Nzm_Chen,gamma_Chen)"
        # Champ interne
Circulateur["Hint"] = "Hk-Nz*Mr"
Circulateur["Hint_apm"] = "1000*Hint/(4*pi)"

    # Dimensions des ports    
Circulateur["hauteur_port"] = str(hauteur_port)+"um"
Circulateur["largeur_port"] = str(largeur_port)+"um"
Circulateur["epaisseur_pec"] = "10um"

#######################
# Ajout des matériaux #
#######################

# Création du dielectrique
Dielectrique = Circulateur.materials.add_material("mon_dielectrique")
Dielectrique.permittivity.value = "$dielectrique_epsilon"
Dielectrique.dielectric_loss_tangent.value = "$dielectrique_tand"

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
GND = Circulateur.modeler.create_cylinder(orientation = 'Z',
                                                    origin = [0, 0, 0],
                                                    radius = "longueur_substrat/cos(30deg)",
                                                    height = "-epaisseur_metallisation",
                                                    num_sides = 6,
                                                    name = "GND",
                                                    material = "gold")
GND.rotate(axis = "Z",
           angle = "30deg")

# Création du substrat
Substrat = Circulateur.modeler.create_cylinder(orientation = 'Z',
                                                         origin = [0, 0, 0],
                                                         radius = "longueur_substrat/cos(30deg)",
                                                         height = "hauteur_substrat",
                                                         num_sides = 6,
                                                         name = "Substrat",
                                                         material = "mon_dielectrique")

Substrat.rotate(axis = "Z",
                angle = "30deg")

# Création de l'insert de ferrite
Ferrite = Circulateur.modeler.create_cylinder(orientation = 'Z',
                                                        origin = [0, 0, 0],
                                                        radius = "rayon_ferrite",
                                                        height = "hauteur_substrat",
                                                        name = "Ferrite",
                                                        material = "mon_ferrite")

# Création de la cavité dans le substrat
Circulateur.modeler.subtract(blank_list = "Substrat",
                                        tool_list = "Ferrite",
                                        keep_originals = True)

# Métallisation supérieure
    # Création du résonateur
Jonction = Circulateur.modeler.create_cylinder(orientation = 'Z',
                                                         origin = [0, 0, "hauteur_substrat"],
                                                         radius = "rayon_jonction",
                                                         height = "epaisseur_metallisation",
                                                         name = "Jonction",
                                                         material = "gold")

    # Création des lignes d'adaptation
Ligne_adaptation_1 = Circulateur.modeler.create_box(origin = [0, "-largeur_adaptation/2", "hauteur_substrat"],
                                                              sizes = ["rayon_jonction+longueur_adaptation","largeur_adaptation","epaisseur_metallisation"],
                                                              name = "Ligne_adaptation_1",
                                                              material = "gold")

Ligne_adaptation_2 = Circulateur.modeler.create_box(origin = [0, "-largeur_adaptation/2", "hauteur_substrat"],
                                                              sizes = ["rayon_jonction+longueur_adaptation","largeur_adaptation","epaisseur_metallisation"],
                                                              name = "Ligne_adaptation_2",
                                                              material = "gold")
Ligne_adaptation_2.rotate(axis = 'Z',
                          angle = "120deg")

Ligne_adaptation_3 = Circulateur.modeler.create_box(origin = [0, "-largeur_adaptation/2", "hauteur_substrat"],
                                                              sizes = ["rayon_jonction+longueur_adaptation","largeur_adaptation","epaisseur_metallisation"],
                                                              name = "Ligne_adaptation_3",
                                                              material = "gold")
Ligne_adaptation_3.rotate(axis = 'Z',
                          angle = "240deg")

    # Création des lignes d'accès
Ligne_50_Ohm_1 = Circulateur.modeler.create_box(origin = ["rayon_jonction+longueur_adaptation", "-largeur_50_Ohm/2", "hauteur_substrat"],
                                                          sizes = ["longueur_substrat-rayon_jonction-longueur_adaptation","largeur_50_Ohm","epaisseur_metallisation"],
                                                          name = "Ligne_50_Ohm_1",
                                                          material = "gold")

Ligne_50_Ohm_2 = Circulateur.modeler.create_box(origin = ["rayon_jonction+longueur_adaptation", "-largeur_50_Ohm/2", "hauteur_substrat"],
                                                          sizes = ["longueur_substrat-rayon_jonction-longueur_adaptation","largeur_50_Ohm","epaisseur_metallisation"],
                                                          name = "Ligne_50_Ohm_2",
                                                          material = "gold")
Ligne_50_Ohm_2.rotate(axis = 'Z',
                      angle = "120deg")

Ligne_50_Ohm_3 = Circulateur.modeler.create_box(origin = ["rayon_jonction+longueur_adaptation", "-largeur_50_Ohm/2", "hauteur_substrat"],
                                                          sizes = ["longueur_substrat-rayon_jonction-longueur_adaptation","largeur_50_Ohm","epaisseur_metallisation"],
                                                          name = "Ligne_50_Ohm_3",
                                                          material = "gold")
Ligne_50_Ohm_3.rotate(axis = 'Z',
                      angle = "240deg")

    # Union des différents éléments
Circulateur.modeler.unite(["Jonction","Ligne_adaptation_1","Ligne_adaptation_2","Ligne_adaptation_3","Ligne_50_Ohm_1","Ligne_50_Ohm_2","Ligne_50_Ohm_3"])

###################
# Ajout des ports #
###################

# Ajout du port 1
PEC_Port_1 = Circulateur.modeler.create_box(origin = ["longueur_substrat", "-largeur_port/2", 0],
                                                      sizes = ["epaisseur_pec","largeur_port","hauteur_port"],
                                                      name = "PEC_Port_1",
                                                      material = "pec")
Circulateur.wave_port(PEC_Port_1.bottom_face_x,
                                name = "1")

# Ajout du port 2
PEC_Port_2 = Circulateur.modeler.create_box(origin = ["longueur_substrat", "-largeur_port/2", 0],
                                                      sizes = ["epaisseur_pec","largeur_port","hauteur_port"],
                                                      name = "PEC_Port_2",
                                                      material = "pec")
Circulateur.wave_port(PEC_Port_2.bottom_face_x,
                                name = "2")
PEC_Port_2.rotate(axis = 'Z',
                  angle = "120deg")

# Ajout du port 3
PEC_Port_3 = Circulateur.modeler.create_box(origin = ["longueur_substrat", "-largeur_port/2", 0],
                                                      sizes = ["epaisseur_pec","largeur_port","hauteur_port"],
                                                      name = "PEC_Port_3",
                                                      material = "pec")
Circulateur.wave_port(PEC_Port_3.bottom_face_x,
                                name = "3")
PEC_Port_3.rotate(axis = 'Z',
                  angle = "240deg")

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




#Circulateur.release_desktop()