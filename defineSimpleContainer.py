# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
import json

containerPath = '/SNS/SNAP/shared/Calibration_dynamic/SimpleContainers/'

containers = ['kapton-capillary','single-toroid-TiZr-gasket','quartz-capillary']

for container in containers:

    if container == 'kapton-capillary':
        containerInfo={'name':'kapton capillary',
            'Manufacturer':'unknown',
            'DateCommissioned':'unknown'}
        containerMaterial={'ChemicalFormula': 'C22-H10-N2-O5',
            'MassDensity': 1.42}
        containerGeometry={'Shape':'HollowCylinder',
            'Height':0.3,
            'InnerRadius':0.1, #2mm diameter
            'OuterRadius':0.1013, #13 um wall thickness
            'Center':[0.,0.,0.]}
    elif container == 'single-toroid-TiZr-gasket':
    #PE Gasket
        containerInfo={'name':'single toroid gasket TiZr',
            'Manufacturer':'unknown',
            'DateCommissioned':'unknown'}
        containerMaterial={'ChemicalFormula': 'Ti0.6765-Zr0.324',
             'NumberDensity': 0.100}
        containerGeometry={'Shape':'HollowCylinder',
             'Height':0.3,
             'InnerRadius':0.3,
             'OuterRadius':0.7,
             'Center':[0.,0.,0.]}
    elif container == 'quartz-capillary':
        containerInfo={'name':'quartz capillary',
            'Manufacturer':'unknown',
            'DateCommissioned':'unknown'}
        containerMaterial={'ChemicalFormula': 'SiO2',
             'MassDensity': 2.65}
        containerGeometry={'Shape':'HollowCylinder',
             'Height':0.3,
             'InnerRadius':0.1,
             'OuterRadius':0.12, #2mm 200um thick container
             'Center':[0.,0.,0.]}
             
    details = [containerInfo,containerMaterial,containerGeometry]
    fname = containerPath + container + '.json'
    with open(fname, "w") as output:
        json.dump(details, output, indent = 2)
        
       
