# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
import json

import importlib, sys
importlib.reload(sys.modules['SampleContainerAssemblyProperties'])
from SampleContainerAssemblyProperties import SampleContainerAssemblyProperties


contFile = '/SNS/SNAP/shared/Calibration_dynamic/SimpleContainers/quartz-capillary.json'
quartz = SampleContainerAssemblyProperties()
quartz.loadContainerFromFile(contFile)
quartz.showContents()

