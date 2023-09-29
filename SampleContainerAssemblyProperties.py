# a class representing a sample inside a simple container
from typing import Any, Dict, List, Tuple
import json
import copy

class SampleContainerAssemblyProperties:
    

    def __init__(self):

        self.sampleMaterial: Dict[str,Any] = {}
        self.sampleGeometry: Dict[str,Any] = {}
        self.containerMaterial: Dict[str,Any] = {}
        self.containerGeometry: Dict[str,Any] = {}
        self.containerInfo: Dict[str, Any] = {}
        self.containerPropertiesFilePath: str
    
    def showContents(self):
        print('Container Properties:')
        print(self.containerInfo)
        print(self.containerMaterial)
        print(self.containerGeometry)
        return

    def loadContainerFromFile(self, filePath):
        self.containerPropertiesFilePath = filePath

        with (open(self.containerPropertiesFilePath)) as json_file:
            data = json.load(json_file)
        
        # print(data[0])
        # print(data[1])
        # print(data[2])
        self.containerInfo = copy.copy(data[0])
        self.containerMaterial = copy.copy(data[1])
        self.containerGeometry = copy.copy(data[2])
        
        return self





