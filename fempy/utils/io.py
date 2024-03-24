import os
import sys
import fempy

from ROOT import TFile, TDirectoryFile, TList


def Load(container, path):
    '''
    Function to extract an object inside a root file.
    Supports nested containers with the following Data Types:
     - TDirectoryFile (TFile)
     - TList

    Parameters
    -----------
    container: TFile of the input file
    path: path of the object inside the root file

    Returns:
    -----------
    obj: target root object
    '''

    # Check that the input file is OK
    path = os.path.normpath(path)
    if container == None:  # pylint: disable=singleton-comparison
        fempy.logger.critical('The container %s is NULL', container.GetName())

    # Start to extract
    for name in path.split(os.sep):
        fempy.logger.debug('Trying to load %s:%s. Available keys:', container.GetName(), name)

        for key in GetKeyNames(container):
            fempy.logger.debug('    %s', key)

        if isinstance(container, TDirectoryFile):
            obj = container.Get(name)
        elif isinstance(container, TList):
            obj = container.FindObject(name)
        else:
            fempy.logger.critical('The container %s of type %s is not valid', container.GetName(), type(container))

        if obj == None:  # pylint: disable=singleton-comparison
            fempy.logger.critical('The container %s does not contain an object named %s', container.GetName(), name)
        container = obj

    fempy.logger.debug('The object %s:%s was succesfully loaded', container.GetName(), path)
    return obj


def GetObjectFromFile(inFile, pathToObj):
    '''
    Function to extract an object inside a root file.
    Supports nested containers with the following Data Types:
     - TFile
     - TDirecotryFile
     - TList

    Parameters
    -----------
    inFile: TFile of the input file
    pathToObj: path of the object inside the root file

    Returns:
    -----------
    outObj: target root object
    '''

    pathToObj = os.path.normpath(pathToObj)
    pathElements = pathToObj.split(os.sep)
    outObj = inFile.Get(pathElements.pop(0))

    for _, containerName in enumerate(pathElements):
        if isinstance(outObj, (TFile, TDirectoryFile)):
            outObj = outObj.Get(containerName)
        elif isinstance(outObj, TList):
            outObj = outObj.FindObject(containerName)
        else:
            print(f'\033[31mError\033[0m: instance of {type(outObj)} not implemented. Exit!')
            sys.exit()

    return outObj


histClasses = [f'TH{n}{t}' for n in [1, 2, 3] for t in ['I', 'F', 'D']]
graphClasses = ['TGraph', 'TGraphErrors', 'TGraphAsymErrors']


def GetHistNamesInDir(rdir):
    return [key.GetName() for key in list(rdir.GetListOfKeys()) if key.GetClassName() in histClasses]


def GetHistsInDir(rdir):
    return [rdir.Get(key.GetName()) for key in list(rdir.GetListOfKeys()) if key.GetClassName() in histClasses]


def GetRegions(rdir):
    return [key.GetName() for key in list(rdir.GetListOfKeys()) if key.GetName() in ['sgn', 'sbl', 'sbr']]


def GetGraphsInDir(rdir):
    return [rdir.Get(key.GetName()) for key in list(rdir.GetListOfKeys()) if key.GetClassName() in graphClasses]


def GetObjNamesInDir(rdir):
    return [key.GetName() for key in list(rdir.GetListOfKeys()) if key.GetClassName() != 'TDirectoryFile']


def GetObjsInDir(rdir):
    return [rdir.Get(key.GetName()) for key in list(rdir.GetListOfKeys()) if key.GetClassName() != 'TDirectoryFile']


def GetSubdirsInDir(rdir):
    if rdir == None:
        return []
    return [key.GetName() for key in list(rdir.GetListOfKeys()) if key.GetClassName() == 'TDirectoryFile']


def GetCombs(rdir):
    if rdir == None:
        return []
    return [key.GetName() for key in list(rdir.GetListOfKeys()) if key.GetName() in ['sc', 'oc', 'pp', 'mm', 'pm', 'mp']]


def GetKeyNames(container):  # pylint: disable=inconsistent-return-statements
    if isinstance(container, TDirectoryFile):
        return [key.GetName() for key in list(container.GetListOfKeys())]

    if isinstance(container, TList):
        it = container.MakeIterator()
        names = []
        while True:
            obj = it.Next()
            if obj == None: # pylint: disable=singleton-comparison
                break
            names.append(obj.GetName())

        return names
    fempy.logger.critical('Unknown container type %s', type(container))


def GetKeyNamesInDir(rdir):
    return [key.GetName() for key in list(rdir.GetListOfKeys())]


def GetPairName(rdir):
    pairs = []
    for keyName in GetKeyNamesInDir(rdir):
        if 'HM_CharmFemto' in keyName and 'Results' in keyName:
            pairs.append(keyName[14:keyName.find('_Results')])
    return list(dict.fromkeys(pairs))
