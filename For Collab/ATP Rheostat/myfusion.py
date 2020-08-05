from biocrnpyler import *
from subsbml import *
import numpy as np
import pandas as pd
import libsbml


def triggeredFusion ( listOfSubsystemTrigger, listOfSubsystemFusion, specietrigger, operator, conctrigger, timepoints, **kwargs):
    '''
    Compartments for listOfSubsystemTrigger should be correctly specified
    specietrigger is the specie_id for which the conctrigger is checked
    Simulates listOfSubsystemTrigger to get the time for conctrigger to be achieved 
    When the point is achieved it fuses the listOfSubsystemFusion
    Returns the sbml document of fusedcell result of simulation for listOfSubsystemFusion and the updated timepoints when the conctrigger is achieved and
    also returns the result of simulation for the fusedcell and the remaining timepoints just after the conctrigger is achieved.
    '''
    
    verbose = kwargs.get('verbose') if 'verbose' in kwargs else None
    
    if type(listOfSubsystemTrigger) is not list:
        raise ValueError('listOfSubsystemTrigger argument should be list of Subsystem objects')
        
    if type(listOfSubsystemFusion) is not list:
        raise ValueError('listOfSubsystemFusion argument should be list of Subsystem objects')
        
    if (type(specietrigger) is not str):
        raise ValueError('The specietrigger argument is not a string.') 
        
    if (type(operator) is not str):
        raise ValueError('The operator argument is not a string.') 

    if type(conctrigger) is not int and type(conctrigger) is not float:
        raise ValueError('The conctrigger argument is not a float.') 

    
    #combining subsystems to get 1 model of listOfSubsystemTrigger (complete)
    subTrig=createNewSubsystem() 
    if len(listOfSubsystemTrigger)==1:
        subTrig= listOfSubsystemTrigger[0]
    else:
        subTrig.combineSubsystems(listOfSubsystemTrigger, mode="virtual", combineNames = True, **kwargs)
        
    model_abc=subTrig.getSBMLDocument().getModel()
    id_abc=subTrig.getAllIds()
    if not specietrigger in id_abc:
        raise ValueError('The specietrigger argument is not as specie id.')
    #simulating the subTrig (complete)
    result1,_= subTrig.simulateWithBioscrape(timepoints)
    
    #getting time (complete)
    triggerTime=0
    indexTirg=0
    i=0
    if operator =='>':
        while i < len(result1.index):
            if result1.loc[i][specietrigger] > conctrigger:
                triggerTime= result1.loc[i]['time']
                indexTirg=i
                i=len(result1.index)+1
            i=i+1    

    elif operator =='<':
        while i < len(result1.index):
            if result1.loc[i][specietrigger] < conctrigger:
                triggerTime= result1.loc[i]['time']
                indexTirg=i
                i=len(result1.index)+1
            i=i+1    

    elif operator =='=':
        while i < len(result1.index):
            if result1.loc[i][specietrigger] == conctrigger:
                triggerTime= result1.loc[i]['time']
                indexTirg=i
                i=len(result1.index)+1
            i=i+1
    else:
        raise ValueError('The operator argument should be ">" or "<" or "=".') 

    initialTimepoints= timepoints[:indexTirg+1]
    finalTimepoints= timepoints[indexTirg:]
    
    result_xyz=[]
    #initial simulation for the listOfSubsystemFusion
    for s in listOfSubsystemFusion:
        result_abc,_= s.simulateWithBioscrape(initialTimepoints)
        result_xyz.append(result_abc)
    result2 = pd.concat(result_xyz, axis='columns', join='outer', copy = 'False')
    
    #get specie conc at the end of initialTimepoints in a dictionary
    result_endtime=result2.loc[indexTirg,:]
    specie_conc=result_endtime.to_dict()
    
    #set specie conc before fusing
    for key in specie_conc:
        for subsys in listOfSubsystemFusion:
            model_2=subsys.getSBMLDocument().getModel()
            id_1=subsys.getAllIds()
            for abc in id_1:
                if key ==abc:
                    specie_id= key
                    specie_value= specie_conc[key]
                    model_2.getElementBySId(specie_id).setInitialConcentration(specie_value)
    
    #fusion and simulating with finalTimepoints with initial condition as at end of initialTimepoints
    fusedsubsystem=createNewSubsystem()
    for k in range(len(listOfSubsystemFusion)):
        subsystem= listOfSubsystemFusion[k]
        model = subsystem.getSBMLDocument().getModel()
        check(model,'retreived model object')
        for j in range (len(model.getListOfCompartments())):
            oldname = model.getCompartment(j).getName()
            if 'external' in oldname:
                newname = 'fusedcell_external'
                subsystem.renameCompartments(oldname, newname)
            elif 'internal' in oldname:
                newname = 'fusedcell_internal'
                subsystem.renameCompartments(oldname, newname)
    newsub= createNewSubsystem()
    newsub.combineSubsystems(listOfSubsystemFusion, mode = "virtual", combineNames = True, **kwargs)

    
    
    result3,_= newsub.simulateWithBioscrape(finalTimepoints)
    
    #joining result 2, 3
    #result = result2.append(result3, ignore_index=True, sort=False)
    return newsub, result2,initialTimepoints , result3, finalTimepoints