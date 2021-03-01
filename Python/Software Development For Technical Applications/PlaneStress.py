# -*- coding: utf-8 -*-

import numpy as np
import json
import calfem.core as cfc
import calfem.geometry as cfg  # <-- Geometry
import calfem.mesh as cfm      # <-- Mesh generating
import calfem.vis as cfv       # <-- Visualisation
import calfem.utils as cfu     # <-- mixed routines 
import pyvtk as vtk          # Paraview module


class InputData(object):
 #   """Class to define inputdata to the model."""
    
    def __init__(self):
        
        # Version 
        self.version = 1
        
        
    def geometry(self):
        """Create a geometry instance based on the 
        previously defined parameters"""

        #Create the geometry instans 

        g = cfg.Geometry()

        w = self.w
        h = self.h
        a = self.a
        b = self.b



        #Create points for the geometry 

        g.point([0, 0])             # 0
        g.point([(w-a)/2, 0])       # 1
        g.point([(w-a)/2, b])       # 2
        g.point([w/2+a/2, b])       # 3
        g.point([w/2+a/2, 0])       # 4
        g.point([w, 0])             # 5
      
        g.point([w, h])             # 6
        g.point([w/2+a/2, h])       # 7
        g.point([w/2+a/2, h-b])     # 8
        g.point([(w-a)/2, h-b])     # 9
        g.point([(w-a)/2, h])       # 10
        g.point([0, h])             # 11
        

        # Link the points with splines 
        # Use a marker to define splines with a load or boundary condition 

        g.spline([0, 1])            # 0
        g.spline([1, 2])            # 1    
        g.spline([2, 3])            # 2 
        g.spline([3, 4])            # 3
        g.spline([4, 5])            # 4
        g.spline([5, 6], marker=10) # 5
        g.spline([6, 7])            # 6
        g.spline([7, 8])            # 7
        g.spline([8, 9])            # 8
        g.spline([9, 10])           # 9
        g.spline([10, 11])          # 10
        g.spline([11, 0],marker=20) # 11
        
   
        # Create the surface
        g.surface([0,1,2,3,4,5,6,7,8,9,10,11])

        # Return the geometry
        return g
        
    def save(self, filename):
        """Save inputdata to file."""
       
        inputData = {}
        inputData["version"] = self.version
        inputData["t"] = self.t
        inputData["E"] = self.E
        inputData["v"] = self.v
        inputData["w"] = self.w
        inputData["h"] = self.h
        inputData["a"] = self.a
        inputData["b"] = self.b
        inputData["q"] = self.q   
        inputData["elSizeFactor"] = self.elSizeFactor
        inputData["aStart"] = self.aStart
        inputData["bStart"] = self.bStart
        inputData["aEnd"] = self.aEnd
        inputData["bEnd"] = self.bEnd
        inputData["Number of steps"] = self.paramSteps

        ofile = open(filename, "w")
        json.dump(inputData, ofile, sort_keys = True, indent = 4)
        ofile.close()

    def load(self, filename):
        """Read inputdata from file."""

        ifile = open(filename, "r")
        inputData = json.load(ifile)
        ifile.close()

        self.version = inputData["version"]
        self.t = inputData["t"]
        self.E = inputData["E"]
        self.v = inputData["v"]
        self.w = inputData["w"]
        self.h = inputData["h"]
        self.a = inputData["a"]
        self.b = inputData["b"]
        self.q = inputData["q"]  
        self.elSizeFactor = inputData["elSizeFactor"]
        self.aStart = inputData["aStart"]
        self.bStart = inputData["bStart"]
        self.aEnd = inputData["aEnd"]
        self.bEnd = inputData["bEnd"]
        self.paramSteps = inputData["Number of steps"]
        
        
        
class Solver(object):
    """Class to manage the solution of the model"""
    def __init__(self, inputData, outputData):
        self.inputData = inputData
        self.outputData = outputData
      
    def execute(self):

        # Transfer to model parameters 
        E = self.inputData.E
        t = self.inputData.t
        v = self.inputData.v
        q = self.inputData.q
        h = self.inputData.h
        ptype = self.inputData.ptype 
        ep=[ptype,t]
        elSizeFactor = self.inputData.elSizeFactor/1000

        
    
       # Call on inputData for geometry
        geometry = self.inputData.geometry()
        
        # Mesh generating
        elType = 3      # Four node element
        dofsPerNode= 2  # stress-strain --> 2 degrees of freedom 

        meshGen = cfm.GmshMeshGenerator(geometry)
        meshGen.elSizeFactor = elSizeFactor
        meshGen.elType = elType
        meshGen.dofsPerNode = dofsPerNode
        
        # Calculating mesh properties
        coords, edof, dofs, bdofs, elementmarkers = meshGen.create()
        nDofs = np.size(dofs) # Number of dofs
        nelm = edof.shape[0] # Number of elements
        
        # Extract ex and ey from coord-matrix
        ex,ey = cfc.coordxtr(edof,coords,dofs)

        # Initiating global matrix
        K=np.matrix(np.zeros((nDofs,nDofs)))
        f=np.matrix(np.zeros((nDofs,1))) 
        stress=np.matrix(np.zeros([nelm,3]))
        
        bc = np.array([],"i")
        bcVal = np.array([],"i")
        
        # Constitutive matrix
        D=cfc.hooke(ptype,E,v)
        
        # Calculating Stiffness matrix
        for elx, ely, eltopo in zip(ex,ey,edof):
            # Element stiffness matrix
            Ke=cfc.planqe(elx,ely,ep,D)
            # Assemble stiffness matrix into global
            cfc.assem(eltopo,K,Ke)
            
        # Force Vector
        cfu.applyforcetotal(bdofs,f,10,q*h,dimension=1)

        # Apply boundary conditions 
        bc, bcVal = cfu.applybc(bdofs,bc,bcVal,20,value=0.0)
       
        # Solve the equation system 
        a,r = cfc.solveq(K,f,bc,bcVal)
        
        # Extracting elemental displacement 
        ed=cfc.extractEldisp(edof,a)
    
        # Calculate the stress and Von Mises Stress
        i=0
        vonMises=[]
        stress_1 = []
        stress_2 = []

        
        for elx, ely, eld in zip(ex,ey,ed):
            [stress[i],_]=cfc.planqs(elx,ely,ep,D,eld)
            vonMises.append(np.sqrt(pow(stress[i,0],2)-stress[i,0]*stress[i,1]+pow(stress[i,1],2)+3*stress[i,2]))

            
            w,v = np.linalg.eig(np.array([[stress[i,0],stress[i,2],0],[stress[i,2],stress[i,1],0],[0,0,0]]))

            stress_1.append(v[0].tolist())
            stress_2.append(v[1].tolist())

                                     
            i=i+1

        
        
        # Extracting the maximum von mises stress
        MaxStress = np.max(vonMises)
        
        # Extracting the minumum von mises stress
        MinStress = np.min(vonMises)

        # Extracting the max displacement node
        MaxDisp = np.max(np.abs(a))
        
        # Extracting the min displacement node
        MinDisp = np.min(np.abs(a))
        
        
        # Outputdata
        self.outputData.coords = coords
        self.outputData.edof = edof
        self.outputData.dofs = dofs 
        self.outputData.a = a
        self.outputData.r = r
        self.outputData.ed = ed
        self.outputData.stress = stress
        self.outputData.geometry = geometry 
        self.outputData.elType = elType      
        self.outputData.dofsPerNode = dofsPerNode
        self.outputData.vonMises = vonMises
        self.outputData.MaxStress = MaxStress
        self.outputData.MinStress = MinStress
       
        self.outputData.topo = meshGen.topo
        self.outputData.stress1 = stress_1
        self.outputData.stress2 = stress_2
        self.outputData.MaxDisp = MaxDisp
        self.outputData.MinDisp = MinDisp
        
    def executeParamStudy(self):
        """Execute parameter study"""   
        # -- Store previous values of a and b
        old_a = self.inputData.a
        old_b = self.inputData.b
    
        i = 1
     
        if self.inputData.paramA:
    
            # --- Create values to perform simulation
    
            aRange = np.linspace(self.inputData.aStart, self.inputData.aEnd,
                self.inputData.paramSteps)
    
            # --- Begin parameterstudy
            i = 1
            for a in aRange:
                print("Executing for a = %g..." % a)
    
                # --- set the desired parameter in the InputData-instance
                self.inputData.a=a
                
                # --- Run simulation
                self.execute()
                
                
                # --- Export vtk-file
                self.exportVtk("paramStudy_0"+str(i)+".vtk")
                i+= 1
    
        elif self.inputData.paramB:
                       # --- Create values to perform simulation
    
            bRange = np.linspace(self.inputData.bStart, self.inputData.bEnd,
                self.inputData.paramSteps)
    
            # --- Begin parameterstudy
            i = 1
            for b in bRange:
                print("Executing for b = %g..." % b)
    
                # --- set the desired parameter in the InputData-instance
                self.inputData.b=b
                
                # --- Run simulation
                self.execute()
                
                # --- Export vtk-file
                self.exportVtk("paramStudy_0"+str(i)+".vtk")
                i+= 1
    

        # --- Reset original values
    
        self.inputData.a = old_a
        self.inputData.b = old_b
        
        
    def exportVtk(self, filename):
        """Export results to VTK"""
    
        print("Exporting results to %s." % filename)
        # ---  Create points and polygon defined from the mesh
    
        points = self.outputData.coords.tolist()

    
        polygons = (self.outputData.topo-1).tolist()

        
        # --- Results from the simulation is created in separable objects. Points in vtk.PointData and
        # ---  elementdata in vtk.CellData. 
        cellData = vtk.CellData(vtk.Scalars(self.outputData.vonMises, name="mises") , vtk.Vectors(self.outputData.stress1, "principal stress 1"), vtk.Vectors(self.outputData.stress2, "principal stress 2"))
    
        # --- Create a structure for the elementmesh.
    
        structure = vtk.PolyData(points = points, polygons = polygons)
    
        # --- Store everything in a vtk.VtkData instance
    
        vtkData = vtk.VtkData(structure, cellData)
    
        # --- Save to a file
        
        vtkData.tofile(filename, "ascii")
        
class OutputData(object):
    """Class to store the results from the simulation"""
    
    def __init__(self):
        self.a = None # none creates variables containing nothing
        self.r = None
        self.ed = None
        self.stress = None
        self.vonMises = None
        self.MaxStress = None
        self.MinStress = None
        self.MaxDisp = None
        self.MinDisp = None
    
        
        self.geometry = None
        self.coords = None
        self.edof = None
        self.dofsPerNode = None
        self.elType = None
        self.dofs = None
       


class Report(object):
    """Class for presentation of inputdata and outpdata in a report "structure" """
    def __init__(self, inputData, outputData):
        self.inputData = inputData
        self.outputData = outputData
        self.report = ""

    def clear(self):
        self.report = ""

    def addText(self, text=""):
        self.report+=str(text)+"\n"

    def __str__(self):
        self.clear()
        self.addText()
        self.addText("-------------- Model input ----------------------------------")
        self.addText()
        
        
        self.addText("Parameter Study, a start:")
        self.addText()
        self.addText(self.inputData.aStart)
        self.addText()
        
        self.addText("Parameter Study, a end:")
        self.addText()
        self.addText(self.inputData.aEnd)
        self.addText()
        
        self.addText("Parameter Study, b start:")
        self.addText()
        self.addText(self.inputData.bStart)
        self.addText()
        
        self.addText("Parameter Study, b end:")
        self.addText()
        self.addText(self.inputData.bEnd)
        self.addText()
        
        self.addText("Parameter Study, parameter step:")
        self.addText()
        self.addText(self.inputData.paramSteps)
        self.addText()
        
        self.addText("Coordinates:")
        self.addText()
        self.addText(self.outputData.coords)
        self.addText()
        
        
        self.addText("Element degree of freedom,edof:")
        self.addText()
        self.addText(self.outputData.edof)
        self.addText()

        self.addText("Degree of freedom, dof:")
        self.addText()
        self.addText(self.outputData.dofs)
        self.addText()
        
    
        self.addText()
        self.addText("-------------- Model output ----------------------------------")
        self.addText()
        
        self.addText("The size of the magnification factor used for the visual results of nodal displacements:")
        self.addText()
        self.addText(self.inputData.magnfac)
        self.addText()
        
        
        self.addText("Nodal displacements, elementwise:")
        self.addText()
        self.addText(self.outputData.ed)
        self.addText()

 
        self.addText("Max Von Mises stress (Pa):")
        self.addText()
        self.addText(self.outputData.MaxStress)
        self.addText()
        
        self.addText("Min Von Mises stress (Pa):")
        self.addText()
        self.addText(self.outputData.MinStress)
        self.addText()
        
        
        self.addText("Max Nodal displacement (m):")
        self.addText()
        self.addText(self.outputData.MaxDisp)
        self.addText()
        
        self.addText("Min Nodal displacement (m):")
        self.addText()
        self.addText(self.outputData.MinDisp)
        self.addText()


        return self.report
    
    def writeToFile(self,filename):
        # Write the report to a text-file
        ofile=open(filename,"w")
        
        for line in self.report:
            ofile.write(line)
        ofile.close()

class Visualisation(object):
    def __init__(self, inputData, outputData):
        self.inputData = inputData
        self.outputData = outputData
        
        # --- Variables which stores references to open figures

        self.geomFig = None
        self.meshFig = None
        self.elValueFig = None
        self.nodalDisplacementFig = None
        
    def showGeometry(self):
        """Show geometry visualization"""

        geometry = self.outputData.geometry

        self.geomFig = cfv.figure(self.geomFig)
        cfv.clf()            
        cfv.drawGeometry(geometry, title="Geometry")
    
    def showMesh(self):
        """Show mesh visualization"""
        
        coords = self.outputData.coords
        edof = self.outputData.edof
        dofsPerNode = self.outputData.dofsPerNode
        elType = self.outputData.elType
    
        self.meshFig = cfv.figure(self.meshFig)
        cfv.clf()            
        cfv.drawMesh(coords=coords,edof=edof,dofsPerNode=dofsPerNode,elType=elType,filled=True,title="Mesh")
        
    def showNodalDisplacement(self):
        """Show nodal displacement visualization"""

        a = self.outputData.a
        coords = self.outputData.coords
        edof = self.outputData.edof
        dofsPerNode = self.outputData.dofsPerNode
        elType = self.outputData.elType
        magnfac = self.inputData.magnfac

        self.nodalDisplacementFig = cfv.figure(self.nodalDisplacementFig)
        cfv.clf()            
        cfv.drawDisplacements(a,coords,edof,dofsPerNode,elType, doDrawUndisplacedMesh=False, magnfac=magnfac, title="Displacement")
        
    def showElementValues(self):
        """Show Von mises visualization"""

        a = self.outputData.a
        coords = self.outputData.coords
        edof = self.outputData.edof
        dofsPerNode = self.outputData.dofsPerNode
        elType = self.outputData.elType
        vonMises = self.outputData.vonMises
        magnfac = self.inputData.magnfac

        self.elValueFig = cfv.figure(self.elValueFig)
        cfv.clf()            
        cfv.drawElementValues(vonMises,coords,edof,dofsPerNode,elType,a,doDrawMesh=True,doDrawUndisplacedMesh=False,magnfac=magnfac, title="Von Mises Effective stress")
        cfv.colorBar().SetLabel("Effective Stress")
        
        
    def closeAll(self):
        cfv.closeAll()
        

    def wait(self):
        """This method ensures that the windows are kept updated and will return 
        when the last window is closed."""

        cfv.showAndWait()