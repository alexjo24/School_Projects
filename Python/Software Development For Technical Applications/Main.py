# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 10:41:44 2017

@author: kajsa
"""

import sys

from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog,QMessageBox
from PyQt5.uic import loadUi

import PlaneStress as ps

from PyQt5 import QtCore

class SolverThread(QtCore.QThread):
    """To manage the solver in the background"""
    
    def __init__(self, solver, paramStudy=False):
        """Class Constructor"""
        QtCore.QThread.__init__(self)
        self.solver = solver
        self.paramStudy = paramStudy

    def __del__(self):
        self.wait()

    def run(self):
        
        if self.paramStudy == True:
            self.solver.executeParamStudy()
        else:
            self.solver.execute()
    
        
class MainWindow(QMainWindow):
    """MainWindow-klass which manages the main window"""

    def __init__(self):
        self.filename = None
        self.filenameReport = None
       
        super(QMainWindow, self).__init__()

        # Store a reference to the applicationinstance in the class
        self.app = app
        
         # Read an interface from file
        self.ui = loadUi('mainWindow.ui', self)
        
                
        # Connect controls to an eventmethod
        # To connect events in the menu: Triggered
        # To connect buttons: Clicked
        
        """Menu"""
        self.ui.actionNew.triggered.connect(self.onActionNew)
        
        self.ui.actionOpen.triggered.connect(self.onActionOpen)
        
        self.ui.actionSave.triggered.connect(self.onActionSave)
        
        self.ui.actionSaveAs.triggered.connect(self.onActionSaveAs)
        
        self.ui.actionExit.triggered.connect(self.onActionExit)
        
        self.ui.actionExecute.triggered.connect(self.onActionExecute)
        
        
        """Buttons"""
        self.ui.showGeometryButton.clicked.connect(self.onShowGeometry)
        
        self.ui.showMeshButton.clicked.connect(self.onShowMesh)
        
        self.ui.showNodalValuesButton.clicked.connect(self.onShowNodalValues)
        
        self.ui.showElementValuesButton.clicked.connect(self.onShowElementValues)
        
        self.ui.WriteOutPutButton.clicked.connect(self.onWriteReport)
        
        self.ui.showParamButton.clicked.connect(self.onExecuteParamStudy)
     
    
        """Slider"""
        self.ui.ElemSlider.setRange(0,100)
        self.ui.ElemSlider.setSingleStep(1)
        self.ui.ElemSlider.valueChanged.connect(self.onValueChanged)
        
        # Ensure to display the window
        self.ui.show()
        self.ui.raise_()



    def onSolverFinished(self):
        """Called when the calculation thread ends"""

        # --- Activate inteface

        self.ui.setEnabled(True)
       
        # --- Generate result report
        self.report = ps.Report(self.inputData,self.outputData)
      
        self.ui.reportEdit.setPlainText(str(self.report))
        
  
    def onActionNew(self):
        """Create a new model"""
      
        
        self.new = QMessageBox.question(self.ui, "Message", "New model?", 
                                         QMessageBox.Yes | QMessageBox.No,QMessageBox.No)
        if self.new == QMessageBox.Yes:
            
            
            self.newSave = QMessageBox.question(self.ui, "Message", "Do you want to save the file?", 
                                QMessageBox.Yes | QMessageBox.No,QMessageBox.No)
            
            if self.newSave == QMessageBox.Yes:
                self.onActionSave()
                
                
            
            vis = ps.Visualisation(self.inputData, self.outputData)
            vis.closeAll()
            
            self.initModel()
            self.updateModel()
            self.filename = None


    def onActionSave(self):
        """Save model"""
       
        self.updateModel()
        
        if self.filename is None :
            self.filename, _  = QFileDialog.getSaveFileName(self.ui,          
                'Save model', '', 'Modell filer (*.json)')   
            
            #If the user selects cancel in the save window, self.filename 
            #will contain self.filename = '', i.e nothing.
            if self.filename == '':
                self.filename = None
                
           
        
        if  self.filename is not None:
            QMessageBox.information(self.ui,"Message", "The model has been saved")
            self.inputData.save(self.filename)  
          
            
            
         
    def onActionSaveAs(self):
        """Save a new model"""
        
        self.updateModel()
        
       
        self.filename, _  = QFileDialog.getSaveFileName(self.ui,          
                'Spara modell', '', 'Modell filer (*.json)')
     
        #If the user selects cancel in the save window, self.filename 
        #will contain self.filename = '', i.e nothing.
        if self.filename == '':
            self.filename = None
        
        if self.filename is not None:
            self.inputData.save(self.filename)
        
    def onActionExit(self):
        """Close"""
       
        self.exit = QMessageBox.question(self.ui, "Meddelande", "Exit?", 
                                         QMessageBox.Yes | QMessageBox.No,QMessageBox.No)
        if self.exit == QMessageBox.Yes:
            sys.exit()


    def onActionOpen(self):
        """Open a existing file"""
        
        
        self.filename, _ = QFileDialog.getOpenFileName(self.ui, 
        "Open model", "", "Modell filer (*.json *.jpg *.bmp)")
        
        if self.filename == '':
            self.filename = None
        
        if self.filename is not None:
            QMessageBox.information(self.ui, "Message", self.filename)
            self.inputData.load(self.filename)
            self.updateControls()
        
            
                
    def onActionExecute(self):
        """Run simulation""" 
        
        #Close all windows 
        ps.Visualisation.closeAll(self)
        
        
        self.exit = QMessageBox.question(self.ui, "Message", "Run simulation?",      
                                         QMessageBox.Yes | QMessageBox.No,QMessageBox.No)

        # --- Deactivate interface/ui during simulation.
        
        if self.exit == QMessageBox.Yes:
            self.ui.setEnabled(False)

        # --- Update values from the controls
            self.updateModel()

        # --- Create a solver
            self.solver = ps.Solver(self.inputData, self.outputData)


        # --- Create a thread to run the simulation in, so the interface
        #     wont freeze.

            self.solverThread = SolverThread(self.solver)        
            self.solverThread.finished.connect(self.onSolverFinished)  
            self.solverThread.start()
            
     
    def calcDone(self):
        
        if self.outputData.coords == None:
            
            QMessageBox.information(self.ui, "Message","Not possible to obtain results before running a simulation.")
            
 
            
    def onShowGeometry(self):
        """Show geometry window"""
        self.calcDone()
        
        vis = ps.Visualisation(self.inputData, self.outputData)
        if self.outputData.coords != None:
            vis.showGeometry()
        
    def onShowMesh(self):
        """Show mesh window"""
        self.calcDone()
        
        vis = ps.Visualisation(self.inputData, self.outputData)
        
        if self.outputData.coords != None:
            vis.showMesh()
        
    def onShowNodalValues(self):
        """Show nodal values window"""
        self.calcDone()
        vis = ps.Visualisation(self.inputData, self.outputData)
        
        if self.outputData.coords != None:
            vis.showNodalDisplacement()
        
        
    def onShowElementValues(self):
        """Show element values window"""
        self.calcDone()
        vis = ps.Visualisation(self.inputData, self.outputData)
        
        if self.outputData.coords != None:
            vis.showElementValues()
            
    def onValueChanged(self):
        self.ui.labelSlider.setText(str(self.ui.ElemSlider.value()/1000))
        
    def onWriteReport(self):
        """ Write report to a text file"""
        
        self.calcDone()
        if self.outputData.coords != None:

            self.writeReport = QMessageBox.question(self.ui, "Message", "Do you want to save the report to file?", 
                                    QMessageBox.Yes | QMessageBox.No,QMessageBox.No)
                
            if self.writeReport == QMessageBox.Yes:
                
                self.filenameReport, _  = QFileDialog.getSaveFileName(self.ui,          
                    'Save model', '', 'Model filer (*.txt)')
         
                if self.filenameReport!= "":
                    self.report.writeToFile(self.filenameReport)   
                
        
    def onExecuteParamStudy(self):
        
        
        """Execute parameter study"""
    
        # --- collect graphical interface
    
        self.inputData.paramA = self.ui.aRadioButton.isChecked()
        self.inputData.paramB = self.ui.bRadioButton.isChecked()
        
       
        if self.inputData.paramA:
            self.inputData.aStart = float(self.ui.aEdit.text())
            self.inputData.aEnd = float(self.ui.aEndEdit.text())
            self.inputData.paramFilename = "paramStudy"
            self.inputData.paramSteps = int(self.ui.paramSpinBox.value())
                    
            # --- Update values from control
            self.updateModel()
                 
            # --- start a solver thread, 
                
            self.solverThread = SolverThread(self.solver, paramStudy = True)        
            self.solverThread.finished.connect(self.onSolverFinished)        
            self.solverThread.start()
            
        elif self.inputData.paramB:
            self.inputData.bStart = float(self.ui.bEdit.text())
            self.inputData.bEnd = float(self.ui.bEndEdit.text())
            self.inputData.paramFilename = "paramStudy"
            self.inputData.paramSteps = int(self.ui.paramSpinBox.value())
                    
            # --- Update values from control
            self.updateModel()
                 
            # --- start a solver thread, 
                
            self.solverThread = SolverThread(self.solver, paramStudy = True)        
            self.solverThread.finished.connect(self.onSolverFinished)        
            self.solverThread.start()
        else:
            QMessageBox.information(self.ui,"Message" , "Please ensure that parameter a or b is checked to perform parameter study.")
        
   
    def updateControls(self):
        """Update the controls from the window"""
    
        self.ui.wEdit.setText(str(self.inputData.w))
        
        self.ui.hEdit.setText(str(self.inputData.h))
        
        self.ui.aEdit.setText(str(self.inputData.a))
        
        self.ui.bEdit.setText(str(self.inputData.b))
        
        self.ui.tEdit.setText(str(self.inputData.t))

        self.ui.eEdit.setText(str(self.inputData.E))

        self.ui.vEdit.setText(str(self.inputData.v))
        
        self.ui.qEdit.setText(str(self.inputData.q))
        
        self.ui.ElemSlider.setValue(int(self.inputData.elSizeFactor))
        
        self.ui.aEndEdit.setText(str(self.inputData.aEnd))
        
        self.ui.bEndEdit.setText(str(self.inputData.bEnd))
        
            
    def updateModel(self):
        """Collect values from the window"""

        self.inputData.w = float(self.ui.wEdit.text())
        
        self.inputData.h = float(self.ui.hEdit.text())
        
        self.inputData.a = float(self.ui.aEdit.text())
        
        self.inputData.b = float(self.ui.bEdit.text())
        
        self.inputData.t = float(self.ui.tEdit.text())
        
        self.inputData.E = float(self.ui.eEdit.text())
        
        self.inputData.v = float(self.ui.vEdit.text()) 
        
        self.inputData.q = float(self.ui.qEdit.text())
        
        self.inputData.elSizeFactor = float(self.ui.ElemSlider.value())
        
        
    def initModel(self):
        """ Create neccesary objects to indata, outdata and solution."""
        
        self.ui.reportEdit.clear()
        self.inputData = ps.InputData()
        self.outputData = ps.OutputData()
        
        self.solver = ps.Solver(self.inputData, self.outputData)
        
        
        self.inputData.ptype = 1
        self.inputData.magnfac = 1000
        
        self.inputData.w = 0.3
        self.inputData.h = 0.1
        self.inputData.a = 0.05
        self.inputData.b = 0.025
        self.inputData.t = 0.15
        self.inputData.E = 2.08e10
        self.inputData.v = 0.2
        self.inputData.q = 100e3
        self.inputData.elSizeFactor = 50
        self.inputData.aStart = self.inputData.a
        self.inputData.bStart = self.inputData.b
        self.inputData.aEnd = self.inputData.a
        self.inputData.bEnd = self.inputData.b
        self.inputData.paramSteps = 1
        self.inputData.paramFileName = None
        self.inputData.paramA = False
        self.inputData.paramB = False

        
        self.updateControls()

if __name__ == '__main__':

    # --- create application

    app = QApplication(sys.argv)
    
    # --- Show main window

    widget = MainWindow()
    widget.show()
    widget.initModel()
    widget.updateModel()
    # --- Start loop
    sys.exit(app.exec_())
