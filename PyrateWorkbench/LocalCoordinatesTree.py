import sys
import os
from PySide import QtGui, QtCore
import FreeCADGui, FreeCAD
from core.coordinates import LocalCoordinates


#-------------------------------------------------------------------------------
# my test data
class MyData(LocalCoordinates):
    def __init__(self, name="", **kwargs):
        super(MyData, self).__init__(name, **kwargs)
        self.icon = None
        self.index = None

    #---------------------------------------------------------------------------
    def position(self):
        position = 0
        if self.parent is not None:
            count = 0
            for child in self.getChildren():
                if child == self:
                    position = count
                    break
                count += 1
        return position

    #---------------------------------------------------------------------------
    # test initialization
#    @staticmethod
#    def init():
#        root = MyData("root")
#        for i in range(0, 2):
#            child1 = MyData("child %i" % (i), root)
#            root.child.append(child1)
#            for x in range(0, 2):
#                child2 = MyData("child %i %i" % (i, x), child1)
#                child1.child.append(child2)
#
#        return root


#-------------------------------------------------------------------------------
class TreeModel(QtCore.QAbstractItemModel):

    #---------------------------------------------------------------------------
    def __init__(self, tree):
        super(TreeModel, self).__init__()
        self.__tree = tree
        self.__current = tree

    #---------------------------------------------------------------------------
    def flags(self, index):
        flag = QtCore.Qt.ItemIsEnabled
        if index.isValid():
            flag |= QtCore.Qt.ItemIsSelectable \
                 | QtCore.Qt.ItemIsUserCheckable \
                 | QtCore.Qt.ItemIsEditable \
                 | QtCore.Qt.ItemIsDragEnabled \
                 | QtCore.Qt.ItemIsDropEnabled
        return flag

    #---------------------------------------------------------------------------
    def index(self, row, column, parent=QtCore.QModelIndex()):
        node = QtCore.QModelIndex()
        if parent.isValid():
            nodeS = parent.internalPointer()
            nodeX = nodeS.children[row]
            node = self.__createIndex(row, column, nodeX)
        else:
            node = self.__createIndex(row, column, self.__tree)
        return node

    #---------------------------------------------------------------------------
    def parent(self, index):
        node = QtCore.QModelIndex()
        if index.isValid():
            nodeS = index.internalPointer()
            parent = nodeS.parent
            if parent is not None:
                node = self.__createIndex(parent.position(), 0, parent)
        return node

    #---------------------------------------------------------------------------
    def rowCount(self, index=QtCore.QModelIndex()):
        count = 1
        node = index.internalPointer()
        if node is not None:
            count = len(node.children)
        return count

    #---------------------------------------------------------------------------
    def columnCount(self, index=QtCore.QModelIndex()):
        return 1

    #---------------------------------------------------------------------------
    def data(self, index, role=QtCore.Qt.DisplayRole):
        data = None
        if role == QtCore.Qt.DisplayRole or role == QtCore.Qt.EditRole:
            node = index.internalPointer()
            data = node.name

        if role == QtCore.Qt.ToolTipRole:
            node = index.internalPointer()
            data = "ToolTip " + node.name

        if role == QtCore.Qt.DecorationRole:
            data = QtGui.QIcon("icon.png")
        return data

    #---------------------------------------------------------------------------
    def setData(self, index, value, role=QtCore.Qt.DisplayRole):
        result = True
        if role == QtCore.Qt.EditRole and value != "":
            node = index.internalPointer()
            node.text = value
            result = True
        return result

    #---------------------------------------------------------------------------
    def __createIndex(self, row, column, node):
        if node.index == None:
            index = self.createIndex(row, column, node)
            node.index = index
            icon = QtGui.QIcon("icon.png")
            b = self.setData(index, icon, QtCore.Qt.DecorationRole)
            b = self.setData(index, "ToolTip "+node.name, QtCore.Qt.ToolTipRole)
        return node.index



#-------------------------------------------------------------------------------
#class TreeView(QtGui.QTreeView):
#    #---------------------------------------------------------------------------
#    def __init__(self, model, parent=None):
#        super(TreeView, self).__init__(parent)
#        self.__model = model
#        self.setModel(model)
#
#
#        self.setCurrentIndex(self.__model.index(0, 0))
#        return


class LocalCoordinatesTaskPanel:
    def __init__(self):
        # this will create a Qt widget from our ui file
        # TODO: relative paths?
        fn = os.path.join(os.path.dirname(__file__), 'Qt/lcdialog.ui')        
        self.form = FreeCADGui.PySideUic.loadUi(fn)

        data = MyData(name="test1")
        lc21 = data.addChild(MyData(name="lc21", decz=40.0))
        lc22 = data.addChild(MyData(name="lc22", decz=50.0))
        lc23 = data.addChild(MyData(name="lc23", decz=60.0))
        lc31 = lc22.addChild(MyData(name="lc31", decz=60.0))
        lc32 = lc22.addChild(MyData(name="lc32", decz=60.0))
        lc33 = lc23.addChild(MyData(name="lc33", decz=60.0))
        lc34 = lc23.addChild(MyData(name="lc34", decz=60.0))
        
        treeModel = TreeModel(data)
        #self.form.treeView = QtGui.QTreeView()
        self.form.treeView.setModel(treeModel)#TreeView(treeModel)
        self.form.treeView.setCurrentIndex(treeModel.index(0, 0))

    def accept(self):
        print("pressed ok")
        FreeCADGui.Control.closeDialog()
        
    def reject(self):
        print("pressed cancel")
        FreeCADGui.Control.closeDialog()
        


class CreateLocalCoordinatesTool:
    "Tool for creating local coordinates"

    def GetResources(self):
        return {"Pixmap"  : ":/icons/pyrate_logo_icon.svg", # resource qrc file needed, and precompile with python-rcc
                "MenuText": "Create local coordinates ...",
                "Accel": "",
                "ToolTip": "Opens dialog for local coordinates"
                }

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
            return False
        else:
            return True

    def Activated(self):

        doc = FreeCAD.ActiveDocument
        panel = LocalCoordinatesTaskPanel()
        FreeCADGui.Control.showDialog(panel)



FreeCADGui.addCommand('CreateLocalCoordinatesCommand', CreateLocalCoordinatesTool())

#-------------------------------------------------------------------------------
#class MyTree(QtGui.QMainWindow):
#    def __init__(self, parent=None):
#        super(MyTree, self).__init__(parent)
#
#        data = MyData.init()
#        treeModel = TreeModel(data)
#        treeView = TreeView(treeModel)
#
#        self.setCentralWidget(treeView)
#
#
##-------------------------------------------------------------------------------
#def main():
#    app = QtGui.QApplication.instance()
#    if app is None:
#        app = QtGui.QApplication(sys.argv)    
#    form = MyTree()
#    form.show()
#    app.exec_()
#
#if __name__ == '__main__':
#    main()
