Ideas
===

FreeCAD Document Structure
---

* Groups for
    * Surfaces,
    * Apertures,
    * Shapes,
    * Local Coordinates,
    * Functions,
    * Materials (perhaps sub groups for Shelf, Book, ... known from refractive-index.info)
    * No further subdivision

Dialogs and UI/UX
---

* Interface classes and dialogs for every ClassWithOptimizableVariables
    * 1:1 correspondence in a generalized manner
    * Changes structures/annotations of ClassWithOptimizableVariables
    * Links to other ClassWithOptimizableVariables if necessary

* Dialogs for system creation
    * Dialogs for the convenience functions (tabular data) => does not need the Mirror and STOP flags
    * in a first step no lens editor (done by tree view in FC)
