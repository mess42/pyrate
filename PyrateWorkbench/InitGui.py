import resources_rc

class PyrateWorkbench ( Workbench ):
    "Pyrate workbench object"
    Icon = ":/icons/pyrate_logo_icon.svg"
#     Icon = """
#             /* XPM */
#             static const char *test_icon[]={
#             "16 16 2 1",
#             "a c #000000",
#             ". c None",
#             "................",
#             ".###............",
#             ".#..#...........",
#             ".#..#...........",
#             ".###..#..#.###..",
#             ".#....#..#.#....",
#             ".#.....###.#....",
#             ".........#......",
#             "......###.......",
#             "................",
#             "........#.......",
#             ".......###..##..",
#             "..###...#..#..#.",
#             ".#..#...#..####.",
#             "..####...#.#....",
#             "............##.."};
#             """
    MenuText = "Pyrate Workbench"
    ToolTip = "Pyrate optical design Workbench"
    def GetClassName(self):
        return "Gui::PythonWorkbench"

    def Initialize(self):
        import CreateSystem
        self.appendToolbar("Pyrate", ["CreateSystemCommand"])
        self.appendMenu("Pyrate", ["CreateSystemCommand"])
        Log ("Loading Create System Module... done\n")

    def Activated(self):
# do something here if needed...
        Msg ("PyrateWorkbench.Activated()\n")

    def Deactivated(self):
# do something here if needed...
        Msg ("PyrateWorkbench.Deactivated()\n")

FreeCADGui.addWorkbench(PyrateWorkbench())
