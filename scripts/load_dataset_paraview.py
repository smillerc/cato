# -*- coding: utf-8 -*-
from paraview.simple import *
import os

paraview.simple._DisableFirstRenderCameraReset()


def load_cato_results(folder):
    animationScene1 = GetAnimationScene()
    timeKeeper1 = GetTimeKeeper()
    name = folder.split("/")[-2]
    results_folder = os.path.join(folder, "results")
    xdmf_files = list(
        sorted(
            [
                os.path.join(results_folder, f)
                for f in os.listdir(results_folder)
                if f.endswith(".xdmf")
            ]
        )
    )
    solution = XDMFReader(FileNames=xdmf_files)
    RenameSource(name, solution)
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # Filter out the ghost layers
    noGhost = Threshold(Input=solution)
    noGhost.Scalars = ["CELLS", "Ghost Cell"]
    noGhost.ThresholdRange = [0.0, 0.1]
    RenameSource("RemoveGhostCells", noGhost)
    Show(noGhost)
    Hide(solution)
