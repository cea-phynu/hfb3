#!/usr/bin/env python3

# Example of an HFB calculation for a deformed 16O nucleus.

# ==============================================================================
# ==============================================================================
# ==============================================================================

import hfb3

# ==============================================================================
# ==============================================================================
# ==============================================================================


def calc16O():

    # construct a DataTree with default values
    defaultDataTree = hfb3.DataTree.getDefault()

    # load a DataTree instance
    dataTree = hfb3.DataTree("16O_deformed.hfb3")

    # complete with default values
    dataTree = defaultDataTree + dataTree

    # create an Action instance
    action = hfb3.Action(dataTree)

    # launch an HFB calculation, starting from a WS solution
    action.calcWSHFB()

# ==============================================================================
# ==============================================================================
# ==============================================================================


if __name__ == "__main__":
    calc16O()
