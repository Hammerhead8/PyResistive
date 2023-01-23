# PyResistive
PyResistive is a rewrite in Python of the C++ tool Resistive (https://github.com/Hammerhead8/Resistive) for calculating the node voltages of in passive DC and AC circuits. The syntax is similar to SPICE, where the user provides the nodes where the element is connected and its value. This means that the user does not need to calculate the conductance matrix themselves. Instead, the conductance matrix is updated automatically as elements are added to the circuit. PyResistive uses nodtal analysis to calculate the node voltages.

Currently, PyResistive supports DC and AC sources, but only if all sources are of the same frequency and all sources need to have one node connected to ground.

# Dependencies
PyResistive depends on Numpy and the linear algebra module from Scipy.

