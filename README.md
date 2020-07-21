# SSSI-Python-Script
Electrical Power Systems are complex dynamical systems that are in constant tug of war between generation and load. In a simplistic way, both ends are constantly being tracked in order to balance the generation side to the load side and reach a stable operating margin. Since load is a stochastic parameter, there is a certain randomness that turns power system stability into a challenging and important task in power system analysis.
There are numerous stability indices that have been proposed throughout the years utilizing a steady state formulation. These formulations require knowledge of the system’s admittance matrix (YM), generator’s active power injection/consumption and the load’s active and reactive power. Alongside these requirements, there are also modeling limitations and simplifications that do not account for the dynamics of the system. 
For real-time power system stability assessment, these static methods are not helpful, with the exception in planning phases. That being said, the proposition in [1] is a step forward into developing an online stability index generator that can be used to help system operator’s in decision making. The development of the small signal stability assessment tool comes from an estimation of the eigenvalues of the system using time-series data from dynamic simulations. 
The stability index is calculated based on three layers of information:

1st layer - A scalar (positive or negative) that grades the overall damping ratio of all the modes.
2nd layer - A vector that points which pre-defined damping ratios were violated.
3rd layer - A matrix that has more information about each mode’s damping condition.
