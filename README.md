## Processes

Structures and functions to support a Kinetic Monte Carlo simulation.

Demo version V0.1.   Most of the data structures and functions to support a full kinetic monte carlo
simulation have been implemented for the case of a simple cubic 001 crystal. 
FCC 111 has been partially implemented in the Crystal module, but only the crystal structure
and the capability to deposit an atom that makes a soft landing on the surface. 
None of the other processes for FCC111 have been developed yet.

First commit on 1/6/2021  to https://github.com/rheadric/Processes

## Dependencies

PyPlot -- for making a plot of Cubic or FCC111 crystals.
StatsBase -- for choosing a process randomly accoring to weights.
Crystal --  a primitive module to store the crystal structure in an array. 
    Crystal also includes the AtomRegistry structure and plotting functions
    to plot the crystal.

Note: Crystal is an unregistered package.  It can be found at https://github.com/rheadric/Crystal

## Help

For help, use the help feature (?) in the Julia REPL:

julia> ?Processes

There is also help for the exported functions listed below:

getprocess, CubicProcessList, getindex, countneighbors!, IrreversibleStickingRates!, CubicProcessRates,
findneighbors, processController!, updateAtomRegistry!

## Here are my raw development notes from the file Processes.jl

#
#   12/25/2020 RLH   I'll add processes such as diffusion to the crystal types in Crystal.
#
# V0.0 RLH 12/26/20 Atom now holds an array of type 'Process'. See the Crystal module.
# V0.1 RLH 12/27/20 CubicProcessList() initializes the main process list, which is an array of the  
#                   type 'Process'. It takes a Cubic crystal as input. 
# V0.1 RLH 12/29/20 Implemented countneighbors!(). Counts nearest neighbors for a given process  
#                   before and after the move. Modifies the Process, which is linked to from the 
#                   corresponding Atom and also from the CubicProcessList.
#                   This version is for Cubic.  The FCC111 version will be impemented seperately.
# V0.1 RLH 12/30/20 The current version of CubicProcessList() doesn't know anything about rates.
#                   It just builds the process list according to which neighboring sites are available,
#                   and sets all rates to 1.0.
#                   New function: IrreversibleStickingRates!() set rates for each possible process, based 
#                   on the starting NN count. If the number of NN's is 0, then the rate will be 1.0.  For
#                   a number of NN's > 0, the rate will be 0 and the Process will be removed from the list.
# V0.1 RLH 12/31/20 Clean up a few things before working on the process controller. (i) CubicProcessList()
#                   attach all 8 processes to each atom whether they are blocked or not. If they are blocked 
#                   or disallowed, the rates will be set to 0.0. (ii) IrreversibleStickingRates!() Processes
#                   with 0 NN's are not being initialized.  Revise to set the rates for these processes to 1.0.
#                   Note: everything is working awesomely before midnight! :)
# V0.1 RLH 12/31/20 CubicProcessList is a struct that hold the process rates for all of the hop directions 
#                   and also the deposition rates.  This brings us a step closer to a deposition controller.
# V0.1 RLH 1/1/21   Implemented getprocess(). It selects a process weighted by the process rates.
#                   findneighbors() returns list of neighbors before and after a given process.
#                   getprocess() selects a process at random from the process lists, weighted by the process rates.
# V0.1 RLH 1/1/21   processController!() takes the process and lists of old neighbors and new neighbors.
#                   updateAtomRegistry!() The neighbors whose environment has changed need to have their 
#                   processes updated.  Do this in Crystal.Cubic first.
# V0.1 RLH 1/3/21   Apply the sticking model.   New version of IrreversibleStickingRates!() to work on the crystal.
#                   updateCubicProcessListandRates!(): Old processes that are no longer allowed will be dropped 
#                   from the process lists.  Any new processes that were not previously on the relevant list
#                   will need to be added. Procees rates in  struct CubicProcessRates are also updated.
# V0.1 RLH 1/5/21   I consider this to be a first working version, good enough for a demo at least. 
#                   First commit on github!
#
# V0.1 RLH 1/5/21   I expect a long dormant period.  A rough roadmap for when I pick it up again -- 
#                   V0.2 : Develop structures and functions for FCC111 processes and rates. 
#                   Re-use code that I've already developed for the Cubic case  whenever possible.
#                   V0.3:  will probably be about implemention rates for specific processes, 
#                   i.e. surface diffusion, edge and corner diffusion, Schwoebel barrier, etc. I've already
#                   implemented IrreversibleStickingRates!() that just looks at the NN count and sets the rates
#                   to zero if NNs > 0.
#                   V0.4 and beyond:  There needs to be a main() function that controls everything.  We'll need 
#                   lots of data archiving and analysis tools.  Maybe we'll use tools that Jeff Ulbrant
#                   has developed for his KMC code.
#                   Or maybe none of that will ever happen.  It's been worth it as a project just to learn Julia!

