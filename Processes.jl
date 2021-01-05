"""
To do: add some description and examples.
 """
module Processes
export getprocess, CubicProcessList, getindex, countneighbors!, IrreversibleStickingRates!, CubicProcessRates
export findneighbors, processController!, updateAtomRegistry!

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
#                   Plan: Old processes that are no longer allowed will be dropped from the process lists.  Any
#                   new processes that were not previously on the relevant list will need to be added. 
#                   

using PyPlot
using Crystal
using StatsBase

"""

### The base type for process lists.

Summary
* abstract type AbstractProcessList <: Any
    Subtypes 
* Cubic
* FCC111

"""
abstract type AbstractProcessList <: Any 
end

"""
Build the process list from scratch. 
The hop directions (North, South, East, West, etc.) are specific to a cubic crystal type.
"""
mutable struct CubicProcessList <: AbstractProcessList
    proclNorth::Array{Process,1}
    proclSouth::Array{Process,1}
    proclEast::Array{Process,1}
    proclWest::Array{Process,1}
    proclDownNorth::Array{Process,1}
    proclDownSouth::Array{Process,1}
    proclDownEast::Array{Process,1}
    proclDownWest::Array{Process,1}
    # Constructor for CubicProcessList
    function CubicProcessList(cubcryst::Cubic)
        atomreg = cubcryst.AtomRegistry 
        numatoms = cubcryst.nextatomID - 1
        # We have 8 process lists for hops in different directions.
        # More process lists can be created later if needed.
        dummyproc = Process(0, [0,0,0], 0, 0, 0.0)
        proclNorth = fill(dummyproc,1); popfirst!(proclNorth)
        proclSouth = fill(dummyproc,1); popfirst!(proclSouth)
        proclEast =  fill(dummyproc,1); popfirst!(proclEast)
        proclWest =  fill(dummyproc,1); popfirst!(proclWest)
        proclDownNorth = fill(dummyproc,1); popfirst!(proclDownNorth)
        proclDownSouth = fill(dummyproc,1); popfirst!(proclDownSouth)
        proclDownEast =  fill(dummyproc,1); popfirst!(proclDownEast)
        proclDownWest =  fill(dummyproc,1); popfirst!(proclDownWest)
        # If numatoms == 1 there is nothing left to do.  Return an empty process list.
        if numatoms != 1
            # Loop over all of the atoms starting with AtomID #2. 
            #(#1 is reserved for 1st layer, which is immobile.)
            # Note that the 'neighbors' array is 3x3x3 with 'thisatom' at [2,2,2].
            # A hop to the North is represented by a move vector [0,1,0] and so on...
            # At this point I will set the rates to a default value of 1.0.
            # We will re-set the rates elsewhere based on a specific model.
            for atomnum = 2:numatoms
                thisatom = atomreg[atomnum]
                thisatomprocl = thisatom.processlist
                neighbors = Neighborhood(cubcryst, thisatom)
                defaultrate = 1.0
                # Initialize all 8 of the atom processes.
                # N, S, E, W
                push!(thisatomprocl, Process(atomnum, [0,1,0], 0, 0, 0.0))  #1
                push!(thisatomprocl, Process(atomnum, [0,-1,0], 0, 0, 0.0)) #2
                push!(thisatomprocl, Process(atomnum, [1,0,0], 0, 0, 0.0))  #3
                push!(thisatomprocl, Process(atomnum, [-1,0,0], 0, 0, 0.0)) #4
                # Down N, S, E, W
                push!(thisatomprocl, Process(atomnum, [0,1,-1], 0, 0, 0.0))  #5
                push!(thisatomprocl, Process(atomnum, [0,-1,-1], 0, 0, 0.0)) #6
                push!(thisatomprocl, Process(atomnum, [1,0,-1], 0, 0, 0.0))  #7
                push!(thisatomprocl, Process(atomnum, [-1,0,-1], 0, 0, 0.0)) #8 
                # Atom will be trapped if there is one above it.  No processes are allowed.
                if neighbors[2,2,3] != 0 continue end
                # Otherwise check the surrounding atoms for unfilled sites to jump to.
                if neighbors[2,3,1] == 0
                    thisatomprocl[5].rate = defaultrate
                    push!(proclDownNorth,thisatomprocl[5])
                    countneighbors!(cubcryst,thisatomprocl[5])
                elseif neighbors[2,3,2] == 0
                    thisatomprocl[1].rate = defaultrate
                    push!(proclNorth,thisatomprocl[1])
                    countneighbors!(cubcryst,thisatomprocl[1])
                end
                if neighbors[2,1,1] == 0
                    thisatomprocl[6].rate = defaultrate
                    push!(proclDownSouth,thisatomprocl[6])
                    countneighbors!(cubcryst,thisatomprocl[6])
                elseif neighbors[2,1,2] == 0
                    thisatomprocl[2].rate = defaultrate
                    push!(proclSouth,thisatomprocl[2])
                    countneighbors!(cubcryst,thisatomprocl[2])
                end
                if neighbors[3,2,1] == 0
                    thisatomprocl[7].rate = defaultrate
                    push!(proclDownEast,thisatomprocl[7])
                    countneighbors!(cubcryst,thisatomprocl[7])
                elseif neighbors[3,2,2] == 0
                    thisatomprocl[3].rate = defaultrate
                    push!(proclEast, thisatomprocl[3])
                    countneighbors!(cubcryst,thisatomprocl[3])
                end
                if neighbors[1,2,1] == 0
                    thisatomprocl[8].rate = defaultrate
                    push!(proclDownWest,thisatomprocl[8])
                    countneighbors!(cubcryst,thisatomprocl[8])
                elseif neighbors[1,2,2] == 0
                    thisatomprocl[4].rate = defaultrate
                    push!(proclWest,thisatomprocl[4])
                    countneighbors!(cubcryst,thisatomprocl[4])
                end
            end
        end
        # The struct is created by new().
        new(proclNorth,proclSouth, proclEast, proclWest, proclDownNorth, proclDownSouth, proclDownEast, proclDownWest)
    end
end

"""

Get the CubicProcessList at the specified index.

# Examples
```jldoctest

x = cryst.Cubic(5,5)
for i=1:10 cryst.Deposit(x) end
y = CubicProcessList(x)

julia> y[1]
6-element Array{Crystal.Process,1}:
 Crystal.Process(3, [0, 1, 0], 0, 0, 1.0)
 Crystal.Process(4, [0, 1, 0], 0, 0, 1.0)
 Crystal.Process(6, [0, 1, 0], 0, 0, 1.0)
 Crystal.Process(9, [0, 1, 0], 0, 0, 1.0)
 Crystal.Process(10, [0, 1, 0], 0, 0, 1.0)
 Crystal.Process(11, [0, 1, 0], 0, 0, 1.0)
```
"""
function Processes.getindex(list::CubicProcessList, index::Int64)
    if index < 1 throw(BoundsError(list,index)) end
    if index == 1 return list.proclNorth end
    if index == 2 return list.proclSouth end
    if index == 3 return list.proclEast end
    if index == 4 return list.proclWest end
    if index == 5 return list.proclDownNorth end
    if index == 6 return list.proclDownSouth end
    if index == 7 return list.proclDownEast end
    if index == 8 return list.proclDownWest end
    if index > 8 throw(BoundsError(list,index)) end
end

"""
Inserts the number of nearest neighbors into the Process.
"""
function countneighbors!(x::Cubic, myprocess::Process)
    wrap(i) = i<1 ? x.size : (i>x.size ? 1 : i)
    thisatom = x.AtomRegistry[myprocess.AtomID]
    neighbors = Neighborhood(x, thisatom)
    oldNN = (neighbors[1,2,2]>0) + (neighbors[3,2,2]>0) + (neighbors[2,1,2]>0) + (neighbors[2,3,2]>0)
    mymovevector = myprocess.movevector
    oldsite = [thisatom.xpos, thisatom.ypos, thisatom.zpos]
    new = broadcast(wrap, oldsite + mymovevector)
    newNN = (x.world[wrap(new[1]+1),wrap(new[2]),new[3]]>0) + (x.world[wrap(new[1]-1),wrap(new[2]),new[3]]>0) + 
        (x.world[wrap(new[1]),wrap(new[2]+1),new[3]]>0) + (x.world[wrap(new[1]),wrap(new[2]-1),new[3]]>0)
    if mymovevector[3]==0 newNN -= 1 end # Counts itself.  Correct for overcounting.
    myprocess.initialNN = oldNN
    myprocess.finalNN = newNN
    return [oldNN, newNN]
end

"""
Return the nearest neighbor arrays before and after a given process.
"""
function findneighbors(x::Cubic, myprocess::Process)
    # Make sure to include the atom directly beneath for both the old and new sites.
    # Also include the atom itself.
    wrap(i) = i<1 ? x.size : (i>x.size ? 1 : i)
    thisatom = x.AtomRegistry[myprocess.AtomID]
    neighbors = Neighborhood(x, thisatom)
    oldNN = [neighbors[1,2,2], neighbors[3,2,2], neighbors[2,1,2], neighbors[2,3,2], neighbors[2,2,1], neighbors[2,2,2]]
    mymovevector = myprocess.movevector
    oldsite = [thisatom.xpos, thisatom.ypos, thisatom.zpos]
    new = broadcast(wrap, oldsite + mymovevector)
    newNN = []
    newatom = x.world[wrap(new[1]+1),wrap(new[2]),new[3]]
    if  newatom != thisatom  push!(newNN, newatom) end
    newatom = x.world[wrap(new[1]-1),wrap(new[2]),new[3]]
    if  newatom != thisatom  push!(newNN, newatom) end
    newatom = x.world[wrap(new[1]),wrap(new[2]+1),new[3]]
    if  newatom != thisatom  push!(newNN, newatom) end
    newatom = x.world[wrap(new[1]),wrap(new[2]-1),new[3]]
    if  newatom != thisatom  push!(newNN, newatom) end
    newatom = x.world[wrap(new[1]),wrap(new[2]),new[3]-1]
    if  newatom != thisatom  push!(newNN, newatom) end
    return [oldNN, newNN]
end

"""
Applies the irreversible sticking criterion to the process rate.
"""
function IrreversibleStickingRates!(x::Crystal.AbstractCrystal, myprocesslist::AbstractProcessList) 
    for i = 1:8
        todelete = []
        for myindex = 1:length(myprocesslist[i])
            myprocess = myprocesslist[i][myindex]
            if myprocess.initialNN > 0
                myprocess.rate = 0.0
                push!(todelete, myindex)
            end
        end
        for myindex2 = length(todelete):-1:1
                splice!(myprocesslist[i], todelete[myindex2])
        end
    end
    return myprocesslist
end

"""
Base class for the process rates.
"""
abstract type AbstractProcessRates <: Any
end

"""
Total process rates for each process list.
"""
mutable struct CubicProcessRates <: AbstractProcessRates
    processtotals::Array{Float64,1}
    ratesNorth::Array{Float64,1}
    ratesSouth::Array{Float64,1}
    ratesEast::Array{Float64,1}
    ratesWest::Array{Float64,1}
    ratesDownNorth::Array{Float64,1}
    ratesDownSouth::Array{Float64,1}
    ratesDownEast::Array{Float64,1}
    ratesDownWest::Array{Float64,1}
    depositionRate::Array{Float64,1}
    indexarray::Array{Int64,1}
    # Constructor
    function CubicProcessRates(myprocl::CubicProcessList)
        indexarray = collect(1:9)
        processtotals = zeros(Float64,9)
        depositionRate = zeros(Float64,1) # It's a vector.  Maybe for future pulsed deposition.
        depositionRate[1] = processtotals[9] = 1.0 
        ratesNorth = []; ratesSouth = []; ratesEast = []; ratesWest = []
        ratesDownNorth = []; ratesDownSouth = []; ratesDownEast = []; ratesDownWest = []
        #
        for myprocess in myprocl[1]
            push!(ratesNorth, myprocess.rate)
            processtotals[1] += myprocess.rate
        end
        for myprocess in myprocl[2]
            push!(ratesSouth, myprocess.rate)
            processtotals[2] += myprocess.rate
        end
        for myprocess in myprocl[3]
            push!(ratesEast, myprocess.rate)
            processtotals[3] += myprocess.rate
        end
        for myprocess in myprocl[4]
            push!(ratesWest, myprocess.rate)
            processtotals[4] += myprocess.rate
        end
        for myprocess in myprocl[5]
            push!(ratesDownNorth, myprocess.rate)
            processtotals[5] += myprocess.rate
        end
        for myprocess in myprocl[6]
            push!(ratesDownSouth, myprocess.rate)
            processtotals[6] += myprocess.rate
        end
        for myprocess in myprocl[7]
            push!(ratesDownEast, myprocess.rate)
            processtotals[7] += myprocess.rate
        end
        for myprocess in myprocl[8]
            push!(ratesDownWest, myprocess.rate)
            processtotals[8] += myprocess.rate
        end
        # The struct is created by new().
        new(processtotals, ratesNorth, ratesSouth, ratesEast, ratesWest, ratesDownNorth,  
            ratesDownSouth, ratesDownEast, ratesDownWest, depositionRate, indexarray)
    end
end

"""
Get the CubicProcessRates at the specified index.
"""
function Processes.getindex(list::CubicProcessRates, index::Int64)
    if index < 1 throw(BoundsError(list,index)) end
    if index == 1 return list.ratesNorth end
    if index == 2 return list.ratesSouth end
    if index == 3 return list.ratesEast end
    if index == 4 return list.ratesWest end
    if index == 5 return list.ratesDownNorth end
    if index == 6 return list.ratesDownSouth end
    if index == 7 return list.ratesDownEast end
    if index == 8 return list.ratesDownWest end
    if index == 9 return list.depositionRate end
    if index > 9 throw(BoundsError(list,index)) end
end

"""
Select a process at random from the process lists, weighted by the process rates.
"""
function getprocess(myprocesslist::AbstractProcessList, myprocessrates::AbstractProcessRates)
    # Select one of the 8 hop directions (1-8) or deposition (9).
    processchoices = myprocessrates.indexarray # this is the index array 1:9.
    select1 = sample(processchoices, weights(myprocessrates.processtotals))
    # If the process selected is deposition, return the integer 9.
    if select1 == 9 return select1 end
    # Otherwise, select one process from the process list selected above.
    listlen = length(myprocesslist[select1])
    select2 = sample(processchoices[1:listlen], weights(myprocessrates[select1]))
    return myprocesslist[select1][select2]
end

"""
Compare two Processes.  The one with the smaller atom number is less.
"""
Processes.isless(A::Process, B::Process) = A.AtomID < B.AtomID

"""
Update the process lists for a given process.  
"""
function processController!(x::Crystal.AbstractCrystal, processtodo::Process, y::AbstractProcessList, 
    z::AbstractProcessRates)
    # make a list of the neighboring atoms before and after the move.
    (oldnn, newnn) = findneighbors(x, processtodo)
    neighbors = filter(a->a>1, [oldnn; newnn])
    sort!(neighbors); unique!(neighbors)
    # loop through neighbors and record the old process lists.
    # the current atom to move is also in the neighbors list.
    processlist = []
    for atomnum in neighbors
        push!(processlist, x.AtomRegistry[atomnum].processlist)
    end
    #oldprocesslist = deepcopy(processlist)
    #print("\n",oldprocesslist,"\n")
    # Move the atom to the new location.
    wrap(i) = i<1 ? x.size : (i>x.size ? 1 : i)
    thisatom = processtodo.AtomID
    atomtomove = x.AtomRegistry[thisatom ]
    mv = processtodo.movevector
    x.world[atomtomove.xpos, atomtomove.ypos, atomtomove.zpos] = 0 # free up the old position.
    newxpos = wrap(atomtomove.xpos+mv[1])
    newypos = wrap(atomtomove.ypos+mv[2])
    newzpos = atomtomove.zpos+mv[3]
    x.world[newxpos, newypos, newzpos] = thisatom  # update the crystal at the new position.
    atomtomove.xpos = newxpos # update the atom itself.
    atomtomove.ypos = newypos
    atomtomove.zpos = newzpos
    # loop through again, update the processes that each atom holds.
    # This has a side effect on the main process lists since NN counts and rates will change. 
    # Manage the main process list -- drop processes with rates of 0.0, and add any new 
    # processes that were not previously allowed or that had rates of 0.0.
    for atomnum in neighbors
        updateAtomRegistry!(x, atomnum)
        updateCubicProcessListandRates!(x, atomnum, y, z)
    end
    # Update the process rate totals.
    for int in 1:8
        z.processtotals[int] = sum(z[int])
    end
end

"""
Update the process list for a given atom in a Cubic crystal.
"""
function updateAtomRegistry!(x::Cubic, atomnum::Int64)
    thisatom = x.AtomRegistry[atomnum]
    neighbors = Neighborhood(x, thisatom)
    thisatomprocl = thisatom.processlist
    defaultrate =  1.0
    # Record the old rates so we can track what changed.
    oldrates = []
    for int = 1:8
        push!(oldrates, thisatomprocl[int].rate)
    end
    print(atomnum, " oldrates = ", oldrates, "\n")
    # Atom will be trapped if there is one above it.  No processes are allowed.
    if neighbors[2,2,3] != 0 
        for int = 1:8
            thisatomprocl[int].initialNN = 0
            thisatomprocl[int].finalNN = 0
            thisatomprocl[int].rate = 0.0
        end
        return thisatomprocl
    end
    # Otherwise check the surroundings for unfilled sites to jump to.
    if neighbors[2,3,1] == 0
        thisatomprocl[5].rate = defaultrate
        thisatomprocl[1].rate = 0.0
        countneighbors!(x,thisatomprocl[5])
        thisatomprocl[1].initialNN = 0
        thisatomprocl[1].finalNN = 0
    elseif neighbors[2,3,2] == 0
        thisatomprocl[5].rate = 0.0
        thisatomprocl[1].rate = defaultrate
        countneighbors!(x,thisatomprocl[1])
        thisatomprocl[5].initialNN = 0
        thisatomprocl[5].finalNN = 0
    else
        thisatomprocl[5].rate = 0.0
        thisatomprocl[1].rate = 0.0
        thisatomprocl[5].initialNN = 0
        thisatomprocl[5].finalNN = 0
        thisatomprocl[1].initialNN = 0
        thisatomprocl[1].finalNN = 0
    end
    if neighbors[2,1,1] == 0
        thisatomprocl[6].rate = defaultrate
        thisatomprocl[2].rate = 0.0
        countneighbors!(x,thisatomprocl[6])
        thisatomprocl[2].initialNN = 0
        thisatomprocl[2].finalNN = 0
    elseif neighbors[2,1,2] == 0
        thisatomprocl[6].rate = 0.0
        thisatomprocl[2].rate = defaultrate
        countneighbors!(x,thisatomprocl[2])
        thisatomprocl[6].initialNN = 0
        thisatomprocl[6].finalNN = 0
    else
        thisatomprocl[6].rate = 0.0
        thisatomprocl[2].rate = 0.0
        thisatomprocl[6].initialNN = 0
        thisatomprocl[6].finalNN = 0
        thisatomprocl[2].initialNN = 0
        thisatomprocl[2].finalNN = 0
    end
    if neighbors[3,2,1] == 0
        thisatomprocl[7].rate = defaultrate
        thisatomprocl[3].rate = 0.0
        countneighbors!(x,thisatomprocl[7])
        thisatomprocl[3].initialNN = 0
        thisatomprocl[3].finalNN = 0
    elseif neighbors[3,2,2] == 0
        thisatomprocl[7].rate = 0.0
        thisatomprocl[3].rate = defaultrate
        countneighbors!(x,thisatomprocl[3])
        thisatomprocl[7].initialNN = 0
        thisatomprocl[7].finalNN = 0
    else
        thisatomprocl[7].rate = 0.0
        thisatomprocl[3].rate = 0.0
        thisatomprocl[7].initialNN = 0
        thisatomprocl[7].finalNN = 0
        thisatomprocl[3].initialNN = 0
        thisatomprocl[3].finalNN = 0
    end
    if neighbors[1,2,1] == 0
        thisatomprocl[8].rate = defaultrate
        thisatomprocl[4].rate = 0.0
        countneighbors!(x,thisatomprocl[8])
        thisatomprocl[4].initialNN = 0
        thisatomprocl[4].finalNN = 0
    elseif neighbors[1,2,2] == 0
        thisatomprocl[8].rate = 0.0
        thisatomprocl[4].rate = defaultrate
        countneighbors!(x,thisatomprocl[4])
        thisatomprocl[8].initialNN = 0
        thisatomprocl[8].finalNN = 0
    else
        thisatomprocl[8].rate = 0.0
        thisatomprocl[4].rate = 0.0
        thisatomprocl[8].initialNN = 0
        thisatomprocl[8].finalNN = 0
        thisatomprocl[4].initialNN = 0
        thisatomprocl[4].finalNN = 0
    end
    #
    return thisatomprocl
end

"""
Drops atoms with zero rates from the process list and adds those that previously had 
rates of 0.0 but are currently non-zero.  Also updates the rates in the ProcessRates arrays.
"""
function updateCubicProcessListandRates!(x::Cubic, atomnum::Int64, y::AbstractProcessList, z::AbstractProcessRates)
    thisatom = x.AtomRegistry[atomnum]
    thisatomprocl = thisatom.processlist
    newrates = []
    for int = 1:8
        push!(newrates, thisatomprocl[int].rate)
    end
    print(atomnum, " newrates = ", newrates, "\n")
    #
    # Drop processes with newrates of 0.0 from the main process list.
    print(thisatomprocl,"\n")
    for i in 1:8
        if newrates[i] == 0  
            p = searchsorted(y[i], thisatomprocl[i])
            if first(p) == last(p) 
                q = deleteat!(y[i], first(p)) 
                r = deleteat!(z[i], first(p)) 
            end
        end
    end
    #
    # Add processes that have changed from newrate of 0.0 to non-zero.
    for i in 1:8
        if newrates[i] > 0  
            p = searchsorted(y[i], thisatomprocl[i])
            pf = first(p)
            pl = last(p)
            if pf > pl 
                if pf > length(y[i])
                    q = splice!(y[i], pl, thisatomprocl[i]) 
                    r = splice!(z[i], pl, newrates[i]) 
                else
                    q = splice!(y[i], pf, thisatomprocl[i]) 
                    r = splice!(z[i], pf, newrates[i]) 
                end
            end
        end
    end
    # The rates will be summed outside of this function.
    return nothing
end

end # end of module Processes