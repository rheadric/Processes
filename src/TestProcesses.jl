##
print("Hello world!")
##
push!(LOAD_PATH, "/Users/randallheadrick/Documents/myjulia/Crystal/src")
using Pkg
Pkg.activate("Crystal")
using Crystal
cryst = Crystal
##
push!(LOAD_PATH, "/Users/randallheadrick/Documents/myjulia/Crystal/src")
using Pkg
Pkg.activate("/Users/randallheadrick/Documents/myjulia/Processes")
import Crystal
cryst = Crystal
using Processes
##
x = cryst.Cubic(5,5)
for i=1:20 cryst.Deposit(x) end
myfig = cryst.PlotCryst(x)
##
close(myfig)
##
myfig = cryst.PlotCryst(x, x.AtomRegistry[12])
##
cryst.PlotCryst(x, cryst.Atom(1,4,1,1,[]))
##
x = cryst.Cubic(5,5)
for i=1:10 cryst.Deposit(x) end
y = CubicProcessList(x)
myfig = cryst.PlotCryst(x)
##
IrreversibleStickingRates!(x,y)
##
z = CubicProcessRates(y)
##
p = getprocess(y,z)
##
n = findneighbors(x,p)
##
processController!(x,p,y,z)
##
close(myfig)
myfig = cryst.PlotCryst(x)
##
