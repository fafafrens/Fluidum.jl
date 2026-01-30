using DelimitedFiles
using Interpolations
using Fluidum

function fastreso_reader(filename)
   a = readdlm(filename, '\t', Float64, '\n', comments = true, comment_char = '#')
   @show size(a)
   @show pbar=a[:,1]
   @show m = a[1,2]
   Ep = sqrt.(pbar.^2 .+ m^2)
   itp = [interpolate((Ep,),a[:,n]./pbar[:],Gridded(Linear())) for n in 3:4] #interpolation of all the columns from 3 on, through linear interpolation
   ext = extrapolate.(itp, Ref(Flat()))   
   return ext, minimum(Ep), maximum(Ep)
end
