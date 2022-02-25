
using ACEbase, BenchmarkTools
using ACEbase: acquire!, release! 

pool = ACEbase.ObjectPools.FlexibleVectorPool()
spool = ACEbase.ObjectPools.VectorPool{Float64}()

x = acquire!(pool, 1000, Float64)
release!(pool, x)

xs = acquire!(spool, 1000, Float64)
release!(spool, xs)

function runn(pool, N, len, T) 
   for n = 1:N 
      x = acquire!(pool, len, T)
      release!(pool, x)
   end
   return nothing 
end

@btime runn($pool, 1_000, 1_000, Float64)

@btime runn(spool, 1_000, 1_000, Float64)


@code_warntype acquire!(pool, 1000, Float64)

@code_warntype release!(pool, x)

@code_warntype runn(pool, 1_000, 1_000, Float64) 

using Profile

Profile.clear()

@profile runn(pool, 1_000_000, 1_000, Float64)

Profile.print()
