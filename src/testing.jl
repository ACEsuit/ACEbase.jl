
module Testing

using Test

using ACEbase.FIO: read_dict, write_dict, save_dict, load_dict

export print_tf, test_fio, h0, h1, h2, h3


function h0(str)
   dashes = "≡"^(length(str)+4)
   printstyled(dashes, color=:magenta); println()
   printstyled("  "*str*"  ", bold=true, color=:magenta); println()
   printstyled(dashes, color=:magenta); println()
end

function h1(str)
   dashes = "="^(length(str)+2)
   printstyled(dashes, color=:magenta); println()
   printstyled(" " * str * " ", bold=true, color=:magenta); println()
   printstyled(dashes, color=:magenta); println()
end

function h2(str)
   dashes = "-"^length(str)
   printstyled(dashes, color=:magenta); println()
   printstyled(str, bold=true, color=:magenta); println()
   printstyled(dashes, color=:magenta); println()
end

h3(str) = (printstyled(str, bold=true, color=:magenta); println())


print_tf(::Test.Pass) = printstyled("+", bold=true, color=:green)
print_tf(::Test.Fail) = printstyled("-", bold=true, color=:red)
print_tf(::Tuple{Test.Error,Bool}) = printstyled("x", bold=true, color=:magenta)


"""
`test_fio(obj): `  performs two tests:

- encodes `obj` as a Dict using `write_dict`, then decodes it using
`read_dict` and tests whether the two objects are equivalent using `==`
- writes `Dict` to file then reads it and decodes it and test the result is
again equivalent to `obj`

The two results are returned as Booleans.
"""
function test_fio(obj)
   D = write_dict(obj)
   test1 = (obj == read_dict(D))
   tmpf = tempname() * ".json"
   save_dict(tmpf, D)
   test2 = (obj == read_dict(load_dict(tmpf)))
   return test1, test2
end



end
