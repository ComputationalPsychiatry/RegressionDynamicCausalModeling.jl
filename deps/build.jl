using Printf
using Scratch
using UUIDs

# get path to scratch folder containing binary file
const euler_integration_bin = get_scratch!(
    UUID("4fb452e6-e726-4794-90dc-adda47b23488"), "euler_int_bin"
)

if !isfile(joinpath(euler_integration_bin, "libdcm_euler_integration.so")) # check if binary already exists
    if !isnothing(Sys.which("gcc")) # check if gcc is installed
        output = joinpath(euler_integration_bin, "libdcm_euler_integration.so") # define output path
        input = joinpath(@__DIR__, "libdcm_euler_integration.c") # path to c file
        s = @sprintf("gcc -fPIC -shared -o %s %s", output, input) # command for compiling c file
        c = Cmd(convert(Vector{String}, split(s)))
        run(c) # run compilation
    else
        @warn "gcc not found. Cannot compile libdcm_euler_integration.c"
    end
end
