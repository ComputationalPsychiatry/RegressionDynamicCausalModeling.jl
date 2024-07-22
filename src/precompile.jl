@setup_workload begin
    try
        dcm = load_example_DCM()

        opt1 = Options(RigidInversionParams(); synthetic=true, verbose=0, testing=true)
        opt2 = Options(
            SparseInversionParams(; reruns=10, restrictInputs=true);
            synthetic=true,
            verbose=0,
            testing=true,
        )

        @compile_workload begin
            invert(RigidRdcm(dcm), opt1)
            invert(SparseRdcm(dcm; p0=0.15), opt2)
        end
    catch
        @info "Precompilation failed."
    end
end
