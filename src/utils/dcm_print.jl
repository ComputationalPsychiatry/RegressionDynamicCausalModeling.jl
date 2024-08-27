
function Base.show(io::IO, dcm::LinearDCM)
    if isnothing(dcm.Y) & !isnothing(dcm.U)
        print(io, "Linear DCM\n",
            "a:     ",size(dcm.a,1),"x",size(dcm.a,2)," matrix\n",
            "c:     ",size(dcm.c,1),"x",size(dcm.c,2), " matrix\n",
            "scans: ",dcm.scans,"\n",
            "nr:    ",dcm.nr,"\n",
            "---------------------------------------------\n",
            dcm.U,
            "---------------------------------------------\n",
            "Y (BOLD signal)\n",
            "   empty\n",
            "---------------------------------------------\n",
            dcm.Ep)
    elseif !isnothing(dcm.Y) & !isnothing(dcm.U)
        print(io, "Linear DCM\n",
            "a:     ",size(dcm.a,1),"x",size(dcm.a,2)," matrix\n",
            "c:     ",size(dcm.c,1),"x",size(dcm.c,2), " matrix\n",
            "scans: ",dcm.scans,"\n",
            "nr:    ",dcm.nr,"\n",
            "---------------------------------------------\n",
            dcm.U,
            "---------------------------------------------\n",
            dcm.Y,
            "---------------------------------------------\n",
            dcm.Ep)
    elseif !isnothing(dcm.Y) & isnothing(dcm.U)
        print(io, "Linear DCM\n",
            "a:     ",size(dcm.a,1),"x",size(dcm.a,2)," matrix\n",
            "c:     ",size(dcm.c,1),"x",size(dcm.c,2), " matrix\n",
            "scans: ",dcm.scans,"\n",
            "nr:    ",dcm.nr,"\n",
            "---------------------------------------------\n",
            "U (Input)\n",
            "   empty\n",
            "---------------------------------------------\n",
            dcm.Y,
            "---------------------------------------------\n",
            dcm.Ep)
    else
        print(io, "Linear DCM\n",
            "a:     ",size(dcm.a,1),"x",size(dcm.a,2)," matrix\n",
            "c:     ",size(dcm.c,1),"x",size(dcm.c,2), " matrix\n",
            "scans: ",dcm.scans,"\n",
            "nr:    ",dcm.nr,"\n",
            "---------------------------------------------\n",
            "U (Input)\n",
            "   empty\n",
            "---------------------------------------------\n",
            "Y (BOLD signal)\n",
            "   empty\n",
            "---------------------------------------------\n",
            dcm.Ep)
    end
end

function Base.show(io::IO, dcm::BiLinearDCM)
    if isnothing(dcm.Y)
        print(io, "Bilinear DCM\n",
            "a:     ",size(dcm.a,1),"x",size(dcm.a,2)," matrix\n",
            "b:     ",size(dcm.b,1),"x",size(dcm.b,2),"x",size(dcm.b,3)," matrix\n",
            "c:     ",size(dcm.c,1),"x",size(dcm.c,2), " matrix\n",
            "scans: ",dcm.scans,"\n",
            "nr:    ",dcm.nr,"\n",
            "---------------------------------------------\n",
            dcm.U,
            "---------------------------------------------\n",
            "Y (BOLD signal)\n",
            "   empty\n",
            "---------------------------------------------\n",
            dcm.Ep)
    else
        print(io, "Bilinear DCM\n",
            "a:     ",size(dcm.a,1),"x",size(dcm.a,2)," matrix\n",
            "b:     ",size(dcm.b,1),"x",size(dcm.b,2),"x",size(dcm.b,3)," matrix\n",
            "c:     ",size(dcm.c,1),"x",size(dcm.c,2), " matrix\n",
            "scans: ",dcm.scans,"\n",
            "nr:    ",dcm.nr,"\n",
            "---------------------------------------------\n",
            dcm.U,
            "---------------------------------------------\n",
            dcm.Y,
            "---------------------------------------------\n",
            dcm.Ep)
    end
end

function Base.show(io::IO, dcm::NonLinearDCM)
    if isnothing(dcm.Y)
        print(io, "Nonlinear DCM\n",
            "a:     ",size(dcm.a,1),"x",size(dcm.a,2)," matrix\n",
            "b:     ",size(dcm.b,1),"x",size(dcm.b,2),"x",size(dcm.b,3)," matrix\n",
            "c:     ",size(dcm.c,1),"x",size(dcm.c,2), " matrix\n",
            "d:     ",size(dcm.d,1),"x",size(dcm.d,2),"x",size(dcm.d,3)," matrix\n",
            "scans: ",dcm.scans,"\n",
            "nr:    ",dcm.nr,"\n",
            "---------------------------------------------\n",
            dcm.U,
            "---------------------------------------------\n",
            "Y (BOLD signal)\n",
            "   empty\n",
            "---------------------------------------------\n",
            dcm.Ep)
    else
        print(io, "Nonlinear DCM\n",
            "a:     ",size(dcm.a,1),"x",size(dcm.a,2)," matrix\n",
            "b:     ",size(dcm.b,1),"x",size(dcm.b,2),"x",size(dcm.b,3)," matrix\n",
            "c:     ",size(dcm.c,1),"x",size(dcm.c,2), " matrix\n",
            "d:     ",size(dcm.d,1),"x",size(dcm.d,2),"x",size(dcm.d,3)," matrix\n",
            "scans: ",dcm.scans,"\n",
            "nr:    ",dcm.nr,"\n",
            "---------------------------------------------\n",
            dcm.U,
            "---------------------------------------------\n",
            dcm.Y,
            "---------------------------------------------\n",
            dcm.Ep)
    end
end

function Base.show(io::IO, Ep::TrueParamLinear)
    @printf(io,"Ep (True parameters)
    A: %d x %d matrix
    C: %d x %d matrix
    transit: %.2f ... %.2f
    decay:   %.2f ... %.2f
    epsilon: %.2f\n",size(Ep.A,1),size(Ep.A,2),size(Ep.C,1),
    size(Ep.C,2),Ep.transit[1],Ep.transit[end],
    Ep.decay[1],Ep.decay[end],Ep.epsilon)
end

function Base.show(io::IO, Ep::TrueParamBiLinear)
    @printf(io,"Ep (True parameters)
    A: %d x %d matrix
    B: %d x %d x %d matrix
    C: %d x %d matrix
    transit: %.2f ... %.2f
    decay:   %.2f ... %.2f
    epsilon: %.2f\n",size(Ep.A,1),size(Ep.A,2),
    size(Ep.B,1),size(Ep.B,2),size(Ep.B,3),size(Ep.C,1),
    size(Ep.C,2),Ep.transit[1],Ep.transit[end],
    Ep.decay[1],Ep.decay[end],Ep.epsilon)
end

function Base.show(io::IO, Ep::TrueParamNonLinear)
    @printf(io,"Ep (True parameters)
    A: %d x %d matrix
    B: %d x %d x %d matrix
    C: %d x %d matrix
    D: %d x %d x %d matrix
    transit: %.2f ... %.2f
    decay:   %.2f ... %.2f
    epsilon: %.2f\n",size(Ep.A,1),size(Ep.A,2),
    size(Ep.B,1),size(Ep.B,2),size(Ep.B,3),size(Ep.C,1),
    size(Ep.C,2),size(Ep.D,1),size(Ep.D,2),size(Ep.D,3),Ep.transit[1],Ep.transit[end],
    Ep.decay[1],Ep.decay[end],Ep.epsilon)
end

function Base.show(io::IO, U::InputU)
    print(io, "U (Input)\n",
                "   u:  ",size(U.u,1),"x",size(U.u,2)," matrix\n",
                "   dt: ",U.dt,"s\n",
                "   names: ",U.name[1],",...,",U.name[end],"\n")
end

function Base.show(io::IO, Y::BoldY)
    if isnothing(Y.y)
        print(io, "Y (BOLD signal)\n",
                    "   y: empty\n",
                    "   dt: ",Y.dt,"s\n",
                    "   names: empty\n")
    else
        y = Matrix{Float64}(Y.y) # needed otherwise JET throws an error
        name = Vector{String}(Y.name)
        print(io, "Y (BOLD signal)\n",
                    "   y:  ",size(y,1),"x",size(y,2)," matrix\n",
                    "   dt: ",Y.dt,"s\n",
                    "   names: ",name[1],",...,",name[end],"\n")
    end
end

function Base.show(io::IO, Conf::Confound)
    print(io, "Confounds\n",
                "   X0:    ",size(Conf.X0,1),"x",size(Conf.X0,2)," matrix\n",
                "   names: ",Conf.name[1],",...,",Conf.name[end],"\n")
end

function Base.show(io::IO, rdcm::RigidRdcm)
    #if isnothing(rdcm.Conf)
    #    print(io, "rigid rDCM\n",
    #                "a:     ",size(rdcm.a,1),"x",size(rdcm.a,2)," matrix\n",
    #                "c:     ",size(rdcm.c,1),"x",size(rdcm.c,2), " matrix\n",
    #                "scans: ",rdcm.scans,"\n",
    #                "nr:    ",rdcm.nr,"\n",
    #                "HRF:   ",length(rdcm.hrf)," element vector\n",
    #                "---------------------------------------------\n",
    #                rdcm.U,
    #                "---------------------------------------------\n",
    #                rdcm.Y,
    #                "---------------------------------------------\n",
    #                rdcm.Ep)
    #else
        print(io, "rigid rDCM\n",
                    "a:     ",size(rdcm.a,1),"x",size(rdcm.a,2)," matrix\n",
                    "c:     ",size(rdcm.c,1),"x",size(rdcm.c,2), " matrix\n",
                    "scans: ",rdcm.scans,"\n",
                    "nr:    ",rdcm.nr,"\n",
                    "HRF:   ",length(rdcm.hrf)," element vector\n",
                    "---------------------------------------------\n",
                    rdcm.U,
                    "---------------------------------------------\n",
                    rdcm.Y,
                    "---------------------------------------------\n",
                    rdcm.Ep,
                    "---------------------------------------------\n",
                    rdcm.Conf)
    #end
end

function Base.show(io::IO, out::RigidOutput)

    @printf(io,"rigid rDCM output
    F:   %.2f
    F_r: %.2f ... %.2f
    iterations until convergence per region: %d ... %d
    Posteriors:
        α: %.2f,...,%.2f
        β: %.2f,...,%.2f
        μ: %d x %d matrix
        Σ: %d element vector of matrices",out.F,out.F_r[1],out.F_r[end],
        out.iter_all[1],out.iter_all[end],
        out.a_all[1],out.a_all[end],out.b_all[1],out.b_all[end],
        size(out.m_all,1),size(out.m_all,2),length(out.Σ_all))
end

function Base.show(io::IO, rdcm::SparseRdcm)
    #if isnothing(rdcm.Conf)
    #    print(io, "sparse rDCM\n",
    #                "a:     ",size(rdcm.a,1),"x",size(rdcm.a,2)," matrix\n",
    #                "c:     ",size(rdcm.c,1),"x",size(rdcm.c,2), " matrix\n",
    #                "scans: ",rdcm.scans,"\n",
    #                "nr:    ",rdcm.nr,"\n",
    #                "HRF:   ",length(rdcm.hrf)," element vector\n",
    #                "---------------------------------------------\n",
    #                rdcm.U,
    #                "---------------------------------------------\n",
    #                rdcm.Y,
    #                "---------------------------------------------\n",
    #                rdcm.Ep,
    #                "---------------------------------------------\n",
    #                "inform_p0:", rdcm.inform_p0,"\n",
    #                "p0:", rdcm.p0,"\n")
    #else
        print(io, "sparse rDCM\n",
                "a:     ",size(rdcm.a,1),"x",size(rdcm.a,2)," matrix\n",
                "c:     ",size(rdcm.c,1),"x",size(rdcm.c,2), " matrix\n",
                "scans: ",rdcm.scans,"\n",
                "nr:    ",rdcm.nr,"\n",
                "HRF:   ",length(rdcm.hrf)," element vector\n",
                "---------------------------------------------\n",
                rdcm.U,
                "---------------------------------------------\n",
                rdcm.Y,
                "---------------------------------------------\n",
                rdcm.Ep,
                "---------------------------------------------\n",
                rdcm.Conf,
                "---------------------------------------------\n",
                "inform_p0:", rdcm.inform_p0,"\n",
                "p0:", rdcm.p0,"\n")
    #end
end

function Base.show(io::IO, out::SparseOutput)
    @printf(io,"sparse rDCM output
    F:   %.2f
    F_r: %.2f ... %.2f
    iterations until convergence per region: %d ... %d
    Posteriors:
        α: %.2f,...,%.2f
        β: %.2f,...,%.2f
        μ: %d x %d matrix
        Σ: %d element vector of matrices
        Z: %d x %d matrix",out.F,out.F_r[1],out.F_r[end],out.iter_all[1],out.iter_all[end],
        out.a_all[1],out.a_all[end],out.b_all[1],out.b_all[end],
        size(out.m_all,1),size(out.m_all,2),length(out.Σ_all),
        size(out.z_all,1),size(out.z_all,2))
end
