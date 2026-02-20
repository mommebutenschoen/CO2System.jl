using CO2System
using Test

@testset "CO2System.jl" begin

    # Inputs to Carbonate System Call:
    inputs = (
        S = 35.0,
        T = 18.0,
        Rho = 1025.711,
        PO4 = 0.0,
        Sil = 0.0,
        DIC = 2000.0,
        TA = 2300.0,
        )
    opt_inputs = (
        patm = 1013.25,
        p0 = 100.0,
        pH0 = -1.0,
        maxit = 50,
        )
    all_inputs = merge(inputs, opt_inputs)

    # Reference output from previous versions:
    prev_o = (CO2=10.325234153494868,HCO3=1781.9086787891802,CO3=207.76639106353602,pCO2=316.96080477176974,OmegaA=3.160654515808703,
        OmegaC=4.882111778502531,fCO2=305.1809935500387,pH=8.143047856257516,Hplus=7.1936970405460915e-9)

    # Reference values from MOCSY test function:
    mocsy_i = (T=18.0,S=35.0,pr_in=100.0,Rho=1025.71105957,PO4=0.0,Si=0.0,DIC=2000.0,
        TA=2300.0,patm=1013.25)
    mocsy_o = (CO2=10.1729711,HCO3=1779.52,CO3=210.31,pCO2=312.28662109,OmegaA=3.19940853,
        OmegaC=4.94189167,fCO2=300.68057251,pH=8.14892578)

    # Reference values from BFM standalone first iteration:
    bfm_i = (Rho=1025.711,DIC=2000.0,PO4=0.0,pr_in=100.0,Si=0.0,T=18.0,TA=2300.0,
        patm=1013.25,S=35.0)
    bfm_o = (HCO3=1781.9086,fCO2 =305.18082,pCO2=316.96063,OmegaC=4.882114f0,CO3=207.7665,
        CO2=10.325228,pH=8.143048,OmegaA=3.160656)

    function test()
        cs = CarbonateSystem(; inputs..., opt_inputs...)

        println("DIC consistency:")
        println("\tDIC from input:\t",inputs.DIC)
        DIC_o = cs.CO2+cs.HCO3+cs.CO3
        println("\tDIC from output components:\t",DIC_o)
        dDIC=inputs.DIC-DIC_o
        println("\tDifference:\t",dDIC)

        for (k,v) in pairs(all_inputs)
            println(k,": ",v,
                "\tMOCSY reference: ", k in keys(mocsy_i) ? mocsy_i[k] : "N/A",
                "\tBFM: ", k in keys(bfm_i) ? bfm_i[k] : "N/A")
        end
        println()
        for (k,v) in pairs(cs)
            println(k,": ",v,
                "\tMOCSY reference: ", k in keys(mocsy_o) ? mocsy_o[k] : "N/A",
                "\tBFM: ", k in keys(bfm_o) ? bfm_o[k] : "N/A")
        end

        return cs,dDIC ;
    end

    CS,dDIC = test()

    @testset "CO2System.jl" begin
        # error checking
        @test CS.errorflag == 0
        # test that sum of DIC components equals total DIC concentation from input
        @test dDIC < 1.e-3
        # test for consistency with previous versions:
        for (k,v) in pairs(prev_o)
            @test isapprox(CS[k],v,atol=1e-6)
        end
        # @test
        # test for consistency with BFM outputs:

    #    @test YourPackageName.greet_your_package_name() != "Hello world!"
    end
end
