using Distributed
using Base.Filesystem
using DataFrames
using CSV
using Query
using Statistics
using UnicodePlots
using ClusterManagers
using Dates
using DelimitedFiles

## load the packages by covid19abm

#using covid19abm

#addprocs(2, exeflags="--project=.")


#@everywhere using covid19abm

addprocs(SlurmManager(500), N=17, topology=:master_worker, exeflags = "--project=.")
@everywhere using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
@everywhere include("covid19abm.jl")
@everywhere const cv=covid19abm


function run(myp::cv.ModelParameters, nsims=1000, folderprefix="./")
    println("starting $nsims simulations...\nsave folder set to $(folderprefix)")
    dump(myp)
   
    # will return 6 dataframes. 1 total, 4 age-specific 
    cdr = pmap(1:nsims) do x                 
            cv.runsim(x, myp)
    end      

    println("simulations finished")
    println("total size of simulation dataframes: $(Base.summarysize(cdr))")
    ## write the infectors     

    ## write contact numbers
    #writedlm("$(folderprefix)/ctnumbers.dat", [cdr[i].ct_numbers for i = 1:nsims])    
    ## stack the sims together
    allag = vcat([cdr[i].a  for i = 1:nsims]...)
    ag1 = vcat([cdr[i].g1 for i = 1:nsims]...)
    ag2 = vcat([cdr[i].g2 for i = 1:nsims]...)
    ag3 = vcat([cdr[i].g3 for i = 1:nsims]...)
    ag4 = vcat([cdr[i].g4 for i = 1:nsims]...)
    ag5 = vcat([cdr[i].g5 for i = 1:nsims]...)
    ag6 = vcat([cdr[i].g6 for i = 1:nsims]...)
    #mydfs = Dict("all" => allag, "ag1" => ag1, "ag2" => ag2, "ag3" => ag3, "ag4" => ag4, "ag5" => ag5, "ag6" => ag6)
    mydfs = Dict("all" => allag)
    
    ## save at the simulation and time level
    ## to ignore for now: miso, iiso, mild 
    #c1 = Symbol.((:LAT, :ASYMP, :INF, :PRE, :MILD,:IISO, :HOS, :ICU, :DED), :_INC)
    #c2 = Symbol.((:LAT, :ASYMP, :INF, :PRE, :MILD,:IISO, :HOS, :ICU, :DED), :_PREV)
    
    c1 = Symbol.((:LAT, :HOS, :ICU, :DED,:LAT2, :HOS2, :ICU2, :DED2,:LAT3, :HOS3, :ICU3, :DED3,:LAT4, :HOS4, :ICU4, :DED4,:LAT5, :HOS5, :ICU5, :DED5,:LAT6, :HOS6, :ICU6, :DED6), :_INC)
    #c2 = Symbol.((:LAT, :HOS, :ICU, :DED,:LAT2, :HOS2, :ICU2, :DED2), :_PREV)
    for (k, df) in mydfs
        println("saving dataframe sim level: $k")
        # simulation level, save file per health status, per age group
        #for c in vcat(c1..., c2...)
        for c in vcat(c1...)
        #for c in vcat(c2...)
            udf = unstack(df, :time, :sim, c) 
            fn = string("$(folderprefix)/simlevel_", lowercase(string(c)), "_", k, ".dat")
            CSV.write(fn, udf)
        end
        println("saving dataframe time level: $k")
        # time level, save file per age group
        #yaf = compute_yearly_average(df)       
        #fn = string("$(folderprefix)/timelevel_", k, ".dat")   
        #CSV.write(fn, yaf)       
    end
    
   
    R01 = [cdr[i].R01 for i=1:nsims]
    R02 = [cdr[i].R02 for i=1:nsims]
    writedlm(string(folderprefix,"/R01.dat"),R01)
    writedlm(string(folderprefix,"/R02.dat"),R02)

    cov1 = [cdr[i].cov1 for i=1:nsims]
    cov2 = [cdr[i].cov2 for i=1:nsims]
    cov12 = [cdr[i].cov12 for i=1:nsims]
    cov22 = [cdr[i].cov22 for i=1:nsims]
    
    writedlm(string(folderprefix,"/init_iso.dat"),[cdr[i].iniiso for i=1:nsims])
    writedlm(string(folderprefix,"/cov.dat"),[cov1 cov2 cov12 cov22])
    return mydfs
end


function create_folder(ip::cv.ModelParameters,vac="none",province="us")
    
    #RF = string("heatmap/results_prob_","$(replace(string(ip.β), "." => "_"))","_vac_","$(replace(string(ip.vaccine_ef), "." => "_"))","_herd_immu_","$(ip.herd)","_$strategy","cov_$(replace(string(ip.cov_val)))") ## 
    main_folder = "/data/thomas-covid/states_us/"
    #main_folder = "."
   
    RF = string(main_folder,"/results_prob_","$(replace(string(ip.β), "." => "_"))","_herd_immu_","$(ip.herd)","_$vac","_$(ip.file_index)_$(province)") ##  
   
    if !Base.Filesystem.isdir(RF)
        Base.Filesystem.mkpath(RF)
    end
    return RF
end



function run_param_scen_cal(b::Float64,province::String="us",h_i::Int64 = 0,ic1::Int64=1,ic2::Int64=1,ic3::Int64=1,ic4::Int64=1,ic5::Int64=1,ic6::Int64=1,when2::Int64=1,when3::Int64 = 1,when4::Int64=1,when5::Int64=1,when6::Int64=1,red::Float64 = 0.0,index::Int64 = 0,dosis::Int64=3,ta::Int64 = 999,rc=[0.0],dc=[0],mt::Int64=500,vac::Bool=true,scen::String="statuscuo",alpha::Float64 = 1.0,alpha2::Float64 = 0.0,alpha3::Float64 = 1.0,nsims::Int64=500)
    
    
    #b = bd[h_i]
    #ic = init_con[h_i]
    @everywhere ip = cv.ModelParameters(β=$b,fsevere = 1.0,fmild = 1.0,vaccinating = $vac,
    herd = $(h_i),start_several_inf=true,initialinf3=$ic3,initialinf6=$ic6,initialinf=$ic1,initialinf2=$ic2,initialinf5=$ic5,initialinf4=$ic4,
    time_sec_strain = $when2,time_third_strain = $when3,time_fourth_strain = $when4,time_fifth_strain = $when5,time_sixth_strain = $when6,
    strain_ef_red3 = $red,strain_ef_red4 = $red,
    status_relax = $dosis, relax_after = $ta,file_index = $index,
    modeltime=$mt, prov = Symbol($province), scenario = Symbol($scen), α = $alpha,
    time_change_contact = $dc,
    change_rate_values = $rc,
    α2 = $alpha2,
    α3 = $alpha3)

    folder = create_folder(ip,"pfizer",province)

    #println("$v_e $(ip.vaccine_ef)")
    run(ip,nsims,folder)
   
end