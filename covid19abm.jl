module covid19abm
using Base
using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
include("matrices_code.jl")
@enum HEALTH SUS LAT PRE ASYMP MILD MISO INF IISO HOS ICU REC DED  LAT2 PRE2 ASYMP2 MILD2 MISO2 INF2 IISO2 HOS2 ICU2 REC2 DED2 LAT3 PRE3 ASYMP3 MILD3 MISO3 INF3 IISO3 HOS3 ICU3 REC3 DED3 LAT4 PRE4 ASYMP4 MILD4 MISO4 INF4 IISO4 HOS4 ICU4 REC4 DED4 LAT5 PRE5 ASYMP5 MILD5 MISO5 INF5 IISO5 HOS5 ICU5 REC5 DED5 LAT6 PRE6 ASYMP6 MILD6 MISO6 INF6 IISO6 HOS6 ICU6 REC6 DED6 UNDEF
Base.@kwdef mutable struct Human
    idx::Int64 = 0 
    health::HEALTH = SUS
    health_status::HEALTH = SUS
    swap::HEALTH = UNDEF
    swap_status::HEALTH = UNDEF
    sickfrom::HEALTH = UNDEF
    wentTo::HEALTH = UNDEF
    sickby::Int64 = -1
    nextday_meetcnt::Int16 = 0 ## how many contacts for a single day
    age::Int16   = 0    # in years. don't really need this but left it incase needed later
    ag::Int16   = 0
    tis::Int16   = 0   # time in state 
    exp::Int16   = 0   # max statetime
    dur::NTuple{4, Int8} = (0, 0, 0, 0)   # Order: (latents, asymps, pres, infs) TURN TO NAMED TUPS LATER
    doi::Int16   = 999   # day of infection.
    iso::Bool = false  ## isolated (limited contacts)
    isovia::Symbol = :null ## isolated via quarantine (:qu), preiso (:pi), intervention measure (:im), or contact tracing (:ct)    
    tracing::Bool = false ## are we tracing contacts for this individual?
    tracestart::Int16 = -1 ## when to start tracing, based on values sampled for x.dur
    traceend::Int16 = -1 ## when to end tracing
    tracedby::UInt32 = 0 ## is the individual traced? property represents the index of the infectious person 
    tracedxp::Int16 = 0 ## the trace is killed after tracedxp amount of days
    comorbidity::Int8 = 0 ##does the individual has any comorbidity?
    vac_status::Int8 = 0 ##

    got_inf::Bool = false
    herd_im::Bool = false
    hospicu::Int8 = -1
    ag_new::Int16 = -1
    hcw::Bool = false
    days_vac::Int64 = -1
    first_one::Bool = false
    strain::Int16 = -1
    index_day::Int64 = 1
    relaxed::Bool = false
    recovered::Bool = false
    vaccine::Symbol = :none
    vaccine_n::Int16 = 0
    protected::Int64 = 0
end

## default system parameters
@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    β = 0.0345       
    seasonal::Bool = false ## seasonal betas or not
    popsize::Int64 = 100000
    prov::Symbol = :mississippi
    calibration::Bool = false
    calibration2::Bool = false 
    start_several_inf::Bool = true
    modeltime::Int64 = 332
    initialinf::Int64 = 20
    initialhi::Int64 = 0 ## initial herd immunity, inserts number of REC individuals
    τmild::Int64 = 0 ## days before they self-isolate for mild cases
    fmild::Float64 = 0.0  ## percent of people practice self-isolation
    fsevere::Float64 = 0.0 #
    eldq::Float64 = 0.0 ## complete isolation of elderly
    eldqag::Int8 = 5 ## default age group, if quarantined(isolated) is ag 5. 
    fpreiso::Float64 = 0.0 ## percent that is isolated at the presymptomatic stage
    tpreiso::Int64 = 0## preiso is only turned on at this time. 
    frelasymp::Float64 = 0.26 ## relative transmission of asymptomatic
    ctstrat::Int8 = 0 ## strategy 
    fctcapture::Float16 = 0.0 ## how many symptomatic people identified
    fcontactst::Float16 = 0.0 ## fraction of contacts being isolated/quarantined
    cidtime::Int8 = 0  ## time to identification (for CT) post symptom onset
    cdaysback::Int8 = 0 ## number of days to go back and collect contacts
    #vaccine_ef::Float16 = 0.0   ## change this to Float32 typemax(Float32) typemax(Float64)
    vac_com_dec_max::Float16 = 0.0 # how much the comorbidity decreases the vac eff
    vac_com_dec_min::Float16 = 0.0 # how much the comorbidity decreases the vac eff
    herd::Int8 = 0 #typemax(Int32) ~ millions
    file_index::Int16 = 0
    nstrains::Int16 = 6
    
    #the cap for coverage should be 90% for 65+; 95% for HCW; 80% for 50-64; 60% for 16-49; and then 50% for 12-15 (starting from June 1).

    hcw_vac_comp::Float64 = 0.95
    hcw_prop::Float64 = 0.05 #prop que é trabalhador da saude
    
    eld_comp::Float64 = 0.95
    old_adults::Float64 = 0.95
    young_adults::Float64 = 0.8
    kid_comp::Float64 = 0.8
    #comor_comp::Float64 = 0.7 #prop comorbidade tomam

    vac_period::Array{Int64,1} = [21]
    
    vaccinating::Bool = true #vaccinating?
    drop_rate::Float64 = 0.0 #probability of not getting second dose
    fixed_cov::Float64 = 0.4 #coverage of population
    pfizer_proportion::Float64 = 1.0
    red_risk_perc::Float64 = 1.0 #relative isolation in vaccinated individuals
    reduction_protection::Float64 = 0.0 #reduction in protection against infection
    extra_dose::Bool = false #if it is turned on, extradoses are given to general population
    extra_dose_n::Int64 = 0 #how many
    extra_dose_day::Int64 = 999 #when extra doses are implemented
    days_Rt::Array{Int64,1} = [100;200;300] #days to get Rt

    ##Alpha - B.1.1.7
    sec_strain_trans::Float64 = 1.5#1.5 #transmissibility of second strain
    ins_sec_strain::Bool = true #insert second strain?
    initialinf2::Int64 = 1 #number of initial infected of second strain
    time_sec_strain::Int64 = 178 #when will the second strain introduced -- Jan 3

    ## Gamma - P.1
    ins_third_strain::Bool = true #insert third strain?
    initialinf3::Int64 = 1 #number of initial infected of third strain
    time_third_strain::Int64 = 206 #when will the third strain introduced - P1 March 20
    third_strain_trans::Float64 = 1.6 #transmissibility of third strain
    reduction_recovered::Float64 = 0.21

    ## Delta - B.1.617.2
    ins_fourth_strain::Bool = true #insert fourth strain?
    initialinf4::Int64 = 1 #number of initial infected of fourth strain
    time_fourth_strain::Int64 = 187 #when will the fourth strain introduced
    fourth_strain_trans::Float64 = 1.3 #transmissibility compared to second strain strain

    ## Iota - B.1.526
    ins_fifth_strain::Bool = true #insert fifth strain?
    initialinf5::Int64 = 1 #number of initial infected of fifth strain
    time_fifth_strain::Int64 = 190 #when will the fifth strain introduced
    fifth_strain_trans::Float64 = 1.35 #transmissibility of fifth strain

    ## Beta - B.1.351
    ins_sixth_strain::Bool = true #insert third strain?
    initialinf6::Int64 = 1 #number of initial infected of sixth strain
    time_sixth_strain::Int64 = 212 #when will the sixth strain introduced
    sixth_strain_trans::Float64 = 1.2 #transmissibility of sixth strain

    reduction_recovered_4::Float64 = 0.21
    strain_ef_red3::Float64 = 0.8 #reduction in efficacy against third strain
    strain_ef_red4::Float64 = 0.8 #reduction in efficacy against third strain
    mortality_inc::Float64 = 1.3 #The mortality increase when infected by strain 2

    #=------------ Vaccine Efficacy ----------------------------=#
    days_to_protection::Array{Array{Int64,1},1} = [[14],[0;7]]
    vac_efficacy_inf::Array{Array{Array{Float64,1},1},1} = [[[0.46],[0.6;0.861]],[[0.295],[0.6;0.895]],[[0.46*(1-strain_ef_red3)],[0.6*(1-strain_ef_red3);0.92*(1-strain_ef_red3)]],[[0.46*(1-strain_ef_red4)],[0.6*(1-strain_ef_red4);0.64]],[[0.46],[0.6;0.92]],[[0.46*(1-strain_ef_red3)],[0.6*(1-strain_ef_red3);0.92*(1-strain_ef_red3)]]] #### 50:5:80
    vac_efficacy_symp::Array{Array{Array{Float64,1},1},1} = [[[0.57],[0.66;0.94]],[[0.536],[0.62;0.937]],[[0.332],[0.66;0.94]],[[0.332],[0.62;0.88]],[[0.57],[0.66;0.94]],[[0.332],[0.66;0.94]]] #### 50:5:80
    vac_efficacy_sev::Array{Array{Array{Float64,1},1},1} = [[[0.62],[0.80;0.92]],[[0.541],[0.8;0.94]],[[0.34],[0.68;0.974]],[[0.34],[0.68;0.94]],[[0.62],[0.80;0.92]],[[0.34],[0.68;0.974]]]#### 50:5:80
   

    time_change::Int64 = 999## used to calibrate the model
    how_long::Int64 = 1## used to calibrate the model
    how_much::Float64 = 0.0## used to calibrate the model
    rate_increase::Float64 = how_much/how_long## used to calibrate the model
    time_change_contact::Array{Int64,1} = [1;map(y-> 95+y,0:3);map(y->134+y,0:9);map(y->166+y,0:13);map(y->199+y,0:35)]
    change_rate_values::Array{Float64,1} = [1.0;map(y-> 1.0-0.01*y,1:4);map(y-> 0.96-(0.055/10)*y,1:10);map(y-> 0.90+(0.1/14)*y,1:14);map(y-> 1.0-(0.34/36)*y,1:36)]
    contact_change_rate::Float64 = 1.0 #the rate that receives the value of change_rate_values
    contact_change_2::Float64 = 0.5 ##baseline number that multiplies the contact rate

    relaxed::Bool = false
    relaxing_time::Int64 = 215 ### relax measures for vaccinated
    status_relax::Int16 = 2
    relax_after::Int64 = 1

    day_inital_vac::Int64 = 107 ###this must match to the matrices in matrice code
    day_final_vac::Int64 = 332
    vac_limiar::Float64 = 0.74
    α::Float64 = 1.0
    α2::Float64 = 0.0
    α3::Float64 = 1.0

    scenario::Symbol = :statuscuo
    time_back_to_normal::Int64 = 999 ###relaxing time of measures for non-vaccinated
    ### after calibration, how much do we want to increase the contact rate... in this case, to reach 70%
    ### 0.5*0.95 = 0.475, so we want to multiply this by 1.473684211
    back_normal_rate::Float64 = 1.5151515#2.118644068=>1 # ####1.483050847 =>0.7
end

Base.@kwdef mutable struct ct_data_collect
    total_symp_id::Int64 = 0  # total symptomatic identified
    totaltrace::Int64 = 0     # total contacts traced
    totalisolated::Int64 = 0  # total number of people isolated
    iso_sus::Int64 = 0        # total susceptible isolated 
    iso_lat::Int64 = 0        # total latent isolated
    iso_asymp::Int64 = 0      # total asymp isolated
    iso_symp::Int64 = 0       # total symp (mild, inf) isolated
end

Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

## constants 
const humans = Array{Human}(undef, 0) 
const p = ModelParameters()  ## setup default parameters
const agebraks = @SVector [0:4, 5:19, 20:49, 50:64, 65:99]
const agebraks_vac = @SVector [12:15, 16:17, 18:24, 25:39, 40:49, 50:64, 65:74, 75:99]
const BETAS = Array{Float64, 1}(undef, 0) ## to hold betas (whether fixed or seasonal), array will get resized
const ct_data = ct_data_collect()
export ModelParameters, HEALTH, Human, humans, BETAS

function runsim(simnum, ip::ModelParameters)
    # function runs the `main` function, and collects the data as dataframes. 
    hmatrix,hh1,hh2,hh3,hh4 = main(ip,simnum)            

    ct_numbers = (ct_data.total_symp_id, ct_data.totaltrace, ct_data.totalisolated, 
                    ct_data.iso_sus, ct_data.iso_lat, ct_data.iso_asymp, ct_data.iso_symp)
    ###use here to create the vector of comorbidity
    # get simulation age groups
    #ags = [x.ag for x in humans] # store a vector of the age group distribution 
    ags = [x.ag_new for x in humans] # store a vector of the age group distribution 
    all = _collectdf(hmatrix)
    spl = _splitstate(hmatrix, ags)
    ag1 = _collectdf(spl[1])
    ag2 = _collectdf(spl[2])
    ag3 = _collectdf(spl[3])
    ag4 = _collectdf(spl[4])
    ag5 = _collectdf(spl[5])
    ag6 = _collectdf(spl[6])
    insertcols!(all, 1, :sim => simnum); insertcols!(ag1, 1, :sim => simnum); insertcols!(ag2, 1, :sim => simnum); 
    insertcols!(ag3, 1, :sim => simnum); insertcols!(ag4, 1, :sim => simnum); insertcols!(ag5, 1, :sim => simnum);
    insertcols!(ag6, 1, :sim => simnum);

    
    R01 = zeros(Float64,size(hh1,1))

    for i = 1:size(hh1,1)
        if length(hh1[i]) > 0
            R01[i] = length(findall(k -> k.sickby in hh1[i],humans))/length(hh1[i])
        end
    end

    R02 = zeros(Float64,size(hh2,1))

    for i = 1:size(hh2,1)
        if length(hh2[i]) > 0
            R02[i] = length(findall(k -> k.sickby in hh2[i],humans))/length(hh2[i])
        end
    end

    R03 = zeros(Float64,size(hh3,1))

    for i = 1:size(hh3,1)
        if length(hh3[i]) > 0
            R03[i] = length(findall(k -> k.sickby in hh3[i],humans))/length(hh3[i])
        end
    end

    R04 = zeros(Float64,size(hh4,1))
    for i = 1:size(hh4,1)
        if length(hh4[i]) > 0
            R04[i] = length(findall(k -> k.sickby in hh4[i],humans))/length(hh4[i])
        end
    end

    coverage1 = length(findall(x-> x.age >= 18 && x.vac_status >= 1,humans))/length(findall(x-> x.age >= 18,humans))
    coverage2 = length(findall(x-> x.age >= 18 && x.vac_status == 2,humans))/length(findall(x-> x.age >= 18,humans))

    coverage12 = length(findall(x-> x.vac_status >= 1,humans))/p.popsize
    coverage22 = length(findall(x-> x.vac_status == 2,humans))/p.popsize

    return (a=all, g1=ag1, g2=ag2, g3=ag3, g4=ag4, g5=ag5,g6=ag6,   
    iniiso = ct_data.totalisolated,
    R01 = R01,
    R02 = R02, cov1 = coverage1,cov2 = coverage2,cov12 = coverage12,cov22 = coverage22)
end
export runsim

function main(ip::ModelParameters,sim::Int64)
    Random.seed!(sim*726)
    ## datacollection            
    # matrix to collect model state for every time step

    # reset the parameters for the simulation scenario
    reset_params(ip)  #logic: outside "ip" parameters are copied to internal "p" which is a global const and available everywhere. 

    p.popsize == 0 && error("no population size given")
    
    hmatrix = zeros(Int16, p.popsize, p.modeltime)
    initialize() # initialize population
    
    vac_rate_1::Matrix{Int64} = vaccination_rate_1(sim)
    vac_rate_2::Matrix{Int64} = vaccination_rate_2(sim)
    vaccination_days::Vector{Int64} = days_vac_f(size(vac_rate_1,1))

    #h_init::Int64 = 0
    # insert initial infected agents into the model
    # and setup the right swap function. 
   
    N = herd_immu_dist_4(sim,1)
    if p.initialinf > 0
        insert_infected(PRE, p.initialinf, 4, 1)[1]
    end
        #findall(x->x.health in (MILD,INF,LAT,PRE,ASYMP),humans)
   
    h_init1 = findall(x->x.health  in (LAT,MILD,MISO,INF,PRE,ASYMP),humans)
    h_init1 = [h_init1]
    h_init2 = []
    h_init3 = []
    h_init4 = []
    ## save the preisolation isolation parameters
    _fpreiso = p.fpreiso
    p.fpreiso = 0

    # split population in agegroups 
    grps = get_ag_dist()
    count_change::Int64 = 1
    
    time_vac::Int64 = 1
    time_pos::Int64 = 0
    if p.vaccinating
        vac_ind = vac_selection(sim)
    else
        time_vac = 9999 #this guarantees that no one will be vaccinated
    end
    # start the time loop
    for st = 1:p.modeltime
        if p.ins_sec_strain && st == p.time_sec_strain ##insert second strain
            insert_infected(PRE, p.initialinf2, 4, 2)[1]
            h_init2 = findall(x->x.health  in (LAT2,MILD2,INF2,PRE2,ASYMP2),humans)
            h_init2 = [h_init2]
        end
        if p.ins_third_strain && st == p.time_third_strain #insert third strain
            insert_infected(PRE, p.initialinf3, 4, 3)[1]
            h_init3 = findall(x->x.health  in (LAT3,MILD3,INF3,PRE3,ASYMP3),humans)
            h_init3 = [h_init3]
        end
        if p.ins_fourth_strain && st == p.time_fourth_strain #insert third strain
            insert_infected(PRE, p.initialinf4, 4, 4)[1]
            h_init4 = findall(x->x.health  in (LAT4,MILD4,INF4,PRE4,ASYMP4),humans)
            h_init4 = [h_init4]
        end
        if p.ins_fifth_strain && st == p.time_fifth_strain #insert third strain
            insert_infected(PRE, p.initialinf5, 4, 5)[1]
        end
        if p.ins_sixth_strain && st == p.time_sixth_strain #insert third strain
            insert_infected(PRE, p.initialinf6, 4, 6)[1]
        end
        if length(p.time_change_contact) >= count_change && p.time_change_contact[count_change] == st ###change contact pattern throughout the time
            setfield!(p, :contact_change_rate, p.change_rate_values[count_change])
            count_change += 1
        end
        # start of day
        #println("$st")

        if st == p.relaxing_time ### time that people vaccinated people is allowed to go back to normal
            setfield!(p, :relaxed, true)
        elseif st == p.time_back_to_normal ##time that non-vaccinated people is allowed to go back to normal
            setfield!(p, :contact_change_2, p.contact_change_2*p.back_normal_rate)
            #setfield!(p, :contact_change_rate, 1.0)
        end

        if time_pos < length(vaccination_days) && time_vac == vaccination_days[time_pos+1]
            time_pos += 1
        end
        time_vac += 1
        time_pos > 0 && vac_time!(sim,vac_ind,time_pos+1,vac_rate_1,vac_rate_2)
        
        println([time_vac length(findall(x-> x.vac_status == 2 && x.age >= 18,humans))])
       
        _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
        dyntrans(st, grps,sim)
        if st in p.days_Rt ### saves individuals that became latent on days_Rt
            aux1 = findall(x->x.swap == LAT,humans)
            h_init1 = vcat(h_init1,[aux1])
            aux2 = findall(x->x.swap == LAT2,humans)
            h_init2 = vcat(h_init2,[aux2])
            aux3 = findall(x->x.swap == LAT3,humans)
            h_init3 = vcat(h_init3,[aux3])
            aux4 = findall(x->x.swap == LAT4,humans)
            h_init4 = vcat(h_init4,[aux4])
        end
        sw = time_update() ###update the system
        # end of day
    end
    
    
    return hmatrix,h_init1,h_init2,h_init3, h_init4 ## return the model state as well as the age groups. 
end
export main

function vac_selection(sim::Int64)
    
    

   
    aux_1 = map(k-> findall(y-> y.age in k && y.age >= 12 && y.comorbidity == 1,humans),agebraks_vac)
    aux_2 = map(k-> findall(y-> y.age in k && y.age >= 12 && y.comorbidity == 0,humans),agebraks_vac)

    v = map(x-> [aux_1[x];aux_2[x]],1:length(aux_1))
    
    return v
end


function vac_time!(sim::Int64,vac_ind::Vector{Vector{Int64}},time_pos::Int64,vac_rate_1::Matrix{Int64},vac_rate_2::Matrix{Int64})
    aux_states = (MILD, MISO, INF, IISO, HOS, ICU, DED)
    ##first dose
    rng = MersenneTwister(123*sim)
    ### lets create distribute the number of doses per age group
    
    remaining_doses::Int64 = 0
    total_given::Int64 = 0
    for i in 1:length(vac_ind)
        pos = findall(y-> humans[y].vac_status == 1 && humans[y].days_vac >= p.vac_period[humans[y].vaccine_n] && !(humans[y].health_status in aux_states),vac_ind[i])
        
        l1 = min(vac_rate_2[time_pos,i],length(pos))
        remaining_doses += (vac_rate_2[time_pos,i] - l1)
        for j = 1:l1
            x = humans[vac_ind[i][pos[j]]]
            x.days_vac = 0
            x.vac_status = 2
            x.index_day = 1
            total_given += 1
        end

        pos = findall(y-> humans[y].vac_status == 0 && !(humans[y].health_status in aux_states),vac_ind[i])
        
        l2 = min(vac_rate_1[time_pos,i],length(pos))

        remaining_doses += (vac_rate_1[time_pos,i] - l2)

        for j = 1:l2
            x = humans[vac_ind[i][pos[j]]]
            x.days_vac = 0
            x.vac_status = 1
            x.index_day = 1
    
            x.vaccine = rand(rng) <= p.pfizer_proportion ? :pfizer : :moderna
            x.vaccine_n = x.vaccine == :pfizer ? 1 : 2
            total_given += 1
        end

    end

    ###remaining_doses are given to any individual within the groups that are vaccinated on that day
    if remaining_doses > 0

        pos = map(k->findall(y-> humans[y].vac_status == 1 && humans[y].days_vac >= p.vac_period[humans[y].vaccine_n] && !(humans[y].health_status in aux_states),vac_ind[k]),1:length(vac_ind))
        pos2 = map(k->findall(y-> humans[y].vac_status == 0 && !(humans[y].health_status in aux_states),vac_ind[k]),1:length(vac_ind))
        
        aux = findall(x-> vac_rate_1[time_pos,x] > 0 || vac_rate_2[time_pos,x] > 0, 1:length(vac_ind))
        position = map(k-> vac_ind[k][pos[k]],aux)
        position2 = map(k-> vac_ind[k][pos2[k]],aux)

        r = vcat(position...,position2...)
        m = min(remaining_doses,length(r))

        rr = sample(rng,r,m,replace=false)
        
        for i in rr
            x = humans[i]
            if x.vac_status == 0
                x.days_vac = 0
                x.vac_status = 1
                x.index_day = 1
                x.vaccine = rand(rng) <= p.pfizer_proportion ? :pfizer : :moderna
                x.vaccine_n = x.vaccine == :pfizer ? 1 : 2
                total_given += 1
            elseif x.vac_status == 1
                x.days_vac = 0
                x.vac_status = 2
                x.index_day = 1
                total_given += 1
            else
                error("error in humans vac status - vac time")
            end
        end
    end

    t = sum(vac_rate_1[time_pos,:]+vac_rate_2[time_pos,:])
    println("Total $time_pos $remaining_doses $total_given $t")
    if total_given > t
        error("vaccination")
    end

end

function vac_update(x::Human)
    
    if x.vac_status == 1
        #x.index_day == 2 && error("saiu com indice 2")
        if x.days_vac == p.days_to_protection[x.vac_status][1]#14
            x.protected = 1
            x.index_day = min(length(p.days_to_protection[x.vac_status]),x.index_day+1)
        elseif x.days_vac == p.days_to_protection[x.vac_status][x.index_day]#14
            x.protected = x.index_day
            x.index_day = min(length(p.days_to_protection[x.vac_status]),x.index_day+1)
        end
        if !x.relaxed
            x.relaxed = p.relaxed &&  x.vac_status >= p.status_relax && x.days_vac >= p.relax_after ? true : false
        end
        x.days_vac += 1

    elseif x.vac_status == 2
        if x.days_vac == p.days_to_protection[x.vac_status][1]#0
            x.protected = 1
            x.index_day = min(length(p.days_to_protection[x.vac_status]),x.index_day+1)

        elseif x.days_vac == p.days_to_protection[x.vac_status][x.index_day]#7
            x.protected = x.index_day
            x.index_day = min(length(p.days_to_protection[x.vac_status]),x.index_day+1)
        end
        if !x.relaxed
            x.relaxed = p.relaxed &&  x.vac_status >= p.status_relax && x.days_vac >= p.relax_after ? true : false
        end
        x.days_vac += 1
    end
   
end
function reset_params(ip::ModelParameters)
    # the p is a global const
    # the ip is an incoming different instance of parameters 
    # copy the values from ip to p. 
    for x in propertynames(p)
        setfield!(p, x, getfield(ip, x))
    end

    # reset the contact tracing data collection structure
    for x in propertynames(ct_data)
        setfield!(ct_data, x, 0)
    end

    # resize and update the BETAS constant array
    #init_betas()

    # resize the human array to change population size
    resize!(humans, p.popsize)
end
export reset_params, reset_params_default

function _model_check() 
    ## checks model parameters before running 
    (p.fctcapture > 0 && p.fpreiso > 0) && error("Can not do contact tracing and ID/ISO of pre at the same time.")
    (p.fctcapture > 0 && p.maxtracedays == 0) && error("maxtracedays can not be zero")
end

## Data Collection/ Model State functions
function _get_model_state(st, hmatrix)
    # collects the model state (i.e. agent status at time st)
    for i=1:length(humans)
        hmatrix[i, st] = Int(humans[i].health)
    end    
end
export _get_model_state

function _collectdf(hmatrix)
    ## takes the output of the humans x time matrix and processes it into a dataframe
    #_names_inci = Symbol.(["lat_inc", "mild_inc", "miso_inc", "inf_inc", "iiso_inc", "hos_inc", "icu_inc", "rec_inc", "ded_inc"])    
    #_names_prev = Symbol.(["sus", "lat", "mild", "miso", "inf", "iiso", "hos", "icu", "rec", "ded"])
    mdf_inc, mdf_prev = _get_incidence_and_prev(hmatrix)
    mdf = hcat(mdf_inc, mdf_prev)    
    _names_inc = Symbol.(string.((Symbol.(instances(HEALTH)[1:end - 1])), "_INC"))
    _names_prev = Symbol.(string.((Symbol.(instances(HEALTH)[1:end - 1])), "_PREV"))
    _names = vcat(_names_inc..., _names_prev...)
    datf = DataFrame(mdf, _names)
    insertcols!(datf, 1, :time => 1:p.modeltime) ## add a time column to the resulting dataframe
    return datf
end

function _splitstate(hmatrix, ags)
    #split the full hmatrix into 4 age groups based on ags (the array of age group of each agent)
    #sizes = [length(findall(x -> x == i, ags)) for i = 1:4]
    matx = []#Array{Array{Int64, 2}, 1}(undef, 4)
    for i = 1:maximum(ags)#length(agebraks)
        idx = findall(x -> x == i, ags)
        push!(matx, view(hmatrix, idx, :))
    end
    return matx
end
export _splitstate

function _get_incidence_and_prev(hmatrix)
    cols = instances(HEALTH)[1:end - 1] ## don't care about the UNDEF health status
    inc = zeros(Int64, p.modeltime, length(cols))
    pre = zeros(Int64, p.modeltime, length(cols))
    for i = 1:length(cols)
        inc[:, i] = _get_column_incidence(hmatrix, cols[i])
        pre[:, i] = _get_column_prevalence(hmatrix, cols[i])
    end
    return inc, pre
end

function _get_column_incidence(hmatrix, hcol)
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for r in eachrow(hmatrix)
        idx = findfirst(x -> x == inth, r)
        if idx !== nothing 
            timevec[idx] += 1
        end
    end
    return timevec
end

function herd_immu_dist_4(sim::Int64,strain::Int64)
    rng = MersenneTwister(200*sim)
    vec_n = zeros(Int32,6)
    N::Int64 = 0
    if p.herd == 5
        vec_n = [9; 148; 262;  68; 4; 9]
        N = 5

    elseif p.herd == 10
        vec_n = [32; 279; 489; 143; 24; 33]

        N = 9

    elseif p.herd == 20
        vec_n = [71; 531; 962; 302; 57; 77]

        N = 14
    elseif p.herd == 30
        vec_n = [105; 757; 1448; 481; 87; 122]

        N = 16
    elseif p.herd == 50
        vec_n = map(y->y*5,[32; 279; 489; 143; 24; 33])

        N = 16
    elseif p.herd == 0
        vec_n = [0;0;0;0;0;0]
       
    else
        vec_n = map(y->Int(round(y*p.herd/10)),[32; 279; 489; 143; 24; 33])
        N = 16
    end

    for g = 1:6
        pos = findall(y->y.ag_new == g && y.health == SUS,humans)
        n_dist = min(length(pos),Int(floor(vec_n[g]*p.popsize/10000)))
        pos2 = sample(rng,pos,n_dist,replace=false)
        for i = pos2
            humans[i].strain = strain
            humans[i].swap = strain == 1 ? REC : REC2
            move_to_recovered(humans[i])
            humans[i].sickfrom = INF
            humans[i].herd_im = true
        end
    end
    return N
end

function _get_column_prevalence(hmatrix, hcol)
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for (i, c) in enumerate(eachcol(hmatrix))
        idx = findall(x -> x == inth, c)
        if idx !== nothing
            ps = length(c[idx])    
            timevec[i] = ps    
        end
    end
    return timevec
end

export _collectdf, _get_incidence_and_prev, _get_column_incidence, _get_column_prevalence

## initialization functions 
function get_province_ag(prov) 
    ret = @match prov begin
        :wisconsin => Distributions.Categorical(@SVector [0.056762515470334,0.187481558399803,0.375762610619545,0.205282361294263,0.174710954216055])
        :newhampshire => Distributions.Categorical(@SVector [0.046790089952939,0.167253923811751,0.370258827059574,0.228992778612514,0.186704380563223])
        :colorado => Distributions.Categorical(@SVector [0.05768644369181,0.186806618674654,0.42620898058185,0.183013772466736,0.146284184584951])
        :california => Distributions.Categorical(@SVector [0.060328572249656,0.190869569651902,0.41753535355376,0.183511846448123,0.147754658096559])
        :michigan => Distributions.Categorical(@SVector [0.056718745447141,0.184367113697533,0.37770471730996,0.204436991537978,0.176772432007387])
        :maine => Distributions.Categorical(@SVector [0.047267097749462,0.161894849919507,0.353223301086436,0.225397481944812,0.212217269299783])
        :connecticut => Distributions.Categorical(@SVector [0.05096644393565,0.181716647215217,0.376013207351891,0.214531396771144,0.176772304726099])
        :oregon => Distributions.Categorical(@SVector [0.054012613873269,0.174896870051404,0.402817435036846,0.186640134271056,0.181632946767425])
        :minnesota => Distributions.Categorical(@SVector [0.06234839436332,0.193935348973124,0.385302622582466,0.195214687766861,0.163198946314228])
        :virginia => Distributions.Categorical(@SVector [0.059220417645371,0.185337177504965,0.400913172356596,0.195323213503479,0.159206018989589])
        :pennsylvania => Distributions.Categorical(@SVector [0.054516841093989,0.178067486232022,0.37460593037535,0.205857386692021,0.186952355606617])
        :illinois => Distributions.Categorical(@SVector [0.058944487931135,0.189268377449461,0.396559894588157,0.193985063393809,0.161242176637438])
        :ohio => Distributions.Categorical(@SVector [0.059100187354031,0.187871179132696,0.377904372449547,0.200062023594631,0.175062237469095])
        :arizona => Distributions.Categorical(@SVector [0.059047219448153,0.19355196801854,0.389072002662008,0.178539844315969,0.179788965555331])
        :kentucky => Distributions.Categorical(@SVector [0.061018342210811,0.189432843451166,0.384296702108682,0.197253245705315,0.167998866524027])
        :tennessee => Distributions.Categorical(@SVector [0.059832272541306,0.185814711998845,0.392468254579544,0.194457045610494,0.167427715269812])
        :southcarolina => Distributions.Categorical(@SVector [0.056803310496563,0.18555604370334,0.379181286822302,0.196467700478217,0.181991658499579])
        :northcarolina => Distributions.Categorical(@SVector [0.05813931314814,0.189039961922502,0.391583438881687,0.194276952778029,0.166960333269642])
        :nevada => Distributions.Categorical(@SVector [0.060248571825583,0.186580484884532,0.404644764745682,0.187504464059613,0.161021714484591])
        :westvirginia => Distributions.Categorical(@SVector [0.05190701432416,0.172668871470923,0.364284291411363,0.206351376310091,0.204788446483464])
        :oklahoma => Distributions.Categorical(@SVector [0.064577930947687,0.202923397720125,0.390831269675719,0.181157759306298,0.160509642350171])
        :maryland => Distributions.Categorical(@SVector [0.059867045559805,0.186543780021437,0.391742698918897,0.203155310899684,0.158691164600177])
        :massachusetts => Distributions.Categorical(@SVector [0.051847928103912,0.174228288330088,0.400024925632967,0.204246120748877,0.169652737184155])
        :newyork => Distributions.Categorical(@SVector [0.057932889510563,0.174466515410726,0.399115102885276,0.199048852803865,0.16943663938957])
        :texas => Distributions.Categorical(@SVector [0.068661166046309,0.214502673672857,0.416443080311993,0.171608270843711,0.128784809125131])
        :alabama => Distributions.Categorical(@SVector [0.06003383515001,0.188057558505339,0.381706584597563,0.196878559548538,0.173323462198551])
        :louisiana => Distributions.Categorical(@SVector [0.064848861876865,0.193984934587336,0.391711484742064,0.190053807503624,0.159400911290111])
        :vermont => Distributions.Categorical(@SVector [0.04654408971953,0.1688683614615,0.366096197208605,0.218104806334727,0.200386545275638])
        :missouri => Distributions.Categorical(@SVector [0.059973004978633,0.188875698419599,0.382552593692341,0.195556021186725,0.173042681722702])
        :georgia => Distributions.Categorical(@SVector [0.061838545944718,0.202038008658033,0.40539884301492,0.18785057353371,0.14287402884862])
        :florida => Distributions.Categorical(@SVector [0.053066205252444,0.166443652792657,0.372408508401048,0.198686342048047,0.209395291505804])
        :mississippi => Distributions.Categorical(@SVector [0.061649467146974,0.200597819531213,0.384003287469814,0.190218298882213,0.163531126969785])
        :arkansas => Distributions.Categorical(@SVector [0.062450709191187,0.195660155530313,0.380966424592187,0.187325618231005,0.173597092455309])
        :usa => Distributions.Categorical(@SVector [0.059444636404977,0.188450296592341,0.396101793107413,0.189694011721906,0.166309262173363])
        :newyorkcity   => Distributions.Categorical(@SVector [0.064000, 0.163000, 0.448000, 0.181000, 0.144000])
        _ => error("shame for not knowing your canadian provinces and territories")
    end       
    return ret  
end
export get_province_ag

function comorbidity(ag::Int16)

    a = [4;19;49;64;79;999]
    g = findfirst(x->x>=ag,a)
    prob = [0.05; 0.1; 0.28; 0.55; 0.74; 0.81]

    com = rand() < prob[g] ? 1 : 0

    return com    
end
export comorbidity


function initialize() 
    agedist = get_province_ag(p.prov)
    for i = 1:p.popsize 
        humans[i] = Human()              ## create an empty human       
        x = humans[i]
        x.idx = i 
        x.ag = rand(agedist)
        x.age = rand(agebraks[x.ag]) 
        a = [4;19;49;64;79;999]
        g = findfirst(y->y>=x.age,a)
        x.ag_new = g
        x.exp = 999  ## susceptible people don't expire.
        x.dur = sample_epi_durations() # sample epi periods   
        if rand() < p.eldq && x.ag == p.eldqag   ## check if elderly need to be quarantined.
            x.iso = true   
            x.isovia = :qu         
        end
        x.comorbidity = comorbidity(x.age)
        # initialize the next day counts (this is important in initialization since dyntrans runs first)
        get_nextday_counts(x)
        
    end
end
export initialize

function init_betas() 
    if p.seasonal  
        tmp = p.β .* td_seasonality()
    else 
        tmp = p.β .* ones(Float64, p.modeltime)
    end
    resize!(BETAS, length(tmp))
    for i = 1:length(tmp)
        BETAS[i] = tmp[i]
    end
end

function td_seasonality()
    ## returns a vector of seasonal oscillations
    t = 1:p.modeltime
    a0 = 6.261
    a1 = -11.81
    b1 = 1.817
    w = 0.022 #0.01815    
    temp = @. a0 + a1*cos((80-t)*w) + b1*sin((80-t)*w)  #100
    #temp = @. a0 + a1*cos((80-t+150)*w) + b1*sin((80-t+150)*w)  #100
    temp = (temp .- 2.5*minimum(temp))./(maximum(temp) .- minimum(temp)); # normalize  @2
    return temp
end

function get_ag_dist() 
    # splits the initialized human pop into its age groups
    grps =  map(x -> findall(y -> y.ag == x, humans), 1:length(agebraks)) 
    return grps
end

function insert_infected(health, num, ag,strain) 
    ## inserts a number of infected people in the population randomly
    ## this function should resemble move_to_inf()
    l = findall(x -> x.health == SUS && x.ag == ag, humans)
    aux_pre = [PRE;PRE2;PRE3;PRE4;PRE5;PRE6]
    aux_lat = [LAT;LAT2;LAT3;LAT4;LAT5;LAT6]
    aux_mild = [MILD;MILD2;MILD3;MILD4;MILD5;MILD6]
    aux_inf = [INF;INF2;INF3;INF4;INF5;INF6]
    aux_asymp = [ASYMP;ASYMP2;ASYMP3;ASYMP4;ASYMP5;ASYMP6]
    aux_rec = [REC;REC2;REC3;REC4;REC5;REC6]
    if length(l) > 0 && num < length(l)
        h = sample(l, num; replace = false)
        @inbounds for i in h 
            x = humans[i]
            x.strain = strain
            x.first_one = true

            if x.strain > 0
                if health == PRE
                    x.swap = aux_pre[x.strain]
                    x.swap_status = PRE
                    move_to_pre(x) ## the swap may be asymp, mild, or severe, but we can force severe in the time_update function
                elseif health == LAT
                    x.swap = aux_lat[x.strain]
                    x.swap_status = LAT
                    move_to_latent(x)
                elseif health == MILD
                    x.swap =  aux_mild[x.strain] 
                    x.swap_status = MILD
                    move_to_mild(x)
                elseif health == INF
                    x.swap = aux_inf[x.strain]
                    x.swap_status = INF
                    move_to_infsimple(x)
                elseif health == ASYMP
                    x.swap = aux_asymp[x.strain]
                    x.swap_status = ASYMP
                    move_to_asymp(x)
                elseif health == REC 
                    x.swap_status = REC
                    x.swap = aux_rec[x.strain]
                    move_to_recovered(x)
                else 
                    error("can not insert human of health $(health)")
                end
            else
                error("no strain insert inf")
            end
            
            x.sickfrom = INF # this will add +1 to the INF count in _count_infectors()... keeps the logic simple in that function.    
            
        end
    end    
    return h
end
export insert_infected

function time_update()
    # counters to calculate incidence

    lat_v = zeros(Int64,p.nstrains)
    pre_v = zeros(Int64,p.nstrains)
    asymp_v = zeros(Int64,p.nstrains)
    mild_v = zeros(Int64,p.nstrains)
    miso_v = zeros(Int64,p.nstrains)
    inf_v = zeros(Int64,p.nstrains)
    infiso_v = zeros(Int64,p.nstrains)
    hos_v = zeros(Int64,p.nstrains)
    icu_v = zeros(Int64,p.nstrains)
    rec_v = zeros(Int64,p.nstrains)
    ded_v = zeros(Int64,p.nstrains)
    
    for x in humans 
        x.tis += 1 
        x.doi += 1 # increase day of infection. variable is garbage until person is latent
        if x.tis >= x.exp             
            @match Symbol(x.swap_status) begin
                :LAT  => begin move_to_latent(x); lat_v[x.strain] += 1; end
                :PRE  => begin move_to_pre(x); pre_v[x.strain] += 1; end
                :ASYMP => begin move_to_asymp(x); asymp_v[x.strain] += 1; end
                :MILD => begin move_to_mild(x); mild_v[x.strain] += 1; end
                :MISO => begin move_to_miso(x); miso_v[x.strain] += 1; end
                :INF  => begin move_to_inf(x); inf_v[x.strain] +=1; end    
                :IISO => begin move_to_iiso(x); infiso_v[x.strain] += 1; end
                :HOS  => begin move_to_hospicu(x); hos_v[x.strain] += 1; end 
                :ICU  => begin move_to_hospicu(x); icu_v[x.strain] += 1; end
                :REC  => begin move_to_recovered(x); rec_v[x.strain] += 1; end
                :DED  => begin move_to_dead(x); ded_v[x.strain] += 1; end
                _    => begin dump(x); error("swap expired, but no swap set."); end
            end
        end
        # run covid-19 functions for other integrated dynamics. 
        #ct_dynamics(x)
        # get the meet counts for the next day 
        get_nextday_counts(x)
        if p.vaccinating
            vac_update(x)
        end
    end


    (lat,lat2,lat3,lat4,lat5,lat6) = lat_v
    (mild,mild2,mild3,mild4,mild5,mild6) = mild_v
    (miso,miso2,miso3,miso4,miso5,miso6) = miso_v
    (inf,inf2,inf3,inf4,inf5,inf6) = inf_v
    (infiso,infiso2,infiso3,infiso4,infiso5,infiso6) = infiso_v
    (hos,hos2,hos3,hos4,hos5,hos6) = hos_v
    (icu,icu2,icu3,icu4,icu5,icu6) = icu_v
    (rec,rec2,rec3,rec4,rec5,rec6) = rec_v
    (ded,ded2,ded3,ded4,ded5,ded6) = ded_v

    return (lat, mild, miso, inf, infiso, hos, icu, rec, ded,lat2, mild2, miso2, inf2, infiso2, hos2, icu2, rec2, ded2,lat3, mild3, miso3, inf3, infiso3, hos3, icu3, rec3, ded3, lat4, mild4, miso4, inf4, infiso4, hos4, icu4, rec4, ded4, lat5, mild5, miso5, inf5, infiso5, hos5, icu5, rec5, ded5, lat6, mild6, miso6, inf6, infiso6, hos6, icu6, rec6, ded6)
end
export time_update

@inline _set_isolation(x::Human, iso) = _set_isolation(x, iso, x.isovia)
@inline function _set_isolation(x::Human, iso, via)
    # a helper setter function to not overwrite the isovia property. 
    # a person could be isolated in susceptible/latent phase through contact tracing
    # --> in which case it will follow through the natural history of disease 
    # --> if the person remains susceptible, then iso = off
    # a person could be isolated in presymptomatic phase through fpreiso
    # --> if x.iso == true from CT and x.isovia == :ct, do not overwrite
    # a person could be isolated in mild/severe phase through fmild, fsevere
    # --> if x.iso == true from CT and x.isovia == :ct, do not overwrite
    # --> if x.iso == true from PRE and x.isovia == :pi, do not overwrite
    x.iso = iso 
    x.isovia == :null && (x.isovia = via)
end

function sample_epi_durations()
    # when a person is sick, samples the 
    lat_dist = Distributions.truncated(LogNormal(1.434, 0.661), 4, 7) # truncated between 4 and 7
    pre_dist = Distributions.truncated(Gamma(1.058, 5/2.3), 0.8, 3)#truncated between 0.8 and 3
    asy_dist = Gamma(5, 1)
    inf_dist = Gamma((3.2)^2/3.7, 3.7/3.2)

    latents = Int.(round.(rand(lat_dist)))
    pres = Int.(round.(rand(pre_dist)))
    latents = latents - pres # ofcourse substract from latents, the presymp periods
    asymps = Int.(ceil.(rand(asy_dist)))
    infs = Int.(ceil.(rand(inf_dist)))
    return (latents, asymps, pres, infs)
end

function move_to_latent(x::Human)
    ## transfers human h to the incubation period and samples the duration
    x.health = x.swap
    x.health_status = x.swap_status
    x.doi = 0 ## day of infection is reset when person becomes latent
    x.tis = 0   # reset time in state 
    x.exp = x.dur[1] # get the latent period
   
    #0-18 31 19 - 59 29 60+ 18 going to asymp
    symp_pcts = [0.7, 0.623, 0.672, 0.672, 0.812, 0.812] #[0.3 0.377 0.328 0.328 0.188 0.188]
    age_thres = [4, 19, 49, 64, 79, 999]
    g = findfirst(y-> y >= x.age, age_thres)

    if x.recovered
        auxiliar = (1-p.vac_efficacy_symp[1][2][end])
    else
        aux = x.vac_status*x.protected > 0 ? p.vac_efficacy_symp[x.strain][x.vac_status][x.protected] : 0.0
        auxiliar = (1-aux)
    end
 
    if rand() < (symp_pcts[g])*auxiliar

        aux_v = [PRE;PRE2;PRE3;PRE4;PRE5;PRE6]
        x.swap = aux_v[x.strain]
        x.swap_status = PRE
        
    else
        aux_v = [ASYMP;ASYMP2;ASYMP3;ASYMP4;ASYMP5;ASYMP6]
        x.swap = aux_v[x.strain]
        x.swap_status = ASYMP
        
    end
    x.wentTo = x.swap
    x.got_inf = true
    ## in calibration mode, latent people never become infectious.
    if p.calibration && !x.first_one
        x.swap = LAT 
        x.exp = 999
    end 
end
export move_to_latent

function move_to_asymp(x::Human)
    ## transfers human h to the asymptomatic stage 
    x.health = x.swap  
    x.health_status = x.swap_status
    x.tis = 0 
    x.exp = x.dur[2] # get the presymptomatic period
   
    aux_v = [REC;REC2;REC3;REC4;REC5;REC6]
    x.swap = aux_v[x.strain]
    x.swap_status = REC
    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, the asymptomatic individual has limited contacts
end
export move_to_asymp

function move_to_pre(x::Human)
    if x.strain == 1 || x.strain == 3 || x.strain == 5 || x.strain == 6
        θ = (0.95, 0.9, 0.85, 0.6, 0.2)  # percentage of sick individuals going to mild infection stage
    elseif x.strain == 2 || x.strain == 4
        θ = (0.89, 0.78, 0.67, 0.48, 0.04)
            if x.strain == 4
                θ = map(y-> max(0,1-(1-y)*1.88),θ)
            end
    else
        error("no strain in move to pre")
    end  # percentage of sick individuals going to mild infection stage
    x.health = x.swap
    x.health_status = x.swap_status
    x.tis = 0   # reset time in state 
    x.exp = x.dur[3] # get the presymptomatic period


    if x.recovered
        auxiliar = (1-p.vac_efficacy_sev[1][2][end])
    else
        aux = x.vac_status*x.protected > 0 ? p.vac_efficacy_sev[x.strain][x.vac_status][x.protected] : 0.0
        auxiliar = (1-aux)
    end

    if rand() < (1-θ[x.ag])*auxiliar
        aux_v = [INF;INF2;INF3;INF4;INF5;INF6]
        x.swap = aux_v[x.strain]
        x.swap_status = INF
    else 
        aux_v = [MILD;MILD2;MILD3;MILD4;MILD5;MILD6]
        x.swap = aux_v[x.strain]
        x.swap_status = MILD
    end
    # calculate whether person is isolated
    rand() < p.fpreiso && _set_isolation(x, true, :pi)
end
export move_to_pre

function move_to_mild(x::Human)
    ## transfers human h to the mild infection stage for γ days
   
    x.health = x.swap 
    x.health_status = x.swap_status
    x.tis = 0 
    x.exp = x.dur[4]
    aux_v = [REC;REC2;REC3;REC4;REC5;REC6]
    x.swap = aux_v[x.strain]
    x.swap_status = REC
    
    #x.swap = x.strain == 1 ? REC : REC2
    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, staying in MILD is same as MISO since contacts will be limited. 
    # we still need the separation of MILD, MISO because if x.iso is false, then here we have to determine 
    # how many days as full contacts before self-isolation
    # NOTE: if need to count non-isolated mild people, this is overestimate as isolated people should really be in MISO all the time
    #   and not go through the mild compartment 
    aux = x.vac_status > 0 ? p.fmild*p.red_risk_perc : p.fmild

    if x.iso || rand() < aux#p.fmild
        aux_v = [MISO;MISO2;MISO3;MISO4;MISO5;MISO6]
        x.swap = aux_v[x.strain]
        x.swap_status = MISO
        #x.swap = x.strain == 1 ? MISO : MISO2  
        x.exp = p.τmild
    end
end
export move_to_mild

function move_to_miso(x::Human)
    ## transfers human h to the mild isolated infection stage for γ days
    x.health = x.swap
    x.health_status = x.swap_status
    aux_v = [REC;REC2;REC3;REC4;REC5;REC6]
    x.swap = aux_v[x.strain]
    x.swap_status = REC
    #x.swap = x.strain == 1 ? REC : REC2
    x.tis = 0 
    x.exp = x.dur[4] - p.τmild  ## since tau amount of days was already spent as infectious
    _set_isolation(x, true, :mi) 
end
export move_to_miso

function move_to_infsimple(x::Human)
    ## transfers human h to the severe infection stage for γ days 
    ## simplified function for calibration/general purposes
    x.health = x.swap
    x.health_status = x.swap_status
    x.tis = 0 
    x.exp = x.dur[4]
    aux_v = [REC;REC2;REC3;REC4;REC5;REC6]
    x.swap = aux_v[x.strain]
    x.swap_status = REC
    #x.swap = x.strain == 1 ? REC : REC2
    _set_isolation(x, false, :null) 
end

function move_to_inf(x::Human)
    ## transfers human h to the severe infection stage for γ days
    ## for swap, check if person will be hospitalized, selfiso, die, or recover
 
    # h = prob of hospital, c = prob of icu AFTER hospital    
    comh = 0.98
    if x.strain == 1 || x.strain == 3 || x.strain == 5 || x.strain == 6
        h = x.comorbidity == 1 ? comh : 0.04 #0.376
        c = x.comorbidity == 1 ? 0.396 : 0.25

    elseif x.strain == 2 || x.strain == 4
        if x.age <  20
            h = x.comorbidity == 1 ? comh : 0.05*1.07*1 #0.376
            c = x.comorbidity == 1 ? 0.396*1.07 : 0.25*1.07
        elseif x.age >= 20 && x.age < 30
            h = x.comorbidity == 1 ? comh : 0.05*1.29*1 #0.376
            c = x.comorbidity == 1 ? 0.396*1.29 : 0.25*1.29
        elseif  x.age >= 30 && x.age < 40
            h = x.comorbidity == 1 ? comh : 0.05*1.45*1 #0.376
            c = x.comorbidity == 1 ? 0.396*1.45 : 0.25*1.45
        elseif  x.age >= 40 && x.age < 50
            h = x.comorbidity == 1 ? comh : 0.05*1.61*1 #0.376
            c = x.comorbidity == 1 ? 0.396*1.61 : 0.25*1.61
        elseif  x.age >= 50 && x.age < 60
            h = x.comorbidity == 1 ? comh : 0.05*1.58*1 #0.376
            c = x.comorbidity == 1 ? 0.396*1.58 : 0.25*1.58
        elseif  x.age >= 60 && x.age < 70
            h = x.comorbidity == 1 ? comh : 0.05*1.65*1 #0.376
            c = x.comorbidity == 1 ? 0.396*1.65 : 0.25*1.65
        elseif  x.age >= 70 && x.age < 80
            h = x.comorbidity == 1 ? comh : 0.05*1.45*1 #0.376
            c = x.comorbidity == 1 ? 0.396*1.45 : 0.25*1.45
        else
            h = x.comorbidity == 1 ? comh : 0.05*1.60*1 #0.376
            c = x.comorbidity == 1 ? 0.396*1.60 : 0.25*1.60
        end
        if x.strain == 4
            h = h*2.26
        end
    else
        error("no strain in movetoinf")
        
    end
    
    groups = [0:34,35:54,55:69,70:84,85:100]
    gg = findfirst(y-> x.age in y,groups)

    mh = [0.0002; 0.0015; 0.011; 0.0802; 0.381] # death rate for severe cases.
   
    ###prop/(prob de sintoma severo)
    if p.calibration && !p.calibration2
        h =  0#, 0, 0, 0)
        c =  0#, 0, 0, 0)
        mh = (0, 0, 0, 0, 0)
    end

    time_to_hospital = Int(round(rand(Uniform(2, 5)))) # duration symptom onset to hospitalization
   	
    x.health = x.swap
    x.health_status = x.swap_status
    x.swap = UNDEF
    
    x.tis = 0 
    if rand() < h     # going to hospital or ICU but will spend delta time transmissing the disease with full contacts 
        x.exp = time_to_hospital
        if rand() < c
            aux_v = [ICU;ICU2;ICU3;ICU4;ICU5;ICU6]
            x.swap = aux_v[x.strain]
            x.swap_status = ICU
            #x.swap = x.strain == 1 ? ICU : ICU2
        else
            aux_v = [HOS;HOS2;HOS3;HOS4;HOS5;HOS6]
            x.swap = aux_v[x.strain]
            x.swap_status = HOS
            #x.swap = x.strain == 1 ? HOS : HOS2
        end
       
    else ## no hospital for this lucky (but severe) individual 
        aux = (p.mortality_inc^Int(x.strain==2 || x.strain == 4))
        aux = x.strain == 4 ? aux*7.0 : aux
        if x.iso || rand() < p.fsevere 
            x.exp = 1  ## 1 day isolation for severe cases 
            aux_v = [IISO;IISO2;IISO3;IISO4;IISO5;IISO6]
            x.swap = aux_v[x.strain]
            x.swap_status = IISO
            #x.swap = x.strain == 1 ? IISO : IISO2
        else
            if rand() < mh[gg]*aux
                x.exp = x.dur[4] 
                aux_v = [DED;DED2;DED3;DED4;DED5;DED6]
                x.swap = aux_v[x.strain]
                x.swap_status = DED
            else 
                x.exp = x.dur[4]  
                aux_v = [REC;REC2;REC3;REC4;REC5;REC6]
                x.swap = aux_v[x.strain]
                x.swap_status = REC
            end

        end  
       
    end
    ## before returning, check if swap is set 
    x.swap == UNDEF && error("agent I -> ?")
end

function move_to_iiso(x::Human)
    ## transfers human h to the sever isolated infection stage for γ days
    x.health = x.swap
    x.health_status = x.swap_status
    groups = [0:34,35:54,55:69,70:84,85:100]
    gg = findfirst(y-> x.age in y,groups)
    
    mh = [0.0002; 0.0015; 0.011; 0.0802; 0.381] # death rate for severe cases.
    aux = (p.mortality_inc^Int(x.strain==2 || x.strain == 4))
    aux = x.strain == 4 ? aux*7.0 : aux

    if rand() < mh[gg]*aux
        x.exp = x.dur[4] 
        aux_v = [DED;DED2;DED3;DED4;DED5;DED6]
        x.swap = aux_v[x.strain]
        x.swap_status = DED
    else 
        x.exp = x.dur[4]  
        aux_v = [REC;REC2;REC3;REC4;REC5;REC6]
        x.swap = aux_v[x.strain]
        x.swap_status = REC
    end
    #x.swap = x.strain == 1 ? REC : REC2
    x.tis = 0     ## reset time in state 
    x.exp = x.dur[4] - 1  ## since 1 day was spent as infectious
    _set_isolation(x, true, :mi)
end 

function move_to_hospicu(x::Human)   
    #death prob taken from https://www.cdc.gov/nchs/nvss/vsrr/covid_weekly/index.htm#Comorbidities
    # on May 31th, 2020
    #= age_thres = [24;34;44;54;64;74;84;999]
    g = findfirst(y-> y >= x.age,age_thres) =#
    aux = [0:4, 5:19, 20:44, 45:54, 55:64, 65:74, 75:84, 85:99]
   
    if x.strain == 1 || x.strain == 3 || x.strain == 5 || x.strain == 6

        mh = [0.001, 0.001, 0.0015, 0.0065, 0.01, 0.02, 0.0735, 0.38]
        mc = [0.002,0.002,0.0022, 0.008, 0.022, 0.04, 0.08, 0.4]

    elseif x.strain == 2  || x.strain == 4
    
        mh = 0.5*[0.0016, 0.0016, 0.0025, 0.0107, 0.02, 0.038, 0.15, 0.66]
        mc = 0.5*[0.0033, 0.0033, 0.0036, 0.0131, 0.022, 0.04, 0.2, 0.70]
        
        if x.strain == 4
            mh = 7*mh
            mc = 7*mc
        end

    else
      
            error("No strain - hospicu")
    end
    
    gg = findfirst(y-> x.age in y,aux)

    psiH = Int(round(rand(Distributions.truncated(Gamma(4.5, 2.75), 8, 17))))
    psiC = Int(round(rand(Distributions.truncated(Gamma(4.5, 2.75), 8, 17)))) + 2
    muH = Int(round(rand(Distributions.truncated(Gamma(5.3, 2.1), 9, 15))))
    muC = Int(round(rand(Distributions.truncated(Gamma(5.3, 2.1), 9, 15)))) + 2

    swaphealth = x.swap_status 
    x.health = x.swap ## swap either to HOS or ICU
    x.health_status = x.swap_status
    x.swap = UNDEF
    x.tis = 0
    _set_isolation(x, true) # do not set the isovia property here.  

    if swaphealth == HOS
        x.hospicu = 1 
        if rand() < mh[gg] ## person will die in the hospital 
            x.exp = muH 
            aux_v = [DED;DED2;DED3;DED4;DED5;DED6]
            x.swap = aux_v[x.strain]
            x.swap_status = DED
            #x.swap = x.strain == 1 ? DED : DED2
        else 
            x.exp = psiH 
            aux_v = [REC;REC2;REC3;REC4;REC5;REC6]
            x.swap = aux_v[x.strain]
            x.swap_status = REC
            #x.swap = x.strain == 1 ? REC : REC2
        end    
    elseif swaphealth == ICU
        x.hospicu = 2 
                
        if rand() < mc[gg] ## person will die in the ICU 
            x.exp = muC
            aux_v = [DED;DED2;DED3;DED4;DED5;DED6]
            x.swap = aux_v[x.strain]
            x.swap_status = DED
            #x.swap = x.strain == 1 ? DED : DED2
        else 
            x.exp = psiC
            aux_v = [REC;REC2;REC3;REC4;REC5;REC6]
            x.swap = aux_v[x.strain]
            x.swap_status = REC
            #x.swap = x.strain == 1 ? REC : REC2
        end
    else
        error("error in hosp")
    end
    
    ## before returning, check if swap is set 
    x.swap == UNDEF && error("agent H -> ?")    
end

function move_to_dead(h::Human)
    # no level of alchemy will bring someone back to life. 
    h.health = h.swap
    h.health_status = h.swap_status
    h.swap = UNDEF
    h.swap_status = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely
    h.iso = true # a dead person is isolated
    _set_isolation(h, true)  # do not set the isovia property here.  
    # isolation property has no effect in contact dynamics anyways (unless x == SUS)
end

function move_to_recovered(h::Human)
    h.health = h.swap
    h.health_status = h.swap_status

    if h.strain in (1,2)
        h.recovered = true
    end

    h.swap = UNDEF
    h.swap_status = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely
    h.iso = false ## a recovered person has ability to meet others
    _set_isolation(h, false)  # do not set the isovia property here.  
    # isolation property has no effect in contact dynamics anyways (unless x == SUS)
end


@inline function _get_betavalue(sys_time, xhealth) 
    #bf = p.β ## baseline PRE
    #length(BETAS) == 0 && return 0
    bf = p.β#BETAS[sys_time]
    # values coming from FRASER Figure 2... relative tranmissibilities of different stages.
    if xhealth == ASYMP
        bf = bf * p.frelasymp #0.11

    elseif xhealth == MILD || xhealth == MISO 
        bf = bf * 0.44

    elseif xhealth == INF || xhealth == IISO 
        bf = bf * 0.89

    elseif xhealth == ASYMP2
        bf = bf*p.frelasymp*p.sec_strain_trans #0.11

    elseif xhealth == MILD2 || xhealth == MISO2
        bf = bf * 0.44*p.sec_strain_trans

    elseif xhealth == INF2 || xhealth == IISO2 
        bf = bf * 0.89*p.sec_strain_trans

    elseif xhealth == PRE2
        bf = bf*p.sec_strain_trans

    elseif xhealth == ASYMP3
        bf = bf*p.frelasymp*p.third_strain_trans #0.11

    elseif xhealth == MILD3 || xhealth == MISO3
        bf = bf * 0.44*p.third_strain_trans

    elseif xhealth == INF3 || xhealth == IISO3 
        bf = bf * 0.89*p.third_strain_trans

    elseif xhealth == PRE3
        bf = bf*p.third_strain_trans

    elseif xhealth == ASYMP4
        bf = bf*p.frelasymp*p.sec_strain_trans*p.fourth_strain_trans #0.11

    elseif xhealth == MILD4 || xhealth == MISO4
        bf = bf * 0.44*p.sec_strain_trans*p.fourth_strain_trans

    elseif xhealth == INF4 || xhealth == IISO4
        bf = bf * 0.89*p.sec_strain_trans*p.fourth_strain_trans

    elseif xhealth == PRE4
        bf = bf*p.sec_strain_trans*p.fourth_strain_trans
    ############### 5 strain    
    elseif xhealth == ASYMP5
        bf = bf*p.frelasymp*p.fifth_strain_trans #0.11

    elseif xhealth == MILD5 || xhealth == MISO5
        bf = bf * 0.44*p.fifth_strain_trans

    elseif xhealth == INF5 || xhealth == IISO5
        bf = bf * 0.89*p.fifth_strain_trans

    elseif xhealth == PRE5
        bf = bf*p.fifth_strain_trans

    elseif xhealth == ASYMP6
        bf = bf*p.frelasymp*p.sixth_strain_trans #0.11

    elseif xhealth == MILD6 || xhealth == MISO6
        bf = bf * 0.44*p.sixth_strain_trans

    elseif xhealth == INF6 || xhealth == IISO6
        bf = bf * 0.89*p.sixth_strain_trans

    elseif xhealth == PRE6
        bf = bf*p.sixth_strain_trans
    end
    return bf
end
export _get_betavalue

@inline function get_nextday_counts(x::Human)
    # get all people to meet and their daily contacts to recieve
    # we can sample this at the start of the simulation to avoid everyday    
    cnt = 0
    ag = x.ag
    #if person is isolated, they can recieve only 3 maximum contacts
    
    if !x.iso 
        #cnt = rand() < 0.5 ? 0 : rand(1:3)
        aux = x.relaxed ? 1.0 : p.contact_change_rate*p.contact_change_2
        cnt = rand(negative_binomials(ag,aux)) ##using the contact average for shelter-in
    else 
        cnt = rand(negative_binomials_shelter(ag,p.contact_change_2))  # expensive operation, try to optimize
    end
    
    if x.health in (DED,DED2,DED3)
        cnt = 0 
    end
    x.nextday_meetcnt = cnt
    return cnt
end

function dyntrans(sys_time, grps,sim)
    totalmet = 0 # count the total number of contacts (total for day, for all INF contacts)
    totalinf = 0 # count number of new infected 
    ## find all the people infectious
    #rng = MersenneTwister(246*sys_time*sim)
    pos = shuffle(1:length(humans))
    # go through every infectious person
    for x in humans[pos]        
        if x.health_status in (PRE, ASYMP, MILD, MISO, INF, IISO)
            
            xhealth = x.health
            cnts = x.nextday_meetcnt
            cnts == 0 && continue # skip person if no contacts
            
            gpw = Int.(round.(cm[x.ag]*cnts)) # split the counts over age groups
            for (i, g) in enumerate(gpw) 
                meet = rand(grps[i], g)   # sample the people from each group
                # go through each person
                for j in meet 
                    y = humans[j]
                    ycnt = y.nextday_meetcnt    
                    ycnt == 0 && continue

                    y.nextday_meetcnt = y.nextday_meetcnt - 1 # remove a contact
                    totalmet += 1
                    
                    beta = _get_betavalue(sys_time, xhealth)
                    adj_beta = 0 # adjusted beta value by strain and vaccine efficacy
                    if y.health == SUS && y.swap == UNDEF
                        aux = y.vac_status*y.protected > 0 ? p.vac_efficacy_inf[x.strain][y.vac_status][y.protected] : 0.0
                        adj_beta = beta*(1-aux)
                    elseif (x.strain in (3,4,6) && y.health in (REC, REC2, REC5) && y.swap == UNDEF)
                        adj_beta = beta*(p.reduction_recovered) #0.21
                    end

                    if rand() < adj_beta
                        totalinf += 1
                        y.exp = y.tis   ## force the move to latent in the next time step.
                        y.sickfrom = xhealth ## stores the infector's status to the infectee's sickfrom
                        y.sickby = x.idx
                        y.strain = x.strain       
                        aux_v = [LAT;LAT2;LAT3;LAT4;LAT5;LAT6]
                        y.swap = aux_v[y.strain]
                        y.swap_status = LAT
                        #y.swap = y.strain == 1 ? LAT : LAT2
                    end  
                end
            end            
        end
    end
    return totalmet, totalinf
end
export dyntrans

function contact_matrix()
    # regular contacts, just with 5 age groups. 
    #  0-4, 5-19, 20-49, 50-64, 65+
    CM = Array{Array{Float64, 1}, 1}(undef, 5)
     CM[1] = [0.2287, 0.1839, 0.4219, 0.1116, 0.0539]
    CM[2] = [0.0276, 0.5964, 0.2878, 0.0591, 0.0291]
    CM[3] = [0.0376, 0.1454, 0.6253, 0.1423, 0.0494]
    CM[4] = [0.0242, 0.1094, 0.4867, 0.2723, 0.1074]
    CM[5] = [0.0207, 0.1083, 0.4071, 0.2193, 0.2446] 
   
    return CM
end
# 
# calibrate for 2.7 r0
# 20% selfisolation, tau 1 and 2.

function negative_binomials(ag,mult) 
    ## the means/sd here are calculated using _calc_avgag
    means = [10.21, 16.793, 13.7950, 11.2669, 8.0027]
    sd = [7.65, 11.7201, 10.5045, 9.5935, 6.9638]
    means = means*mult
    totalbraks = length(means)
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, totalbraks)
    for i = 1:totalbraks
        p = 1 - (sd[i]^2-means[i])/(sd[i]^2)
        r = means[i]^2/(sd[i]^2-means[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms[ag]
end
#const nbs = negative_binomials()
const cm = contact_matrix()
#export negative_binomials, contact_matrix, nbs, cm

export negative_binomials


function negative_binomials_shelter(ag,mult) 
    ## the means/sd here are calculated using _calc_avgag
    means = [2.86, 4.7, 3.86, 3.15, 2.24]
    sd = [2.14, 3.28, 2.94, 2.66, 1.95]
    means = means*mult
    #sd = sd*mult
    totalbraks = length(means)
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, totalbraks)
    for i = 1:totalbraks
        p = 1 - (sd[i]^2-means[i])/(sd[i]^2)
        r = means[i]^2/(sd[i]^2-means[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms[ag]   
end
## internal functions to do intermediate calculations
function _calc_avgag(lb, hb) 
    ## internal function to calculate the mean/sd of the negative binomials
    ## returns a vector of sampled number of contacts between age group lb to age group hb
    dists = _negative_binomials_15ag()[lb:hb]
    totalcon = Vector{Int64}(undef, 0)
    for d in dists 
        append!(totalcon, rand(d, 10000))
    end    
    return totalcon
end
export _calc_avgag

function _negative_binomials_15ag()
    ## negative binomials 15 agegroups
    AgeMean = Vector{Float64}(undef, 15)
    AgeSD = Vector{Float64}(undef, 15)
    #0-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55-59, 60-64, 65-69, 70+
    #= AgeMean = [10.21, 14.81, 18.22, 17.58, 13.57, 13.57, 14.14, 14.14, 13.83, 13.83, 12.3, 12.3, 9.21, 9.21, 6.89]
    AgeSD = [7.65, 10.09, 12.27, 12.03, 10.6, 10.6, 10.15, 10.15, 10.86, 10.86, 10.23, 10.23, 7.96, 7.96, 5.83]
     =#
     AgeMean = repeat([14.14],15)#[10.21, 14.81, 18.22, 17.58, 13.57, 13.57, 14.14, 14.14, 13.83, 13.83, 12.3, 12.3, 9.21, 9.21, 6.89]
    AgeSD = repeat([10.86],15)#[7.65, 10.09, 12.27, 12.03, 10.6, 10.6, 10.15, 10.15, 10.86, 10.86, 10.23, 10.23, 7.96, 7.96, 5.83]
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, 15)
    for i = 1:15
        p = 1 - (AgeSD[i]^2-AgeMean[i])/(AgeSD[i]^2)
        r = AgeMean[i]^2/(AgeSD[i]^2-AgeMean[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms    
end

#const vaccination_days = days_vac_f()
#const vac_rate_1 = vaccination_rate_1()
#const vac_rate_2 = vaccination_rate_2()
## references: 
# critical care capacity in Canada https://www.ncbi.nlm.nih.gov/pubmed/25888116
end # module end
