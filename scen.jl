#### not re-calibrated yet
dc = [1]#;map(y->125+y,0:23);map(y->170+y,0:19);map(y->222+y,0:29);map(y->301+y,0:10)]
rc = [1.0]#;map(y->1.0-(0.1/24)*y,1:24);map(y->0.9+(0.12/20)*y,1:20);map(y->1.02-(0.2/30)*y,1:30);map(y->0.82+(0.43/11)*y,1:11)]#;map(y->0.901+(0.03/5)*y,1:5)]
run_param_scen_cal(0.089,"newyork",20,7,1,999,0.2,1,0.5,0.02,0.75,rc,dc,350,true)

#proportion,dw,efrate,rc=[1.0],dc=[1],mt::Int64=500,vac::Bool=true,scen::String="statuscuo",alpha::Float64 = 1.0,alpha2::Float64 = 0.0,alpha3::Float64 = 1.0,nsims::Int64=500)
    