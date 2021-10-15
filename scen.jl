



dc = [1]
rc = [1.0]
run_param_scen_cal(0.15,"newyork",20,7,1,999,0.2,1,0.5,0.02,0.75,rc,dc,350,true)






rw = 0.0:0.004:0.04
re = 0.5:0.05:1.0

for i in rw, j in re
    run_param_scen_cal(0.15,"newyork",20,7,1,999,0.2,1,0.5,i,j,[1.0],[1],350,true)
end
