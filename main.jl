using JuMP,PowerModels,DataFrames,CSV,CPLEX,Random
case = PowerModels.parse_file("testcase.m")
pm = build_generic_model(case, ACPPowerModel,
PowerModels.post_opf)
m = JuMP.Model(solver = CplexSolver())
ref = pm.ref[:nw][0]
for gen in keys(ref[:gen])
    if ref[:gen][gen]["pmin"] == 0
        ref[:gen][gen]["pmin"] = 0.1 * ref[:gen][gen]["pmax"]
    end
    if ref[:gen][gen]["pmax"] == 0
        delete!(ref[:gen],gen)
    end
end
T = 24
aggr = [gen for gen in keys(ref[:gen]) if ref[:gen][gen]["pmin"] < 0]
aggr = [5,6,7]
aggr_v = @variable(m,[x in aggr])
gens = setdiff([gen for gen in keys(ref[:gen])],aggr)
pmax = maximum([ref[:gen][gen]["pmax"] for gen in keys(ref[:gen])])
hydro = Dict("pmax_pump"=>pmax/6,"pmin_pump"=>0.1,"pmax_gen"=>pmax/6,"pmin_gen"=>0.1,"e2w"=>2.5,"w2e"=>1.6,"C_max"=>154.900,"C_min"=>37.820,"C_int"=>97.820)
ref[:hydro] = Dict([(i,hydro) for i in range(1,6)])
hydros = [gen for gen in keys(ref[:hydro])]
price = @variable(m,[1:T],lowerbound=0.5,upperbound=1.5)
pw = @variable(m,[x in gens,1:T])
hydro_gen = @variable(m,[x in hydros,1:T])
hydro_pump = @variable(m,[x in hydros,1:T])
hydro_v = @variable(m,[x in hydros])
hydro_W = @variable(m,[x in hydros,1:T])
u_gen = @variable(m,[x in hydros,1:T])
u_pump = @variable(m,[x in hydros,1:T])
pv = @variable(m,[x in gens],lowerbound=0)
pub = @variable(m,[x in gens,1:T])
plb = @variable(m,[x in gens,1:T])
st = @variable(m,[x in gens,1:T],Bin)
st_z =@variable(m,[x in gens,1:T],Bin)
st_y =@variable(m,[x in gens,1:T],Bin)
load = sum([ref[:load][x]["pd"] for x in keys(ref[:load])])
Random.seed!(1234)
# load_real = load + 0.2*rand(Float64,T)*load
load_data = CSV.read("load data.csv")[1:24,:]
load_real = load_data[:load]/100

for confidence in [50,70,90] #风机出力置信区间循环
    wind_data = CSV.read(string(confidence,"wind.csv"))
    ramp = 0.5
    r = 0.6 * load/maximum(wind_data[:max])
    windmax = r * wind_data[:max]
    windmid = r * wind_data[:mid]
    windmin = r * wind_data[:min]
    pmax = [ref[:gen][x]["pmax"] for x in aggr]
    @constraint(m,sum(pv[gen] for gen in gens) + sum(hydro_v[gen] for gen in hydros) + sum(aggr_v[gen] for gen in aggr)== 1)
    for t in 1:T
        for gen in hydros
            @constraint(m,hydro_gen[gen,t] >= u_gen[gen,t] * ref[:hydro][gen]["pmin_gen"])
            @constraint(m,hydro_gen[gen,t] <= u_gen[gen,t] * ref[:hydro][gen]["pmax_gen"])
            @constraint(m,hydro_pump[gen,t] >= u_pump[gen,t] * ref[:hydro][gen]["pmin_pump"])
            @constraint(m,hydro_pump[gen,t] <= u_pump[gen,t] * ref[:hydro][gen]["pmax_pump"])
            @constraint(m,u_pump[gen,t] + u_gen[gen,t] <= 1)
            @constraint(m,hydro_gen[gen,t] - hydro_pump[gen,t] - (windmax[t] - windmid[t])*hydro_v[gen] >=  - u_pump[gen,t] * ref[:hydro][gen]["pmax_pump"])
            @constraint(m,hydro_gen[gen,t] - hydro_pump[gen,t] - (windmax[t] - windmid[t])*hydro_v[gen] <=  u_gen[gen,t] * ref[:hydro][gen]["pmax_gen"])
            if t >=2
                @constraint(m,(windmax[t] - windmid[t])*hydro_v[gen] - (windmin[t-1] - windmid[t-1])*hydro_v[gen] <= 0.2 * ref[:hydro][gen]["pmax_gen"])
                @constraint(m,(windmax[t] - windmid[t])*hydro_v[gen] - (windmax[t-1] - windmid[t-1])*hydro_v[gen] >= -0.2 * ref[:hydro][gen]["pmax_gen"])
            end
            if t == 1
                @constraint(m,hydro_W[gen,t] == ref[:hydro][gen]["C_int"])
            else
                @constraint(m,hydro_W[gen,t] == hydro_W[gen,t-1] + ref[:hydro][gen]["e2w"]*hydro_pump[gen,t] - ref[:hydro][gen]["e2w"]*hydro_gen[gen,t])
            end
            @constraint(m,hydro_W[gen,t]<=ref[:hydro][gen]["C_max"])
            @constraint(m,hydro_W[gen,t]>=ref[:hydro][gen]["C_min"])
            if t == 24
                @constraint(m,hydro_W[gen,t]==ref[:hydro][gen]["C_int"])
            end
        end
        for gen in gens
            @constraint(m,pub[gen,t] <= st[gen,t]*ref[:gen][gen]["pmax"])
            @constraint(m,plb[gen,t] >= st[gen,t]*ref[:gen][gen]["pmin"])
            @constraint(m,pw[gen,t] - (windmin[t] - windmid[t])*pv[gen]  <= pub[gen,t])
            @constraint(m,pw[gen,t] - (windmax[t] - windmid[t])*pv[gen] >= plb[gen,t])
            if t >=2
                @constraint(m,st[gen,t] == st[gen,t-1] + st_y[gen,t] - st_z[gen,t])
                @constraint(m,st_y[gen,t] + st_z[gen,t] <= 1)
                @constraint(m,pub[gen,t] - plb[gen,t-1] <= ramp*ref[:gen][gen]["pmax"]) # add adjustable var
                @constraint(m,plb[gen,t-1] - pub[gen,t] <= ramp*ref[:gen][gen]["pmax"]) # add adjustable var
                @constraint(m,sum(st_y[gen,k] for k in max(t-4,1):t)<=st[gen,t])
                @constraint(m,sum(st_z[gen,k] for k in max(t-4,1):t)<=1 - st[gen,t])
            else
                @constraint(m,st[gen,t] == st_y[gen,t] - st_z[gen,t])
                @constraint(m,st_y[gen,t] + st_z[gen,t] <= 1)
            end
        end
        for x in aggr
            @constraint(m, ref[:gen][x]["pmax"] - (windmin[t] - windmid[t])*aggr_v[x]  <= 1.2 * ref[:gen][x]["pmax"])
            @constraint(m, ref[:gen][x]["pmax"] - (windmax[t] - windmid[t])*aggr_v[x] >= 0.5 * ref[:gen][x]["pmax"])
        end
        if length(aggr) > 0
            @constraint(m,sum(pw[gen,t] for gen in gens)
            - (2 - price[t]) * sum([ref[:gen][x]["pmax"] for x in aggr])
            + windmid[t]
            + sum(hydro_gen[gen,t] for gen in hydros) == load_real[t] + sum(hydro_pump[gen,t] for gen in hydros))
        else
            @constraint(m,sum(pw[gen,t] for gen in gens)
            + windmid[t]
            + sum(hydro_gen[gen,t] for gen in hydros)== load_real[t] + sum(hydro_pump[gen,t] for gen in hydros))
        end
    end
    startup_cost = sum(ref[:gen][gen]["startup"]*st_y[gen,t] for gen in gens,t in 1:T-1)
    fuel_cost = sum(ref[:gen][gen]["cost"][2]*pw[gen,t] + ref[:gen][gen]["cost"][1]*pw[gen,t]^2 for gen in gens,t in 1:T)
    run_cost = sum(ref[:gen][gen]["cost"][3] * st[gen,t] for gen in gens,t in 1:T) # May slow down the program
    interval = 100*sum(pub[gen,t] - plb[gen,t] for gen in gens,t in 1:T)
    CO2Emission = sum(900+(100000/ref[:gen][gen]["cost"][1])*pw[gen,t] for gen in gens,t in 1:T)
    # if length(aggr) > 0
    #     sell_cost =   1000 * sum(price[t]*(2 - price[t]) * sum([ref[:gen][x]["pg"] for x in aggr]) for t in 1:T)
    #     # @objective(m,Min,startup_cost + fuel_cost  + interval + run_cost)
    #     # @objective(m,Min,startup_cost + fuel_cost  + interval + run_cost + CO2Emission)
    #     # @objective(m,Min,startup_cost + fuel_cost  + interval + run_cost + 2*CO2Emission)
    #     @objective(m,Min,startup_cost + fuel_cost  + interval + run_cost + 5*CO2Emission)
    #     # @objective(m,Min,CO2Emission)
    # else
    #     # @objective(m,Min,startup_cost + fuel_cost + interval + run_cost)
    #     # @objective(m,Min,startup_cost + fuel_cost  + interval + run_cost + CO2Emission)
    #     # @objective(m,Min,startup_cost + fuel_cost  + interval + run_cost + 2*CO2Emission)
    #     @objective(m,Min,startup_cost + fuel_cost  + interval + run_cost +  5*CO2Emission)
    #     # @objective(m,Min,CO2Emission)
    # end

    for k in [0,1,2,3,4,6,8,12,16,20]#多目标权重循环
        @objective(m,Min,startup_cost + fuel_cost  + interval + run_cost +  k*CO2Emission)
        status = solve(m)
        result = Dict()
        result["power"] = getvalue(pw)
        result["on_off"] = getvalue(st)
        result["power_ub"] = getvalue(pub)
        result["power_lb"] = getvalue(plb)
        result["factor"] = getvalue(pv)
        nonagc = DataFrame()
        agc = DataFrame()
        wind = DataFrame()
        for gen in gens
            if getvalue(pv[gen]) == 0
                nonagc[Symbol("gen",gen," state")] = getvalue(st[gen,:])
                # nonagc[Symbol("gen",gen," on")] = getvalue(st_y[gen,:])
                # nonagc[Symbol("gen",gen," off")] = getvalue(st_z[gen,:])
                nonagc[Symbol("gen",gen," power")] = getvalue(pw[gen,:])
            else
                agc[Symbol("gen",gen," state")] = getvalue(st[gen,:])
                agc[Symbol("gen",gen," power")] = getvalue(pw[gen,:])
                agc[Symbol("gen",gen," upperbound")] = getvalue(pub[gen,:])
                agc[Symbol("gen",gen," lowerbound")] = getvalue(plb[gen,:])
            end
        end
        for gen in hydros
            agc[Symbol("hydro",gen,"power")] = getvalue(hydro_gen[gen,:]) -  getvalue(hydro_pump[gen,:])
            agc[Symbol("hydro",gen," upperbound")] = getvalue(hydro_gen[gen,:]) -  getvalue(hydro_pump[gen,:]) - (windmin - windmid)*getvalue(hydro_v[gen])
            agc[Symbol("hydro",gen," lowerbound")] = getvalue(hydro_gen[gen,:]) -  getvalue(hydro_pump[gen,:]) - (windmax - windmid)*getvalue(hydro_v[gen])
        end
        wind[:min] = windmin
        wind[:mean] = windmid
        wind[:max] = windmax
        wind[:price] = getvalue(price)
        dir = mkdir(string(confidence,"wind_","weight",k))
        CSV.write(string(dir,"/","agc.csv"),agc)
        CSV.write(string(dir,"/","nonagc.csv"),nonagc)
        CSV.write(string(dir,"/","wind.csv"),wind)
        cost = getvalue(startup_cost + fuel_cost + run_cost)
        CO2 = getvalue(CO2Emission)
        print(string("运行成本：","$cost"))
        print(string("CO2排放量：","$CO2"))
