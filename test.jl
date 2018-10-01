using JuMP,PowerModels,DataFrames,CSV,CPLEX
case = PowerModels.parse_file("case2869pegase.m")
pm = build_generic_model(case, ACPPowerModel,
PowerModels.post_opf)
m = JuMP.Model(solver = CplexSolver())
ref = pm.ref[:nw][0]
T = 24
gens = [gen for gen in keys(ref[:gen]) if ref[:gen][gen]["pmin"] >= 0]
aggr = [gen for gen in keys(ref[:gen]) if ref[:gen][gen]["pmin"] < 0]
price = @variable(m,[1:T],lowerbound=0.5,upperbound=1.5)
pw = @variable(m,[x in gens,1:T])
pv = @variable(m,[x in gens],lowerbound=0)
pub = @variable(m,[x in gens,1:T])
plb = @variable(m,[x in gens,1:T])
st = @variable(m,[x in gens,1:T],Bin)
st_z =@variable(m,[x in gens,1:T],Bin)
st_y =@variable(m,[x in gens,1:T],Bin)
load = sum([ref[:load][x]["pd"] for x in keys(ref[:load])])
srand(1234)
# load_real = load + 0.2*rand(Float64,T)*load
load_data = CSV.read("load data.csv")[1:24,:]
load_real = (load / maximum([get(x) for x in load_data[:load]]))*[get(x) for x in load_data[:load]]
wind_data = CSV.read("1.csv")
ramp = 0.5
r = 0.25 * load/maximum([get(x) for x in wind_data[:max]])
windmax = r * [get(x) for x in wind_data[:max]]
windmid = r * [get(x) for x in wind_data[:mid]]
windmin = r * [get(x) for x in wind_data[:min]]
pmax = [ref[:gen][x]["pmax"] for x in aggr]
@constraint(m,sum(pv[gen] for gen in gens) == 1)
for t in 1:T
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
    @constraint(m,sum(pw[gen,t] for gen in gens) + (2 - price[t]) * sum([ref[:gen][x]["pg"] for x in aggr]) + windmid[t]== load_real[t])
end
startup_cost = sum(ref[:gen][gen]["startup"]*st_y[gen,t] for gen in gens,t in 1:T-1)
fuel_cost = sum(ref[:gen][gen]["cost"][2]*pw[gen,t] + ref[:gen][gen]["cost"][1]*pw[gen,t]^2 for gen in gens,t in 1:T)
run_cost = sum((ref[:gen][gen]["cost"][3]+10)*st[gen,t] for gen in gens,t in 1:T)
sell_cost =   sum(price[t]*(2 - price[t]) * sum([ref[:gen][x]["pg"]*ref[:gen][x]["cost"][2] for x in aggr]) for t in 1:T)
interval = sum(pub[gen,t] - plb[gen,t] for gen in gens,t in 1:T)
@objective(m,Min,startup_cost + fuel_cost + run_cost + interval)
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
wind[:min] = windmin
wind[:mean] = windmid
wind[:max] = windmax
wind[:price] = getvalue(price)
CSV.write("agc.csv",agc)
CSV.write("nonagc.csv",nonagc)
CSV.write("wind.csv",wind)
